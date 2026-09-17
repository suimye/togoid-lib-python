"""
TogoID CLI - Unified command-line interface

Provides convert, annotate, and config subcommands.
"""

import argparse
import json
import sys
import os
from typing import Optional, List

from .converter import TogoIDConverter
from .annotations import AnnotationsConverter, parse_filters, load_ids, ensure_fields, output_table, output_json
from .label_converter import LabelConverter, parse_labels, output_results
from .enrichment.cli import add_enrich_parser, handle_enrich


def create_parser() -> argparse.ArgumentParser:
    """Create main argument parser with subcommands"""
    parser = argparse.ArgumentParser(
        prog='togoid',
        description='TogoID - Biological database ID conversion and annotation tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert IDs
  togoid convert --ids 1,9 --route ncbigene,ensembl_gene

  # Convert labels to IDs (SPARQList)
  togoid label2id --labels "BRCA1,TP53" --dataset ncbigene --taxonomy 9606

  # Convert labels to IDs (PubDictionaries)
  togoid label2id --labels "caffeine" --dataset chebi --label_types "togoid_chebi_label"

  # Get annotations
  togoid annotate --dataset ncbigene --ids 672,7157 --field gene_synonym

  # Get configuration
  togoid config dataset ncbigene

  # Enrichment analysis of a gene list
  togoid enrich --genes "CD3D,CD3E,LCK,ZAP70" --preset reactome
        """
    )

    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    parser.add_argument('--api-url', type=str, help='TogoID API base URL')

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # ========== CONVERT subcommand ==========
    convert_parser = subparsers.add_parser('convert', help='Convert IDs between databases')
    convert_parser.add_argument('--route', required=True, help='Comma-separated route (e.g., "refseq_rna,uniprot")')
    convert_parser.add_argument('--ids', required=True, help='Comma-separated IDs to convert')
    convert_parser.add_argument('--report', choices=['target', 'all', 'pair', 'full'], default='target',
                                help='Report format (default: target)')
    convert_parser.add_argument('--format', choices=['json', 'csv', 'tsv', 'dict', 'table', 'dataframe'],
                                default='json', help='Output format (default: json)')
    convert_parser.add_argument('--limit', type=int, help='Limit number of results')
    convert_parser.add_argument('--offset', type=int, help='Offset for results')
    convert_parser.add_argument('--annotate', action='append', nargs=2, metavar=('DATASET', 'FIELD'),
                                help='Add annotation column (dataset field). Can be used multiple times.')
    convert_parser.add_argument('--filter', action='append', nargs=3, metavar=('DATASET', 'FIELD', 'VALUES'),
                                help='Filter by annotation value (dataset field values). Values should be comma-separated.')
    convert_parser.add_argument('--output', help='Output file path')
    convert_parser.add_argument('--raw', action='store_true',
                                help='Return raw IDs without dataset prefixes (sends prefix=no). '
                                     'By default the API formats IDs as CURIEs, e.g. GO:0005634.')

    # ========== LABEL2ID subcommand ==========
    label2id_parser = subparsers.add_parser('label2id', help='Convert labels to IDs')
    label_input = label2id_parser.add_mutually_exclusive_group(required=True)
    label_input.add_argument('--labels', help='Comma-separated labels')
    label_input.add_argument('--label-file', help='File containing labels (one per line or comma-separated)')

    # Dataset argument (required)
    label2id_parser.add_argument('--dataset', required=True, help='Dataset name (e.g., ncbigene)')

    # Label types argument (used for both SPARQList and PubDictionaries)
    label2id_parser.add_argument('--label_types',
                                 help='Label types or dictionary names (comma-separated, uses dataset config if not specified)')

    # SPARQList specific options
    label2id_parser.add_argument('--taxonomy', help='Taxonomy ID for SPARQList (e.g., 9606 for human)')

    # PubDictionaries specific options
    label2id_parser.add_argument('--tags', help='Taxonomy tags for PubDictionaries (e.g., 9606)')
    label2id_parser.add_argument('--threshold', type=float, default=0.5,
                                 help='Matching score threshold for PubDictionaries (0-1, default: 0.5)')
    label2id_parser.add_argument('--preferred-dictionary', help='Preferred dictionary for synonym resolution')

    # Output options
    label2id_parser.add_argument('--format', choices=['json', 'csv', 'tsv'], default='json',
                                 help='Output format (default: json)')
    label2id_parser.add_argument('--output', help='Output file path')
    label2id_parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')

    # ========== ANNOTATE subcommand ==========
    annotate_parser = subparsers.add_parser('annotate', help='Get annotations for IDs')
    annotate_parser.add_argument('--dataset', required=True, help='Dataset key (e.g. ncbigene)')
    annotate_parser.add_argument('--ids', nargs='*', default=[], help='IDs to annotate')
    annotate_parser.add_argument('--ids-file', help='File containing IDs (one per line)')
    annotate_parser.add_argument('--field', dest='fields', action='append',
                                 help='Annotation field to include (repeat for multiple)')
    annotate_parser.add_argument('--include-label', action='store_true',
                                 help='Include the GraphQL label field')
    annotate_parser.add_argument('--filter', dest='filters', action='append',
                                 help='Filter by annotation values (format: field=value1,value2)')
    annotate_parser.add_argument('--list-fields', action='store_true',
                                 help='List available fields and exit')
    annotate_parser.add_argument('--format', choices=['table', 'csv', 'json'], default='table',
                                 help='Output format (default: table)')
    annotate_parser.add_argument('--delimiter', default='\t', help='Delimiter for table/CSV')
    annotate_parser.add_argument('--compact', action='store_true',
                                 help='Keep list-valued annotations in single cell')
    annotate_parser.add_argument('--no-header', action='store_true', help='Suppress header row')
    annotate_parser.add_argument('--output', '-o', help='Output file path')
    annotate_parser.add_argument('--graphql-endpoint', help='GRASP GraphQL endpoint')
    annotate_parser.add_argument('--timeout', type=float, default=30.0, help='Request timeout (seconds)')
    annotate_parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')

    # ========== SEARCH subcommand ==========
    search_parser = subparsers.add_parser('search', help='Search databases or IDs')
    search_subparsers = search_parser.add_subparsers(dest='search_type', help='Search type')

    search_db_parser = search_subparsers.add_parser('databases', help='Search databases by name')
    search_db_parser.add_argument('name', help='Database name to search for')

    search_id_parser = search_subparsers.add_parser('id', help='Search databases by ID pattern')
    search_id_parser.add_argument('id_string', help='ID string to search for')

    # ========== LOOKUP subcommand ==========
    lookup_parser = subparsers.add_parser('lookup', help='Lookup ID in tables')
    lookup_subparsers = lookup_parser.add_subparsers(dest='lookup_type', help='Lookup type')

    lookup_id_parser = lookup_subparsers.add_parser('id', help='Lookup which tables contain an ID')
    lookup_id_parser.add_argument('id_string', help='ID to lookup')

    # ========== ROUTE subcommand ==========
    route_parser = subparsers.add_parser('route', help='Find routes between databases')
    route_parser.add_argument('src', help='Source database')
    route_parser.add_argument('dst', help='Destination database')
    route_parser.add_argument('--max-hops', type=int, default=3, help='Maximum hops (default: 3)')

    # ========== COUNT subcommand ==========
    count_parser = subparsers.add_parser('count', help='Count mappings between databases')
    count_parser.add_argument('src', help='Source database')
    count_parser.add_argument('dst', help='Destination database')
    count_parser.add_argument('--ids', required=True, help='Comma-separated IDs')
    count_parser.add_argument('--link', help='Link name')

    # ========== CONFIG subcommand ==========
    config_parser = subparsers.add_parser('config', help='Get configuration')
    config_subparsers = config_parser.add_subparsers(dest='config_type', help='Config type')

    config_dataset_parser = config_subparsers.add_parser('dataset', help='Get dataset configuration')
    config_dataset_parser.add_argument('name', nargs='?', help='Dataset name (optional)')

    config_relation_parser = config_subparsers.add_parser('relation', help='Get relation configuration')
    config_relation_parser.add_argument('src', nargs='?', help='Source database (optional)')
    config_relation_parser.add_argument('dst', nargs='?', help='Destination database (optional)')

    config_subparsers.add_parser('descriptions', help='Get database descriptions')
    config_subparsers.add_parser('statistics', help='Get database statistics')
    config_subparsers.add_parser('taxonomy', help='Get taxonomy list')

    # ========== GET-ORTHOLOG subcommand ==========
    ortholog_parser = subparsers.add_parser('get-ortholog', help='Get orthologs via round-trip conversion')
    ortholog_parser.add_argument('--ids', required=True, help='Comma-separated source IDs')
    ortholog_parser.add_argument('--route', required=True, help='Comma-separated route (e.g., "ncbigene,homologene")')
    ortholog_parser.add_argument('--target-taxids', required=True, help='Comma-separated target taxonomy IDs (e.g., "10090,10116")')
    ortholog_parser.add_argument('--format', choices=['json', 'csv', 'tsv', 'dict', 'table', 'dataframe'],
                                 default='table', help='Output format (default: table)')
    ortholog_parser.add_argument('--output', help='Output file path')

    # ========== ENRICH subcommand ==========
    # Defined in togoid.enrichment.cli so the feature stays self-contained.
    add_enrich_parser(subparsers)

    return parser


def handle_convert(args, converter: TogoIDConverter):
    """Handle convert command"""
    route = [r.strip() for r in args.route.split(',')]
    ids = [i.strip() for i in args.ids.split(',')]

    kwargs = {'format': args.format}
    if args.limit is not None:
        kwargs['limit'] = args.limit
    if args.offset is not None:
        kwargs['offset'] = args.offset

    # Handle annotate option
    if args.annotate:
        kwargs['annotate'] = [(dataset, field) for dataset, field in args.annotate]

    # Handle filter option
    if args.filter:
        filter_list = []
        for dataset, field, values in args.filter:
            value_list = [v.strip() for v in values.split(',')]
            filter_list.append((dataset, field, value_list))
        kwargs['filter'] = filter_list

    # Set report parameter (automatically set to 'full' if annotate/filter used)
    if not (args.annotate or args.filter):
        kwargs['report'] = args.report

    # --raw opts out of the API's CURIE prefixing (PR #149)
    if getattr(args, 'raw', False):
        kwargs['prefix'] = False

    result = converter.convert(route, ids, **kwargs)

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            if isinstance(result, str):
                f.write(result)
            elif args.format == 'table':
                for row in result:
                    f.write('\t'.join(str(cell) for cell in row) + '\n')
            else:
                json.dump(result, f, indent=2, ensure_ascii=False)
    else:
        if isinstance(result, str):
            print(result)
        elif args.format == 'table':
            for row in result:
                print('\t'.join(str(cell) for cell in row))
        else:
            print(json.dumps(result, indent=2, ensure_ascii=False))


def handle_label2id(args, api_url: str):
    """Handle label2id command"""
    # Parse labels
    if args.labels:
        labels = parse_labels(args.labels)
    else:
        with open(args.label_file, 'r', encoding='utf-8') as f:
            labels = parse_labels(f.read())

    # Create converter
    converter = LabelConverter(verbose=args.verbose, api_base_url=api_url)

    # Convert labels to IDs
    try:
        # Convert label_types from comma-separated string to list if provided
        label_types = None
        if args.label_types:
            label_types = [t.strip() for t in args.label_types.split(',')]

        results = converter.convert(
            labels=labels,
            dataset=args.dataset,
            label_types=label_types,
            tags=args.tags,
            threshold=args.threshold,
            preferred_dictionary=args.preferred_dictionary,
            taxonomy=args.taxonomy,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Output results
    output_results(results, args.format, args.output)
    return 0


def handle_get_ortholog(args, converter: TogoIDConverter):
    """Handle get-ortholog command"""
    ids = [i.strip() for i in args.ids.split(',')]
    route = [r.strip() for r in args.route.split(',')]
    target_taxids = [t.strip() for t in args.target_taxids.split(',')]

    result = converter.get_ortholog(
        ids=ids,
        route=route,
        target_taxids=target_taxids,
        format=args.format
    )

    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            if isinstance(result, str):
                f.write(result)
            elif args.format == 'table':
                for row in result:
                    f.write('\t'.join(str(cell) for cell in row) + '\n')
            else:
                json.dump(result, f, indent=2, ensure_ascii=False)
    else:
        if isinstance(result, str):
            print(result)
        elif args.format == 'table':
            for row in result:
                print('\t'.join(str(cell) for cell in row))
        else:
            print(json.dumps(result, indent=2, ensure_ascii=False))


def handle_annotate(args, api_url: str):
    """Handle annotate command"""
    grasp_endpoint = args.graphql_endpoint or os.environ.get(
        "TOGOID_GRASP_ENDPOINT", "https://dx.dbcls.jp/grasp-dev-togoid"
    )

    converter = AnnotationsConverter(
        api_endpoint=api_url,
        grasp_endpoint=grasp_endpoint,
        verbose=args.verbose,
        timeout=args.timeout,
    )

    try:
        field_entries = converter.list_fields(args.dataset)
        field_catalog = {name: meta for name, meta in field_entries}
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if args.list_fields:
        print(f"Available annotation fields for dataset '{args.dataset}':")
        for name, meta in field_entries:
            label = meta.get("label", name)
            type_hint = "list" if meta.get("is_list") else "scalar"
            descriptor = f"{label} ({name}, {type_hint})"
            if meta.get("items"):
                descriptor += f" - allowed values: {', '.join(meta['items'])}"
            print(f"  - {descriptor}")
        return 0

    filters = parse_filters(args.filters)

    try:
        fields = ensure_fields(args.fields, field_catalog, args.include_label)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    ids = load_ids(args.ids, args.ids_file)
    if not ids:
        print("Error: No identifiers provided. Use --ids or --ids-file.", file=sys.stderr)
        return 1

    try:
        dataset = converter.get_dataset(args.dataset)
        records = converter.execute_query(args.dataset, ids, fields, filters)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if args.format == 'json':
        output_json(fields, records, output_path=args.output)
        return 0

    delimiter = args.delimiter
    if args.format == 'csv' and delimiter == '\t':
        delimiter = ','

    header, rows = converter.build_rows(
        dataset_label=dataset.get("label", args.dataset),
        fields=fields,
        field_meta=field_catalog,
        records=records,
        filters=filters,
        compact=args.compact,
    )

    if not rows:
        print("No results found for the supplied identifiers.", file=sys.stderr)
        return 0

    output_table(
        header, rows,
        delimiter=delimiter,
        output_path=args.output,
        include_header=not args.no_header,
    )

    return 0


def main():
    """Main CLI entry point"""
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    api_url = args.api_url or os.getenv('TOGOID_API_ENDPOINT') or 'https://api.togoid.dbcls.jp'

    try:
        if args.command == 'convert':
            converter = TogoIDConverter(api_base_url=api_url)
            handle_convert(args, converter)

        elif args.command == 'label2id':
            return handle_label2id(args, api_url)

        elif args.command == 'annotate':
            return handle_annotate(args, api_url)

        elif args.command == 'search':
            converter = TogoIDConverter(api_base_url=api_url)
            if args.search_type == 'databases':
                result = converter.search_databases(args.name)
            elif args.search_type == 'id':
                result = converter.search_id(args.id_string)
            else:
                parser.parse_args(['search', '--help'])
                return 1
            print(json.dumps(result, indent=2, ensure_ascii=False))

        elif args.command == 'lookup':
            converter = TogoIDConverter(api_base_url=api_url)
            if args.lookup_type == 'id':
                result = converter.lookup_id(args.id_string)
            else:
                parser.parse_args(['lookup', '--help'])
                return 1
            print(json.dumps(result, indent=2, ensure_ascii=False))

        elif args.command == 'route':
            converter = TogoIDConverter(api_base_url=api_url)
            result = converter.route(args.src, args.dst, args.max_hops)
            print(json.dumps(result, indent=2, ensure_ascii=False))

        elif args.command == 'count':
            converter = TogoIDConverter(api_base_url=api_url)
            ids = [i.strip() for i in args.ids.split(',')]
            result = converter.count(args.src, args.dst, ids, args.link)
            print(json.dumps(result, indent=2, ensure_ascii=False))

        elif args.command == 'get-ortholog':
            converter = TogoIDConverter(api_base_url=api_url)
            return handle_get_ortholog(args, converter)

        elif args.command == 'enrich':
            return handle_enrich(args, api_url)

        elif args.command == 'config':
            converter = TogoIDConverter(api_base_url=api_url)
            if args.config_type == 'dataset':
                result = converter.config_dataset(args.name)
            elif args.config_type == 'relation':
                if hasattr(args, 'src') and hasattr(args, 'dst') and args.src and args.dst:
                    result = converter.config_relation(args.src, args.dst)
                else:
                    result = converter.config_relation()
            elif args.config_type == 'descriptions':
                result = converter.config_descriptions()
            elif args.config_type == 'statistics':
                result = converter.config_statistics()
            elif args.config_type == 'taxonomy':
                result = converter.config_taxonomy()
            else:
                parser.parse_args(['config', '--help'])
                return 1
            print(json.dumps(result, indent=2, ensure_ascii=False))

    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
