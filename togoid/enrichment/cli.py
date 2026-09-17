"""
``togoid enrich`` subcommand.

Kept in its own module so the enrichment feature stays self-contained: the main
CLI only registers the parser and dispatches to :func:`handle_enrich`.
"""
import csv
import json
import sys
from typing import Dict, List, Optional, Sequence

from .analysis import enrich_clusters
from .genesets import GeneSetLibrary, build_gene_sets
from .presets import GO_ASPECTS, ROUTES

__all__ = ["add_enrich_parser", "handle_enrich"]


def add_enrich_parser(subparsers) -> None:
    """
    Register the ``enrich`` subcommand on an argparse subparser group.

    Args:
        subparsers: The result of ``parser.add_subparsers(...)``.
    """
    parser = subparsers.add_parser(
        "enrich",
        help="Enrichment analysis of gene lists via a TogoID route",
        description=(
            "Convert a gene list along a TogoID route into gene sets, then test "
            "for over-representation with a hypergeometric test and BH-FDR."
        ),
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--genes", help='Comma-separated gene symbols (e.g. "CD3D,MS4A1,LYZ")'
    )
    source.add_argument(
        "--genes-file", help="File with one gene symbol per line (a single query)"
    )
    source.add_argument(
        "--clusters-file",
        help=(
            "CSV/TSV with a cluster column and a gene column, giving one query "
            "per cluster"
        ),
    )

    parser.add_argument(
        "--route",
        help=(
            'Comma-separated TogoID route, e.g. "ncbigene,uniprot,reactome_pathway". '
            "Mutually exclusive with --preset."
        ),
    )
    parser.add_argument(
        "--preset",
        choices=sorted(ROUTES),
        help=f"Named route instead of --route. Available: {sorted(ROUTES)}",
    )
    parser.add_argument(
        "--go-aspect",
        choices=list(GO_ASPECTS),
        default="biological_process",
        help="GO aspect to keep when --preset go is used (default: biological_process)",
    )
    parser.add_argument(
        "--taxonomy",
        default="9606",
        help="Taxonomy ID used to resolve gene symbols (default: 9606)",
    )
    parser.add_argument(
        "--id-source",
        choices=["symbol", "id"],
        default="symbol",
        help=(
            "Whether the input genes are labels to resolve ('symbol', the "
            "default) or identifiers already in the route's first dataset ('id')"
        ),
    )

    parser.add_argument("--cluster-column", default="cluster", help="Column holding the cluster label")
    parser.add_argument("--gene-column", default="gene", help="Column holding the gene symbol")

    parser.add_argument("--min-set-size", type=int, default=5, help="Minimum term size (default: 5)")
    parser.add_argument("--max-set-size", type=int, default=500, help="Maximum term size (default: 500)")
    parser.add_argument("--fdr", type=float, default=None, help="Keep only terms below this FDR")
    parser.add_argument("--top", type=int, default=None, help="Keep only the top N terms per cluster")

    parser.add_argument(
        "--save-genesets", help="Write the gene-set library to this JSON file for reuse"
    )
    parser.add_argument(
        "--load-genesets",
        help="Reuse a gene-set library written by --save-genesets instead of querying the API",
    )
    parser.add_argument("--output", help="Write results to this CSV file (default: stdout)")
    parser.add_argument(
        "--format",
        choices=["csv", "json", "summary"],
        default="csv",
        help="Output format (default: csv)",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages")


def _read_gene_file(path: str) -> List[str]:
    """Read a plain-text gene list, one symbol per line."""
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def _read_clusters_file(
    path: str, cluster_column: str, gene_column: str
) -> Dict[str, List[str]]:
    """
    Read a per-cluster gene table.

    Args:
        path: CSV or TSV file. The delimiter is sniffed from the extension.
        cluster_column: Column holding the cluster label.
        gene_column: Column holding the gene symbol.

    Returns:
        Mapping of cluster label to its gene list, preserving file order.

    Raises:
        KeyError: If a requested column is missing.
    """
    delimiter = "\t" if path.lower().endswith((".tsv", ".tab", ".txt")) else ","
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))

    if not rows:
        return {}

    for column in (cluster_column, gene_column):
        if column not in rows[0]:
            raise KeyError(
                f"{path} has no column {column!r}; found {sorted(rows[0])}. "
                "Use --cluster-column / --gene-column to override."
            )

    grouped: Dict[str, List[str]] = {}
    for row in rows:
        cluster = str(row[cluster_column]).strip()
        gene = str(row[gene_column]).strip()
        if cluster and gene:
            grouped.setdefault(cluster, []).append(gene)

    # Keep order but drop duplicate symbols within a cluster.
    return {c: list(dict.fromkeys(genes)) for c, genes in grouped.items()}


def _resolve_route(args) -> Sequence[str]:
    """Turn ``--route``/``--preset`` into an explicit route."""
    if args.route and args.preset:
        raise ValueError("--route and --preset are mutually exclusive")
    if args.route:
        return [part.strip() for part in args.route.split(",") if part.strip()]
    if args.preset:
        return ROUTES[args.preset]
    raise ValueError("One of --route or --preset is required")


def handle_enrich(args, api_url: Optional[str] = None) -> int:
    """
    Execute the ``enrich`` subcommand.

    Args:
        args: Parsed arguments from :func:`add_enrich_parser`.
        api_url: TogoID API base URL (currently only used for future wiring;
            the converters read the same environment variable).

    Returns:
        Process exit code.
    """
    verbose = not args.quiet

    # Step 1: collect the query gene lists.
    if args.clusters_file:
        cluster_genes = _read_clusters_file(
            args.clusters_file, args.cluster_column, args.gene_column
        )
    elif args.genes_file:
        cluster_genes = {"query": _read_gene_file(args.genes_file)}
    else:
        cluster_genes = {
            "query": [g.strip() for g in args.genes.split(",") if g.strip()]
        }

    if not cluster_genes:
        print("No genes to analyse.", file=sys.stderr)
        return 1

    all_genes = sorted({gene for genes in cluster_genes.values() for gene in genes})

    # Step 2: build or load the gene-set library.
    if args.load_genesets:
        library = GeneSetLibrary.load_json(args.load_genesets)
        if verbose:
            print(f"Loaded gene sets from {args.load_genesets}: {library!r}")
    else:
        route = _resolve_route(args)
        term_filters = None
        if args.preset == "go" and args.go_aspect:
            term_filters = {"go_aspect": [args.go_aspect]}

        library = build_gene_sets(
            all_genes,
            route=route,
            id_source="symbol" if args.id_source == "symbol" else None,
            taxonomy=args.taxonomy,
            term_filters=term_filters,
            verbose=verbose,
        )

    if args.save_genesets:
        library.save_json(args.save_genesets)
        if verbose:
            print(f"Saved gene sets to {args.save_genesets}")

    if not library.sets:
        print("No gene sets could be built for this input.", file=sys.stderr)
        return 1

    # Step 3: enrichment.
    results = enrich_clusters(
        cluster_genes,
        library,
        min_set_size=args.min_set_size,
        max_set_size=args.max_set_size,
        verbose=verbose,
    )

    # --fdr / --top narrow the exported table. The summary keeps every tested
    # term so its "terms tested" count stays meaningful, and reports the same
    # threshold as the significance cut-off.
    filtered = results
    if args.fdr is not None:
        filtered = filtered.significant(args.fdr)
    if args.top is not None:
        filtered = filtered.top(args.top)

    # Step 4: output.
    if args.format == "summary":
        text = results.summary(alpha=args.fdr if args.fdr is not None else 0.05)
        if args.output:
            with open(args.output, "w", encoding="utf-8") as handle:
                handle.write(text)
        else:
            print(text)
    elif args.format == "json":
        text = json.dumps(filtered.to_rows(), indent=2, ensure_ascii=False)
        if args.output:
            with open(args.output, "w", encoding="utf-8") as handle:
                handle.write(text)
        else:
            print(text)
    else:  # csv
        rows = filtered.to_rows()
        if args.output:
            filtered.to_csv(args.output)
        elif rows:
            writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    if args.output and verbose:
        count = len(results) if args.format == "summary" else len(filtered)
        print(f"Wrote {count} rows to {args.output}")

    return 0
