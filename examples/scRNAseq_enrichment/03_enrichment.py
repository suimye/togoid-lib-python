#!/usr/bin/env python3
"""
Step 3: convert marker genes with TogoID and test them for enrichment.

This is the step the TogoID library exists for. One call to ``build_gene_sets``
takes gene symbols all the way to an annotation database — resolving symbols to
NCBI Gene IDs, walking the route, and fetching term labels — and the same code
serves Reactome, GO and MONDO by changing nothing but the route.
"""
import argparse
import os
from typing import Dict, List

from togoid.enrichment import GeneSetLibrary, build_gene_sets, enrich_clusters
from togoid.enrichment.presets import GO_ASPECTS, ROUTES

#: Databases analysed by default, with the route and term-size bounds each one
#: wants. GO and Reactome have large sets; MONDO's disease sets are much smaller,
#: so a lower minimum keeps them testable.
TARGETS = {
    "reactome": {
        "route": ROUTES["reactome"],
        "term_filters": None,
        "min_set_size": 5,
        "max_set_size": 500,
    },
    "go": {
        "route": ROUTES["go"],
        "term_filters": {"go_aspect": ["biological_process"]},
        "min_set_size": 5,
        "max_set_size": 500,
    },
    "mondo": {
        "route": ROUTES["mondo"],
        "term_filters": None,
        "min_set_size": 3,
        "max_set_size": 500,
    },
}


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="results", help="Directory holding step 2 outputs")
    parser.add_argument(
        "--targets",
        default=",".join(TARGETS),
        help=f"Comma-separated databases to analyse. Available: {sorted(TARGETS)}",
    )
    parser.add_argument(
        "--go-aspect",
        default="biological_process",
        choices=list(GO_ASPECTS),
        help="GO aspect to keep (default: biological_process)",
    )
    parser.add_argument("--taxonomy", default="9606", help="Taxonomy ID (default: 9606, human)")
    parser.add_argument(
        "--reuse-genesets",
        action="store_true",
        help="Reuse cached gene-set JSON files instead of querying the API again",
    )
    return parser.parse_args()


def load_markers(results_dir: str) -> Dict[str, List[str]]:
    """
    Read the per-cluster marker table written by step 2.

    Args:
        results_dir: Directory holding ``02_markers.csv``.

    Returns:
        Mapping of cluster label to its marker gene symbols.
    """
    import csv

    path = os.path.join(results_dir, "02_markers.csv")
    markers: Dict[str, List[str]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            markers.setdefault(row["cluster"], []).append(row["gene"])

    total = sum(len(v) for v in markers.values())
    print(f"Loaded {total} marker genes across {len(markers)} clusters")
    return markers


def run_target(
    name: str,
    markers: Dict[str, List[str]],
    all_genes: List[str],
    args: argparse.Namespace,
) -> None:
    """
    Build gene sets and run enrichment for one annotation database.

    Args:
        name: Target key in :data:`TARGETS`.
        markers: Per-cluster marker genes.
        all_genes: Every marker gene, deduplicated; this is the background.
        args: Parsed command-line arguments.
    """
    config = TARGETS[name]
    term_filters = config["term_filters"]
    if name == "go":
        term_filters = {"go_aspect": [args.go_aspect]}

    print("\n" + "=" * 60)
    print(f"Target: {name}  ({' -> '.join(config['route'])})")
    print("=" * 60)

    cache_path = os.path.join(args.results_dir, f"03_genesets_{name}.json")

    if args.reuse_genesets and os.path.exists(cache_path):
        library = GeneSetLibrary.load_json(cache_path)
        print(f"Reusing cached gene sets from {cache_path}: {library!r}")
    else:
        library = build_gene_sets(
            all_genes,
            route=config["route"],
            taxonomy=args.taxonomy,
            term_filters=term_filters,
            verbose=True,
        )
        library.save_json(cache_path)
        print(f"Saved gene sets to {cache_path}")

    if not library.sets:
        print(f"No gene sets built for {name}; skipping.")
        return

    results = enrich_clusters(
        markers,
        library,
        min_set_size=config["min_set_size"],
        max_set_size=config["max_set_size"],
        verbose=True,
    )

    all_path = os.path.join(args.results_dir, f"03_enrichment_{name}_all.csv")
    results.to_csv(all_path)
    print(f"Wrote {len(results)} rows to {all_path}")

    significant = results.significant(0.05)
    sig_path = os.path.join(args.results_dir, f"03_enrichment_{name}_significant.csv")
    significant.to_csv(sig_path)
    print(f"Wrote {len(significant)} significant rows to {sig_path}")

    summary_path = os.path.join(args.results_dir, f"03_enrichment_{name}_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(results.summary())
    print(f"Wrote {summary_path}")


def main() -> int:
    """Run step 3."""
    args = parse_arguments()

    print("=" * 60)
    print("Step 3: TogoID conversion and enrichment analysis")
    print("=" * 60)

    markers = load_markers(args.results_dir)
    all_genes = sorted({gene for genes in markers.values() for gene in genes})
    print(f"Unique marker genes (background): {len(all_genes)}")

    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    unknown = [t for t in targets if t not in TARGETS]
    if unknown:
        raise SystemExit(f"Unknown target(s) {unknown}; available: {sorted(TARGETS)}")

    for name in targets:
        run_target(name, markers, all_genes, args)

    print("\nNext: 04_visualize_umap.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
