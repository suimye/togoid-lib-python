#!/usr/bin/env python3
"""
Step 4: draw the enrichment results on the UMAP embedding.

Each figure has two panels: the usual cluster UMAP on the left, and the same
embedding on the right with each cluster's enriched terms written around its
centroid, sized by significance.

This step needs neither scanpy nor the AnnData object — only the plain UMAP
table from step 1 and the enrichment CSV from step 3.
"""
import argparse
import csv
import os
from typing import Any, Dict, List

import matplotlib

# A non-interactive backend, so the script also works over SSH and in CI.
matplotlib.use("Agg")

from togoid.enrichment import plot_umap_centroids, plot_umap_enrichment


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="results", help="Directory holding earlier outputs")
    parser.add_argument(
        "--targets",
        default="reactome,go,mondo",
        help="Comma-separated databases to plot (must have been run in step 3)",
    )
    parser.add_argument("--top-n", type=int, default=3, help="Terms drawn per cluster (default: 3)")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05, help="FDR cut-off (default: 0.05)")
    parser.add_argument("--pval-cutoff", type=float, default=None, help="Optional p-value cut-off")
    parser.add_argument(
        "--clusters", default=None, help="Comma-separated cluster IDs to label (default: all)"
    )
    parser.add_argument("--max-label-chars", type=int, default=40, help="Truncate labels beyond this")
    parser.add_argument(
        "--no-centroids",
        action="store_true",
        help=(
            "Hide the cluster centroid markers. They are still drawn "
            "transparently and still reserve their space, so the labels do not move"
        ),
    )
    parser.add_argument(
        "--centroid-marker",
        default="o",
        help="Matplotlib marker for the centroids (default: o, a filled circle)",
    )
    parser.add_argument(
        "--centroid-size", type=float, default=26.0, help="Centroid marker size (default: 26)"
    )
    parser.add_argument("--dpi", type=int, default=200, help="Raster resolution for the PNG output")
    parser.add_argument(
        "--formats", default="pdf,png", help="Comma-separated output formats (default: pdf,png)"
    )
    return parser.parse_args()


def read_umap(results_dir: str) -> Dict[str, List[Any]]:
    """
    Read the UMAP table written by step 1.

    Returns plain lists rather than a DataFrame, to show that the plotting API
    does not require pandas.

    Args:
        results_dir: Directory holding ``01_umap.csv``.

    Returns:
        Mapping with ``umap_1``, ``umap_2`` and ``cluster`` keys.
    """
    path = os.path.join(results_dir, "01_umap.csv")
    embedding: Dict[str, List[Any]] = {"umap_1": [], "umap_2": [], "cluster": []}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            embedding["umap_1"].append(float(row["umap_1"]))
            embedding["umap_2"].append(float(row["umap_2"]))
            embedding["cluster"].append(str(row["cluster"]))
    print(f"Loaded {len(embedding['cluster'])} cells from {path}")
    return embedding


def read_enrichment(results_dir: str, target: str) -> List[Dict[str, Any]]:
    """
    Read one enrichment CSV written by step 3.

    Args:
        results_dir: Directory holding the CSV.
        target: Database name used in the file name.

    Returns:
        List of row dictionaries with numeric p-value and FDR.
    """
    path = os.path.join(results_dir, f"03_enrichment_{target}_all.csv")
    rows: List[Dict[str, Any]] = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row["pvalue"] = float(row["pvalue"])
            row["fdr"] = float(row["fdr"])
            rows.append(row)
    print(f"  loaded {len(rows)} enrichment rows from {path}")
    return rows


def save_figure(fig, results_dir: str, stem: str, formats: List[str], dpi: int) -> None:
    """
    Save a figure in every requested format.

    Args:
        fig: Matplotlib figure.
        results_dir: Output directory.
        stem: File name without extension.
        formats: Extensions to write.
        dpi: Raster resolution.
    """
    for extension in formats:
        path = os.path.join(results_dir, f"{stem}.{extension}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"  wrote {path}")


def main() -> int:
    """Run step 4."""
    args = parse_arguments()

    print("=" * 60)
    print("Step 4: UMAP visualisation of enrichment results")
    print("=" * 60)

    embedding = read_umap(args.results_dir)
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]
    clusters = (
        [c.strip() for c in args.clusters.split(",") if c.strip()]
        if args.clusters
        else None
    )

    # A reference figure showing where the labels will be anchored.
    fig = plot_umap_centroids(
        embedding,
        title="PBMC clusters and centroids",
        centroid_marker=args.centroid_marker,
    )
    save_figure(fig, args.results_dir, "04_umap_centroids", formats, args.dpi)

    titles = {
        "reactome": "Enriched Reactome pathways",
        "go": "Enriched GO biological processes",
        "mondo": "Enriched MONDO diseases",
    }

    for target in [t.strip() for t in args.targets.split(",") if t.strip()]:
        print(f"\nTarget: {target}")
        rows = read_enrichment(args.results_dir, target)

        fig = plot_umap_enrichment(
            embedding,
            rows,
            top_n=args.top_n,
            fdr_cutoff=args.fdr_cutoff,
            pval_cutoff=args.pval_cutoff,
            clusters=clusters,
            max_label_chars=args.max_label_chars,
            title_left="PBMC clusters (UMAP)",
            title_right=f"{titles.get(target, target)} (top {args.top_n} per cluster)",
            show_centroids=not args.no_centroids,
            centroid_marker=args.centroid_marker,
            centroid_size=args.centroid_size,
            verbose=True,
        )
        save_figure(
            fig,
            args.results_dir,
            f"04_umap_enrichment_{target}_top{args.top_n}",
            formats,
            args.dpi,
        )

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
