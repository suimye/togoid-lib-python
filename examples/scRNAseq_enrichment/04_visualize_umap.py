#!/usr/bin/env python3
"""
Step 4: draw the enrichment results on the UMAP embedding.

Each figure has two panels: the usual cluster UMAP on the left, and the same
embedding on the right with each cluster's enriched terms written around its
centroid, sized by significance.

Beside every figure this writes the same terms as a TSV, so the table and the
picture always agree: a long table with one row per term, and a wide one with one
row per cluster.

This step needs neither scanpy nor the AnnData object — only the plain UMAP
table from step 1 and the enrichment CSV from step 3.
"""
import argparse
import csv
import glob
import os
from datetime import date
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

# A non-interactive backend, so the script also works over SSH and in CI.
matplotlib.use("Agg")

from togoid.enrichment import (
    plot_umap_centroids,
    plot_umap_enrichment,
    read_metadata,
    select_terms,
    selected_terms_table,
    timestamp,
)


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
    parser.add_argument(
        "--no-tables",
        action="store_true",
        help="Skip the TSV tables written beside each figure",
    )
    parser.add_argument(
        "--date-suffix",
        default=date.today().strftime("%Y%m%d"),
        help=(
            "Date stamp added to output file names, so a re-run against an "
            "updated API does not overwrite an earlier result "
            "(default: today, YYYYMMDD). Pass an empty string to omit it."
        ),
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


def find_enrichment_file(results_dir: str, target: str) -> Optional[str]:
    """
    Locate step 3's output for one target.

    Step 3 date-stamps its files, so pick the most recent one rather than
    guessing today's date; a plain unstamped name still works.

    Args:
        results_dir: Directory holding the file.
        target: Database name used in the file name.

    Returns:
        The path, or ``None`` when nothing matches.
    """
    patterns = [
        f"03_enrichment_{target}_all*.tsv",
        f"03_enrichment_{target}_all*.csv",   # files from an older run
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.join(results_dir, pattern)))
        if matches:
            return matches[-1]
    return None


def read_enrichment(path: str) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    """
    Read one enrichment table written by step 3, header included.

    Args:
        path: Path to the table.

    Returns:
        Tuple of ``(rows, metadata)``. The metadata carries the date the TogoID
        API was queried, which this step passes on to its own outputs.
    """
    metadata = read_metadata(path)
    delimiter = "\t" if path.endswith(".tsv") else ","

    rows: List[Dict[str, Any]] = []
    with open(path, newline="", encoding="utf-8") as handle:
        # Skip the "#" provenance header; csv has no notion of comment lines.
        data = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(data, delimiter=delimiter):
            row["pvalue"] = float(row["pvalue"])
            row["fdr"] = float(row["fdr"])
            rows.append(row)

    print(f"  loaded {len(rows)} enrichment rows from {path}")
    if metadata.get("api_retrieved_at"):
        print(f"  TogoID API retrieved on {metadata['api_retrieved_at']}")
    return rows, metadata


def write_tsv(
    rows: List[Dict[str, Any]],
    path: str,
    columns: List[str],
    header: List[str],
) -> None:
    """
    Write rows to a tab-separated file, behind a ``#`` provenance header.

    Tabs rather than commas: term labels routinely contain commas, which a CSV
    has to quote and some spreadsheet imports then mis-parse.

    Args:
        rows: Rows to write.
        path: Destination file path.
        columns: Column order.
        header: Comment lines written above the table.
    """
    with open(path, "w", newline="", encoding="utf-8") as handle:
        for line in header:
            handle.write(f"# {line}\n")
        writer = csv.DictWriter(
            handle, fieldnames=columns, extrasaction="ignore", delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"  wrote {path}")


def build_header(metadata: Dict[str, str], args: argparse.Namespace, top_n: int) -> List[str]:
    """
    Build the provenance header for the tables this step writes.

    It carries step 3's header forward — above all the date the TogoID API was
    queried — and adds the options used here, so the file explains itself.

    Args:
        metadata: Header read back from step 3's output.
        args: Parsed command-line arguments.
        top_n: Terms per cluster drawn on the figure.

    Returns:
        Comment lines, without the leading "#".
    """
    lines = ["togoid enrichment: terms drawn on the UMAP figure"]
    lines.append(f"figure_generated_at: {timestamp()}")
    for key in (
        "api_retrieved_at", "api_base_url", "route", "target_dataset",
        "taxonomy", "term_filters", "min_set_size", "max_set_size",
        "background_size", "background", "togoid_version",
    ):
        if metadata.get(key):
            lines.append(f"{key}: {metadata[key]}")

    lines.append(f"selection: top {top_n} terms per cluster")
    lines.append(f"fdr_cutoff: {args.fdr_cutoff}")
    if args.pval_cutoff is not None:
        lines.append(f"pval_cutoff: {args.pval_cutoff}")
    if args.clusters:
        lines.append(f"clusters: {args.clusters}")
    lines.append(f"show_centroids: {not args.no_centroids}")
    return lines


def write_tables(
    rows: List[Dict[str, Any]],
    results_dir: str,
    stem: str,
    top_n: int,
    header: List[str],
    **filters,
) -> None:
    """
    Write the terms drawn on a figure as long- and wide-format TSV tables.

    The filters are the same ones the figure used, so the tables and the figure
    cannot disagree.

    Args:
        rows: The full enrichment table for this target.
        results_dir: Output directory.
        stem: File name stem shared with the figure.
        top_n: Terms per cluster, used for the wide table's columns.
        header: Provenance comment lines written above each table.
        **filters: Passed through to ``select_terms``.
    """
    selected = select_terms(rows, top_n=top_n, **filters)
    long_rows = selected_terms_table(selected)

    if not long_rows:
        print("  no terms passed the filters; skipping the tables")
        return

    long_columns = [
        "cluster", "term_id", "term_label", "overlap_count", "term_size",
        "query_size", "background_size", "pvalue", "fdr", "fold_enrichment", "genes",
    ]
    long_columns = [c for c in long_columns if c in long_rows[0]]
    write_tsv(long_rows, os.path.join(results_dir, f"{stem}.tsv"), long_columns, header)

    # Wide: one row per cluster, its best terms side by side. This is the shape
    # you want when labelling clusters or reading the figure as a table.
    by_cluster: Dict[str, List[Dict[str, Any]]] = {}
    for row in long_rows:
        by_cluster.setdefault(str(row["cluster"]), []).append(row)

    wide_columns = ["cluster", "n_terms"]
    for index in range(1, top_n + 1):
        wide_columns += [f"top{index}_term_id", f"top{index}_term_label", f"top{index}_fdr"]

    wide_rows = []
    for cluster, terms in by_cluster.items():
        wide: Dict[str, Any] = {"cluster": cluster, "n_terms": len(terms)}
        for index in range(1, top_n + 1):
            term = terms[index - 1] if index <= len(terms) else None
            wide[f"top{index}_term_id"] = term["term_id"] if term else ""
            wide[f"top{index}_term_label"] = term["term_label"] if term else ""
            wide[f"top{index}_fdr"] = f"{float(term['fdr']):.3e}" if term else ""
        wide_rows.append(wide)

    write_tsv(
        wide_rows,
        os.path.join(results_dir, f"{stem}_by_cluster.tsv"),
        wide_columns,
        header + [f"table: one row per cluster, top {top_n} terms"],
    )


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

    suffix = f"_{args.date_suffix}" if args.date_suffix else ""

    # A reference figure showing where the labels will be anchored.
    fig = plot_umap_centroids(
        embedding,
        title="PBMC clusters and centroids",
        centroid_marker=args.centroid_marker,
    )
    save_figure(fig, args.results_dir, f"04_umap_centroids{suffix}", formats, args.dpi)

    titles = {
        "reactome": "Enriched Reactome pathways",
        "go": "Enriched GO biological processes",
        "mondo": "Enriched MONDO diseases",
    }

    for target in [t.strip() for t in args.targets.split(",") if t.strip()]:
        print(f"\nTarget: {target}")

        path = find_enrichment_file(args.results_dir, target)
        if path is None:
            print(f"  no step 3 output found for {target}; skipping")
            continue
        rows, metadata = read_enrichment(path)

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
        stem = f"04_umap_enrichment_{target}_top{args.top_n}{suffix}"
        save_figure(fig, args.results_dir, stem, formats, args.dpi)

        if not args.no_tables:
            write_tables(
                rows,
                args.results_dir,
                stem,
                top_n=args.top_n,
                header=build_header(metadata, args, args.top_n),
                fdr_cutoff=args.fdr_cutoff,
                pval_cutoff=args.pval_cutoff,
                clusters=clusters,
                max_label_chars=None,  # keep full labels in the table
            )

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
