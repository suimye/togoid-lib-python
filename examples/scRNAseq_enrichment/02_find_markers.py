#!/usr/bin/env python3
"""
Step 2: find marker genes for every cluster.

A Wilcoxon rank-sum test per cluster gives the gene symbols that feed TogoID in
step 3. The gene lists are also written as plain text, one file per cluster, so
they can be inspected or reused outside this pipeline.
"""
import argparse
import os

import scanpy as sc

from togoid.enrichment.adapters import marker_genes_from_anndata, write_marker_gene_lists


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir", default="results", help="Directory holding step 1 outputs"
    )
    parser.add_argument("--cluster-key", default="leiden", help="Cluster column in adata.obs")
    parser.add_argument("--top-n", type=int, default=100, help="Marker genes kept per cluster")
    parser.add_argument("--pval-cutoff", type=float, default=0.05, help="Adjusted p-value cut-off")
    parser.add_argument("--logfc-min", type=float, default=0.25, help="Minimum log2 fold change")
    return parser.parse_args()


def main() -> int:
    """Run step 2."""
    args = parse_arguments()

    print("=" * 60)
    print("Step 2: marker genes per cluster")
    print("=" * 60)

    h5ad_path = os.path.join(args.results_dir, "01_clustered.h5ad")
    print(f"Loading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    print(f"Testing differential expression ({args.cluster_key}, Wilcoxon)...")
    sc.tl.rank_genes_groups(
        adata, args.cluster_key, method="wilcoxon", use_raw=True, pts=True
    )

    markers = marker_genes_from_anndata(
        adata,
        top_n=args.top_n,
        pval_cutoff=args.pval_cutoff,
        logfc_min=args.logfc_min,
    )

    for cluster in sorted(markers, key=lambda c: (len(c), c)):
        print(f"  cluster {cluster}: {len(markers[cluster])} marker genes")

    # One plain-text list per cluster, for inspection and reuse.
    lists_dir = os.path.join(args.results_dir, "02_marker_gene_lists")
    write_marker_gene_lists(markers, lists_dir)
    print(f"  wrote gene lists to {lists_dir}")

    # A single long-format table, which is what `togoid enrich --clusters-file`
    # and step 3 both read.
    table_path = os.path.join(args.results_dir, "02_markers.csv")
    with open(table_path, "w", encoding="utf-8") as handle:
        handle.write("cluster,gene\n")
        for cluster in sorted(markers, key=lambda c: (len(c), c)):
            for gene in markers[cluster]:
                handle.write(f"{cluster},{gene}\n")
    print(f"  wrote {table_path}")

    adata.write_h5ad(os.path.join(args.results_dir, "02_with_markers.h5ad"))

    print("\nNext: 03_enrichment.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
