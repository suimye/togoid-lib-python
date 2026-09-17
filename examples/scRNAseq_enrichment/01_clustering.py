#!/usr/bin/env python3
"""
Step 1: cluster a 10x Genomics scRNA-seq dataset and compute a UMAP embedding.

This step uses scanpy only; TogoID enters the pipeline in step 3. The clustering
produced here is reused unchanged by every later step, so that the UMAP shown
next to the enrichment results is exactly the one the gene lists came from.
"""
import argparse
import os

import scanpy as sc

from togoid.enrichment.adapters import umap_dataframe_from_anndata


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory holding the 10x matrix (barcodes/features/matrix files)",
    )
    parser.add_argument(
        "--results-dir", default="results", help="Directory for outputs (default: results)"
    )
    parser.add_argument("--min-genes", type=int, default=200, help="Minimum genes per cell")
    parser.add_argument("--max-genes", type=int, default=6000, help="Maximum genes per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum cells per gene")
    parser.add_argument("--max-mt-pct", type=float, default=15.0, help="Maximum mitochondrial percentage")
    parser.add_argument("--n-hvg", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--n-pcs", type=int, default=30, help="Number of principal components")
    parser.add_argument("--n-neighbors", type=int, default=15, help="Neighbours for the kNN graph")
    parser.add_argument("--resolution", type=float, default=0.5, help="Leiden resolution")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    return parser.parse_args()


def load_and_filter(args: argparse.Namespace):
    """
    Load the 10x matrix and apply standard quality-control filters.

    Args:
        args: Parsed command-line arguments.

    Returns:
        The filtered AnnData object.
    """
    print(f"Loading 10x data from {args.data_dir}")
    adata = sc.read_10x_mtx(args.data_dir, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    print(f"  loaded {adata.n_obs} cells x {adata.n_vars} genes")

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs.n_genes_by_counts < args.max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :].copy()
    print(f"  after QC: {adata.n_obs} cells x {adata.n_vars} genes")

    return adata


def normalize_and_cluster(adata, args: argparse.Namespace):
    """
    Normalise, embed and cluster the cells.

    Args:
        adata: Filtered AnnData object.
        args: Parsed command-line arguments.

    Returns:
        The AnnData object with ``X_umap`` and a ``leiden`` column.
    """
    # Keep the raw counts around; differential expression in step 2 needs the
    # full gene set, not just the highly variable genes.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
    adata = adata[:, adata.var.highly_variable].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs, random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=args.seed,
    )

    counts = adata.obs["leiden"].value_counts().sort_index()
    print(f"  found {len(counts)} clusters")
    for cluster, n_cells in counts.items():
        print(f"    cluster {cluster}: {n_cells} cells")

    return adata


def save_outputs(adata, results_dir: str) -> None:
    """
    Write the clustered object and a plain UMAP table.

    Args:
        adata: Clustered AnnData object.
        results_dir: Output directory.
    """
    os.makedirs(results_dir, exist_ok=True)

    h5ad_path = os.path.join(results_dir, "01_clustered.h5ad")
    adata.write_h5ad(h5ad_path)
    print(f"  wrote {h5ad_path}")

    # The plain table is what the plotting API consumes, and it keeps the later
    # steps usable without scanpy.
    umap_df = umap_dataframe_from_anndata(adata, cluster_key="leiden")
    umap_path = os.path.join(results_dir, "01_umap.csv")
    umap_df.to_csv(umap_path, index_label="cell")
    print(f"  wrote {umap_path}")


def main() -> int:
    """Run step 1."""
    args = parse_arguments()

    print("=" * 60)
    print("Step 1: clustering and UMAP")
    print("=" * 60)

    adata = load_and_filter(args)
    adata = normalize_and_cluster(adata, args)
    save_outputs(adata, args.results_dir)

    print("\nNext: 02_find_markers.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
