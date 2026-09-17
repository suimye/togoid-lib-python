"""
Adapters between single-cell toolkits and the plain inputs the core expects.

``togoid.enrichment`` deliberately knows nothing about AnnData or Seurat. These
helpers do the translation, and import their heavy dependencies lazily so that
importing ``togoid.enrichment`` never pulls in scanpy.
"""
import csv
from typing import Any, Dict, List, Mapping, Optional, Sequence

__all__ = [
    "umap_dataframe_from_anndata",
    "marker_genes_from_anndata",
    "umap_dataframe_from_csv",
    "marker_genes_from_csv",
    "write_marker_gene_lists",
]


def _require_pandas():
    """Import pandas, with an actionable error when it is missing."""
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "pandas is required for this adapter. "
            'Install with: pip install "togoid[enrichment]"'
        ) from exc
    return pd


def umap_dataframe_from_anndata(
    adata: Any,
    basis: str = "X_umap",
    cluster_key: str = "leiden",
    x_key: str = "umap_1",
    y_key: str = "umap_2",
):
    """
    Extract a plotting-ready embedding table from an AnnData object.

    Args:
        adata: AnnData object carrying the embedding and cluster assignments.
        basis: Key in ``adata.obsm`` holding the 2-D embedding.
        cluster_key: Column in ``adata.obs`` holding cluster labels.
        x_key: Name given to the first embedding column in the output.
        y_key: Name given to the second embedding column in the output.

    Returns:
        A pandas DataFrame indexed by cell barcode with ``x_key``, ``y_key`` and
        ``cluster`` columns — exactly what :func:`~togoid.enrichment.plot.plot_umap_enrichment`
        expects.

    Raises:
        KeyError: If ``basis`` or ``cluster_key`` is absent.
        ValueError: If the embedding has fewer than two dimensions.
    """
    pd = _require_pandas()

    if basis not in adata.obsm:
        raise KeyError(
            f"adata.obsm has no {basis!r}; available: {list(adata.obsm.keys())}"
        )
    if cluster_key not in adata.obs:
        raise KeyError(
            f"adata.obs has no {cluster_key!r}; available: {list(adata.obs.columns)}"
        )

    coords = adata.obsm[basis]
    if getattr(coords, "shape", (0, 0))[1] < 2:
        raise ValueError(f"{basis!r} must have at least 2 dimensions")

    return pd.DataFrame(
        {
            x_key: [float(v) for v in coords[:, 0]],
            y_key: [float(v) for v in coords[:, 1]],
            "cluster": [str(v) for v in adata.obs[cluster_key]],
        },
        index=adata.obs_names,
    )


def marker_genes_from_anndata(
    adata: Any,
    key: str = "rank_genes_groups",
    top_n: Optional[int] = 100,
    pval_cutoff: Optional[float] = 0.05,
    logfc_min: Optional[float] = 0.25,
    use_adjusted: bool = True,
) -> Dict[str, List[str]]:
    """
    Read per-cluster marker genes out of a scanpy differential expression run.

    Call ``scanpy.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon")``
    first; this helper only formats what that stored in ``adata.uns``.

    Args:
        adata: AnnData object with differential expression results.
        key: Key in ``adata.uns`` holding the results.
        top_n: Keep at most this many genes per cluster; ``None`` keeps all.
        pval_cutoff: Drop genes at or above this (adjusted) p-value.
        logfc_min: Drop genes with a smaller log fold change; ``None`` disables.
        use_adjusted: Filter on ``pvals_adj`` rather than the raw p-value.

    Returns:
        Mapping of cluster label to its marker gene symbols, most significant
        first.

    Raises:
        KeyError: If ``key`` is not present in ``adata.uns``.
    """
    if key not in adata.uns:
        raise KeyError(
            f"adata.uns has no {key!r}. Run scanpy.tl.rank_genes_groups() first."
        )

    result = adata.uns[key]
    groups = list(result["names"].dtype.names)
    pval_field = "pvals_adj" if use_adjusted else "pvals"

    markers: Dict[str, List[str]] = {}
    for group in groups:
        names = list(result["names"][group])
        pvals = list(result[pval_field][group]) if pval_field in result else None
        logfcs = list(result["logfoldchanges"][group]) if "logfoldchanges" in result else None

        genes: List[str] = []
        for index, name in enumerate(names):
            if pval_cutoff is not None and pvals is not None:
                value = pvals[index]
                # scanpy writes NaN for genes it could not test.
                if value != value or value >= pval_cutoff:
                    continue
            if logfc_min is not None and logfcs is not None:
                value = logfcs[index]
                if value != value or value < logfc_min:
                    continue
            genes.append(str(name))
            if top_n is not None and len(genes) >= top_n:
                break

        markers[str(group)] = genes

    return markers


def umap_dataframe_from_csv(
    path: str,
    x_column: str = "UMAP_1",
    y_column: str = "UMAP_2",
    cluster_column: str = "cluster",
    x_key: str = "umap_1",
    y_key: str = "umap_2",
):
    """
    Load an embedding exported from another toolkit, such as Seurat.

    In Seurat, write the table with::

        df <- cbind(Embeddings(obj, "umap"), cluster = as.character(Idents(obj)))
        write.csv(df, "umap.csv")

    Args:
        path: Path to the CSV file.
        x_column: Column in the file holding the first dimension.
        y_column: Column in the file holding the second dimension.
        cluster_column: Column in the file holding cluster labels.
        x_key: Name given to the first embedding column in the output.
        y_key: Name given to the second embedding column in the output.

    Returns:
        A pandas DataFrame with ``x_key``, ``y_key`` and ``cluster`` columns.

    Raises:
        KeyError: If a requested column is missing from the file.
    """
    pd = _require_pandas()
    frame = pd.read_csv(path)

    missing = [
        column
        for column in (x_column, y_column, cluster_column)
        if column not in frame.columns
    ]
    if missing:
        raise KeyError(
            f"{path} is missing column(s) {missing}; found {list(frame.columns)}"
        )

    return pd.DataFrame(
        {
            x_key: frame[x_column].astype(float),
            y_key: frame[y_column].astype(float),
            "cluster": frame[cluster_column].astype(str),
        }
    )


def marker_genes_from_csv(
    path: str,
    cluster_column: str = "cluster",
    gene_column: str = "gene",
    pval_column: Optional[str] = "p_val_adj",
    logfc_column: Optional[str] = "avg_log2FC",
    pval_cutoff: Optional[float] = 0.05,
    logfc_min: Optional[float] = 0.25,
    top_n: Optional[int] = 100,
) -> Dict[str, List[str]]:
    """
    Load per-cluster marker genes from a long-format CSV.

    The default column names match ``Seurat::FindAllMarkers()`` output.

    Args:
        path: Path to the CSV file.
        cluster_column: Column holding the cluster label.
        gene_column: Column holding the gene symbol.
        pval_column: Column holding the adjusted p-value, or ``None`` to skip.
        logfc_column: Column holding the log fold change, or ``None`` to skip.
        pval_cutoff: Drop rows at or above this p-value.
        logfc_min: Drop rows below this log fold change.
        top_n: Keep at most this many genes per cluster.

    Returns:
        Mapping of cluster label to marker gene symbols, most significant first.

    Raises:
        KeyError: If a required column is missing.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        return {}

    available = set(rows[0])
    for column in (cluster_column, gene_column):
        if column not in available:
            raise KeyError(
                f"{path} is missing column {column!r}; found {sorted(available)}"
            )

    def _sort_key(row: Mapping[str, str]) -> float:
        if pval_column and pval_column in row:
            try:
                return float(row[pval_column])
            except (TypeError, ValueError):
                return 1.0
        return 1.0

    grouped: Dict[str, List[Mapping[str, str]]] = {}
    for row in rows:
        if pval_column and pval_cutoff is not None and pval_column in available:
            try:
                if float(row[pval_column]) >= pval_cutoff:
                    continue
            except (TypeError, ValueError):
                continue
        if logfc_column and logfc_min is not None and logfc_column in available:
            try:
                if float(row[logfc_column]) < logfc_min:
                    continue
            except (TypeError, ValueError):
                continue
        grouped.setdefault(str(row[cluster_column]), []).append(row)

    markers: Dict[str, List[str]] = {}
    for cluster, cluster_rows in grouped.items():
        cluster_rows.sort(key=_sort_key)
        genes = [str(row[gene_column]) for row in cluster_rows]
        # Preserve order while removing duplicate symbols.
        genes = list(dict.fromkeys(genes))
        markers[cluster] = genes[:top_n] if top_n is not None else genes

    return markers


def write_marker_gene_lists(
    markers: Mapping[str, Sequence[str]], directory: str, prefix: str = "cluster_"
) -> Dict[str, str]:
    """
    Write one plain-text gene list per cluster.

    Args:
        markers: Mapping of cluster label to gene symbols.
        directory: Destination directory, created if needed.
        prefix: File name prefix; files are ``<prefix><cluster>_markers.txt``.

    Returns:
        Mapping of cluster label to the path written.
    """
    import os

    os.makedirs(directory, exist_ok=True)
    written: Dict[str, str] = {}
    for cluster, genes in markers.items():
        path = os.path.join(directory, f"{prefix}{cluster}_markers.txt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(genes))
            handle.write("\n")
        written[str(cluster)] = path
    return written
