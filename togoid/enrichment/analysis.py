"""
Over-representation analysis over a :class:`~togoid.enrichment.genesets.GeneSetLibrary`.

Results are plain dataclasses so the analysis runs without pandas; a DataFrame
view is available for anyone who has it installed.
"""
import csv
import os
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set

from .genesets import GeneSetLibrary
from .stats import benjamini_hochberg, fold_enrichment, hypergeometric_sf

__all__ = [
    "EnrichmentRow",
    "EnrichmentResult",
    "ClusterEnrichmentResult",
    "enrich",
    "enrich_clusters",
    "RESULT_COLUMNS",
]

#: Column order used for every table this module produces.
RESULT_COLUMNS = [
    "cluster",
    "term_id",
    "term_label",
    "overlap_count",
    "term_size",
    "query_size",
    "background_size",
    "pvalue",
    "fdr",
    "fold_enrichment",
    "genes",
]


@dataclass
class EnrichmentRow:
    """A single enriched term."""

    term_id: str
    term_label: str
    overlap_count: int
    term_size: int
    query_size: int
    background_size: int
    pvalue: float
    fdr: float
    fold_enrichment: float
    genes: List[str] = field(default_factory=list)
    cluster: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Return the row as a dict with genes joined into a single string."""
        data = asdict(self)
        data["genes"] = ",".join(self.genes)
        return {column: data.get(column) for column in RESULT_COLUMNS}


class _BaseResult:
    """Shared output behaviour for single-query and per-cluster results."""

    rows: List[EnrichmentRow]

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def to_rows(self) -> List[Dict[str, Any]]:
        """Return the result as a list of plain dictionaries."""
        return [row.to_dict() for row in self.rows]

    def to_dataframe(self):
        """
        Return the result as a pandas DataFrame.

        Raises:
            ImportError: If pandas is not installed.
        """
        try:
            import pandas as pd
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise ImportError(
                "pandas is required for to_dataframe(). "
                'Install with: pip install "togoid[enrichment]"'
            ) from exc
        frame = pd.DataFrame(self.to_rows(), columns=RESULT_COLUMNS)
        # A single-query result has no cluster; drop the column rather than
        # filling the table with None.
        if frame.empty or frame["cluster"].isna().all():
            frame = frame.drop(columns=["cluster"])
        return frame

    def to_csv(self, path: str) -> str:
        """
        Write the result to a CSV file.

        Args:
            path: Destination file path.

        Returns:
            The path that was written.
        """
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)

        rows = self.to_rows()
        columns = list(RESULT_COLUMNS)
        if not any(row.get("cluster") is not None for row in rows):
            columns.remove("cluster")

        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return path


class EnrichmentResult(_BaseResult):
    """Enrichment results for a single query gene list."""

    def __init__(self, rows: Sequence[EnrichmentRow], query_genes: Sequence[str]):
        self.rows = list(rows)
        self.query_genes = list(query_genes)

    def __repr__(self) -> str:
        return f"<EnrichmentResult terms={len(self.rows)} query={len(self.query_genes)}>"

    def significant(self, alpha: float = 0.05, use: str = "fdr") -> "EnrichmentResult":
        """
        Keep only terms passing a significance threshold.

        Args:
            alpha: Cut-off value.
            use: ``"fdr"`` (default) or ``"pvalue"``.

        Returns:
            A new ``EnrichmentResult``.
        """
        kept = [row for row in self.rows if getattr(row, use) < alpha]
        return EnrichmentResult(kept, self.query_genes)

    def top(self, n: int) -> "EnrichmentResult":
        """Return the ``n`` most significant terms."""
        return EnrichmentResult(self.rows[:n], self.query_genes)


class ClusterEnrichmentResult(_BaseResult):
    """Enrichment results for several clusters, concatenated."""

    def __init__(self, per_cluster: Mapping[str, EnrichmentResult]):
        self.per_cluster: Dict[str, EnrichmentResult] = dict(per_cluster)
        self.rows = [row for result in self.per_cluster.values() for row in result.rows]

    def __repr__(self) -> str:
        return (
            f"<ClusterEnrichmentResult clusters={len(self.per_cluster)} "
            f"terms={len(self.rows)}>"
        )

    def __getitem__(self, cluster: str) -> EnrichmentResult:
        return self.per_cluster[str(cluster)]

    @property
    def clusters(self) -> List[str]:
        """Cluster identifiers, in insertion order."""
        return list(self.per_cluster)

    def significant(
        self, alpha: float = 0.05, use: str = "fdr"
    ) -> "ClusterEnrichmentResult":
        """Keep only terms passing a significance threshold, per cluster."""
        return ClusterEnrichmentResult(
            {
                cluster: result.significant(alpha, use)
                for cluster, result in self.per_cluster.items()
            }
        )

    def top(self, n: int) -> "ClusterEnrichmentResult":
        """Keep the ``n`` most significant terms of each cluster."""
        return ClusterEnrichmentResult(
            {cluster: result.top(n) for cluster, result in self.per_cluster.items()}
        )

    def summary(self, alpha: float = 0.05) -> str:
        """
        Build a human-readable per-cluster summary.

        Args:
            alpha: FDR threshold used to count significant terms.

        Returns:
            A multi-line report suitable for printing or writing to a file.
        """
        lines = ["Enrichment summary", "=" * 60, ""]
        for cluster, result in self.per_cluster.items():
            significant = result.significant(alpha)
            lines.append(f"Cluster {cluster}:")
            lines.append(f"  terms tested: {len(result)}")
            lines.append(f"  significant (FDR < {alpha}): {len(significant)}")
            for row in significant.rows[:5]:
                lines.append(f"    - {row.term_label} [{row.term_id}]")
                lines.append(
                    f"      FDR={row.fdr:.2e}  genes={row.overlap_count}/{row.term_size}"
                    f"  fold={row.fold_enrichment:.2f}"
                )
            lines.append("")
        return "\n".join(lines)


def enrich(
    query_genes: Sequence[str],
    library: GeneSetLibrary,
    background: Optional[Iterable[str]] = None,
    min_set_size: int = 5,
    max_set_size: Optional[int] = 500,
    min_overlap: int = 1,
    cluster: Optional[str] = None,
) -> EnrichmentResult:
    """
    Test a gene list for over-representation of the library's terms.

    Args:
        query_genes: Genes of interest, in the same spelling as the library's
            gene sets (gene symbols when the library was built from symbols).
        library: Gene sets to test against.
        background: Universe of genes. Defaults to every gene that appears in
            the library, i.e. all genes that TogoID could annotate.
        min_set_size: Skip terms with fewer background genes than this.
        max_set_size: Skip terms with more background genes than this; ``None``
            disables the upper bound.
        min_overlap: Skip terms with fewer overlapping genes than this.
        cluster: Optional cluster label recorded on every row.

    Returns:
        An ``EnrichmentResult`` sorted by FDR, then p-value.
    """
    background_set: Set[str] = set(background) if background is not None else library.genes
    # Testing genes that are not in the universe would inflate the query size
    # relative to the background and bias every p-value downwards.
    query_set = {gene for gene in query_genes if gene in background_set}

    background_size = len(background_set)
    query_size = len(query_set)

    if background_size == 0 or query_size == 0:
        return EnrichmentResult([], sorted(query_set))

    rows: List[EnrichmentRow] = []
    for term_id, members in library.sets.items():
        # Restrict each gene set to the background, so term sizes and the
        # background agree with one another.
        term_genes = members & background_set
        term_size = len(term_genes)
        if term_size < min_set_size:
            continue
        if max_set_size is not None and term_size > max_set_size:
            continue

        overlap = query_set & term_genes
        overlap_count = len(overlap)
        if overlap_count < min_overlap:
            continue

        pvalue = hypergeometric_sf(overlap_count, background_size, term_size, query_size)
        rows.append(
            EnrichmentRow(
                term_id=term_id,
                term_label=library.label(term_id),
                overlap_count=overlap_count,
                term_size=term_size,
                query_size=query_size,
                background_size=background_size,
                pvalue=pvalue,
                fdr=1.0,  # replaced below, once every p-value is known
                fold_enrichment=fold_enrichment(
                    overlap_count, background_size, term_size, query_size
                ),
                genes=sorted(overlap),
                cluster=str(cluster) if cluster is not None else None,
            )
        )

    if rows:
        for row, adjusted in zip(rows, benjamini_hochberg([r.pvalue for r in rows])):
            row.fdr = adjusted
        rows.sort(key=lambda r: (r.fdr, r.pvalue, -r.fold_enrichment))

    return EnrichmentResult(rows, sorted(query_set))


def enrich_clusters(
    cluster_genes: Mapping[str, Sequence[str]],
    library: GeneSetLibrary,
    background: Optional[Iterable[str]] = None,
    min_set_size: int = 5,
    max_set_size: Optional[int] = 500,
    min_overlap: int = 1,
    verbose: bool = False,
) -> ClusterEnrichmentResult:
    """
    Run :func:`enrich` for every cluster against a shared background.

    Args:
        cluster_genes: Mapping of cluster label to its gene list.
        library: Gene sets to test against.
        background: Shared universe of genes; defaults to the library's genes.
            Sharing one background across clusters is what makes the FDRs
            comparable between clusters.
        min_set_size: Minimum term size.
        max_set_size: Maximum term size, or ``None``.
        min_overlap: Minimum overlap required to report a term.
        verbose: Print a one-line progress report per cluster.

    Returns:
        A ``ClusterEnrichmentResult`` keyed by cluster label.
    """
    background_set: Set[str] = set(background) if background is not None else library.genes

    per_cluster: Dict[str, EnrichmentResult] = {}
    for cluster in sorted(cluster_genes, key=_cluster_sort_key):
        result = enrich(
            cluster_genes[cluster],
            library,
            background=background_set,
            min_set_size=min_set_size,
            max_set_size=max_set_size,
            min_overlap=min_overlap,
            cluster=str(cluster),
        )
        per_cluster[str(cluster)] = result
        if verbose:
            significant = len(result.significant(0.05))
            print(
                f"  cluster {cluster}: {len(result)} terms tested, "
                f"{significant} significant (FDR < 0.05)"
            )

    return ClusterEnrichmentResult(per_cluster)


def _cluster_sort_key(cluster: Any):
    """Sort clusters numerically when they look like numbers, else as text."""
    text = str(cluster)
    try:
        return (0, float(text), "")
    except ValueError:
        return (1, 0.0, text)
