"""
Enrichment analysis and UMAP visualisation built on TogoID ID conversion.

A TogoID route that ends in an annotation dataset turns a gene list into a
gene-set library, which can then be tested for over-representation and drawn on
a single-cell embedding::

    from togoid.enrichment import build_gene_sets, enrich_clusters, plot_umap_enrichment

    library = build_gene_sets(
        genes=all_marker_genes,
        route=["ncbigene", "uniprot", "reactome_pathway"],
    )
    results = enrich_clusters(markers_per_cluster, library)
    fig = plot_umap_enrichment(umap_df, results, top_n=3)
    fig.savefig("umap_enrichment.pdf")

The core (``stats``, ``genesets``, ``analysis``) needs nothing beyond the
standard library and ``requests``. pandas unlocks the ``to_dataframe()`` views,
matplotlib the plots, and scanpy only the AnnData adapters.
"""
from .analysis import (
    RESULT_COLUMNS,
    ClusterEnrichmentResult,
    EnrichmentResult,
    EnrichmentRow,
    enrich,
    enrich_clusters,
)
from .genesets import GeneSetLibrary, build_gene_sets, map_labels_to_ids
from .plot import (
    cluster_centroids,
    plot_umap_centroids,
    plot_umap_enrichment,
    select_terms,
)
from .presets import (
    GO_ASPECTS,
    ROUTES,
    gene_sets_from_preset,
    go_gene_sets,
    mondo_gene_sets,
    reactome_gene_sets,
)
from .stats import benjamini_hochberg, fold_enrichment, hypergeometric_sf

__all__ = [
    # gene sets
    "GeneSetLibrary",
    "build_gene_sets",
    "map_labels_to_ids",
    # presets
    "ROUTES",
    "GO_ASPECTS",
    "reactome_gene_sets",
    "go_gene_sets",
    "mondo_gene_sets",
    "gene_sets_from_preset",
    # analysis
    "enrich",
    "enrich_clusters",
    "EnrichmentResult",
    "ClusterEnrichmentResult",
    "EnrichmentRow",
    "RESULT_COLUMNS",
    # statistics
    "hypergeometric_sf",
    "benjamini_hochberg",
    "fold_enrichment",
    # visualisation
    "plot_umap_enrichment",
    "plot_umap_centroids",
    "cluster_centroids",
    "select_terms",
]
