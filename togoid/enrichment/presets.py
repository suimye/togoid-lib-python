"""
Ready-made routes for the annotation databases used most often.

These are thin wrappers around :func:`~togoid.enrichment.genesets.build_gene_sets`.
They exist for convenience and as worked examples; any other TogoID route works
exactly the same way by passing ``route=[...]`` directly.
"""
from typing import Any, Optional, Sequence

from .genesets import GeneSetLibrary, build_gene_sets

__all__ = [
    "ROUTES",
    "GO_ASPECTS",
    "reactome_gene_sets",
    "go_gene_sets",
    "mondo_gene_sets",
    "gene_sets_from_preset",
]

#: Named routes. ``ncbigene`` is the source in each case because gene symbols
#: resolve to NCBI Gene IDs through TogoID's label resolver.
ROUTES = {
    "reactome": ["ncbigene", "uniprot", "reactome_pathway"],
    "go": ["ncbigene", "uniprot", "go"],
    "mondo": ["ncbigene", "medgen", "mondo"],
}

#: Accepted values of the ``go_aspect`` annotation field.
GO_ASPECTS = ("biological_process", "molecular_function", "cellular_component")


def reactome_gene_sets(
    genes: Sequence[str],
    taxonomy: Optional[str] = "9606",
    **kwargs: Any,
) -> GeneSetLibrary:
    """
    Build Reactome pathway gene sets (``ncbigene -> uniprot -> reactome_pathway``).

    Args:
        genes: Gene symbols, or ``route[0]`` IDs when ``id_source=None`` is passed.
        taxonomy: Taxonomy ID for label resolution.
        **kwargs: Forwarded to :func:`build_gene_sets`.

    Returns:
        A ``GeneSetLibrary`` of Reactome pathways.
    """
    return build_gene_sets(
        genes, route=ROUTES["reactome"], taxonomy=taxonomy, **kwargs
    )


def go_gene_sets(
    genes: Sequence[str],
    aspect: Optional[str] = "biological_process",
    taxonomy: Optional[str] = "9606",
    **kwargs: Any,
) -> GeneSetLibrary:
    """
    Build Gene Ontology gene sets (``ncbigene -> uniprot -> go``).

    Args:
        genes: Gene symbols, or ``route[0]`` IDs when ``id_source=None`` is passed.
        aspect: Restrict to one GO aspect — ``"biological_process"`` (default),
            ``"molecular_function"`` or ``"cellular_component"``. Pass ``None``
            to keep all three.
        taxonomy: Taxonomy ID for label resolution.
        **kwargs: Forwarded to :func:`build_gene_sets`.

    Returns:
        A ``GeneSetLibrary`` of GO terms.

    Raises:
        ValueError: If ``aspect`` is not a recognised GO aspect.
    """
    if aspect is not None and aspect not in GO_ASPECTS:
        raise ValueError(
            f"aspect must be one of {GO_ASPECTS} or None, got {aspect!r}"
        )

    term_filters = kwargs.pop("term_filters", None)
    if aspect is not None and term_filters is None:
        term_filters = {"go_aspect": [aspect]}

    return build_gene_sets(
        genes,
        route=ROUTES["go"],
        taxonomy=taxonomy,
        term_filters=term_filters,
        **kwargs,
    )


def mondo_gene_sets(
    genes: Sequence[str],
    taxonomy: Optional[str] = "9606",
    **kwargs: Any,
) -> GeneSetLibrary:
    """
    Build MONDO disease gene sets (``ncbigene -> medgen -> mondo``).

    Args:
        genes: Gene symbols, or ``route[0]`` IDs when ``id_source=None`` is passed.
        taxonomy: Taxonomy ID for label resolution.
        **kwargs: Forwarded to :func:`build_gene_sets`.

    Returns:
        A ``GeneSetLibrary`` of MONDO diseases.
    """
    return build_gene_sets(genes, route=ROUTES["mondo"], taxonomy=taxonomy, **kwargs)


def gene_sets_from_preset(
    preset: str, genes: Sequence[str], **kwargs: Any
) -> GeneSetLibrary:
    """
    Dispatch to a preset by name.

    Args:
        preset: ``"reactome"``, ``"go"`` or ``"mondo"``.
        genes: Gene symbols, or ``route[0]`` IDs when ``id_source=None`` is passed.
        **kwargs: Forwarded to the preset function.

    Returns:
        A ``GeneSetLibrary``.

    Raises:
        ValueError: If the preset name is unknown.
    """
    builders = {
        "reactome": reactome_gene_sets,
        "go": go_gene_sets,
        "mondo": mondo_gene_sets,
    }
    if preset not in builders:
        raise ValueError(
            f"Unknown preset {preset!r}. Available: {sorted(builders)}. "
            "For any other database, pass route=[...] to build_gene_sets()."
        )
    return builders[preset](genes, **kwargs)
