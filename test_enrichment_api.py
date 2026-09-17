#!/usr/bin/env python3
"""
Tests for togoid.enrichment that call the TogoID API.

These check that the documented routes still work and that the results are
biologically sensible: the classic T cell, B cell and myeloid marker genes must
come out enriched for the matching pathways, GO terms and diseases.

Run:
    python3 test_enrichment_api.py

Requires network access to https://api.togoid.dbcls.jp/.
"""
import sys
import warnings

from togoid.enrichment import (
    GeneSetLibrary,
    build_gene_sets,
    enrich_clusters,
    go_gene_sets,
    map_labels_to_ids,
    mondo_gene_sets,
    reactome_gene_sets,
)
from togoid.enrichment.presets import ROUTES

passed = 0
failed = 0
warned = 0

#: Well-known PBMC lineage markers, used as the query throughout.
CLUSTERS = {
    "T": ["CD3D", "CD3E", "CD3G", "IL7R", "LCK", "ZAP70", "CD2", "CD28", "LAT", "TRAC"],
    "B": [
        "MS4A1", "CD79A", "CD79B", "CD19", "BLNK",
        "BANK1", "TNFRSF13C", "IGHM", "CR2", "PAX5",
    ],
    "Myeloid": [
        "LYZ", "CD14", "FCGR3A", "CSF1R", "ITGAM",
        "TLR2", "TLR4", "FCN1", "S100A8", "S100A9",
    ],
}
ALL_GENES = sorted({gene for genes in CLUSTERS.values() for gene in genes})


def check(condition: bool, message: str) -> None:
    """Record one assertion."""
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✓ {message}")
    else:
        failed += 1
        print(f"  ✗ {message}")


def warn(message: str) -> None:
    """Record a non-fatal problem, usually an API hiccup."""
    global warned
    warned += 1
    print(f"  ⚠ {message}")


def section(title: str) -> None:
    """Print a test-section header."""
    print(f"\n=== {title} ===")


def top_labels(results, cluster: str, n: int = 10):
    """Lower-cased labels of a cluster's most significant terms."""
    return [row.term_label.lower() for row in results[cluster].rows[:n]]


# --------------------------------------------------------------------- #


def test_label_resolution() -> None:
    """Gene symbols must resolve to NCBI Gene IDs."""
    section("Label resolution")

    resolved = map_labels_to_ids(["CD3D", "MS4A1", "LYZ"], verbose=False)
    check(len(resolved) == 3, "resolves three known gene symbols")
    check(resolved.get("CD3D") == "915", "CD3D resolves to NCBI Gene 915")
    check(resolved.get("MS4A1") == "931", "MS4A1 resolves to NCBI Gene 931")

    mixed = map_labels_to_ids(["CD3D", "NOT_A_REAL_GENE_XYZ"], verbose=False)
    check("CD3D" in mixed, "a valid symbol still resolves alongside an invalid one")
    check(
        "NOT_A_REAL_GENE_XYZ" not in mixed,
        "an unknown symbol is simply absent from the mapping",
    )


def test_reactome() -> None:
    """The Reactome route must recover TCR, BCR and neutrophil biology."""
    section("Reactome pathways")

    library = reactome_gene_sets(ALL_GENES, verbose=False)
    check(len(library) > 20, f"builds a gene-set library ({len(library)} pathways)")
    check(library.route == ROUTES["reactome"], "records the route used")
    check(library.target_dataset == "reactome_pathway", "records the target dataset")
    check(len(library.unmapped) == 0, "every marker gene resolved")
    check(
        all(term.startswith("R-") for term in list(library.sets)[:10]),
        "term IDs look like Reactome stable IDs",
    )
    check(
        sum(1 for term in library.sets if term in library.labels) > 0.9 * len(library),
        "almost every term has a label",
    )

    results = enrich_clusters(CLUSTERS, library, min_set_size=3, max_set_size=500)
    check(len(results) > 0, "produces enrichment results")

    t_labels = " | ".join(top_labels(results, "T"))
    check("tcr" in t_labels or "t cell" in t_labels, f"T cluster finds TCR biology")

    b_labels = " | ".join(top_labels(results, "B"))
    check(
        "b cell" in b_labels or "bcr" in b_labels, "B cluster finds B cell receptor biology"
    )

    m_labels = " | ".join(top_labels(results, "Myeloid"))
    check(
        "neutrophil" in m_labels or "toll" in m_labels or "myd88" in m_labels,
        "myeloid cluster finds innate immune biology",
    )

    significant = results.significant(0.05)
    check(len(significant) > 0, f"{len(significant)} terms reach FDR < 0.05")
    check(
        all(row.fdr >= row.pvalue - 1e-15 for row in results),
        "every FDR is at least its raw p-value",
    )
    check(
        all(row.overlap_count <= row.term_size for row in results),
        "no overlap exceeds its term size",
    )
    check(
        all(row.overlap_count <= row.query_size for row in results),
        "no overlap exceeds its query size",
    )


def test_go() -> None:
    """The GO route must work, and the aspect filter must bite."""
    section("GO terms")

    biological = go_gene_sets(ALL_GENES, aspect="biological_process", verbose=False)
    check(len(biological) > 50, f"builds BP gene sets ({len(biological)} terms)")

    results = enrich_clusters(CLUSTERS, biological, min_set_size=3)
    t_labels = " | ".join(top_labels(results, "T"))
    check("t cell" in t_labels, "T cluster finds T cell processes")

    b_labels = " | ".join(top_labels(results, "B"))
    check("b cell" in b_labels, "B cluster finds B cell processes")

    cellular = go_gene_sets(ALL_GENES, aspect="cellular_component", verbose=False)
    check(len(cellular) > 0, f"builds CC gene sets ({len(cellular)} terms)")
    overlap = set(biological.sets) & set(cellular.sets)
    check(
        len(overlap) == 0,
        f"the aspect filter partitions the terms (overlap: {len(overlap)})",
    )

    unfiltered = go_gene_sets(ALL_GENES, aspect=None, verbose=False)
    check(
        len(unfiltered) > len(biological),
        "aspect=None keeps more terms than a single aspect",
    )


def test_mondo() -> None:
    """The MONDO route must return disease terms."""
    section("MONDO diseases")

    library = mondo_gene_sets(ALL_GENES, verbose=False)
    if len(library) == 0:
        warn("MONDO route returned no gene sets for these genes")
        return

    check(len(library) > 0, f"builds disease gene sets ({len(library)} terms)")
    check(library.route == ROUTES["mondo"], "records the MONDO route")
    check(len(library.labels) > 0, "disease terms carry labels")

    results = enrich_clusters(CLUSTERS, library, min_set_size=2)
    check(len(results) >= 0, "enrichment runs against the disease sets")

    # Disease annotation is sparse, so this is informational rather than a check.
    labels = [row.term_label for row in results.rows[:3]]
    if labels:
        print(f"  - top disease terms: {labels}")
    else:
        warn("no MONDO terms passed the size filter (disease coverage is sparse)")


def test_generic_route() -> None:
    """A route passed explicitly must behave like the presets."""
    section("Generic routes")

    library = build_gene_sets(
        ALL_GENES,
        route=["ncbigene", "uniprot", "reactome_pathway"],
        verbose=False,
    )
    preset = reactome_gene_sets(ALL_GENES, verbose=False)
    check(
        set(library.sets) == set(preset.sets),
        "an explicit route matches the equivalent preset",
    )

    # A route that does not exist must fail loudly, not return an empty library.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            build_gene_sets(["CD3D"], route=["ncbigene", "hp_phenotype"], verbose=False)
            check(False, "a broken route raises")
        except RuntimeError as exc:
            check(
                "failed for all" in str(exc),
                "a broken route raises RuntimeError rather than returning nothing",
            )

    try:
        build_gene_sets(ALL_GENES, route=["ncbigene"], verbose=False)
        check(False, "a one-dataset route raises")
    except ValueError:
        check(True, "a one-dataset route raises ValueError")


def test_ids_as_input() -> None:
    """id_source=None must accept identifiers directly."""
    section("Identifier input")

    library = build_gene_sets(
        ["915", "931", "3320"],  # CD3D, MS4A1, HSP90AA1
        route=ROUTES["reactome"],
        id_source=None,
        verbose=False,
    )
    check(len(library) > 0, f"builds sets straight from NCBI Gene IDs ({len(library)} terms)")
    check(len(library.unmapped) == 0, "nothing is reported unmapped")
    check("915" in library.genes, "the input IDs are the set members")


def test_prefix_tolerance() -> None:
    """Results must be identical whether or not the API prefixes IDs."""
    section("CURIE prefix tolerance")

    from togoid import TogoIDConverter

    converter = TogoIDConverter()
    ids = ["915", "931"]

    try:
        prefixed = converter.convert(
            ids=ids, route=ROUTES["reactome"], format="table", prefix=True
        )
        raw = converter.convert(
            ids=ids, route=ROUTES["reactome"], format="table", prefix=False
        )
    except Exception as exc:  # noqa: BLE001 - API availability, not our logic
        warn(f"could not compare prefix modes: {exc}")
        return

    check(len(prefixed) == len(raw), "both prefix modes return the same row count")

    from togoid._ids import local_id

    normalised_prefixed = {(local_id(r[0]), local_id(r[-1])) for r in prefixed}
    normalised_raw = {(local_id(r[0]), local_id(r[-1])) for r in raw}
    check(
        normalised_prefixed == normalised_raw,
        "local_id() normalisation makes both modes identical",
    )


def test_caching() -> None:
    """A saved library must reproduce the same enrichment offline."""
    section("Gene-set caching")

    import os
    import tempfile

    library = reactome_gene_sets(ALL_GENES, verbose=False)
    live = enrich_clusters(CLUSTERS, library, min_set_size=3)

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "reactome.json")
        library.save_json(path)
        cached = enrich_clusters(
            CLUSTERS, GeneSetLibrary.load_json(path), min_set_size=3
        )

    check(len(cached) == len(live), "the cached library yields the same row count")
    check(
        [(r.cluster, r.term_id) for r in cached] == [(r.cluster, r.term_id) for r in live],
        "the cached library yields identical terms in the same order",
    )
    check(
        all(
            abs(a.pvalue - b.pvalue) < 1e-15
            for a, b in zip(cached.rows, live.rows)
        ),
        "the cached library yields identical p-values",
    )


def test_plot_with_real_results() -> None:
    """The figure must build from real API results."""
    section("Plotting real results")

    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        print("  - matplotlib not installed, skipping")
        return

    from togoid.enrichment import plot_umap_enrichment

    library = reactome_gene_sets(ALL_GENES, verbose=False)
    results = enrich_clusters(CLUSTERS, library, min_set_size=3)

    # A stand-in embedding: three well-separated blobs, one per cluster.
    embedding = {"umap_1": [], "umap_2": [], "cluster": []}
    for cluster, (cx, cy) in zip(CLUSTERS, [(0, 0), (10, 2), (5, -9)]):
        for index in range(50):
            embedding["umap_1"].append(cx + (index % 10) * 0.2)
            embedding["umap_2"].append(cy + (index // 10) * 0.2)
            embedding["cluster"].append(cluster)

    figure = plot_umap_enrichment(embedding, results, top_n=3, fdr_cutoff=0.05)
    check(figure is not None, "builds a figure from real enrichment results")
    check(len(figure.axes[1].texts) > 0, "draws pathway labels")


def main() -> int:
    """Run every API-backed test."""
    print("=" * 60)
    print("togoid.enrichment - TogoID API tests")
    print("=" * 60)
    print("These tests require network access to https://api.togoid.dbcls.jp/")

    tests = [
        test_label_resolution,
        test_reactome,
        test_go,
        test_mondo,
        test_generic_route,
        test_ids_as_input,
        test_prefix_tolerance,
        test_caching,
        test_plot_with_real_results,
    ]

    for test in tests:
        try:
            test()
        except Exception as exc:  # noqa: BLE001 - report, then run the rest
            global failed
            failed += 1
            print(f"  ✗ {test.__name__} raised: {exc}")

    print("\n" + "=" * 60)
    print(f"Test Results: {passed} passed, {failed} failed, {warned} warnings")
    print("=" * 60)
    if failed:
        print("Some tests failed. Check the TogoID API status before assuming a bug.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
