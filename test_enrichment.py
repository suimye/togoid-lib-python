#!/usr/bin/env python3
"""
Tests for togoid.enrichment that need no network access.

Covers the statistics, the gene-set container, the enrichment logic and the
word-cloud layout. Where SciPy and statsmodels are installed the statistics are
compared against them; otherwise they are compared against values computed by
hand, so the suite still means something on a bare install.

Run:
    python3 test_enrichment.py
"""
import json
import math
import os
import sys
import tempfile
import warnings

from togoid.enrichment import (
    GeneSetLibrary,
    benjamini_hochberg,
    cluster_centroids,
    enrich,
    enrich_clusters,
    fold_enrichment,
    hypergeometric_sf,
    select_terms,
)
from togoid.enrichment.plot import _overlaps, spiral_positions
from togoid.enrichment.presets import ROUTES, gene_sets_from_preset
from togoid.enrichment.stats import hypergeometric_pmf

passed = 0
failed = 0


def check(condition: bool, message: str) -> None:
    """Record one assertion."""
    global passed, failed
    if condition:
        passed += 1
        print(f"  ✓ {message}")
    else:
        failed += 1
        print(f"  ✗ {message}")


def section(title: str) -> None:
    """Print a test-section header."""
    print(f"\n=== {title} ===")


# --------------------------------------------------------------------- #


def test_hypergeometric() -> None:
    """The hypergeometric survival function must be exact."""
    section("Hypergeometric test")

    # P(X >= 1) with N=4, M=2, n=2 is 1 - C(2,0)C(2,2)/C(4,2) = 1 - 1/6 = 5/6.
    check(
        abs(hypergeometric_sf(1, 4, 2, 2) - 5 / 6) < 1e-12,
        "matches a hand-computed value",
    )
    check(hypergeometric_sf(0, 100, 10, 10) == 1.0, "P(X >= 0) is 1")
    check(hypergeometric_sf(-3, 100, 10, 10) == 1.0, "negative overlap gives 1")
    check(
        hypergeometric_sf(11, 100, 10, 10) == 0.0,
        "overlap beyond the support gives 0",
    )
    check(
        abs(hypergeometric_sf(10, 100, 10, 10) - hypergeometric_pmf(10, 100, 10, 10))
        < 1e-18,
        "the last point of the support equals its pmf",
    )

    # The pmf must sum to 1 over the whole support.
    total = sum(hypergeometric_pmf(k, 50, 20, 10) for k in range(0, 11))
    check(abs(total - 1.0) < 1e-12, "the pmf sums to 1 over its support")

    # Monotonicity: a larger overlap can never be more probable.
    values = [hypergeometric_sf(k, 2000, 100, 200) for k in range(1, 30)]
    check(
        all(a >= b for a, b in zip(values, values[1:])),
        "the survival function decreases with overlap",
    )

    check(
        0.0 <= hypergeometric_sf(300, 20000, 500, 800) <= 1.0,
        "stays within [0, 1] at genome scale",
    )

    try:
        from scipy.stats import hypergeom

        cases = [
            (5, 20000, 150, 300),
            (3, 50, 20, 10),
            (25, 20000, 400, 1200),
            (40, 20000, 500, 800),
            (7, 300, 25, 30),
            (2, 1000, 3, 4),
        ]
        worst = 0.0
        for k, N, M, n in cases:
            mine = hypergeometric_sf(k, N, M, n)
            reference = float(hypergeom.sf(k - 1, N, M, n))
            if reference > 0:
                worst = max(worst, abs(mine - reference) / reference)
        check(worst < 1e-9, f"agrees with scipy (worst relative error {worst:.2e})")
    except ImportError:
        print("  - scipy not installed, skipping the comparison")


def test_benjamini_hochberg() -> None:
    """FDR adjustment must be monotone and match statsmodels."""
    section("Benjamini-Hochberg correction")

    check(benjamini_hochberg([]) == [], "handles an empty input")
    check(benjamini_hochberg([0.03]) == [0.03], "a single p-value is unchanged")

    pvalues = [0.01, 0.04, 0.03, 0.005, 0.2, 0.9, 1e-8]
    adjusted = benjamini_hochberg(pvalues)
    check(len(adjusted) == len(pvalues), "returns one value per input")
    check(all(0.0 <= a <= 1.0 for a in adjusted), "all values stay within [0, 1]")
    check(
        all(a >= p - 1e-12 for a, p in zip(adjusted, pvalues)),
        "never reports an FDR below the raw p-value",
    )

    # Monotonicity: sorting by raw p-value must not un-sort the adjusted values.
    by_p = sorted(zip(pvalues, adjusted))
    check(
        all(a <= b + 1e-12 for (_, a), (_, b) in zip(by_p, by_p[1:])),
        "is monotone in the raw p-value",
    )

    try:
        from statsmodels.stats.multitest import multipletests

        reference = list(multipletests(pvalues, method="fdr_bh")[1])
        worst = max(abs(a - b) for a, b in zip(adjusted, reference))
        check(worst < 1e-12, f"agrees with statsmodels (worst error {worst:.2e})")
    except ImportError:
        print("  - statsmodels not installed, skipping the comparison")


def test_fold_enrichment() -> None:
    """Fold enrichment is the observed overlap over the expected one."""
    section("Fold enrichment")
    # Expected overlap = n*M/N = 100*100/1000 = 10, so 20 observed is 2-fold.
    check(abs(fold_enrichment(20, 1000, 100, 100) - 2.0) < 1e-12, "computes 2-fold")
    check(fold_enrichment(5, 0, 10, 10) == 0.0, "an empty background gives 0")
    check(fold_enrichment(5, 100, 0, 10) == 0.0, "an empty term gives 0")


def build_library() -> GeneSetLibrary:
    """A small hand-made library used by several tests."""
    return GeneSetLibrary(
        sets={
            "T:1": {"CD3D", "CD3E", "CD3G", "LCK", "ZAP70"},
            "T:2": {"CD3D", "CD3E", "LAT"},
            "B:1": {"MS4A1", "CD79A", "CD79B", "CD19", "BLNK"},
            "BIG": {f"G{i}" for i in range(600)},
            "TINY": {"CD3D"},
        },
        labels={"T:1": "TCR signalling", "B:1": "BCR signalling", "BIG": "Huge set"},
        route=["ncbigene", "uniprot", "demo"],
    )


def test_gene_set_library() -> None:
    """The container behaves and round-trips through JSON."""
    section("GeneSetLibrary")

    library = build_library()
    check(len(library) == 5, "reports the number of sets")
    check("T:1" in library, "supports the `in` operator")
    check(library.label("T:1") == "TCR signalling", "returns a known label")
    check(library.label("T:2") == "T:2", "falls back to the ID for unknown labels")
    check("CD3D" in library.genes, "collects every gene")
    check(
        library.gene_to_terms["CD3D"] == {"T:1", "T:2", "TINY"},
        "builds the reverse gene-to-term mapping",
    )

    filtered = library.filter_by_size(min_size=3, max_size=500)
    check(
        set(filtered.sets) == {"T:1", "T:2", "B:1"},
        "filter_by_size drops sets outside the range",
    )
    check(len(library) == 5, "filter_by_size leaves the original untouched")
    check("TINY" not in filtered.labels, "filter_by_size prunes labels too")

    rows = library.to_rows()
    check(rows[0]["term_id"] == "BIG", "to_rows sorts by set size")
    check(rows[0]["n_genes"] == 600, "to_rows reports the set size")

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "nested", "library.json")
        library.save_json(path)
        check(os.path.exists(path), "save_json creates missing directories")

        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
        check(payload["route"] == library.route, "save_json records the route")

        restored = GeneSetLibrary.load_json(path)
        check(restored.sets == library.sets, "load_json restores the sets")
        check(restored.labels == library.labels, "load_json restores the labels")
        check(restored.route == library.route, "load_json restores the route")


def test_enrich() -> None:
    """Over-representation analysis on a hand-made library."""
    section("enrich()")

    library = build_library()
    query = ["CD3D", "CD3E", "CD3G", "LCK", "ZAP70"]

    result = enrich(query, library, min_set_size=3, max_set_size=500)
    terms = [row.term_id for row in result]
    check("T:1" in terms, "finds the fully overlapping term")
    check("BIG" not in terms, "excludes sets above max_set_size")
    check("TINY" not in terms, "excludes sets below min_set_size")

    top = result.rows[0]
    check(top.term_id == "T:1", "ranks the best term first")
    check(top.overlap_count == 5, "reports the overlap")
    check(top.term_size == 5, "reports the term size")
    check(top.genes == sorted(query), "lists the overlapping genes")
    check(top.term_label == "TCR signalling", "attaches the term label")
    check(top.fold_enrichment > 1.0, "reports enrichment above 1")
    check(
        all(a.fdr <= b.fdr + 1e-15 for a, b in zip(result.rows, result.rows[1:])),
        "rows are sorted by FDR",
    )

    # Restricting the background changes the size the p-value is computed against.
    background = library.genes
    wider = enrich(query, library, background=list(background) + ["X1", "X2"], min_set_size=3)
    check(
        wider.rows[0].background_size == len(background) + 2,
        "an explicit background is used as given",
    )
    check(
        wider.rows[0].pvalue < result.rows[0].pvalue,
        "a larger background makes the same overlap more significant",
    )

    # Genes outside the background must not inflate the query size.
    padded = enrich(list(query) + ["NOT_A_GENE"], library, min_set_size=3)
    check(
        padded.rows[0].query_size == result.rows[0].query_size,
        "genes outside the background are dropped from the query",
    )

    check(len(enrich([], library)) == 0, "an empty query gives no rows")
    check(
        len(enrich(query, GeneSetLibrary())) == 0, "an empty library gives no rows"
    )

    filtered = result.significant(alpha=1e-9)
    check(len(filtered) < len(result), "significant() filters on FDR")
    check(len(result.top(1)) == 1, "top() limits the row count")


def test_enrich_clusters() -> None:
    """Per-cluster analysis shares one background."""
    section("enrich_clusters()")

    library = build_library()
    clusters = {
        "1": ["MS4A1", "CD79A", "CD79B", "CD19", "BLNK"],
        "0": ["CD3D", "CD3E", "CD3G", "LCK", "ZAP70"],
    }

    results = enrich_clusters(clusters, library, min_set_size=3)
    check(results.clusters == ["0", "1"], "clusters are ordered numerically")
    check(results["0"].rows[0].term_id == "T:1", "cluster 0 finds the T term")
    check(results["1"].rows[0].term_id == "B:1", "cluster 1 finds the B term")
    check(
        all(row.cluster is not None for row in results),
        "every row records its cluster",
    )
    check(
        len({row.background_size for row in results}) == 1,
        "all clusters share one background size",
    )

    rows = results.to_rows()
    check("cluster" in rows[0], "to_rows includes the cluster column")

    summary = results.summary()
    check("Cluster 0:" in summary and "Cluster 1:" in summary, "summary lists clusters")

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "results.csv")
        results.to_csv(path)
        with open(path, encoding="utf-8") as handle:
            header = handle.readline().strip()
        check(header.startswith("cluster,term_id"), "to_csv writes the cluster column")

        single = os.path.join(directory, "single.csv")
        enrich(clusters["0"], library, min_set_size=3).to_csv(single)
        with open(single, encoding="utf-8") as handle:
            header = handle.readline().strip()
        check(
            header.startswith("term_id"),
            "to_csv omits the cluster column for a single query",
        )


def test_tables() -> None:
    """TSV export and the human-readable table views."""
    section("Table output")

    library = build_library()
    clusters = {
        "0": ["CD3D", "CD3E", "CD3G", "LCK", "ZAP70"],
        "1": ["MS4A1", "CD79A", "CD79B", "CD19", "BLNK"],
    }
    results = enrich_clusters(clusters, library, min_set_size=3)

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "out.tsv")
        results.to_tsv(path)
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().splitlines()

        check(lines[0].split("\t")[0] == "cluster", "to_tsv writes a tab-separated header")
        check(len(lines) == len(results) + 1, "to_tsv writes one line per term")
        check(
            all(len(line.split("\t")) == len(lines[0].split("\t")) for line in lines),
            "every to_tsv row has the same number of fields",
        )
        # A label containing a comma must not need quoting in a TSV.
        check('"' not in "\n".join(lines), "to_tsv needs no quoting")

        wide = os.path.join(directory, "wide.tsv")
        results.write_cluster_table(wide, top_n=2)
        with open(wide, encoding="utf-8") as handle:
            wide_lines = handle.read().splitlines()
        check(
            wide_lines[0].startswith("cluster\tn_tested\tn_significant"),
            "write_cluster_table writes the cluster columns first",
        )
        check(
            len(wide_lines) == len(results.clusters) + 1,
            "write_cluster_table writes one row per cluster",
        )
        check("top2_fdr" in wide_lines[0], "write_cluster_table honours top_n")

    table = results.to_cluster_table(top_n=2)
    check(len(table) == len(results.clusters), "to_cluster_table has one row per cluster")
    check(table[0]["cluster"] == "0", "to_cluster_table keeps the cluster order")
    check(
        table[0]["top1_term_id"] == results["0"].significant().rows[0].term_id,
        "to_cluster_table reports the most significant term first",
    )

    text = results.to_text()
    check("cluster" in text and "term_label" in text, "to_text has a header")
    check(len(text.splitlines()) == len(results) + 2, "to_text writes one line per term")
    check("overlap" in text, "to_text folds the overlap and term size into one column")
    check(
        "background_size" not in text,
        "to_text omits the constant background_size column",
    )
    check(str(results) == text, "printing a result shows the table")

    html = results._repr_html_()
    check(html.startswith("<table"), "_repr_html_ returns a table")
    # One <tr> per term, plus the header row.
    check(
        html.count("<tr>") == len(results) + 1,
        "_repr_html_ has one row per term plus a header",
    )
    check("</table>" in html, "_repr_html_ is closed properly")

    # HTML special characters in a label must be escaped, not injected.
    risky = GeneSetLibrary(
        sets={"X": {"A", "B", "C"}}, labels={"X": "<b>bold</b> & risky"}
    )
    risky_html = enrich(["A", "B", "C"], risky, min_set_size=3)._repr_html_()
    check("&lt;b&gt;" in risky_html, "_repr_html_ escapes HTML in labels")
    check("<b>bold</b>" not in risky_html, "_repr_html_ does not inject raw HTML")

    empty = enrich([], library)
    check(empty.to_text() == "(no enriched terms)", "to_text handles an empty result")
    check("no enriched terms" in empty._repr_html_(), "_repr_html_ handles an empty result")


def test_selected_terms_table() -> None:
    """The exported table must match what the figure draws."""
    section("selected_terms_table()")

    from togoid.enrichment import selected_terms_table

    rows = [
        {"cluster": "1", "term_id": "B", "term_label": "beta", "pvalue": 1e-3, "fdr": 2e-3},
        {"cluster": "0", "term_id": "A", "term_label": "alpha", "pvalue": 1e-6, "fdr": 1e-5},
        {"cluster": "0", "term_id": "C", "term_label": "gamma", "pvalue": 0.4, "fdr": 0.5},
    ]

    selected = select_terms(rows, top_n=None, fdr_cutoff=0.05)
    table = selected_terms_table(selected)

    check(len(table) == 2, "keeps only the terms that pass the filters")
    check([r["cluster"] for r in table] == ["0", "1"], "orders rows by cluster")
    check(
        "label" not in table[0] and "weight" not in table[0],
        "drops the layout-only columns",
    )
    check(table[0]["term_id"] == "A", "keeps the original fields")
    check(selected_terms_table({}) == [], "handles an empty selection")


def test_dataframe_views() -> None:
    """The pandas views mirror the plain rows."""
    section("DataFrame views")

    try:
        import pandas  # noqa: F401
    except ImportError:
        print("  - pandas not installed, skipping")
        return

    library = build_library()
    results = enrich_clusters(
        {"0": ["CD3D", "CD3E", "CD3G", "LCK", "ZAP70"]}, library, min_set_size=3
    )

    frame = results.to_dataframe()
    check("cluster" in frame.columns, "per-cluster frames keep the cluster column")
    check(len(frame) == len(results.rows), "the frame has one row per result")

    single = enrich(["CD3D", "CD3E", "CD3G"], library, min_set_size=3).to_dataframe()
    check(
        "cluster" not in single.columns,
        "single-query frames drop the empty cluster column",
    )

    library_frame = library.to_dataframe()
    check(list(library_frame.columns)[0] == "term_id", "the library frame starts with term_id")


def test_presets() -> None:
    """Preset routes are well formed; no network access involved."""
    section("Presets")

    check(ROUTES["reactome"][-1] == "reactome_pathway", "the Reactome route ends in Reactome")
    check(ROUTES["go"][-1] == "go", "the GO route ends in GO")
    check(ROUTES["mondo"][-1] == "mondo", "the MONDO route ends in MONDO")
    check(
        all(route[0] == "ncbigene" for route in ROUTES.values()),
        "every preset starts from ncbigene",
    )

    try:
        gene_sets_from_preset("nonexistent", ["CD3D"])
        check(False, "an unknown preset raises")
    except ValueError as exc:
        check("Unknown preset" in str(exc), "an unknown preset raises ValueError")

    from togoid.enrichment.presets import go_gene_sets

    try:
        go_gene_sets(["CD3D"], aspect="not_an_aspect")
        check(False, "an invalid GO aspect raises")
    except ValueError as exc:
        check("aspect must be one of" in str(exc), "an invalid GO aspect raises ValueError")


def test_centroids_and_layout() -> None:
    """Centroids, spirals and the overlap test."""
    section("Layout geometry")

    embedding = {
        "umap_1": [0.0, 2.0, 10.0, 12.0],
        "umap_2": [0.0, 2.0, 10.0, 12.0],
        "cluster": ["0", "0", "1", "1"],
    }
    centroids = cluster_centroids(embedding)
    check(centroids["0"] == {"x": 1.0, "y": 1.0, "n_cells": 2}, "computes cluster 0's centroid")
    check(centroids["1"]["x"] == 11.0, "computes cluster 1's centroid")

    try:
        cluster_centroids({"umap_1": [0.0], "umap_2": [0.0]})
        check(False, "a missing column raises")
    except KeyError:
        check(True, "a missing column raises KeyError")

    positions = spiral_positions(0.0, 0.0, 16, radius_start=1.0, radius_step=1.0)
    check(len(positions) == 16, "generates the requested number of positions")
    first = math.hypot(*positions[0])
    last = math.hypot(*positions[-1])
    check(abs(first - 1.0) < 1e-9, "starts at radius_start")
    check(last > first, "spirals outwards")
    check(spiral_positions(0.0, 0.0, 0, 1.0, 1.0) == [], "handles zero positions")

    check(_overlaps((0, 0, 1, 1), (0.5, 0.5, 1.5, 1.5)), "detects overlapping boxes")
    check(not _overlaps((0, 0, 1, 1), (2, 2, 3, 3)), "detects disjoint boxes")
    check(
        _overlaps((0, 0, 1, 1), (1.05, 0, 2, 1), padding=0.1),
        "padding widens the box",
    )
    check(
        not _overlaps((0, 0, 1, 1), (1.05, 0, 2, 1), padding=0.0),
        "no padding leaves a gap intact",
    )


def test_select_terms() -> None:
    """Term selection applies every filter."""
    section("select_terms()")

    rows = [
        {"cluster": "0", "term_id": "A", "term_label": "alpha", "pvalue": 1e-6, "fdr": 1e-5},
        {"cluster": "0", "term_id": "B", "term_label": "beta", "pvalue": 1e-3, "fdr": 2e-3},
        {"cluster": "0", "term_id": "C", "term_label": "gamma", "pvalue": 0.4, "fdr": 0.5},
        {"cluster": "1", "term_id": "D", "term_label": "d" * 60, "pvalue": 1e-4, "fdr": 1e-3},
    ]

    selected = select_terms(rows, top_n=None, fdr_cutoff=0.05)
    check(set(selected) == {"0", "1"}, "groups terms by cluster")
    check(len(selected["0"]) == 2, "drops terms above the FDR cut-off")

    check(
        len(select_terms(rows, top_n=1, fdr_cutoff=0.05)["0"]) == 1,
        "top_n limits the terms per cluster",
    )
    check(
        selected["0"][0]["term_id"] == "A",
        "terms are ordered by significance",
    )
    check(
        len(select_terms(rows, top_n=None, fdr_cutoff=None, pval_cutoff=1e-5)["0"]) == 1,
        "pval_cutoff filters independently",
    )
    check(
        set(select_terms(rows, top_n=None, fdr_cutoff=0.05, clusters=["1"])) == {"1"},
        "clusters restricts the output",
    )

    truncated = select_terms(rows, top_n=None, fdr_cutoff=0.05, max_label_chars=20)["1"][0]
    check(len(truncated["label"]) == 20, "labels are truncated to max_label_chars")
    check(truncated["label"].endswith("..."), "truncated labels end with an ellipsis")

    # weight drives the font size and must survive p == 0.
    zero = select_terms(
        [{"cluster": "0", "term_id": "Z", "term_label": "z", "pvalue": 0.0, "fdr": 0.0}],
        top_n=None,
    )
    check(math.isfinite(zero["0"][0]["weight"]), "a zero p-value gives a finite weight")

    check(select_terms([], top_n=3) == {}, "handles empty input")


def test_plot_smoke() -> None:
    """The figure builds, and labels never overlap."""
    section("plot_umap_enrichment()")

    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        print("  - matplotlib not installed, skipping")
        return

    from togoid.enrichment import plot_umap_centroids, plot_umap_enrichment
    from togoid.enrichment.plot import _overlaps as overlaps

    embedding = {"umap_1": [], "umap_2": [], "cluster": []}
    for index, (cx, cy) in enumerate([(0.0, 0.0), (10.0, 0.0), (5.0, 9.0)]):
        for step in range(60):
            angle = step / 60 * 2 * math.pi
            embedding["umap_1"].append(cx + math.cos(angle))
            embedding["umap_2"].append(cy + math.sin(angle))
            embedding["cluster"].append(str(index))

    rows = []
    for cluster in ("0", "1", "2"):
        for term in range(4):
            rows.append(
                {
                    "cluster": cluster,
                    "term_id": f"{cluster}:{term}",
                    "term_label": f"term {cluster}-{term}",
                    "pvalue": 10 ** -(6 - term),
                    "fdr": 10 ** -(5 - term),
                }
            )

    figure = plot_umap_enrichment(embedding, rows, top_n=4, fdr_cutoff=0.05)
    check(figure is not None, "returns a figure")
    check(len(figure.axes) == 2, "draws two panels")

    texts = [artist for artist in figure.axes[1].texts]
    check(len(texts) > 0, "places term labels on the right panel")

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    boxes = []
    for artist in texts:
        extent = artist.get_window_extent(renderer=renderer)
        boxes.append((extent.x0, extent.y0, extent.x1, extent.y1))
    collisions = sum(
        1
        for i in range(len(boxes))
        for j in range(i + 1, len(boxes))
        if overlaps(boxes[i], boxes[j])
    )
    check(collisions == 0, f"no two labels overlap ({len(boxes)} labels placed)")

    # show_centroids only sets the marker opacity: the marker is still drawn and
    # still reserves its space, so the labels must not move.
    with_markers = plot_umap_enrichment(embedding, rows, top_n=4, fdr_cutoff=0.05)
    without = plot_umap_enrichment(
        embedding, rows, top_n=4, fdr_cutoff=0.05, show_centroids=False
    )

    def centroid_alphas(figure):
        # The centroid scatters are the ones drawn at zorder 5.
        return [
            c.get_alpha()
            for c in figure.axes[1].collections
            if c.get_zorder() == 5
        ]

    check(
        len(centroid_alphas(without)) == len(centroid_alphas(with_markers)) > 0,
        "show_centroids=False still draws the centroid markers",
    )
    check(
        all(a == 0.0 for a in centroid_alphas(without)),
        "show_centroids=False makes the centroid markers transparent",
    )
    check(
        all(a and a > 0.0 for a in centroid_alphas(with_markers)),
        "show_centroids=True makes the centroid markers visible",
    )

    def label_positions(figure):
        return sorted(
            (t.get_text(), round(t.get_position()[0], 6), round(t.get_position()[1], 6))
            for t in figure.axes[1].texts
        )

    check(
        label_positions(without) == label_positions(with_markers),
        "hiding the centroids leaves every label in exactly the same place",
    )

    square = plot_umap_enrichment(
        embedding, rows, top_n=4, fdr_cutoff=0.05, centroid_marker="s"
    )
    check(square is not None, "centroid_marker accepts another marker code")

    empty = plot_umap_enrichment(embedding, rows, fdr_cutoff=1e-300)
    check(len(empty.axes) == 2, "handles the case where nothing is significant")

    check(plot_umap_centroids(embedding) is not None, "plot_umap_centroids returns a figure")
    check(
        len(plot_umap_centroids(embedding, show_labels=False).axes[0].texts)
        < len(plot_umap_centroids(embedding).axes[0].texts),
        "plot_umap_centroids(show_labels=False) omits the cluster labels",
    )

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "figure.pdf")
        figure.savefig(path)
        check(os.path.getsize(path) > 0, "the figure saves to PDF")

    try:
        plot_umap_enrichment({"umap_1": [], "umap_2": [], "cluster": []}, rows)
        check(False, "an empty embedding raises")
    except ValueError:
        check(True, "an empty embedding raises ValueError")


def test_input_flexibility() -> None:
    """The plotting API accepts results in several shapes."""
    section("Input flexibility")

    library = build_library()
    results = enrich_clusters(
        {"0": ["CD3D", "CD3E", "CD3G", "LCK", "ZAP70"]}, library, min_set_size=3
    )

    from_object = select_terms(results, top_n=None, fdr_cutoff=None)
    from_rows = select_terms(results.to_rows(), top_n=None, fdr_cutoff=None)
    check(
        {c: len(v) for c, v in from_object.items()} == {c: len(v) for c, v in from_rows.items()},
        "a result object and its rows give the same selection",
    )

    try:
        import pandas  # noqa: F401

        from_frame = select_terms(results.to_dataframe(), top_n=None, fdr_cutoff=None)
        check(
            {c: len(v) for c, v in from_frame.items()}
            == {c: len(v) for c, v in from_rows.items()},
            "a DataFrame gives the same selection",
        )
    except ImportError:
        print("  - pandas not installed, skipping the DataFrame case")


def main() -> int:
    """Run every offline test."""
    print("=" * 60)
    print("togoid.enrichment - offline tests")
    print("=" * 60)

    # Turn the library's batch-failure warnings into output, not silence.
    warnings.simplefilter("always", RuntimeWarning)

    test_hypergeometric()
    test_benjamini_hochberg()
    test_fold_enrichment()
    test_gene_set_library()
    test_enrich()
    test_enrich_clusters()
    test_tables()
    test_selected_terms_table()
    test_dataframe_views()
    test_presets()
    test_centroids_and_layout()
    test_select_terms()
    test_plot_smoke()
    test_input_flexibility()

    print("\n" + "=" * 60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("=" * 60)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
