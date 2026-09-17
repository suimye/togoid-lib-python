"""
Word-cloud style visualisation of enrichment results on a UMAP embedding.

The layout places each cluster's enriched terms around that cluster's centroid,
sized by significance, and rejects candidate positions that would overlap text
already placed. Nothing here is specific to scanpy or Seurat: the embedding is
passed in as plain coordinates plus cluster labels.
"""
import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

__all__ = [
    "cluster_centroids",
    "select_terms",
    "plot_umap_enrichment",
    "plot_umap_centroids",
    "spiral_positions",
]

#: Bounding box type: (x_min, y_min, x_max, y_max) in data coordinates.
BBox = Tuple[float, float, float, float]


# ---------------------------------------------------------------------- #
# Input normalisation
# ---------------------------------------------------------------------- #


def _as_columns(
    embedding: Any,
    x_key: str = "umap_1",
    y_key: str = "umap_2",
    cluster_key: str = "cluster",
) -> Tuple[List[float], List[float], List[str]]:
    """
    Normalise an embedding into three parallel lists.

    Args:
        embedding: A pandas DataFrame or any mapping of column name to sequence.
        x_key: Column holding the first embedding dimension.
        y_key: Column holding the second embedding dimension.
        cluster_key: Column holding the cluster label.

    Returns:
        Tuple of ``(x, y, cluster)`` lists.

    Raises:
        KeyError: If a required column is missing.
    """
    for key in (x_key, y_key, cluster_key):
        try:
            _ = embedding[key]
        except Exception as exc:  # noqa: BLE001 - surface a clear message
            raise KeyError(
                f"embedding is missing the {key!r} column; "
                f"expected columns {x_key!r}, {y_key!r} and {cluster_key!r}"
            ) from exc

    xs = [float(v) for v in list(embedding[x_key])]
    ys = [float(v) for v in list(embedding[y_key])]
    clusters = [str(v) for v in list(embedding[cluster_key])]
    return xs, ys, clusters


def _result_rows(enrichment: Any) -> List[Dict[str, Any]]:
    """
    Normalise enrichment results into a list of dictionaries.

    Accepts a ``ClusterEnrichmentResult``, an ``EnrichmentResult``, a pandas
    DataFrame, or an already-normalised list of dicts.
    """
    if hasattr(enrichment, "to_rows"):
        return enrichment.to_rows()
    if hasattr(enrichment, "to_dict"):  # pandas DataFrame
        return enrichment.to_dict(orient="records")
    return [dict(row) for row in enrichment]


# ---------------------------------------------------------------------- #
# Geometry
# ---------------------------------------------------------------------- #


def cluster_centroids(
    embedding: Any,
    x_key: str = "umap_1",
    y_key: str = "umap_2",
    cluster_key: str = "cluster",
) -> Dict[str, Dict[str, float]]:
    """
    Compute the centre of mass of every cluster in the embedding.

    Args:
        embedding: DataFrame or mapping with embedding coordinates.
        x_key: Column holding the first embedding dimension.
        y_key: Column holding the second embedding dimension.
        cluster_key: Column holding the cluster label.

    Returns:
        Mapping of cluster label to ``{"x", "y", "n_cells"}``.
    """
    xs, ys, clusters = _as_columns(embedding, x_key, y_key, cluster_key)

    sums: Dict[str, List[float]] = {}
    for x, y, cluster in zip(xs, ys, clusters):
        acc = sums.setdefault(cluster, [0.0, 0.0, 0.0])
        acc[0] += x
        acc[1] += y
        acc[2] += 1

    return {
        cluster: {"x": sx / n, "y": sy / n, "n_cells": int(n)}
        for cluster, (sx, sy, n) in sums.items()
        if n > 0
    }


def spiral_positions(
    center_x: float,
    center_y: float,
    n_positions: int,
    radius_start: float,
    radius_step: float,
    angle_step: float = math.pi / 4,
) -> List[Tuple[float, float]]:
    """
    Generate candidate label positions spiralling out from a centre.

    Args:
        center_x: Centre x coordinate.
        center_y: Centre y coordinate.
        n_positions: Number of candidate positions to generate.
        radius_start: Radius of the first candidate.
        radius_step: Radius added per full turn.
        angle_step: Angle between consecutive candidates, in radians.

    Returns:
        List of ``(x, y)`` candidate positions, closest to the centre first.
    """
    positions: List[Tuple[float, float]] = []
    angle = 0.0
    radius = radius_start
    per_step = radius_step / max(1.0, (2 * math.pi / angle_step))

    for _ in range(max(0, n_positions)):
        positions.append(
            (center_x + radius * math.cos(angle), center_y + radius * math.sin(angle))
        )
        angle += angle_step
        radius += per_step

    return positions


def _contains(outer: Optional[BBox], inner: BBox) -> bool:
    """Return True when ``inner`` lies entirely within ``outer``."""
    if outer is None:
        return True
    ox0, oy0, ox1, oy1 = outer
    ix0, iy0, ix1, iy1 = inner
    return ix0 >= ox0 and iy0 >= oy0 and ix1 <= ox1 and iy1 <= oy1


def _overlaps(a: BBox, b: BBox, padding: float = 0.0) -> bool:
    """Return True when two bounding boxes intersect, after padding ``a``."""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    ax0 -= padding
    ay0 -= padding
    ax1 += padding
    ay1 += padding
    return not (ax1 < bx0 or ax0 > bx1 or ay1 < by0 or ay0 > by1)


# ---------------------------------------------------------------------- #
# Term selection
# ---------------------------------------------------------------------- #


def select_terms(
    enrichment: Any,
    top_n: Optional[int] = 3,
    pval_cutoff: Optional[float] = None,
    fdr_cutoff: Optional[float] = 0.05,
    clusters: Optional[Sequence[str]] = None,
    max_label_chars: Optional[int] = 40,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Pick the terms to draw for each cluster.

    Args:
        enrichment: Enrichment results in any form accepted by this module.
        top_n: Keep at most this many terms per cluster; ``None`` keeps all.
        pval_cutoff: Drop terms with a p-value at or above this value.
        fdr_cutoff: Drop terms with an FDR at or above this value.
        clusters: Restrict to these clusters; ``None`` keeps all.
        max_label_chars: Truncate labels longer than this, adding an ellipsis.

    Returns:
        Mapping of cluster label to its selected terms. Each term carries the
        original fields plus ``label`` (possibly truncated) and ``weight``
        (``-log10(pvalue)``), sorted most significant first.
    """
    wanted = {str(c) for c in clusters} if clusters is not None else None

    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in _result_rows(enrichment):
        cluster = str(row.get("cluster", ""))
        if wanted is not None and cluster not in wanted:
            continue

        pvalue = float(row.get("pvalue", 1.0))
        fdr = float(row.get("fdr", 1.0))
        if pval_cutoff is not None and pvalue >= pval_cutoff:
            continue
        if fdr_cutoff is not None and fdr >= fdr_cutoff:
            continue

        label = str(row.get("term_label") or row.get("term_id") or "")
        if max_label_chars is not None and len(label) > max_label_chars:
            label = label[: max_label_chars - 3] + "..."

        entry = dict(row)
        entry["label"] = label
        # Guard against p == 0 from extreme enrichment, which would be -inf.
        entry["weight"] = -math.log10(pvalue) if pvalue > 0 else 300.0
        grouped.setdefault(cluster, []).append(entry)

    for cluster, entries in grouped.items():
        entries.sort(key=lambda e: (e.get("fdr", 1.0), e.get("pvalue", 1.0)))
        if top_n is not None:
            grouped[cluster] = entries[:top_n]

    return grouped


# ---------------------------------------------------------------------- #
# Plotting
# ---------------------------------------------------------------------- #


def _require_matplotlib():
    """Import matplotlib, with an actionable error when it is missing."""
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "matplotlib is required for plotting. "
            'Install with: pip install "togoid[plot]"'
        ) from exc
    return plt


def _cluster_palette(
    clusters: Sequence[str], palette: Optional[str] = None
) -> Dict[str, Any]:
    """
    Assign a stable colour to every cluster label.

    Args:
        clusters: Cluster label of every cell (duplicates are fine).
        palette: Matplotlib colormap name, or ``None`` to choose automatically.
            ``tab20`` alternates dark and light shades, so its light half is hard
            to read as text; for ten clusters or fewer ``tab10`` is used instead.

    Returns:
        Mapping of cluster label to an RGBA colour.
    """
    from matplotlib import colormaps

    ordered = sorted(set(clusters), key=_cluster_sort_key)
    if palette is None:
        palette = "tab10" if len(ordered) <= 10 else "tab20"
    cmap = colormaps[palette]
    # Qualitative colormaps are indexed by position; continuous ones need the
    # index spread across [0, 1].
    n_colors = getattr(cmap, "N", 256)
    if n_colors <= 32:
        return {c: cmap(i % n_colors) for i, c in enumerate(ordered)}
    span = max(1, len(ordered) - 1)
    return {c: cmap(i / span) for i, c in enumerate(ordered)}


def _cluster_sort_key(cluster: str):
    """Sort cluster labels numerically when possible, otherwise as text."""
    try:
        return (0, float(cluster), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(cluster))


def _scatter_clusters(
    ax,
    xs: Sequence[float],
    ys: Sequence[float],
    clusters: Sequence[str],
    colors: Mapping[str, Any],
    size: float,
    alpha: float,
    label: bool,
) -> None:
    """Draw the embedding, one scatter call per cluster so colours group."""
    by_cluster: Dict[str, Tuple[List[float], List[float]]] = {}
    for x, y, cluster in zip(xs, ys, clusters):
        acc = by_cluster.setdefault(cluster, ([], []))
        acc[0].append(x)
        acc[1].append(y)

    for cluster in sorted(by_cluster, key=_cluster_sort_key):
        cx, cy = by_cluster[cluster]
        ax.scatter(
            cx,
            cy,
            c=[colors[cluster]],
            s=size,
            alpha=alpha,
            linewidths=0,
            label=f"Cluster {cluster}" if label else None,
        )


def _measure_text(ax, renderer, text: str, fontsize: float) -> Tuple[float, float]:
    """
    Return the width and height of a text label in data coordinates.

    The size of a label does not depend on where it is drawn, so it is measured
    once and reused for every candidate position. Measuring per candidate (which
    forces a full canvas draw each time) is what made naive layouts unusable on
    real datasets.
    """
    artist = ax.text(0, 0, text, fontsize=fontsize, alpha=0.0)
    extent = artist.get_window_extent(renderer=renderer)
    data_box = extent.transformed(ax.transData.inverted())
    artist.remove()
    return abs(data_box.x1 - data_box.x0), abs(data_box.y1 - data_box.y0)


def _place_labels(
    ax,
    renderer,
    selected: Mapping[str, Sequence[Mapping[str, Any]]],
    centroids: Mapping[str, Mapping[str, float]],
    colors: Mapping[str, Any],
    fontsize_range: Tuple[float, float],
    weight_scale: float,
    candidates_per_term: int,
    padding: float,
    span: float,
    reserved: Optional[Sequence[BBox]] = None,
    bounds: Optional[BBox] = None,
) -> Tuple[List[BBox], int]:
    """
    Place every selected term around its cluster centroid without overlaps.

    Args:
        reserved: Boxes that are already occupied (e.g. centroid markers) and
            which labels must therefore avoid.
        bounds: Axis limits as ``(x_min, y_min, x_max, y_max)``. Positions that
            keep a label fully inside are tried first, so labels only leave the
            visible area when there is genuinely nowhere else to put them.

    Returns:
        Tuple of ``(label boxes placed, skipped count)``.
    """
    min_size, max_size = fontsize_range
    obstacles: List[BBox] = list(reserved or [])
    label_boxes: List[BBox] = []
    skipped = 0

    # Draw the most significant terms first so they win the space nearest their
    # centroid; later terms settle further out or are dropped.
    ordered_clusters = sorted(
        selected,
        key=lambda c: -max((e.get("weight", 0.0) for e in selected[c]), default=0.0),
    )

    for cluster in ordered_clusters:
        entries = selected[cluster]
        centroid = centroids.get(cluster)
        if not entries or centroid is None:
            continue

        positions = spiral_positions(
            centroid["x"],
            centroid["y"],
            n_positions=max(1, len(entries)) * candidates_per_term,
            radius_start=span * 0.02,
            radius_step=span * 0.06,
        )

        for entry in entries:
            fontsize = min(
                max_size, min_size + float(entry.get("weight", 0.0)) * weight_scale
            )
            width, height = _measure_text(ax, renderer, entry["label"], fontsize)
            half_w, half_h = width / 2.0, height / 2.0

            # Pass 1 keeps the label inside the axes; pass 2 drops that
            # requirement so an edge cluster still gets labelled.
            chosen = None
            for require_inside in (True, False):
                if require_inside and bounds is None:
                    continue
                for px, py in positions:
                    box: BBox = (px - half_w, py - half_h, px + half_w, py + half_h)
                    if require_inside and not _contains(bounds, box):
                        continue
                    if any(
                        _overlaps(box, other, padding)
                        for other in obstacles + label_boxes
                    ):
                        continue
                    chosen = (px, py, box)
                    break
                if chosen is not None:
                    break

            if chosen is None:
                # No free position anywhere on the spiral; dropping the label is
                # better than drawing unreadable overlapping text.
                skipped += 1
                continue

            px, py, box = chosen
            ax.text(
                px,
                py,
                entry["label"],
                fontsize=fontsize,
                color=colors.get(cluster, "black"),
                ha="center",
                va="center",
                alpha=0.9,
                zorder=4,
            )
            label_boxes.append(box)

    return label_boxes, skipped


def plot_umap_enrichment(
    embedding: Any,
    enrichment: Any,
    top_n: Optional[int] = 3,
    pval_cutoff: Optional[float] = None,
    fdr_cutoff: Optional[float] = 0.05,
    clusters: Optional[Sequence[str]] = None,
    max_label_chars: Optional[int] = 40,
    x_key: str = "umap_1",
    y_key: str = "umap_2",
    cluster_key: str = "cluster",
    palette: Optional[str] = None,
    figsize: Tuple[float, float] = (20.0, 8.0),
    point_size: float = 8.0,
    fontsize_range: Tuple[float, float] = (6.0, 14.0),
    weight_scale: float = 0.8,
    candidates_per_term: int = 40,
    padding_fraction: float = 0.004,
    title_left: str = "UMAP clustering",
    title_right: Optional[str] = None,
    show_centroids: bool = True,
    centroid_marker: str = "o",
    centroid_size: float = 26.0,
    centroid_color: str = "black",
    legend: bool = True,
    verbose: bool = False,
):
    """
    Draw a two-panel figure: clusters on the left, enriched terms on the right.

    Args:
        embedding: DataFrame or mapping with the UMAP coordinates and cluster
            labels of every cell.
        enrichment: Enrichment results (``ClusterEnrichmentResult``, DataFrame or
            list of dicts) carrying ``cluster``, ``term_label``, ``pvalue`` and
            ``fdr``.
        top_n: Terms drawn per cluster; ``None`` draws all that pass the cut-offs.
        pval_cutoff: Optional p-value cut-off.
        fdr_cutoff: FDR cut-off; pass ``None`` to disable.
        clusters: Restrict the labels to these clusters.
        max_label_chars: Truncate term labels longer than this.
        x_key: Embedding column for the first dimension.
        y_key: Embedding column for the second dimension.
        cluster_key: Embedding column holding cluster labels.
        palette: Matplotlib colormap name; ``None`` picks tab10 or tab20
            depending on the number of clusters.
        figsize: Figure size in inches.
        point_size: Marker size in the left panel; the right panel uses half.
        fontsize_range: Smallest and largest term font size, in points.
        weight_scale: Points of font size added per unit of ``-log10(p)``.
        candidates_per_term: Spiral positions tried per term before giving up.
        padding_fraction: Gap kept between labels, as a fraction of the plot span.
        title_left: Title of the left panel.
        title_right: Title of the right panel; defaults to a generated one.
        show_centroids: Whether the cluster centroid markers are visible. The
            markers are drawn either way — transparently when this is False —
            and always reserve their space, so the label positions are identical
            with and without them.
        centroid_marker: Matplotlib marker for the centroids; ``"o"`` (a filled
            circle) by default.
        centroid_size: Centroid marker area in points squared.
        centroid_color: Centroid marker colour.
        legend: Draw the cluster legend on the left panel.
        verbose: Print how many labels were placed.

    Returns:
        The matplotlib ``Figure``. Save it with ``fig.savefig(path)``.
    """
    plt = _require_matplotlib()

    xs, ys, cell_clusters = _as_columns(embedding, x_key, y_key, cluster_key)
    if not xs:
        raise ValueError("embedding contains no cells")

    colors = _cluster_palette(cell_clusters, palette)
    centroids = cluster_centroids(embedding, x_key, y_key, cluster_key)
    selected = select_terms(
        enrichment,
        top_n=top_n,
        pval_cutoff=pval_cutoff,
        fdr_cutoff=fdr_cutoff,
        clusters=clusters,
        max_label_chars=max_label_chars,
    )

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=figsize)

    # Left panel: the embedding as usual.
    _scatter_clusters(
        ax_left, xs, ys, cell_clusters, colors, point_size, 0.6, label=legend
    )
    ax_left.set_xlabel("UMAP 1", fontsize=12)
    ax_left.set_ylabel("UMAP 2", fontsize=12)
    ax_left.set_title(title_left, fontsize=14, fontweight="bold")
    if legend:
        ax_left.legend(
            bbox_to_anchor=(1.0, 1.0), loc="upper left", fontsize=8, frameon=False
        )

    # Right panel: the same embedding, faded, with the enriched terms on top.
    _scatter_clusters(
        ax_right, xs, ys, cell_clusters, colors, point_size / 2, 0.2, label=False
    )
    ax_right.set_xlabel("UMAP 1", fontsize=12)
    ax_right.set_ylabel("UMAP 2", fontsize=12)

    if title_right is None:
        title_right = "Enriched terms"
        if top_n is not None:
            title_right += f" (top {top_n} per cluster)"
    ax_right.set_title(title_right, fontsize=14, fontweight="bold")

    if not any(selected.values()):
        ax_right.text(
            0.5,
            0.5,
            "No enrichment results match the selected criteria",
            transform=ax_right.transAxes,
            fontsize=13,
            ha="center",
            va="center",
            color="gray",
        )
        fig.tight_layout()
        return fig

    # Freeze the axis limits before measuring text: placing labels must not
    # rescale the axes, or every measurement taken so far becomes wrong.
    ax_right.set_xlim(ax_right.get_xlim())
    ax_right.set_ylim(ax_right.get_ylim())

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    x_lo, x_hi = ax_right.get_xlim()
    y_lo, y_hi = ax_right.get_ylim()
    span = max(abs(x_hi - x_lo), abs(y_hi - y_lo))

    # The centroid markers occupy space in the layout whether or not they are
    # visible: the marker is always drawn and its box always reserved, and
    # show_centroids only sets the opacity. That way toggling it changes what you
    # see without moving a single label, so the two figures stay comparable.
    reserved: List[BBox] = []
    # Reserve a box a little larger than the marker itself, scaled with the
    # requested size so labels keep clear of it.
    marker_half = span * 0.010 * max(1.0, (centroid_size / 26.0) ** 0.5)
    for cluster, centroid in centroids.items():
        if not selected.get(cluster):
            continue
        ax_right.scatter(
            centroid["x"],
            centroid["y"],
            marker=centroid_marker,
            c=centroid_color,
            s=centroid_size,
            alpha=0.8 if show_centroids else 0.0,
            linewidths=0,
            zorder=5,
        )
        reserved.append(
            (
                centroid["x"] - marker_half,
                centroid["y"] - marker_half,
                centroid["x"] + marker_half,
                centroid["y"] + marker_half,
            )
        )

    label_boxes, skipped = _place_labels(
        ax_right,
        renderer,
        selected,
        centroids,
        colors,
        fontsize_range=fontsize_range,
        weight_scale=weight_scale,
        candidates_per_term=candidates_per_term,
        padding=span * padding_fraction,
        span=span,
        reserved=reserved,
        bounds=(min(x_lo, x_hi), min(y_lo, y_hi), max(x_lo, x_hi), max(y_lo, y_hi)),
    )

    # A label that had to go outside the axes would be cut off on save, so grow
    # the limits to cover everything placed. Widening the axes only shrinks each
    # label in data units, so it cannot introduce a new overlap.
    if label_boxes:
        margin = span * 0.01
        ax_right.set_xlim(
            min([x_lo] + [b[0] for b in label_boxes]) - margin,
            max([x_hi] + [b[2] for b in label_boxes]) + margin,
        )
        ax_right.set_ylim(
            min([y_lo] + [b[1] for b in label_boxes]) - margin,
            max([y_hi] + [b[3] for b in label_boxes]) + margin,
        )

    if verbose:
        print(f"  placed {len(label_boxes)} term labels, skipped {skipped} (no free space)")

    fig.tight_layout()
    return fig


def plot_umap_centroids(
    embedding: Any,
    x_key: str = "umap_1",
    y_key: str = "umap_2",
    cluster_key: str = "cluster",
    palette: Optional[str] = None,
    figsize: Tuple[float, float] = (9.0, 8.0),
    point_size: float = 8.0,
    title: str = "UMAP with cluster centroids",
    centroid_marker: str = "o",
    centroid_size: float = 80.0,
    centroid_color: str = "black",
    show_labels: bool = True,
):
    """
    Draw the embedding with each cluster's centroid marked and labelled.

    Useful as a reference figure when checking where the enrichment labels of
    :func:`plot_umap_enrichment` are anchored.

    Args:
        embedding: DataFrame or mapping with coordinates and cluster labels.
        x_key: Embedding column for the first dimension.
        y_key: Embedding column for the second dimension.
        cluster_key: Embedding column holding cluster labels.
        palette: Matplotlib colormap name; ``None`` picks it automatically.
        figsize: Figure size in inches.
        point_size: Marker size.
        title: Figure title.
        centroid_marker: Matplotlib marker for the centroids; ``"o"`` (a filled
            circle) by default.
        centroid_size: Centroid marker area in points squared.
        centroid_color: Centroid marker colour.
        show_labels: Write the cluster label beside each centroid.

    Returns:
        The matplotlib ``Figure``.
    """
    plt = _require_matplotlib()

    xs, ys, cell_clusters = _as_columns(embedding, x_key, y_key, cluster_key)
    colors = _cluster_palette(cell_clusters, palette)
    centroids = cluster_centroids(embedding, x_key, y_key, cluster_key)

    fig, ax = plt.subplots(figsize=figsize)
    _scatter_clusters(ax, xs, ys, cell_clusters, colors, point_size, 0.5, label=False)

    for cluster, centroid in centroids.items():
        ax.scatter(
            centroid["x"],
            centroid["y"],
            marker=centroid_marker,
            c=centroid_color,
            s=centroid_size,
            zorder=5,
            linewidths=0.8,
            edgecolors="white",
        )
        if show_labels:
            ax.text(
                centroid["x"],
                centroid["y"],
                f" {cluster}",
                fontsize=11,
                fontweight="bold",
                va="center",
                zorder=6,
            )

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig
