# Enrichment analysis and UMAP visualisation

`togoid.enrichment` turns TogoID ID conversion into a gene-set analysis toolkit.
A conversion route that ends in an annotation dataset *is* a gene-set library:
every term reached by the route becomes a set containing the input genes that
map to it. From there it is a standard over-representation analysis, and the
results can be drawn straight onto a single-cell embedding.

- [Why routes make gene sets](#why-routes-make-gene-sets)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Building gene sets](#building-gene-sets)
- [Running the enrichment](#running-the-enrichment)
- [Visualising on a UMAP](#visualising-on-a-umap)
- [Single-cell adapters](#single-cell-adapters)
- [Command line](#command-line)
- [API reference](#api-reference)
- [Notes on the statistics](#notes-on-the-statistics)

## Why routes make gene sets

Most enrichment tools ship a fixed set of annotation files, which go stale and
cover only the databases the author chose. TogoID resolves the annotations live,
so the *route* is the only thing that changes between databases:

```python
build_gene_sets(genes, route=["ncbigene", "uniprot", "reactome_pathway"])  # pathways
build_gene_sets(genes, route=["ncbigene", "uniprot", "go"])                # GO terms
build_gene_sets(genes, route=["ncbigene", "medgen", "mondo"])              # diseases
```

Anything TogoID can reach works the same way. Use `togoid route --src ... --dst ...`
or `converter.config_list_targets("ncbigene")` to find a route, and if you pass
one that does not exist the error names the alternatives:

```
RuntimeError: Conversion along route ncbigene -> hp_phenotype failed for all 1 batch(es).
First error: No direct connection between 'ncbigene' and 'hp_phenotype'.
Try one of these routes instead:
  - ncbigene -> medgen -> hp_phenotype
  - ncbigene -> nando -> hp_phenotype
```

## Installation

The analysis core needs nothing beyond what `togoid` already requires:

```bash
pip install togoid                  # gene sets, statistics, CSV/JSON output
pip install "togoid[enrichment]"    # + pandas, for the to_dataframe() views
pip install "togoid[plot]"          # + matplotlib, for the UMAP figures
pip install "togoid[singlecell]"    # + scanpy and leidenalg, for the full example
```

The hypergeometric test and the FDR correction are implemented in the standard
library, so an enrichment analysis runs on a bare install. pandas, matplotlib and
scanpy are each imported lazily and only when you use the feature that needs them.

## Quick start

```python
from togoid.enrichment import build_gene_sets, enrich_clusters

clusters = {
    "T cells":  ["CD3D", "CD3E", "CD3G", "IL7R", "LCK", "ZAP70", "CD2", "CD28", "LAT"],
    "B cells":  ["MS4A1", "CD79A", "CD79B", "CD19", "BLNK", "BANK1", "PAX5"],
    "Myeloid":  ["LYZ", "CD14", "FCGR3A", "CSF1R", "ITGAM", "TLR2", "TLR4", "S100A8"],
}
all_genes = sorted({gene for genes in clusters.values() for gene in genes})

library = build_gene_sets(all_genes, route=["ncbigene", "uniprot", "reactome_pathway"])
results = enrich_clusters(clusters, library, min_set_size=3)

print(results.summary())
results.to_csv("enrichment.csv")
```

```
Cluster B cells:
  terms tested: 6
  significant (FDR < 0.05): 2
    - Antigen activates B Cell Receptor (BCR) leading to generation of second messengers
      FDR=1.40e-02  genes=4/4  fold=3.71

Cluster T cells:
  terms tested: 11
  significant (FDR < 0.05): 5
    - Generation of second messenger molecules
      FDR=4.01e-03  genes=6/6  fold=2.89
    - Translocation of ZAP-70 to Immunological synapse
      FDR=1.05e-02  genes=5/5  fold=2.89
```

## Building gene sets

```python
library = build_gene_sets(
    genes,                      # gene symbols, or route[0] IDs with id_source=None
    route=["ncbigene", "uniprot", "reactome_pathway"],
    id_source="symbol",         # or None when genes are already route[0] IDs
    taxonomy="9606",            # used when resolving symbols
    label_dataset=None,         # defaults to route[-1]
    label_field="label",
    term_filters=None,          # e.g. {"go_aspect": ["biological_process"]}
    batch_size=100,
    convert_batch_size=50,
    pause=0.0,                  # seconds between API calls
    verbose=True,
)
```

What it does, in three API stages:

1. **Resolve** the gene symbols to `route[0]` identifiers with `LabelConverter`.
   Skipped when `id_source=None`.
2. **Convert** those identifiers along the route with `TogoIDConverter`, in
   batches, collecting `(source ID, term ID)` pairs. IDs are normalised through
   `togoid._ids.local_id()`, so the result is identical whether the API returns
   raw IDs or CURIEs such as `GO:0005634`.
3. **Annotate** the resulting terms with `AnnotationsConverter` to get their
   labels, plus any field named in `term_filters`.

Genes are stored under their *input* spelling, so the `genes` column of the
results shows the symbols you started with.

### The `GeneSetLibrary` object

| Member | Description |
|---|---|
| `sets` | `{term_id: set(genes)}` |
| `labels` | `{term_id: label}` |
| `genes` | Every gene in at least one set |
| `gene_to_terms` | Reverse mapping |
| `id_map` | `{input gene: route[0] ID}` |
| `unmapped` | Genes that did not resolve |
| `route`, `target_dataset` | Provenance |
| `filter_by_size(min, max)` | New library, filtered by set size |
| `to_rows()` / `to_dataframe()` | Tabular views |
| `save_json(path)` / `GeneSetLibrary.load_json(path)` | Cache to disk |

Caching matters: building a library for a few thousand genes is many API calls,
and the saved JSON reproduces the analysis exactly with no network access.

```python
library.save_json("reactome.json")
library = GeneSetLibrary.load_json("reactome.json")   # offline from here on
```

The format is shared with the R library's `togoid_save_gene_sets()`, so a
library built in either language can be read by the other and gives identical
p-values.

### Presets

```python
from togoid.enrichment import reactome_gene_sets, go_gene_sets, mondo_gene_sets

reactome_gene_sets(genes, taxonomy="9606")
go_gene_sets(genes, aspect="biological_process")   # or molecular_function,
                                                   # cellular_component, or None
mondo_gene_sets(genes)
```

| Preset | Route |
|---|---|
| `reactome` | `ncbigene → uniprot → reactome_pathway` |
| `go` | `ncbigene → uniprot → go` |
| `mondo` | `ncbigene → medgen → mondo` |

They are thin wrappers over `build_gene_sets` and accept the same keyword
arguments. They exist for convenience, not because these three databases are
special.

## Running the enrichment

```python
from togoid.enrichment import enrich, enrich_clusters

result  = enrich(query_genes, library, background=None,
                 min_set_size=5, max_set_size=500, min_overlap=1)

results = enrich_clusters(cluster_genes, library, background=None,
                          min_set_size=5, max_set_size=500, verbose=False)
```

Columns, sorted by FDR:

| Column | Meaning |
|---|---|
| `cluster` | Cluster label (per-cluster results only) |
| `term_id`, `term_label` | The annotation term |
| `overlap_count` | Query genes in the term |
| `term_size` | Term genes present in the background |
| `query_size` | Query genes present in the background |
| `background_size` | Universe size |
| `pvalue` | One-sided hypergeometric p-value |
| `fdr` | Benjamini-Hochberg adjusted p-value |
| `fold_enrichment` | Observed overlap / expected overlap |
| `genes` | The overlapping genes |

### Choosing a background

This is the decision that most affects the results.

- **Default** (`background=None`): every gene in the library, i.e. every gene in
  your experiment that TogoID could annotate. This asks *which terms distinguish
  this cluster from the rest of the experiment*, which is usually what you want
  for cell types.
- **Explicit**: pass a larger universe (say all expressed genes, or a
  genome-wide list) for the more conventional *which terms are over-represented
  relative to the genome* question. Expect many more significant hits.

`enrich_clusters` shares one background across all clusters, which is what makes
the FDR values comparable between them. Query genes outside the background are
dropped rather than counted, so the query size and the background stay consistent.

### Working with the results

```python
results.significant(0.05)          # filter by FDR
results.significant(0.01, "pvalue")# ... or by raw p-value
results.top(3)                     # best N terms per cluster
results["0"]                       # one cluster's EnrichmentResult
results.clusters                   # cluster labels, numerically ordered
results.rows                       # list of EnrichmentRow dataclasses
results.to_rows()                  # list of dicts
results.to_dataframe()             # pandas (needs [enrichment])
results.to_csv("out.csv")
results.summary(alpha=0.05)        # readable per-cluster report
```

## Tables

### Reading the result

A result prints as a formatted table in a console, and renders as an HTML table in
Jupyter and similar notebooks:

```python
print(results)     # or just `results` in a notebook cell
results.to_text(max_rows=None, max_label=60, max_genes=10)
```

The display folds `overlap_count` and `term_size` into one `k/M` column and drops
`query_size` and `background_size`, which are constant within a query. Nothing is
lost — `to_rows()`, `to_dataframe()` and the file exports all carry every column.

### Writing the result

```python
results.to_tsv("enrichment.tsv")     # one row per term, tab-separated
results.to_csv("enrichment.csv")     # ... comma-separated
results.to_csv("enrichment.txt", sep="|")
```

Tabs are the default for a reason: term labels routinely contain commas, which a
CSV has to quote and some spreadsheet imports then mis-parse.

### One row per cluster

```python
results.to_cluster_table(top_n=3, alpha=0.05)
results.write_cluster_table("by_cluster.tsv", top_n=3)
```

| Column | Meaning |
|---|---|
| `cluster` | Cluster label |
| `n_tested` | Terms tested for that cluster |
| `n_significant` | Terms below `alpha` |
| `top1_term_id`, `top1_term_label`, `top1_fdr` | The best term |
| `top2_…`, `top3_…` | Runners-up, up to `top_n` |

This is the shape you want when labelling clusters or reading a figure as a table.

### Matching a figure exactly

`plot_umap_enrichment` picks the terms it draws with `select_terms`. Call it
yourself with the same filters and the table cannot disagree with the picture:

```python
from togoid.enrichment import select_terms, selected_terms_table

selected = select_terms(results, top_n=3, fdr_cutoff=0.05, max_label_chars=None)
rows = selected_terms_table(selected)   # flat, ordered by cluster then FDR
```

`examples/scRNAseq_enrichment/04_visualize_umap.py` does exactly this, writing
`<figure-name>.tsv` and `<figure-name>_by_cluster.tsv` next to every figure.
Pass `--no-tables` to skip them.

## Visualising on a UMAP

```python
from togoid.enrichment import plot_umap_enrichment

fig = plot_umap_enrichment(
    embedding,              # umap_1, umap_2, cluster
    results,
    top_n=3,
    fdr_cutoff=0.05,
    pval_cutoff=None,
    clusters=None,          # e.g. ["0", "3", "7"]
    max_label_chars=40,
    palette=None,           # tab10 for <=10 clusters, tab20 above
    figsize=(20, 8),
    fontsize_range=(6, 14),
    title_left="UMAP clustering",
    title_right=None,
    show_centroids=True,     # mark each cluster centroid
    centroid_marker="o",     # any matplotlib marker; "o" is a filled circle
    centroid_size=26,
    centroid_color="black",
)
fig.savefig("umap_enrichment.pdf")
```

The left panel is the usual cluster UMAP. The right panel repeats it, faded, with
each cluster's enriched terms written around that cluster's centroid in the
cluster's colour, sized by `-log10(p)`.

`embedding` is a pandas DataFrame or *any* mapping with `umap_1`, `umap_2` and
`cluster` keys — plain lists are fine, so this step does not require pandas.
`results` may be a `ClusterEnrichmentResult`, a DataFrame, or a list of dicts.

How the layout works:

- Each label is measured once, then tested against candidate positions spiralling
  outwards from the centroid.
- A candidate is rejected if the label would overlap another label or a centroid
  marker.
- Positions that keep the label inside the axes are tried first; the axes are
  then widened to cover anything that had to go outside, so nothing is clipped.
- Labels that cannot be placed anywhere are dropped rather than drawn on top of
  each other. `verbose=True` reports how many.
- The centroid markers count as obstacles, and they do so whether or not they
  are visible: `show_centroids=False` draws them transparently rather than
  skipping them, so hiding the markers never moves a label. The reserved space
  scales with `centroid_size`.

If the figure is crowded, lower `top_n`, raise `figsize`, or shorten labels with
`max_label_chars`.

A companion figure shows where the labels are anchored:

```python
from togoid.enrichment import plot_umap_centroids

plot_umap_centroids(
    embedding,
    centroid_marker="o",   # filled circle, as in the enrichment panel
    centroid_size=80,
    show_labels=True,      # write the cluster label beside each centroid
).savefig("umap_centroids.pdf")
```

## Single-cell adapters

The core knows nothing about scanpy or Seurat. These helpers do the conversion,
and import their dependencies only when called.

### scanpy / AnnData

```python
from togoid.enrichment.adapters import (
    umap_dataframe_from_anndata,
    marker_genes_from_anndata,
)

import scanpy as sc
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

embedding = umap_dataframe_from_anndata(adata, basis="X_umap", cluster_key="leiden")
markers   = marker_genes_from_anndata(adata, top_n=100,
                                      pval_cutoff=0.05, logfc_min=0.25)
```

### Seurat, via CSV

```r
# In R
write.csv(cbind(Embeddings(obj, "umap"), cluster = as.character(Idents(obj))), "umap.csv")
write.csv(FindAllMarkers(obj), "markers.csv")
```

```python
from togoid.enrichment.adapters import umap_dataframe_from_csv, marker_genes_from_csv

embedding = umap_dataframe_from_csv("umap.csv")          # defaults match Seurat
markers   = marker_genes_from_csv("markers.csv")         # p_val_adj, avg_log2FC
```

`write_marker_gene_lists(markers, directory)` writes one plain-text gene list per
cluster, which is handy for sharing or for feeding other tools.

## Command line

```bash
# A single gene list
togoid enrich --genes "CD3D,CD3E,LCK,ZAP70,LAT" --preset reactome --format summary

# Per-cluster, from a CSV with cluster and gene columns
togoid enrich --clusters-file markers.csv --preset reactome \
              --fdr 0.05 --output enrichment.csv

# Any route, not just the presets
togoid enrich --genes-file genes.txt \
              --route ncbigene,uniprot,reactome_pathway --output out.csv

# GO, restricted to one aspect
togoid enrich --genes-file genes.txt --preset go \
              --go-aspect molecular_function --output go_mf.csv

# Build the library once, reuse it offline
togoid enrich --genes-file all_genes.txt --preset reactome --save-genesets sets.json
togoid enrich --clusters-file markers.csv --load-genesets sets.json --output out.csv

# Mouse
togoid enrich --genes-file genes.txt --preset reactome --taxonomy 10090
```

| Option | Description |
|---|---|
| `--genes`, `--genes-file`, `--clusters-file` | Input (mutually exclusive) |
| `--route`, `--preset` | Which annotation database to use |
| `--go-aspect` | GO aspect when `--preset go` |
| `--taxonomy` | Taxonomy for symbol resolution (default 9606) |
| `--id-source` | `symbol` (default) or `id` |
| `--cluster-column`, `--gene-column` | Column names in `--clusters-file` |
| `--min-set-size`, `--max-set-size` | Term size bounds |
| `--fdr`, `--top` | Narrow the exported table |
| `--save-genesets`, `--load-genesets` | Cache the library |
| `--output`, `--format` | `csv` (default), `json` or `summary` |

## API reference

### `togoid.enrichment`

| Name | Purpose |
|---|---|
| `build_gene_sets(genes, route, ...)` | Build a library from any route |
| `map_labels_to_ids(labels, dataset, taxonomy)` | Resolve labels to IDs |
| `GeneSetLibrary` | The gene-set container |
| `reactome_gene_sets`, `go_gene_sets`, `mondo_gene_sets` | Preset routes |
| `gene_sets_from_preset(name, genes)` | Dispatch by preset name |
| `ROUTES`, `GO_ASPECTS` | The preset routes and valid GO aspects |
| `enrich(genes, library, ...)` | Over-representation for one gene list |
| `enrich_clusters(mapping, library, ...)` | ... for several clusters |
| `EnrichmentResult`, `ClusterEnrichmentResult`, `EnrichmentRow` | Result types |
| `RESULT_COLUMNS` | The output column order |
| `plot_umap_enrichment(embedding, results, ...)` | The two-panel figure |
| `plot_umap_centroids(embedding, ...)` | Centroid reference figure |
| `cluster_centroids(embedding)` | Centroid coordinates |
| `select_terms(results, ...)` | The term selection used by the plot |
| `hypergeometric_sf`, `benjamini_hochberg`, `fold_enrichment` | Statistics |

### `togoid.enrichment.adapters`

`umap_dataframe_from_anndata`, `marker_genes_from_anndata`,
`umap_dataframe_from_csv`, `marker_genes_from_csv`, `write_marker_gene_lists`.

## Notes on the statistics

**The test.** Over-representation is the one-sided hypergeometric upper tail,
`P(X >= k)` for `X ~ Hypergeometric(N, M, n)`, where `k` is the observed overlap,
`N` the background size, `M` the term size and `n` the query size. This is the
same test as Fisher's exact test one-sided.

**The implementation.** `hypergeometric_sf` sums the probability mass function
from the far tail inwards, evaluating each term in log space with `math.lgamma`.
Log space avoids overflow at genome scale (`N ≈ 20000` would otherwise overflow
the factorials), and summing the smallest terms first avoids losing them to
rounding when the p-value is very small. The results agree with
`scipy.stats.hypergeom.sf` to better than 1e-9 relative error across the range
the package is tested on; `test_enrichment.py` checks this whenever SciPy is
installed.

**Multiple testing.** `benjamini_hochberg` applies the standard BH step-up
procedure with the cumulative-minimum monotonicity correction, matching
`statsmodels.stats.multitest.multipletests(method="fdr_bh")` to within floating
point. Correction is applied per cluster, over the terms actually tested for that
cluster.

**Term size filters.** `min_set_size` and `max_set_size` are applied *after*
intersecting each term with the background, so the sizes in the output are the
sizes the test used. Very small sets have no power; very large ones are usually
too general to be informative. Disease annotations (MONDO) cover far fewer genes
than pathways, so a lower `min_set_size` is often needed there.

**Annotation bias.** Only genes that TogoID could annotate enter the default
background. Well-studied genes are annotated more thoroughly than others, which
biases any over-representation analysis; this is a property of the annotation
databases, not of this implementation. The `unmapped` list on the library tells
you how many genes dropped out.
