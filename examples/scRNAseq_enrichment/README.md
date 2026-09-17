# scRNA-seq enrichment analysis with TogoID

A complete worked example: cluster a single-cell RNA-seq dataset, convert each
cluster's marker genes with TogoID, test them for enrichment, and draw the
result on the UMAP embedding.

The point of the example is step 3. One call takes gene symbols all the way to
an annotation database, and the *only* thing that changes between Reactome, GO
and MONDO is the route:

```python
build_gene_sets(genes, route=["ncbigene", "uniprot", "reactome_pathway"])
build_gene_sets(genes, route=["ncbigene", "uniprot", "go"])
build_gene_sets(genes, route=["ncbigene", "medgen", "mondo"])
```

## Data

The example uses the 10x Genomics public PBMC dataset. It is **not** included in
this repository; download it from 10x Genomics and note the licence terms on
their site.

```bash
mkdir -p data && cd data
# Pick any "Filtered feature-barcode matrix (MTX)" PBMC dataset from
# https://www.10xgenomics.com/datasets and unpack it here, so that you end up
# with a directory containing barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz
tar -xzf filtered_feature_bc_matrix.tar.gz
```

Any 10x-format directory works, as does any other dataset scanpy can read — only
step 1 touches the raw data.

## Requirements

```bash
pip install "togoid[singlecell]"
```

That pulls in pandas, matplotlib, scanpy and leidenalg. Steps 3 and 4 alone need
only `togoid[plot]`; the enrichment core itself needs nothing beyond `requests`.

## Running the pipeline

```bash
./run_pipeline.sh ../path/to/filtered_feature_bc_matrix results
```

or step by step:

```bash
python 01_clustering.py     --data-dir path/to/filtered_feature_bc_matrix --results-dir results
python 02_find_markers.py   --results-dir results
python 03_enrichment.py     --results-dir results
python 04_visualize_umap.py --results-dir results --top-n 3
```

Steps 1 and 2 need scanpy; steps 3 and 4 do not, and read only the CSV files the
earlier steps wrote.

## What each step does

### 1. `01_clustering.py` — clustering and UMAP

Loads the 10x matrix, applies standard QC (200–6000 genes per cell, <15%
mitochondrial reads), normalises, runs PCA, builds the neighbour graph, computes
the UMAP and clusters with Leiden.

**Outputs**

| File | Contents |
|---|---|
| `01_clustered.h5ad` | The AnnData object, for steps 2 and any further analysis |
| `01_umap.csv` | `cell, umap_1, umap_2, cluster` — all step 4 needs |

The clustering is computed once here and reused unchanged by every later step,
so the UMAP shown beside the enrichment results is exactly the one the gene
lists came from.

### 2. `02_find_markers.py` — marker genes

Wilcoxon rank-sum test per cluster, filtered by adjusted p-value and log fold
change.

**Outputs**

| File | Contents |
|---|---|
| `02_markers.csv` | Long format `cluster, gene` — the input to step 3 |
| `02_marker_gene_lists/cluster_*_markers.txt` | One plain gene list per cluster |
| `02_with_markers.h5ad` | AnnData with the DE results attached |

Useful options: `--top-n` (genes kept per cluster, default 100), `--pval-cutoff`,
`--logfc-min`.

### 3. `03_enrichment.py` — TogoID conversion and enrichment

For each target database, `build_gene_sets()` resolves the gene symbols to NCBI
Gene IDs, walks the route, fetches the term labels, and returns a gene-set
library. `enrich_clusters()` then runs a hypergeometric test per cluster with
BH-FDR correction, against a background of every marker gene in the experiment.

| Target | Route | Term sizes |
|---|---|---|
| `reactome` | `ncbigene → uniprot → reactome_pathway` | 5–500 |
| `go` | `ncbigene → uniprot → go` (biological process) | 5–500 |
| `mondo` | `ncbigene → medgen → mondo` | 3–500 |

**Outputs** (per target)

| File | Contents |
|---|---|
| `03_genesets_<target>.json` | The gene-set library, cached |
| `03_enrichment_<target>_all.csv` | Every tested term |
| `03_enrichment_<target>_significant.csv` | FDR < 0.05 only |
| `03_enrichment_<target>_summary.txt` | Readable per-cluster report |

Useful options: `--targets reactome,go`, `--go-aspect molecular_function`,
`--taxonomy 10090` (mouse), and `--reuse-genesets` to skip the API calls on a
re-run.

To add another database, add an entry to `TARGETS` at the top of the script with
its route — nothing else changes.

### 4. `04_visualize_umap.py` — UMAP figures

Draws a two-panel figure per target: cluster UMAP on the left, and the same
embedding on the right with each cluster's enriched terms written around its
centroid. Font size scales with `-log10(p)`, and labels that cannot be placed
without overlapping are dropped rather than drawn illegibly.

**Outputs**

| File | Contents |
|---|---|
| `04_umap_centroids.pdf/.png` | Reference figure showing where labels are anchored |
| `04_umap_enrichment_<target>_top<N>.pdf/.png` | The two-panel enrichment figure |
| `04_umap_enrichment_<target>_top<N>.tsv` | The terms drawn on that figure, one row per term |
| `04_umap_enrichment_<target>_top<N>_by_cluster.tsv` | The same terms, one row per cluster |

The TSV tables are written with the same filters the figure used, so the two can
never disagree. Tabs rather than commas, because term labels contain commas.

Useful options: `--top-n 5`, `--fdr-cutoff 0.01`, `--clusters 0,3,7`,
`--max-label-chars 30`, `--formats pdf` (skip the PNG), `--no-tables` (skip the
TSVs), `--no-centroids` (hide the centroid markers).

## Reading the results

A few things worth knowing when you look at the output:

- **Background choice matters.** The default background is every marker gene in
  the experiment, not the whole genome. This asks "which terms distinguish this
  cluster from the other clusters?", which is usually the question you want for
  cell types. Pass `background=` to `enrich_clusters()` for a genome-wide
  universe instead.
- **Ribosomal genes dominate some clusters.** In PBMC data, highly expressed
  ribosomal protein genes often appear in several clusters' marker lists, so
  translation pathways come out strongly. Either raise `--logfc-min` in step 2 or
  filter `RPL*`/`RPS*` out of the gene lists if you want cell-type-specific
  biology only.
- **MONDO sets are small.** Disease annotations cover far fewer genes than
  pathways, so step 3 uses a lower minimum term size for MONDO, and fewer terms
  reach significance.

## Doing this without the scripts

The whole analysis is four calls:

```python
from togoid.enrichment import build_gene_sets, enrich_clusters, plot_umap_enrichment
from togoid.enrichment.adapters import umap_dataframe_from_anndata, marker_genes_from_anndata

markers = marker_genes_from_anndata(adata)
library = build_gene_sets(
    sorted({g for genes in markers.values() for g in genes}),
    route=["ncbigene", "uniprot", "reactome_pathway"],
)
results = enrich_clusters(markers, library)

fig = plot_umap_enrichment(umap_dataframe_from_anndata(adata), results, top_n=3)
fig.savefig("umap_enrichment.pdf")
```

Or from the command line, without writing any Python:

```bash
togoid enrich --clusters-file results/02_markers.csv \
              --preset reactome --fdr 0.05 \
              --output results/enrichment.csv
```

## Licence note

The 10x Genomics PBMC data is distributed by 10x Genomics under their own terms.
Check the licence on the dataset page before redistributing it or the results
derived from it. Nothing from that dataset is included in this repository.
