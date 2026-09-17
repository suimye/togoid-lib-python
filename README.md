# TogoID Python Library

Python library and CLI tool for biological database ID conversion and annotation using [TogoID](https://togoid.dbcls.jp/).

## Features

- **ID Conversion**: Convert IDs between biological databases
- **ID Conversion with Annotations**: Add annotation columns during conversion
- **ID Conversion with Filtering**: Filter conversion results by annotation values
- **Ortholog Retrieval**: Get orthologs through round-trip conversion and taxonomy filtering
- **Label to ID**: Convert biological labels (gene names, etc.) to database IDs with dataset-based API selection
- **Annotations**: Get labels and annotations for database IDs
- **Multiple Formats**: Support for JSON, CSV, TSV, dict, table, and pandas DataFrame
- **Dual Interface**: Use as Python library or command-line tool
- **Comprehensive**: Search databases, find routes, get configurations
- **Enrichment Analysis**: Turn any conversion route into gene sets, test them for over-representation, and draw the result on a single-cell UMAP

## Installation

### Using uv (recommended - faster)

[uv](https://github.com/astral-sh/uv) is a blazingly fast Python package installer and resolver (10-100x faster than pip).

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup
git clone https://github.com/togoid/togoid-lib-python.git
cd togoid-lib-python

# Create virtual environment and install
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install the package
uv pip install -e .

# With pandas support
uv pip install -e ".[pandas]"

# With development tools
uv pip install -e ".[dev]"
```

**💡 Tip:** See [QUICKSTART_UV.md](QUICKSTART_UV.md) for a detailed uv quick start guide.

### Using pip (traditional)

```bash
# From source
git clone https://github.com/togoid/togoid-lib-python.git
cd togoid-lib-python
pip install -e .

# With pandas support
pip install -e ".[pandas]"
```

## Quick Start

### As a Python Library

```python
from togoid import TogoIDConverter, AnnotationsConverter, LabelConverter

# ID Conversion
converter = TogoIDConverter()

# JSON format (default)
result = converter.convert(ids=["1", "9"], route=["ncbigene", "ensembl_gene"])

# Dict format
result_dict = converter.convert(ids=["1", "9"], route=["ncbigene", "ensembl_gene"], format="dict")
# Output: {'ids': ['1', '9'], 'route': ['ncbigene', 'ensembl_gene'], 'results': {'1': ['ENSG00000121410'], '9': ['ENSG00000171428']}}

# Table format
result_table = converter.convert(ids=["1", "9"], route=["ncbigene", "ensembl_gene"], format="table")
# Output: [["1", "ENSG00000121410"], ["9", "ENSG00000075624"]]

# DataFrame format (requires pandas)
result_df = converter.convert(ids=["1", "9"], route=["ncbigene", "ensembl_gene"], format="dataframe")

# ID Conversion with Annotations
result_with_annotations = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene", "ensembl_transcript"],
    format="table",
    annotate=[("ncbigene", "label")]  # Add ncbigene label as annotation column
)

# ID Conversion with Filtering
result_filtered = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene", "ensembl_transcript"],
    format="table",
    annotate=[("ncbigene", "label")],
    filter=[("ensembl_transcript", "transcript_flag", ["MANE Select"])]  # Only MANE Select transcripts
)

# Get Orthologs
orthologs = converter.get_ortholog(
    ids=["1", "9"],
    route=["ncbigene", "homologene"],
    target_taxids=["10090", "10116"]  # Mouse and Rat
)

# Label to ID Conversion
label_converter = LabelConverter()

# Convert labels with dataset specification
results = label_converter.convert(
    labels=["BRCA1", "TP53"],
    dataset="ncbigene",
    taxonomy="9606"  # Human
)

# Convert labels for other datasets
results = label_converter.convert(
    labels=["caffeine"],
    dataset="chebi",
    label_types=["togoid_chebi_label"]
)

# Get Annotations
annotator = AnnotationsConverter()
annotations = annotator.execute_query(
    dataset_name="ncbigene",
    ids=["672", "7157"],
    fields=["label", "gene_synonym"],
    filters={}
)
```

### As a Command-Line Tool

```bash
# Basic ID Conversion
togoid convert --ids 1,9 --route ncbigene,ensembl_gene

# Convert with different output formats
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format dict
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format table

# ID Conversion with Annotations
togoid convert --ids 1,9 --route ncbigene,ensembl_gene,ensembl_transcript \
  --format table \
  --annotate ncbigene label \
  --annotate ncbigene full_name

# ID Conversion with Filtering
togoid convert --ids 1,9 --route ncbigene,ensembl_gene,ensembl_transcript \
  --format table \
  --annotate ncbigene label \
  --filter ensembl_transcript transcript_flag "MANE Select"

# Get Orthologs
togoid get-ortholog --ids 672,7157 \
  --route ncbigene,homologene \
  --target-taxids 10090,10116 \
  --format table

# Label to ID Conversion
togoid label2id --labels "BRCA1,TP53,EGFR" --dataset ncbigene --taxonomy 9606

# Get Annotations
togoid annotate --dataset ncbigene --ids 672,7157 \
  --field gene_synonym \
  --field full_name

# List available annotation fields
togoid annotate --dataset ncbigene --list-fields

# Configuration
togoid config dataset ncbigene
togoid config descriptions
togoid count ncbigene ensembl_gene --ids 1,9
```

## ID Prefixes (CURIE format)

The TogoID API can format converted IDs with their dataset prefix (CURIE), e.g.
`0005634` → `GO:0005634`, `217124` → `ORPHA:217124`. This is the **default** output
of `convert()`:

```python
converter.convert(route=["orphanet_gene", "uniprot", "go"], ids=["217124"], report="full")
# -> [['ORPHA:217124', 'P23284', 'GO:0005737'], ...]
```

Pass `prefix=False` (CLI: `--raw`) to get raw IDs without prefixes:

```python
converter.convert(route=["orphanet_gene", "uniprot", "go"], ids=["217124"], report="full", prefix=False)
# -> [['217124', 'P23284', '0005737'], ...]
```

```bash
togoid convert --route orphanet_gene,uniprot,go --ids 217124 --report full --raw
```

Notes:

- Prefixing requires the server-side `?prefix` support ([togoid-api PR #149](https://github.com/togoid/togoid-api/pull/149)).
  Against API deployments that predate it, output is raw regardless of `prefix`.
- `annotate` / `filter` and `get_ortholog` transparently handle prefixed IDs: the
  library matches on the raw local ID internally, so joins keep working while the
  displayed IDs stay prefixed.
- `label2id` (`LabelConverter`) returns raw identifiers — it resolves labels via
  PubDictionaries / SPARQList, which are not covered by the API's `?prefix`.

## Breaking Changes

### Version 0.2.0+

**1. label_types parameter now requires list format**

The `label_types` parameter in `LabelConverter.convert()` has been changed from string to list type.

```python
# ❌ Old (will not work)
label_converter.convert(
    labels=["BRCA1"],
    dataset="ncbigene",
    label_types="symbol,synonym"  # String format
)

# ✅ New (correct)
label_converter.convert(
    labels=["BRCA1"],
    dataset="ncbigene",
    label_types=["symbol", "synonym"]  # List format
)
```

**2. format="dict" deprecated for routes with 3+ datasets**

When using routes with 3 or more datasets, `format="dict"` is no longer supported. Use `format="table"` or `format="dataframe"` instead.

```python
# ❌ Old (will raise error)
converter.convert(
    ids=["1"],
    route=["ncbigene", "ensembl_gene", "ensembl_transcript"],
    format="dict"
)

# ✅ New (correct)
converter.convert(
    ids=["1"],
    route=["ncbigene", "ensembl_gene", "ensembl_transcript"],
    format="table"  # or "dataframe"
)
```

**3. annotator.execute_query filters parameter is now optional**

The `filters` parameter in `AnnotationsConverter.execute_query()` is now optional and defaults to an empty dictionary.

```python
# Both work now
annotations = annotator.execute_query(
    dataset_name="ncbigene",
    ids=["672"],
    fields=["label"],
    filters={}  # Can be omitted
)

annotations = annotator.execute_query(
    dataset_name="ncbigene",
    ids=["672"],
    fields=["label"]  # No filters parameter needed
)
```

## Usage Examples

### ID Conversion

#### Different Output Formats

```python
from togoid import TogoIDConverter

converter = TogoIDConverter()

# JSON (default) - raw API response
json_result = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene"]
)

# Dict - Includes ids, route, and results mapping {source_id: [target_ids]}
dict_result = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene"],
    format="dict"
)

# Table - [[source_id, target_id], ...] 2D array
table_result = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene"],
    format="table"
)

# DataFrame - pandas DataFrame with source_id and target_id columns
df_result = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene"],
    format="dataframe"
)
```

#### ID Conversion with Annotations

```python
# Add annotation columns to conversion results
result = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene", "ensembl_transcript"],
    format="table",
    annotate=[
        ("ncbigene", "label"),           # Add gene label from ncbigene
        ("ncbigene", "full_name"),       # Add full gene name from ncbigene
        ("ensembl_gene", "label")        # Add gene label from ensembl_gene
    ]
)
# Result includes original conversion + 3 annotation columns
```

#### ID Conversion with Filtering

```python
# Filter conversion results by annotation values
result = converter.convert(
    ids=["1", "9"],
    route=["ncbigene", "ensembl_gene", "ensembl_transcript"],
    format="table",
    annotate=[("ncbigene", "label")],
    filter=[
        ("ensembl_transcript", "transcript_flag", ["MANE Select"])
    ]
)
# Only returns transcripts with "MANE Select" flag
# 15 transcripts → 2 transcripts (filtered)
```

#### Get Orthologs

```python
# Get orthologs through round-trip conversion and taxonomy filtering
# Process: ncbigene -> homologene -> ncbigene -> taxonomy -> filter by taxid
result = converter.get_ortholog(
    ids=["1", "9"],                      # Human genes
    route=["ncbigene", "homologene"],    # Via homologene
    target_taxids=["10090", "10116"]     # Mouse and Rat
)
# Returns: [
#   ['1', '11167', '117586', '10090'],   # source_id, homologene_id, mouse_gene_id, taxid
#   ['1', '11167', '140656', '10116'],   # same source via same homologene group
#   ['9', '37329', '116632', '10116'],
#   ['9', '37329', '17961', '10090']
# ]
# Rows are ordered as: [source_id, homologene_id, target_gene_id, taxonomy_id]
```

#### Search and Route

```python
# Search databases by name
databases = converter.search_databases("uniprot")

# Find routes between databases
routes = converter.route(src="ncbigene", dst="uniprot", max_hops=3)

# Lookup which tables contain an ID
tables = converter.lookup_id("672")
```

### Label to ID Conversion

```python
from togoid import LabelConverter

converter = LabelConverter(verbose=True)

# Convert gene symbols (uses SPARQList API based on dataset config)
results = converter.convert(
    labels=["BRCA1", "TP53", "EGFR"],
    dataset="ncbigene",
    taxonomy="9606"  # Human
)
# Returns: [{"input": "BRCA1", "match_type": "symbol", "symbol": "BRCA1", "identifier": "672"}, ...]

# Convert chemical names (uses PubDictionaries API based on dataset config)
results = converter.convert(
    labels=["caffeine"],
    dataset="chebi",
    label_types=["togoid_chebi_label"]  # Optional: override dataset config (list format)
)

# Convert disease names
results = converter.convert(
    labels=["breast cancer"],
    dataset="mondo",
    threshold=0.5  # PubDictionaries matching threshold
)

# Label types are auto-configured from dataset, or can be manually specified (as list)
results = converter.convert(
    labels=["BRCA1"],
    dataset="ncbigene",
    label_types=["symbol"],  # Override: only search by symbol (list format)
    taxonomy="9606"
)
```

### Annotations

```python
from togoid import AnnotationsConverter

annotator = AnnotationsConverter()

# List available fields for a dataset
fields = annotator.list_fields("ncbigene")
for field_name, field_meta in fields:
    print(f"{field_name}: {field_meta['label']}")

# Get annotations for IDs
result = annotator.execute_query(
    dataset_name="ncbigene",
    ids=["672", "7157"],
    fields=["label", "gene_synonym", "type_of_gene"],
    filters={"type_of_gene": ["protein-coding"]}
)

for id, annotations in result.items():
    print(f"{id}: {annotations}")
```

## Command-Line Interface

### Convert Command

```bash
# Basic conversion
togoid convert --ids 1,9 --route ncbigene,ensembl_gene

# With output format
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format dict

# Save to file
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format csv --output results.csv

# With additional parameters
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --report pair --limit 100
```

### Label2ID Command

```bash
# Basic conversion
togoid label2id --dataset ncbigene --labels "BRCA1,TP53,EGFR" --taxonomy 9606
togoid label2id --dataset chebi --labels 'caffeine' --label_types 'togoid_chebi_label'

# From file
echo -e "BRCA1\nTP53\nEGFR" > genes.txt
togoid label2id --dataset ncbigene --label-file genes.txt --taxonomy 9606

# CSV output
togoid label2id --dataset ncbigene --labels "BRCA1,TP53" --taxonomy 9606 --format csv --output results.csv

# With PubDictionaries (for non-gene labels)
togoid label2id --dataset chebi --labels "breast cancer" --label_types "togoid_mondo_label"

# Verbose mode
togoid label2id --dataset ncbigene --labels "BRCA1,TP53" --taxonomy 9606 --verbose
```

### Annotate Command

```bash
# Get annotations
togoid annotate --dataset ncbigene --ids 672,7157 --field gene_synonym --field full_name

# List available fields
togoid annotate --dataset ncbigene --list-fields

# With filters
togoid annotate --dataset ncbigene --ids 672,7157 \
    --field type_of_gene --field gene_synonym \
    --filter type_of_gene=protein-coding

# CSV output
togoid annotate --dataset ncbigene --ids 672,7157 \
    --field gene_synonym --format csv --output genes.csv

# From file
togoid annotate --dataset ncbigene --ids-file gene_ids.txt --field gene_synonym
```

### Other Commands

```bash
# Search databases
togoid search databases uniprot
togoid search id NM_001110

# Lookup ID
togoid lookup id 672

# Find routes
togoid route ncbigene uniprot --max-hops 3

# Count mappings
togoid count ncbigene ensembl_gene --ids 1,9

# Get configuration
togoid config dataset ncbigene
togoid config relation ncbigene-ensembl_gene
togoid config descriptions
togoid config statistics
togoid config taxonomy
```

## API Documentation

### TogoIDConverter

Main class for ID conversion operations.

**Methods:**
- `convert(route, ids, format='json', **kwargs)` - Convert IDs between databases
- `count(src, dst, ids, link=None)` - Count mappings
- `search_databases(name)` - Search databases by name
- `search_id(id_string)` - Search databases by ID pattern
- `lookup_id(id_string)` - Lookup which tables contain an ID
- `route(src, dst, max_hops=3)` - Find routes between databases
- `config_dataset(name=None)` - Get dataset configuration
- `config_relation(src=None, dst=None)` - Get relation configuration
- `config_descriptions()` - Get database descriptions
- `config_statistics()` - Get database statistics
- `config_taxonomy()` - Get taxonomy list

### LabelConverter

Main class for converting biological labels to database IDs. The upstream API is selected from the dataset configuration.

**Methods:**
- `convert(labels, dataset, label_types=None, tags=None, threshold=0.5, preferred_dictionary=None, taxonomy=None, format='json')` - Convert labels to IDs (auto-selects API based on dataset config)
- `convert_pubdictionaries(labels, dictionaries, tags=None, threshold=0.5, preferred_dictionary=None)` - Convert using PubDictionaries API
- `convert_sparqlist(labels, sparqlist, label_types, taxonomy=None)` - Convert using SPARQList API

**API selection:**
- The dataset's `label_resolver` configuration is fetched from the TogoID API
  (`/config/dataset`) and decides which upstream to use — it is not inferred from
  the labels themselves.
- `label_resolver.sparqlist` present → SPARQList (e.g. `ncbigene`)
- otherwise → PubDictionaries (e.g. `chebi`, `mondo`)

### AnnotationsConverter

Main class for getting annotations and labels for IDs.

**Methods:**
- `list_fields(dataset_name)` - List available annotation fields
- `execute_query(dataset_name, ids, fields, filters)` - Execute GraphQL query to get annotations
- `build_rows(dataset_label, fields, field_meta, records, filters, compact)` - Build table rows from query results

## CLI Command Reference

### Basic Commands

```bash
# Convert IDs between databases
togoid convert --ids 1,9 --route ncbigene,ensembl_gene
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format dict

# Label to ID conversion
togoid label2id --labels "BRCA1,TP53" --dataset ncbigene --taxonomy 9606

# Get annotations
togoid annotate --dataset ncbigene --ids 672,7157 --field label --field gene_synonym
togoid annotate --dataset ncbigene --list-fields

# Utilities
togoid count ncbigene ensembl_gene --ids 1,9

# Configuration
togoid config dataset ncbigene
togoid config descriptions
```

### Advanced Features

#### ID Conversion with Annotations

Add annotation columns to your conversion results:

```bash
# Add single annotation
togoid convert --ids 1,9 \
  --route ncbigene,ensembl_gene,ensembl_transcript \
  --format table \
  --annotate ncbigene label

# Add multiple annotations
togoid convert --ids 1,9 \
  --route ncbigene,ensembl_gene,ensembl_transcript \
  --format table \
  --annotate ncbigene label \
  --annotate ncbigene full_name \
  --annotate ensembl_transcript transcript_flag
```

#### ID Conversion with Filtering

Filter results by annotation values:

```bash
# Filter by single value
togoid convert --ids 1,9 \
  --route ncbigene,ensembl_gene,ensembl_transcript \
  --format table \
  --filter ensembl_transcript transcript_flag "MANE Select"

# Combine annotations and filtering
togoid convert --ids 1,9 \
  --route ncbigene,ensembl_gene,ensembl_transcript \
  --format table \
  --annotate ncbigene label \
  --filter ensembl_transcript transcript_flag "MANE Select"
```

#### Ortholog Retrieval

Get orthologs using round-trip conversion:

```bash
# Get mouse and rat orthologs for human genes
togoid get-ortholog \
  --ids 672,7157 \
  --route ncbigene,homologene \
  --target-taxids 10090,10116 \
  --format table

# Output as JSON
togoid get-ortholog \
  --ids 672,7157 \
  --route ncbigene,homologene \
  --target-taxids 10090 \
  --format json
```

Table output columns are `[source_id, homologene_id, target_id, taxonomy_id]`.

### Input/Output Options

```bash
# Read IDs from file
echo "1\n9\n672" > ids.txt
togoid convert --ids-file ids.txt --route ncbigene,ensembl_gene

# Save output to file
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --output result.json

# Different output formats
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format json
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format dict
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format table
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --format csv
```

### Report Options

Control what information is returned:

```bash
# Only target IDs (default)
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --report target

# Source-target pairs
togoid convert --ids 1,9 --route ncbigene,ensembl_gene --report pair

# Full path including intermediate IDs
togoid convert --ids 1,9 --route ncbigene,ensembl_gene,ensembl_transcript --report full
```

**Note:** When using routes with 3+ datasets or annotations, the library automatically uses `report=full` to include all intermediate IDs.

### Finding Reachable Datasets

Get a list of datasets that are reachable from a source dataset in one hop:

```bash
# CLI
togoid config list-targets ncbigene

# Python
converter = TogoIDConverter()
targets = converter.config_list_targets("ncbigene")
print(targets)  # ['ensembl_gene', 'hgnc', 'mgi', ...]
```

### Route Suggestions

When datasets are not directly connected, the library automatically suggests alternative routes:

```python
# If ncbigene → chembl_compound isn't directly connected
converter.convert(ids=["1"], route=["ncbigene", "chembl_compound"])

# Error message will suggest alternatives:
# RuntimeError: No direct connection between 'ncbigene' and 'chembl_compound'.
#
# Suggested routes (2 hops):
# - ncbigene → ensembl_gene → chembl_compound
# - ncbigene → uniprot → chembl_compound
#
# Suggested routes (3 hops):
# - ncbigene → ensembl_gene → pdb → chembl_compound
```

## Use Case: Enrichment Analysis and UMAP Visualization

`togoid.enrichment` turns ID conversion into gene-set analysis. A TogoID route
that ends in an annotation dataset *is* a gene-set library: every term reached by
the route becomes a set containing the input genes that map to it. From there it
is a standard over-representation analysis, and the results can be drawn directly
onto a single-cell embedding.

The key property is that **only the route changes** between annotation databases:

```python
from togoid.enrichment import build_gene_sets

build_gene_sets(genes, route=["ncbigene", "uniprot", "reactome_pathway"])  # pathways
build_gene_sets(genes, route=["ncbigene", "uniprot", "go"])                # GO terms
build_gene_sets(genes, route=["ncbigene", "medgen", "mondo"])              # diseases
```

Anything TogoID can reach works the same way. The enrichment core needs **no
dependencies beyond those `togoid` already requires** — the hypergeometric test
and the BH-FDR correction are implemented with `math.lgamma`. pandas, matplotlib
and scanpy are optional extras, imported lazily.

See [docs/enrichment.md](docs/enrichment.md) for the full manual and
[examples/scRNAseq_enrichment/](examples/scRNAseq_enrichment/) for a runnable
pipeline.

### Installation

```bash
pip install togoid                  # gene sets, statistics, CSV/JSON output
pip install "togoid[enrichment]"    # + pandas, for the to_dataframe() views
pip install "togoid[plot]"          # + matplotlib, for the UMAP figures
pip install "togoid[singlecell]"    # + scanpy and leidenalg, for the full example
```

### Quick Start

```python
from togoid.enrichment import build_gene_sets, enrich_clusters

clusters = {
    "T cells": ["CD3D", "CD3E", "CD3G", "IL7R", "LCK", "ZAP70", "CD2", "CD28", "LAT"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19", "BLNK", "BANK1", "PAX5"],
    "Myeloid": ["LYZ", "CD14", "FCGR3A", "CSF1R", "ITGAM", "TLR2", "TLR4", "S100A8"],
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

### Three Annotation Databases

```python
from togoid.enrichment import reactome_gene_sets, go_gene_sets, mondo_gene_sets

# Reactome pathways: ncbigene -> uniprot -> reactome_pathway
pathways = reactome_gene_sets(all_genes, taxonomy="9606")

# GO terms: ncbigene -> uniprot -> go, filtered by aspect
processes = go_gene_sets(all_genes, aspect="biological_process")
functions = go_gene_sets(all_genes, aspect="molecular_function")

# MONDO diseases: ncbigene -> medgen -> mondo
diseases = mondo_gene_sets(all_genes)
```

These presets are thin wrappers over `build_gene_sets`. If you pass a route that
does not exist, the error names the working alternatives:

```
RuntimeError: Conversion along route ncbigene -> hp_phenotype failed for all 1 batch(es).
First error: No direct connection between 'ncbigene' and 'hp_phenotype'.
Try one of these routes instead:
  - ncbigene -> medgen -> hp_phenotype
  - ncbigene -> nando -> hp_phenotype
```

### Caching Gene Sets

Building a library for a few thousand genes is many API calls. Save it once and
the analysis reproduces exactly, offline:

```python
from togoid.enrichment import GeneSetLibrary

library.save_json("reactome.json")
library = GeneSetLibrary.load_json("reactome.json")   # no network access
```

The format is shared with the R library's `togoid_save_gene_sets()`, so a
library built in either language can be read by the other.

### UMAP Visualization

```python
from togoid.enrichment import plot_umap_enrichment

fig = plot_umap_enrichment(
    embedding,          # DataFrame or dict with umap_1, umap_2, cluster
    results,
    top_n=3,
    fdr_cutoff=0.05,
    show_centroids=True,     # mark each cluster centroid
    centroid_marker="o",     # a black filled circle
    centroid_size=26,
)
fig.savefig("umap_enrichment.pdf")
```

`show_centroids=False` hides the markers by drawing them transparently. They are
still drawn and still reserve their space, so **the labels stay exactly where
they were** — the two figures differ only in whether you can see the dots. The
marker itself is any matplotlib marker code — `"o"` (filled circle, the
default), `"x"`, `"s"`, and so on — with `centroid_size` and `centroid_color` to
match.

The left panel is the usual cluster UMAP; the right repeats it with each
cluster's enriched terms written around its centroid, sized by `-log10(p)`.
Labels that cannot be placed without overlapping are dropped rather than drawn on
top of each other, and the axes are widened so nothing is clipped.

`embedding` is any mapping with `umap_1`, `umap_2` and `cluster` keys — plain
lists work, so this step does not require pandas.

#### Reactome pathways

10x Genomics public PBMC data, 3,733 cells after QC, 13 Leiden clusters.

![UMAP with enriched Reactome pathways](https://raw.githubusercontent.com/suimye/togoid-lib-python/docs-figures/umap_enrichment_reactome.png)

Platelet degranulation marks the platelet cluster, MHC class II antigen
presentation the monocyte/DC clusters, and B cell receptor signalling the B cell
cluster.

#### GO biological process

![UMAP with enriched GO biological processes](https://raw.githubusercontent.com/suimye/togoid-lib-python/docs-figures/umap_enrichment_go.png)

An independent confirmation of the Reactome result via a different route:
platelet aggregation on the same platelet cluster, T cell activation on the T
cell clusters.

#### MONDO diseases

![UMAP with enriched MONDO diseases](https://raw.githubusercontent.com/suimye/togoid-lib-python/docs-figures/umap_enrichment_mondo.png)

Disease annotation is far sparser than pathway annotation, so only a couple of
terms pass the size and significance filters. Shown as-is, because the contrast
with the two figures above illustrates how coverage differs between TogoID
targets.

Full-resolution figures and the analysis notes are collected in
[issue #1](https://github.com/suimye/togoid-lib-python/issues/1).

### Tables and Notebook Display

Results print as a formatted table in a console, and render as an HTML table in
Jupyter:

```python
print(results)          # or just `results` in a notebook
```

```
cluster  term_id        term_label                            overlap    pvalue       fdr  fold_enrichment  genes
-------  -------------  ------------------------------------  -------  --------  --------  ---------------  ---------------------------------
0        R-HSA-6798695  Neutrophil degranulation                32/60  8.15e-20  6.85e-18             5.39  ANPEP, ASAH1, CD14, CD36 (+28)
0        R-HSA-166058   MyD88:MAL(TIRAP) cascade initiated o…      6/7  5.02e-06  2.11e-04             8.65  CD14, CD36, IRAK3, S100A8, TLR2
1        R-HSA-156902   Peptide chain elongation                32/70  3.75e-17  2.47e-15             4.62  EEF1A1, RPL10, RPL11 (+29)
```

The display folds `overlap_count` and `term_size` into one `k/M` column and drops
`query_size` and `background_size`, which are constant. Every column is still
there in the data — only the view is shortened.

Writing tables out:

```python
results.to_tsv("enrichment.tsv")            # one row per term
results.to_csv("enrichment.csv")            # same, comma-separated
results.write_cluster_table("by_cluster.tsv", top_n=3)   # one row per cluster
results.to_cluster_table(top_n=3)           # ... as a list of dicts
results.to_dataframe()                      # ... as a pandas DataFrame
```

Tabs are the default for a reason: term labels routinely contain commas, which a
CSV has to quote and some spreadsheet imports then mis-parse.

To export exactly what a figure shows, pass the same filters to `select_terms()`:

```python
from togoid.enrichment import select_terms, selected_terms_table

selected = select_terms(results, top_n=3, fdr_cutoff=0.05)
rows = selected_terms_table(selected)   # flat, ordered by cluster then FDR
```

Step 4 of the example pipeline does this automatically, writing
`<figure-name>.tsv` and `<figure-name>_by_cluster.tsv` beside every figure, so
the table and the picture can never disagree.

### Single-Cell Adapters

The core knows nothing about scanpy or Seurat; these adapters do the translation
and import their dependencies only when called.

```python
from togoid.enrichment.adapters import (
    umap_dataframe_from_anndata, marker_genes_from_anndata,   # scanpy / AnnData
    umap_dataframe_from_csv, marker_genes_from_csv,           # Seurat, via CSV
)

import scanpy as sc
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

embedding = umap_dataframe_from_anndata(adata, cluster_key="leiden")
markers   = marker_genes_from_anndata(adata, top_n=100, pval_cutoff=0.05)
```

### Enrich Command

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

# Mouse instead of human
togoid enrich --genes-file genes.txt --preset reactome --taxonomy 10090
```

### Choosing a Background

This is the decision that most affects the results.

- **Default** (`background=None`): every gene in the library, i.e. every gene in
  your experiment that TogoID could annotate. This asks *which terms distinguish
  this cluster from the rest of the experiment* — usually the right question for
  cell types.
- **Explicit**: pass a larger universe for the conventional *over-represented
  relative to the genome* question. Expect many more significant hits.

`enrich_clusters` shares one background across clusters, which is what makes the
FDR values comparable between them.

## Configuration

### Environment Variables

- `TOGOID_API_ENDPOINT` - TogoID API base URL (default: https://api.togoid.dbcls.jp)
- `TOGOID_GRASP_ENDPOINT` - GRASP GraphQL endpoint (default: https://dx.dbcls.jp/grasp-dev-togoid)

### Custom API Endpoints

```python
# Python
converter = TogoIDConverter(api_base_url="http://localhost:5000")

# CLI
togoid --api-url http://localhost:5000 convert --ids 1,9 --route ncbigene,ensembl_gene
```

## Requirements

- Python 3.7+
- requests >= 2.20.0
- pandas >= 1.0.0 (optional, for DataFrame format and `togoid[enrichment]`)
- matplotlib >= 3.5 (optional, for the UMAP figures; `togoid[plot]`)
- scanpy >= 1.9 and leidenalg >= 0.9 (optional, for the scRNA-seq example;
  `togoid[singlecell]`)

The enrichment analysis itself runs on a bare install: the hypergeometric test
and the BH-FDR correction use only the standard library.

## Testing

This package includes comprehensive test scripts to verify all functionality:

```bash
# Test Python library examples
python3 test_readme_examples.py

# Test CLI examples
bash test_cli_examples.sh

# Test the enrichment analysis (no network access required)
python3 test_enrichment.py

# Test the enrichment analysis against the live TogoID API
python3 test_enrichment_api.py
```

See [TESTING.md](TESTING.md) for what each script covers.

## License

MIT License

## Links

- [TogoID Website](https://togoid.dbcls.jp/)
- [TogoID API](https://api.togoid.dbcls.jp/)
- [GitHub Repository](https://github.com/togoid/togoid-lib-python)

## Credits

Developed by DBCLS (Database Center for Life Science)
