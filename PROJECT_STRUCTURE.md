# TogoID Python Library - Project Structure

```
togoid-lib-python/
├── togoid/                     # Main package
│   ├── __init__.py            # Package initialization, exports main classes
│   ├── __main__.py            # Enable 'python -m togoid' execution
│   ├── converter.py           # TogoIDConverter - ID conversion
│   ├── annotations.py         # AnnotationsConverter - Annotations retrieval
│   ├── label_converter.py     # LabelConverter - Label to ID conversion
│   ├── cli.py                 # Unified CLI with all subcommands
│   ├── _ids.py                # Shared CURIE/local-ID helper
│   └── enrichment/            # Enrichment analysis and UMAP visualization
│       ├── __init__.py        # Public API
│       ├── genesets.py        # TogoID route -> GeneSetLibrary
│       ├── stats.py           # Hypergeometric test, BH-FDR (no dependencies)
│       ├── analysis.py        # enrich() / enrich_clusters() and result types
│       ├── presets.py         # Reactome / GO / MONDO route presets
│       ├── adapters.py        # AnnData and CSV adapters
│       ├── plot.py            # UMAP word-cloud figures
│       └── cli.py             # `togoid enrich` subcommand
│
├── examples/
│   └── scRNAseq_enrichment/   # End-to-end scRNA-seq pipeline
│
├── docs/
│   └── enrichment.md          # Enrichment analysis manual
│
├── test_readme_examples.py    # Python library tests
├── test_cli_examples.sh       # CLI tests
├── test_enrichment.py         # Enrichment tests (offline)
├── test_enrichment_api.py     # Enrichment tests (TogoID API)
├── TESTING.md                 # Testing documentation
├── README.md                  # Main documentation
├── pyproject.toml             # Package configuration
└── requirements.txt           # Dependencies

```

## Package Components

### togoid/converter.py
- **Class**: `TogoIDConverter`
- **Source**: Integrated from togoid-api-wrapper
- **Features**:
  - ID conversion between biological databases
  - Multiple output formats (JSON, dict, table, DataFrame)
  - Search, route, count, and config operations

### togoid/annotations.py
- **Class**: `AnnotationsConverter`
- **Source**: Integrated from togoid-id-to-annotations
- **Features**:
  - GraphQL-based annotation retrieval
  - Field listing and filtering
  - Multiple output formats (table, CSV, JSON)

### togoid/label_converter.py
- **Class**: `LabelConverter`
- **Source**: Integrated from togoid-label-to-id
- **Features**:
  - Label to ID conversion
  - Automatic API detection (SPARQList vs PubDictionaries)
  - Support for gene symbols and disease labels

### togoid/enrichment/
- **Modules**: `genesets`, `stats`, `analysis`, `presets`, `adapters`, `plot`, `cli`
- **Features**:
  - Any TogoID route becomes a gene-set library (`build_gene_sets`)
  - Over-representation analysis with hypergeometric test and BH-FDR
  - Word-cloud visualization of results on a UMAP embedding
  - Presets for Reactome, GO and MONDO; adapters for scanpy and Seurat
- **Dependencies**: the core needs nothing beyond `requests`; pandas, matplotlib
  and scanpy are optional extras, imported lazily
- **Documentation**: [docs/enrichment.md](docs/enrichment.md)

### togoid/cli.py
- **Unified CLI** with subcommands:
  - `convert` - ID conversion
  - `label2id` - Label to ID conversion
  - `annotate` - Get annotations
  - `search` - Search databases
  - `lookup` - Lookup IDs
  - `route` - Find conversion routes
  - `count` - Count mappings
  - `config` - Get configuration
  - `enrich` - Enrichment analysis of gene lists (defined in `togoid/enrichment/cli.py`)

## Installation

### Using uv (recommended)

```bash
# Create virtual environment
uv venv
source .venv/bin/activate

# Install in development mode
uv pip install -e .

# With pandas support
uv pip install -e ".[pandas]"

# With the enrichment extras
uv pip install -e ".[plot]"          # + matplotlib, for the UMAP figures
uv pip install -e ".[singlecell]"    # + scanpy, for the scRNA-seq example

# With development tools
uv pip install -e ".[dev]"
```

### Using pip

```bash
# Install in development mode
pip install -e .

# With pandas support
pip install -e ".[pandas]"
```

## Usage

### As Python Library

```python
from togoid import TogoIDConverter, AnnotationsConverter, LabelConverter

# ID conversion
converter = TogoIDConverter()
result = converter.convert(ids=["1", "9"], route=["ncbigene", "ensembl_gene"])

# Label conversion
label_converter = LabelConverter()
results = label_converter.convert(labels=["BRCA1", "TP53"], dataset="ncbigene", taxonomy="9606")

# Annotations
annotator = AnnotationsConverter()
annotations = annotator.execute_query(
    dataset_name="ncbigene",
    ids=["672", "7157"],
    fields=["label", "gene_synonym"],
    filters={}
)
```

### As CLI Tool

```bash
# ID conversion
togoid convert --ids 1,9 --route ncbigene,ensembl_gene

# Label to ID conversion
togoid label2id --labels "BRCA1,TP53" --dataset ncbigene --taxonomy 9606

# Get annotations
togoid annotate --dataset ncbigene --ids 672,7157 --field gene_synonym

# Configuration
togoid config dataset ncbigene
```

## Testing

```bash
# Test Python library
python3 test_readme_examples.py

# Test CLI
bash test_cli_examples.sh

# Test enrichment (offline, then against the live API)
python3 test_enrichment.py
python3 test_enrichment_api.py
```

See [TESTING.md](TESTING.md) for detailed testing documentation.

## Migration from Previous Structure

This package consolidates three previously separate repositories:
- `togoid-api-wrapper` → `togoid/converter.py`
- `togoid-id-to-annotations` → `togoid/annotations.py`
- `togoid-label-to-id` → `togoid/label_converter.py`

All functionality is preserved and enhanced with:
- Unified package interface
- Consistent API design
- Comprehensive CLI
- Complete test coverage
