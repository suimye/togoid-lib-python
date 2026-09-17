# Testing Guide

This document describes how to test the TogoID Python library to ensure all functionality works as expected.

## Test Files

### Core Test Scripts

### 1. `test_readme_examples.py`

Python script that tests all the Python library examples from README.md.

**Run the tests:**

```bash
python3 test_readme_examples.py
```

**What it tests:**
- Quick Start examples (ID conversion, label conversion, annotations)
- Different output formats (JSON, dict, table, DataFrame)
- Label to ID conversion with automatic API detection
- Annotations with filters
- Various API methods (count, config, etc.)

**Expected output:**

```
============================================================
Testing README.md Examples
============================================================

=== Testing Quick Start - ID Conversion ===
✓ JSON format
✓ Dict format
✓ Table format
✓ DataFrame format

... (more tests) ...

============================================================
Test Results: 7 passed, 0 failed
============================================================
```

### 2. `test_cli_examples.sh`

Bash script that tests all CLI examples from README.md.

**Run the tests:**

```bash
bash test_cli_examples.sh
# or
./test_cli_examples.sh
```

**What it tests:**
- All CLI commands: convert, label2id, annotate, config, count
- File input/output
- Different output formats
- Various command-line options

**Expected output:**

```
======================================================================
Testing CLI Examples from README.md
======================================================================

=== Quick Start CLI Examples ===

Testing: Convert IDs (basic)
Command: timeout 15 python3 -m togoid convert --ids 1,9 --route ncbigene,ensembl_gene
✓ PASSED

... (more tests) ...

======================================================================
Test Results Summary
======================================================================
Passed:   XX
Failed:   0
Warnings: X
======================================================================
All critical tests passed!
```

### Feature-Specific Test Files

The following test files validate specific features and bug fixes:

### 3. `test_customer_feedback.py`

Tests format output improvements from initial customer feedback.

**What it tests:**
- Dict format returns proper structure with ids, route, results
- Table format includes source-target pairs
- DataFrame format uses dataset names as column headers
- Annotation columns appear in correct positions

**Run:**
```bash
python3 test_customer_feedback.py
```

### 4. `test_config_list_targets.py`

Tests the `config_list_targets()` method for finding reachable datasets.

**What it tests:**
- Retrieving list of datasets reachable from a source in one hop
- Return type is a list of strings
- Works for various datasets (ncbigene, ensembl_gene, etc.)

**Run:**
```bash
python3 test_config_list_targets.py
```

### 5. `test_feedback_fixes.py`

Tests customer feedback improvements including route handling and parameter changes.

**What it tests:**
- Routes with 3+ datasets use `report=full` to show all intermediate IDs
- `format="dict"` raises error for routes with 3+ datasets
- `label_types` parameter accepts list format
- `annotator.execute_query()` filters parameter is optional

**Run:**
```bash
python3 test_feedback_fixes.py
```

### 6. `test_route_suggestion.py`

Tests the route suggestion feature for error handling.

**What it tests:**
- Error messages suggest alternative routes when datasets aren't connected
- Suggestions include 2-hop and 3-hop routes
- Error handling for invalid dataset pairs

**Run:**
```bash
python3 test_route_suggestion.py
```

### 7. `test_enrichment.py`

Tests `togoid.enrichment` without touching the network. Safe to run in CI.

**What it tests:**
- The hypergeometric survival function, against hand-computed values and, when
  SciPy is installed, against `scipy.stats.hypergeom.sf`
- Benjamini-Hochberg correction, including monotonicity, and against
  `statsmodels.stats.multitest.multipletests` when it is installed
- `GeneSetLibrary`: size filtering, reverse mapping, JSON save/load round trip
- `enrich()` and `enrich_clusters()`: term-size filters, background handling,
  ranking, CSV output, and the pandas views
- Term selection: `top_n`, FDR and p-value cut-offs, label truncation
- The word-cloud layout, including a check that **no two placed labels overlap**

**Run:**
```bash
python3 test_enrichment.py
```

**Expected output:**

```
============================================================
togoid.enrichment - offline tests
============================================================

=== Hypergeometric test ===
  ✓ matches a hand-computed value
  ✓ agrees with scipy (worst relative error 1.81e-11)
... (more tests) ...

============================================================
Test Results: 102 passed, 0 failed
============================================================
```

### 8. `test_enrichment_api.py`

Tests `togoid.enrichment` against the live TogoID API, checking that the
documented routes still work and that the results stay biologically sensible.

**What it tests:**
- Gene symbols resolve to the expected NCBI Gene IDs
- All three preset routes build gene sets:
  `ncbigene → uniprot → reactome_pathway`, `ncbigene → uniprot → go`,
  `ncbigene → medgen → mondo`
- T cell, B cell and myeloid marker genes come out enriched for TCR, BCR and
  innate immune terms respectively — in both Reactome and GO
- The GO aspect filter partitions terms with no overlap between aspects
- An explicit `route=` matches the equivalent preset
- A broken route raises `RuntimeError` instead of silently returning nothing
- CURIE-prefixed and raw API output give identical results after `local_id()`
  normalisation, so the code works before and after the API's `?prefix` ships
- A cached gene-set library reproduces identical p-values offline

**Run:**
```bash
python3 test_enrichment_api.py
```

**Expected output:**

```
============================================================
togoid.enrichment - TogoID API tests
============================================================

=== Reactome pathways ===
  ✓ builds a gene-set library (86 pathways)
  ✓ T cluster finds TCR biology
  ✓ B cluster finds B cell receptor biology
  ✓ myeloid cluster finds innate immune biology
... (more tests) ...

============================================================
Test Results: 42 passed, 0 failed, 0 warnings
============================================================
```

### 9. End-to-end example pipeline

`examples/scRNAseq_enrichment/` runs the whole use case on real data. It is not
part of the automated suite — it needs a downloaded 10x dataset and takes a few
minutes — but it is the best end-to-end check of the enrichment feature.

```bash
pip install "togoid[singlecell]"
cd examples/scRNAseq_enrichment
./run_pipeline.sh path/to/filtered_feature_bc_matrix results
```

On the 10x Genomics public PBMC dataset (3,733 cells after QC, 13 Leiden
clusters) this produces, among other outputs:

```
results/03_enrichment_reactome_significant.csv   131 significant pathway-cluster pairs
results/03_enrichment_go_significant.csv          45 significant GO BP terms
results/04_umap_enrichment_reactome_top3.pdf
results/04_umap_enrichment_go_top3.pdf
results/04_umap_enrichment_mondo_top3.pdf
```

The Reactome figure — cluster UMAP on the left, enriched pathways on the right:

![UMAP with enriched Reactome pathways](https://raw.githubusercontent.com/suimye/togoid-lib-python/docs-figures/umap_enrichment_reactome.png)

The same analysis through a different route, GO biological process, which
independently recovers the platelet and T cell assignments:

![UMAP with enriched GO biological processes](https://raw.githubusercontent.com/suimye/togoid-lib-python/docs-figures/umap_enrichment_go.png)

Sanity checks worth making on the output:

- Cluster 2 should be enriched for platelet biology in both Reactome
  (*Platelet degranulation*) and GO (*platelet aggregation*)
- Clusters expressing MHC class II genes should show *MHC class II antigen
  presentation*
- Every `fdr` must be at least its `pvalue`, and `overlap_count` must never
  exceed `term_size` or `query_size` — `test_enrichment_api.py` asserts this
- Translation pathways recurring across several clusters are expected, not a
  bug: ribosomal protein genes dominate those marker lists

## Running Specific Tests

### Test only Python library functionality

```bash
python3 test_readme_examples.py
```

### Test only CLI functionality

```bash
bash test_cli_examples.sh
```

### Test both core scripts

```bash
python3 test_readme_examples.py && bash test_cli_examples.sh
```

### Test all feature-specific tests

```bash
python3 test_customer_feedback.py && \
python3 test_config_list_targets.py && \
python3 test_feedback_fixes.py && \
python3 test_route_suggestion.py
```

### Test the enrichment analysis

```bash
# Offline only (no network access needed)
python3 test_enrichment.py

# Including the live API tests
python3 test_enrichment.py && python3 test_enrichment_api.py
```

### Test everything

```bash
python3 test_readme_examples.py && \
bash test_cli_examples.sh && \
python3 test_customer_feedback.py && \
python3 test_config_list_targets.py && \
python3 test_feedback_fixes.py && \
python3 test_route_suggestion.py && \
python3 test_enrichment.py && \
python3 test_enrichment_api.py
```

## Understanding Test Results

### ✓ PASSED (Green)
The test completed successfully with expected output.

### ✗ FAILED (Red)
The test failed. This indicates a potential issue with the implementation.

### ⚠ WARNING (Yellow)
The test completed but may have issues. Often occurs when:
- External APIs are unavailable (e.g., TogoID API endpoints)
- Network timeouts
- Optional features not installed (e.g., pandas)

## Continuous Testing

You can run these tests after:
- Making changes to the codebase
- Updating dependencies
- Installing the package in a new environment
- Before creating a release

## Test Coverage

The test scripts cover:

### Python Library Tests
- [x] TogoIDConverter class
  - [x] convert() with different formats
  - [x] count() method
  - [x] config methods
- [x] LabelConverter class
  - [x] Automatic API detection
  - [x] SPARQList conversion
  - [x] PubDictionaries conversion
- [x] AnnotationsConverter class
  - [x] list_fields() method
  - [x] execute_query() method
  - [x] Filters support
- [x] togoid.enrichment
  - [x] Hypergeometric test and BH-FDR, verified against SciPy/statsmodels
  - [x] build_gene_sets() over arbitrary routes
  - [x] Reactome, GO and MONDO presets
  - [x] GeneSetLibrary JSON caching
  - [x] enrich() and enrich_clusters()
  - [x] UMAP word-cloud layout (no overlapping labels)
  - [x] CURIE prefix tolerance

### CLI Tests
- [x] convert subcommand
  - [x] Basic conversion
  - [x] Different output formats
  - [x] File output
- [x] label2id subcommand
  - [x] Basic conversion
  - [x] File input
  - [x] CSV output
  - [x] PubDictionaries usage
- [x] annotate subcommand
  - [x] Get annotations
  - [x] List fields
  - [x] Filters
  - [x] CSV output
  - [x] File input
- [x] config subcommand
- [x] count subcommand
- [x] enrich subcommand

## Known Limitations

Some API endpoints may be unavailable or return 404 errors:
- `search databases`
- `route` finding
- `lookup id`

These are API-side issues and not implementation problems. The tests will show warnings for these cases.

## Troubleshooting

### Import Errors

If you get import errors:

```bash
# Install the package in development mode
pip install -e .

# Or add to PYTHONPATH
export PYTHONPATH=/home/souta/projects/togoid-lib-python:$PYTHONPATH
```

### Timeout Errors

If tests timeout:
- Check internet connection
- Verify TogoID API is accessible: https://api.togoid.dbcls.jp/
- Increase timeout values in test scripts

### Pandas Not Found

If DataFrame tests are skipped:

```bash
# Install pandas
pip install pandas

# Or install with pandas support
pip install -e ".[pandas]"
```

## Contributing

When adding new features:
1. Add examples to README.md
2. Add corresponding tests to test scripts
3. Run tests to ensure they pass
4. Update this document if needed
