# DownloadExpressionData

A reusable DSH skill for turning public bulk expression data into a clean, documented RData object.

Repository: <https://github.com/DK-Zhao1999/DownloadExpressionData>

## Install as a DSH skill

Clone the repository as a skill bundle under your DSH home:

```bash
git clone https://github.com/DK-Zhao1999/DownloadExpressionData.git \
  "$HOME/.dsh/skills/geo-expression-matrix"
```

Then load it in a DSH session:

```text
geo-expression-matrix
```

It covers the common case behind most requests:

1. GEO-like database → processed expression matrix (GEO is automated; ArrayExpress, GDC, recount3, ARCHS4 and local files use the same downstream path).
2. Bulk RNA-seq → prefer author TPM; never silently relabel counts or FPKM as TPM.
3. Missing values → `impute::impute.knn`.
4. Negative values → `rt <- rt - min(rt)`.
5. Values outside `[0, 30]` → `rt <- log2(rt + 1)`.
6. Save `rt` (a data.frame; rows = clean gene symbols, columns = sample names), `metadata`, and `dataset_info` as `<accession>_expression.RData`. Duplicate gene symbols are merged by mean; whitespace is removed; `"A /// B"` becomes `"A"`. Excel date corruption of `MARCH1`–`MARCH11`, `SEPT1`–`SEPT15`, and `DEC1`–`DEC2` (for example `1-Mar`, `1-Sep`, `1-Dec`) is restored from a fixed table.
7. Ask the user whether to remove normal samples.
8. Ask whether clinical/prognosis data should be prepared, and provide a schema.

## Requirements

- Python 3.10+ with `pandas` and `numpy`.
- R 4.x with:
  - `impute` (Bioconductor)
  - optional `openxlsx` for `.xlsx` input
  - optional `SummarizedExperiment` for RDS objects

Install `impute` if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("impute")
```

## Quick start

```bash
# 1. Download a GEO series and inspect candidate matrices
python scripts/geo_fetch.py fetch GSE50760 --outdir data/GSE50760

# 2. Prepare RData (FPKM because this series provides FPKM; TPM preferred when available)
Rscript scripts/prepare_expression.R \
  --input data/GSE50760/per_sample_matrix_FPKM.tsv.gz \
  --metadata data/GSE50760/samples.tsv \
  --accession GSE50760 --platform GPL11154 --species human \
  --value-type fpkm \
  --outdir processed/GSE50760

# 3. Only after the user confirms normal-tissue removal:
Rscript scripts/subset_expression.R \
  --rdata processed/GSE50760/GSE50760_expression.RData \
  --outdir processed/GSE50760/tumor_only \
  --column tissue \
  --drop-regex "normal|adjacent|non-tumor|healthy|benign" \
  --label tumor_only \
  --dry-run

# 4. Clinical table validation (optional)
Rscript scripts/validate_clinical.R \
  --clinical clinical.csv \
  --rdata processed/GSE50760/GSE50760_expression.RData \
  --outdir processed/GSE50760/clinical \
  --time-unit months
```

If per-sample files need manual assembly, use:

```bash
python scripts/geo_fetch.py assemble \
  --input data/GSE50760/per_sample \
  --out data/GSE50760/per_sample_matrix.tsv.gz \
  --duplicate mean
```

Inspect `matrix_candidates.tsv` before choosing an assembled file. `--duplicate error` is the safe default; use `mean` for the final gene-level matrix.

## Excel date-corrupted gene symbols

Excel sometimes converts valid gene names into dates. The pipeline restores the fixed set:

```text
1-Mar ... 11-Mar  -> MARCH1 ... MARCH11
1-Sep ... 15-Sep  -> SEPT1 ... SEPT15
1-Dec ... 2-Dec   -> DEC1 ... DEC2
```

The full 28-row map is in `assets/excel_date_gene_map.csv`. Matching is case-insensitive and accepts leading-zero forms such as `01-Mar`. The restoration happens before duplicate symbols are merged by mean, so `1-Mar` and `MARCH1` become one gene. Other date-like strings are kept, counted in `n_unresolved_date_like`, and reported as a warning rather than guessed.

## Output contract

```r
load("processed/GSE50760/GSE50760_expression.RData")
# rt           data.frame: rows = gene symbols, columns = sample names only
#              no NA/Inf; every column numeric; no duplicate genes
#              gene symbols have no whitespace; "A /// B" becomes "A"
#              duplicate gene symbols are merged by mean
# metadata     data.frame: one row per sample, sample_id first
# dataset_info list: source, platform, value types, transformations, KNN params,
#              gene_id_rule, date_gene_fix, collapse (mean + n_collapsed), QC
```

`dataset_info$output_value_type` is the source of truth, for example:

- `TPM`
- `log2(TPM + 1)`
- `log2(FPKM + 1)`
- `log2(intensity)`

## Files

```text
SKILL.md                         operational instructions and gates
scripts/geo_fetch.py             GEO fetch, supplement assembly, candidate ranking
scripts/prepare_expression.R     KNN → shift → log2 → RData + QC
scripts/subset_expression.R      explicit normal/adjacent sample removal
scripts/validate_clinical.R      clinical/prognosis validation and RData export
references/database-adapters.md  GEO/ArrayExpress/GDC/recount3 notes
references/clinical-schema.md    clinical field definitions and pitfalls
assets/clinical_template.csv     broad clinical template
assets/prognosis_template.csv    survival-only template
assets/samples_template.csv      sample-to-patient mapping template
assets/excel_date_gene_map.csv   fixed Excel-date -> MARCH/SEPT/DEC gene map
```

## Scope and limits

- Bulk expression only. Single-cell, proteomics, and metabolomics are out of scope.
- GEO download is automated. Other repositories are downloaded with their own tools/API and then passed via `--input`.
- Raw FASTQ/BAM/CEL processing is not automated. If only raw data exist, ask the user before starting a quantification workflow.
- The pipeline does not perform batch correction, quantile normalization, or differential analysis. Those are separate analysis decisions.
- The clinical script validates structure and coding; it does not invent or infer missing clinical values.

## Test data

The helpers were exercised on:

- A synthetic TPM-like matrix with NA, negative values, and a max >30 (KNN + shift + log2).
- A synthetic matrix with `1-Mar`, `1-Sep`, `1-Dec`, spaces, `///`, and duplicate symbols, verifying fixed-date restoration and mean merging.
- A synthetic array matrix with a GEO-style annotation preamble (probe mapping + no log2 because range was inside `[0, 30]`).
- GEO `GSE50760` (54 samples; per-sample FPKM supplements; normal/tumor/liver-metastasis metadata).

Test artifacts are intentionally excluded from the skill bundle.
