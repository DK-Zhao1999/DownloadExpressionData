# Database adapters

The downstream contract is always the same: a numeric **genes × samples** matrix plus a sample table with `sample_id`. The adapters below differ only in how those two objects are obtained.

## GEO

### Series files

For `GSExxxxx`, the series directory is:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE<prefix>nnn/GSExxxxx/
```

Examples:

```text
GSE50760 -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50760/
GSE1000  -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/
GSE999   -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE999/
```

Useful subdirectories:

```text
matrix/GSExxxxx_series_matrix.txt.gz   sample metadata + optional expression table
suppl/                                 series-level supplementary files
```

- The series matrix is the first choice for array data and for RNA-seq series where the authors deposited a processed table.
- Many RNA-seq series have an **empty** series-matrix table and put per-sample FPKM/TPM/count files in each GSM's `suppl/` directory. These are linked from `!Sample_supplementary_file_*`; `geo_fetch.py` follows those links.
- `GSE..._RAW.tar` often contains CEL files, FASTQ, or BAM, not a clean matrix. Inspect `filelist.txt` before downloading hundreds of MB.
- A series can contain multiple platforms (`GPL...`). Do not combine them without a platform/batch plan.

### Platform annotation (arrays)

Platform directory:

```text
https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz
https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL11nnn/GPL11154/annot/
```

For 1–3 digit GPL numbers, the group is `GPLnnn`; for longer numbers it is `GPL<prefix>nnn`. Annotation availability varies; RNA-seq platforms often have no `annot/` directory. `geo_fetch.py fetch --annot` tries the common layout and skips missing files.

Annotation format:

```text
^Annotation
!Annotation_platform = GPL96
#ID = ...
!platform_table_begin
ID    Gene symbol    Gene title    ...
1007_s_at    MIR4640///DDR1    ...
!platform_table_end
```

The preparation script can read this preamble and map probe IDs to `Gene symbol`. Composite mappings (`A /// B`) always keep the first symbol `A`; whitespace is removed; duplicate gene symbols are merged by mean in the final data.frame.

Excel can silently convert a few valid gene symbols into dates. The fixed map in `assets/excel_date_gene_map.csv` restores:

```text
1-Mar ... 11-Mar -> MARCH1 ... MARCH11
1-Sep ... 15-Sep -> SEPT1 ... SEPT15
1-Dec ... 2-Dec  -> DEC1 ... DEC2
```

This happens before duplicate merging, so a corrupted `1-Mar` and a correct `MARCH1` collapse to one row. Other date-like values are flagged, not guessed. If the gene column was read as a real Excel date/serial, re-export it as text first.

If automatic fetching fails, use Bioconductor:

```r
library(GEOquery)
gpl <- getGEO("GPL96", destdir = "data/GSE1000")
annotation <- Table(gpl)
write.table(annotation, "data/GSE1000/platform_annotation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

Then pass `--annotation data/GSE1000/platform_annotation.tsv` and, if needed, `--annotation-probe-col` / `--annotation-gene-col`.

### GEO pitfalls

- Series matrix values are **not guaranteed** to be normalized or log2. Read the `!Sample_data_processing` lines.
- RNA-seq series may provide FPKM, RPKM, or counts even when the title says “expression profiling”. Check the file names and headers.
- Sample titles are free text; `normal` may appear as “normal-looking surrounding colonic epithelium”. Classify with the metadata column after showing the user the values.
- Some series have no processed matrix at all (only SRA FASTQ). Stop and ask before quantification.
- NCBI may return 403/429 under aggressive access. Keep `--delay` (default 0.2 s), avoid parallel downloads, and retry later.

## ArrayExpress / BioStudies

ArrayExpress is served through BioStudies. Use the study page/API to locate the author-processed matrix and sample annotation:

```text
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-xxxx
https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-xxxx
```

Prefer:

- `*.processed.*` matrix files (often `.txt`/`.tsv`/`.zip`)
- `*.sdrf.txt` sample annotation
- the study's `IDF`/`ADF` only if needed to understand processing

Then run `prepare_expression.R --input ... --metadata ...`. Record the exact file names and accession in `dataset_info`.

Do not scrape the HTML table into a matrix; download the deposited file.

## TCGA / GDC

GDC files are per-sample, not one ready-made matrix. Use the GDC API or portal:

```text
https://api.gdc.cancer.gov/files
https://portal.gdc.cancer.gov/
```

For bulk expression:

- `data_type = Gene Expression Quantification`
- `experimental_strategy = RNA-Seq`
- Prefer `STAR - Counts` files; each file contains `gene_id`, `gene_name`, `unstranded`, `stranded_first`, `stranded_second`, `tpm_unstranded`, `fpkm_unstranded`, etc.
- Use the **`tpm_unstranded`** column for TPM. Do not use `unstranded` as TPM.
- Merge files by `gene_id` (or `gene_name` after checking duplicates); sample IDs come from the file/case metadata.
- Clinical data: use the GDC clinical TSV or XML, map `case_submitter_id` to `patient_id` and `sample_submitter_id` to `sample_id`.

The same preparation script applies; pass `--value-type tpm`.

## recount3

```text
https://rna.recount.bio/
```

recount3 provides uniformly processed human/mouse RNA-seq. Use the provider's documented TPM/processed layer; confirm whether the matrix is TPM, counts, or coverage. Do not assume a matrix is TPM from the project name.

## ARCHS4

```text
https://maayanlab.cloud/archs4/
```

ARCHS4 provides processed human/mouse expression data. Check the download page for the exact value type and gene ID mapping; pass the chosen matrix via `--input`.

## Local files

Any tabular matrix can enter the pipeline:

```bash
Rscript scripts/prepare_expression.R \
  --input matrix.tsv.gz \
  --metadata samples.tsv \
  --value-type tpm \
  --outdir processed/local
```

Before running, verify:

- first column or row names are gene/probe IDs;
- all remaining columns are numeric;
- columns are samples and match metadata;
- the value type is known, not guessed.

Use `python scripts/geo_fetch.py inspect <file-or-dir>` to preview and rank candidates. Use `--transpose yes` only when the input is samples × genes.

## Raw data (not automated)

If only raw data exist:

- **Array CEL**: normalize with `affy`/`oligo` (RMA) first, then pass the normalized matrix.
- **RNA-seq FASTQ**: quantify with salmon/kallisto/featureCounts/STAR, build TPM, then pass the matrix.
- **BAM**: quantify counts first; do not treat BAM depth as TPM.

These are separate, compute-heavy workflows. Ask the user before starting them; they are outside the default skill path.
