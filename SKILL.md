---
name: geo-expression-matrix
description: 从 GEO 等公共数据库获取 bulk RNA-seq/array 表达矩阵，按 TPM 优先、KNN 填补 NA、负值平移、必要时 log2(rt+1) 整理为 genes×samples 的 data.frame RData；随后追问是否删除 normal 组织、是否整理预后/临床数据。
whenToUse: 当用户要求下载 GEO/GSE、ArrayExpress、TCGA/GDC 等公共表达数据，或要求把表达矩阵整理成 RData、做 KNN 填补/对数转换、删除正常样本、整理生存临床表时使用。
---

# 公共数据库 bulk 表达矩阵 → RData

把“下载 → 选矩阵 → 统一方向 → KNN 补 NA → 负值平移 → 条件 log2 → RData → 追问正常样本/临床数据”做成可复现流程。默认只处理 **bulk RNA-seq 或 microarray**；单细胞、蛋白组、代谢组不适用。

**最终 `RData` 里的 `rt` 是 `data.frame`，不是 `matrix`**：行名是 cleaning 后的 gene symbols，列名是样本名；重复 gene 按 mean 合并，gene symbol 不含空格，`"A /// B"` 只保留 `A`。

## 0. 不可违反的规则

1. **RNA-seq 优先用 TPM，但不能把 counts/FPKM 冒充 TPM。** 如果作者只提供 raw counts，先停下来问用户：找作者 TPM、用 recount3/ARCHS4，还是先定量 raw reads。只有用户明确接受时才可 `--allow-counts`，并在结果里标注“log2(counts+1)，不是 TPM”。
2. **不静默删除样本。** normal 组织是否删除必须由用户确认；删除前列出将删样本、按组织分组的计数，并保留 removed-samples 表。
3. **不伪造临床数据。** 没有的字段就留空；不要从标题猜生存时间、事件或分期。临床表必须能追溯到原始文件/补充表。
4. **保留原始文件。** 原始下载、矩阵候选、QC、处理日志、`sessionInfo()` 全部保留；不要覆盖 raw 文件。
5. **每个转换都要记录。** `dataset_info$transformations` 必须能解释 `rt` 是 TPM 还是 log2(TPM+1)，以及是否做过 shift/KNN/subset。
6. **最终对象固定为 `data.frame`，方向为 genes × samples。** `rt` 的行名是 **gene symbols**，列名是 **样本名称**；`metadata$sample_id` 必须与 `colnames(rt)` 完全一致且同序。
7. **最终 gene symbol 必须干净。** 不允许重复；重复用**平均值**合并；不允许空格（包括 `"gene 1"` → `"gene1"`）；`"gene1 /// gene2"` 只保留第一个 `"gene1"`。不允许额外加 `gene_symbol` 列，gene symbols 放在 row names。

## 1. 开始前只问缺少的关键信息

如果用户没有给全，最多问这些，不要一次问一堆：

- 数据集：GEO accession（如 `GSE50760`）、ArrayExpress（`E-MTAB-...`）、GDC 项目，或本地矩阵路径。
- 物种与组织：human/mouse；肿瘤类型。
- 数据类型：RNA-seq / microarray；用户是否已知作者提供的是 TPM、FPKM 还是 counts。
- 分组变量：哪个 metadata 列代表 tumor/normal（如 `tissue_type`、`tissue`、`group`）。
- 输出目录：默认 `data/<accession>/` 与 `processed/<accession>/`。

如果用户只给 accession，先自动下载并检查，不要先要求用户回答所有问题。

## 2. 获取数据

### 2.1 GEO（首选自动路径）

先确认数据是 processed matrix，不是 raw reads：

```bash
python scripts/geo_fetch.py fetch GSE50760 \
  --outdir data/GSE50760 \
  --annot
```

- `--annot` 会尝试下载 GPL 平台注释，供 array 探针映射使用；RNA-seq 平台通常没有，脚本会提示。
- 默认只下载 series matrix 和 table-like supplementary；FASTQ/BAM/SRA/CEL/RAW.tar 默认跳过。没有 processed matrix 时才考虑 `--suppl all`，并先读 `fetch_manifest.json`。
- 下载完成后先看 `data/GSE50760/matrix_candidates.tsv`，不要凭文件名直接选。

关键文件：

```text
data/GSE50760/
  GSE50760_series_matrix.txt.gz      原始 series matrix
  series_matrix.tsv.gz               若 series matrix 自带表
  samples.tsv                        每个样本一行，含 tissue/characteristics
  series_meta.json                   标题、平台、设计、PMID
  suppl/                             系列级补充文件
  per_sample/<GSM>/                  每样本补充文件
  per_sample_matrix_<TYPE>.tsv.gz    脚本自动合并的每样本矩阵
  matrix_candidates.tsv              候选矩阵排名
  fetch_manifest.json                URL、大小、sha256、警告
```

如果 series matrix 没有表达表，脚本会把 series matrix 中 `!Sample_supplementary_file_*` 指向的每样本文件下载下来并尝试合并。合并失败时不要硬拼；先看 `matrix_candidates.tsv` 和 `fetch_manifest.json` 的警告。需要手动合并两列（gene + value）的每样本表时：

```bash
python scripts/geo_fetch.py assemble \
  --input data/GSE50760/per_sample \
  --out data/GSE50760/per_sample_matrix.tsv.gz \
  --duplicate mean
```

`--duplicate error` 是默认值；最终交付需要 `mean`。只有明确知道重复 gene ID 的含义时才用 `max`/`first` 做探索。

### 2.2 选择矩阵：先判断“值是什么”

| 情况 | 处理 |
|---|---|
| RNA-seq，作者提供 TPM | 首选。`--value-type tpm` |
| RNA-seq，只有 FPKM/RPKM | 可用，但必须标注 FPKM/RPKM，不能叫 TPM。`--value-type fpkm` |
| RNA-seq，只有 raw counts | **暂停**。先找 TPM/作者 processed matrix，或用 recount3/ARCHS4；确实要跑 count 流程才 `--allow-counts` |
| array，series matrix 已是 log2 强度 | 可用。`--value-type array_log2` |
| array，线性强度 | 可用。`--value-type array_linear` |
| 值范围/含义不明 | 先看分布和文件说明，问用户；不要猜 |

判断技巧（不是替代文件说明）：

- 文件名/表头出现 `TPM`、`FPKM`、`RPKM`、`counts`、`log2` 是最强证据。
- TPM/FPKM 通常有大量 0，max 可到 10^4–10^6；log2 array 通常 2–16。
- 有负值通常是 log2/中心化或已做 batch correction，仍按第 4 步规则处理。
- 表达矩阵的 `min/max` 只是线索，不能单独决定 value type。

### 2.3 其他数据库

- **ArrayExpress/BioStudies**：从 BioStudies/ArrayExpress 下载作者 processed matrix 和 sample annotation，然后走 `--input` 路径；不要手工抄表。
- **TCGA/GDC**：下载 `Gene Expression Quantification` 的 STAR counts 文件，文件内的 `tpm` 列可作为 TPM；不要把 `unstranded` counts 当 TPM。临床数据用 GDC clinical TSV/XML，按本 skill 的临床 schema 整理。
- **recount3/ARCHS4**：使用其文档明确的 TPM/ processed 层；先确认单位，再 `--input`。
- **本地矩阵**：直接 `--input`，仍需确认 genes × samples、ID 类型、value type。

详细数据库路径和注意事项见 `references/database-adapters.md`。

## 3. 统一矩阵并生成 RData

### 3.1 先检查输入

```bash
python scripts/geo_fetch.py inspect data/GSE50760
```

确认：

- 行是基因/探针，列是样本；如果是 samples × genes，准备 `--transpose yes`。
- 第一列是 gene/probe ID，或 row names 是 ID。
- 样本列名与 `samples.tsv$sample_id` 能对上。
- 没有把 annotation 列（`Gene title`、`Gene assignment`）误当表达值。

### 3.2 运行主处理

```bash
Rscript scripts/prepare_expression.R \
  --input data/GSE50760/per_sample_matrix_FPKM.tsv.gz \
  --metadata data/GSE50760/samples.tsv \
  --accession GSE50760 \
  --platform GPL11154 \
  --species human \
  --value-type fpkm \
  --outdir processed/GSE50760
```

array + 平台注释：

```bash
Rscript scripts/prepare_expression.R \
  --input data/GSE1000/series_matrix.tsv.gz \
  --metadata data/GSE1000/samples.tsv \
  --annotation data/GSE1000/platform_annotation/GPL96.annot.gz \
  --accession GSE1000 --platform GPL96 --species human \
  --value-type array_linear \
  --outdir processed/GSE1000
```

常用参数：

| 参数 | 用途 |
|---|---|
| `--gene-column` | 指定 gene/probe ID 列（列名或 1-based 序号） |
| `--transpose auto/yes/no` | 方向；默认 auto 根据 metadata 判断 |
| `--annotation` | GPL annotation 文件（array 探针映射） |
| `--annotation-probe-col` / `--annotation-gene-col` | 注释列名；不自动识别时显式指定 |
| `--value-type tpm/fpkm/rpkm/counts/array_linear/array_log2/normalized/unknown` | **必填**，决定输出标签和 counts guard |
| `--allow-counts` | 仅在明确接受 count 流程时使用 |
| `--knn-k` / `--knn-rowmax` / `--knn-colmax` / `--knn-maxp` | KNN 参数；默认 10/0.5/0.8/1500 |
| `--drop-all-na-samples` | 默认遇到全 NA 样本直接报错；明确要删才加 |
| `--export-tsv` / `--rds` | 额外导出 TSV/RDS；RData 始终生成 |
| `--force` | 覆盖已有输出 |

### 3.3 处理顺序（必须按此顺序）

1. **统一方向/ID**：genes × samples；探针映射到 gene symbol；`"A /// B"` 只取 `A`；删除所有空格；重复 symbol 按 `mean` 合并。
2. **NA 检查**：只要 `anyNA(rt)` 为真，就用 `impute::impute.knn`：
   - 默认 `k=10, rowmax=0.5, colmax=0.8, maxp=1500`，固定 `rng.seed=1`。
   - 全 NA 的基因行无法 KNN，脚本会删除并记录；全 NA 样本默认报错，不静默删。
   - 某行缺失 >50% 时，`impute.knn` 会按包的设计用 row mean 兜底；脚本在 `dataset_info$knn$mean_fallback_rows` 和 QC 中报告。
3. **负值**：如果 `min(rt) < 0`，执行 `rt <- rt - min(rt)`（即全体加 `-min(rt)`），记录 `negative_shift`。
4. **范围与 log2**：如果处理后的范围不在 `[0, 30]`（`min < 0` 或 `max > 30`），执行 `rt <- log2(rt + 1)`。在 `[0, 30]` 内则**不**再 log。
5. **保存**：`<accession>_expression.RData`，包含 `rt`、`metadata`、`dataset_info`。

注意：

- KNN 在原始尺度上做，然后才 shift/log2；这是本 SOP 的固定顺序，不要擅自把 log 提前。
- `log2(x+1)` 不是 TPM 归一化；它只是稳定方差/压缩范围。
- 如果 `value_type=array_log2` 且 max 仍 >30，脚本会警告“可能重复 log”；此时先确认输入是否真是 log2。

### 3.3.1 Excel 日期型 gene symbol 固定回补

Excel 会把部分合法 gene symbol 自动转成日期。脚本内置固定映射，在去空格、去 `///` 之后、重复合并之前执行：

| Excel 错误格式 | 正确 gene symbol |
|---|---|
| `1-Mar` … `11-Mar` | `MARCH1` … `MARCH11` |
| `1-Sep` … `15-Sep` | `SEPT1` … `SEPT15` |
| `1-Dec` … `2-Dec` | `DEC1` … `DEC2` |

完整表在 `assets/excel_date_gene_map.csv`，共 28 条：

```text
1-Mar -> MARCH1      2-Mar -> MARCH2      3-Mar -> MARCH3
4-Mar -> MARCH4      5-Mar -> MARCH5      6-Mar -> MARCH6
7-Mar -> MARCH7      8-Mar -> MARCH8      9-Mar -> MARCH9
10-Mar -> MARCH10    11-Mar -> MARCH11

1-Sep -> SEPT1       2-Sep -> SEPT2       3-Sep -> SEPT3
4-Sep -> SEPT4       5-Sep -> SEPT5       6-Sep -> SEPT6
7-Sep -> SEPT7       8-Sep -> SEPT8       9-Sep -> SEPT9
10-Sep -> SEPT10     11-Sep -> SEPT11     12-Sep -> SEPT12
13-Sep -> SEPT13     14-Sep -> SEPT14     15-Sep -> SEPT15

1-Dec -> DEC1        2-Dec -> DEC2
```

规则：

- 只回补这 28 个已知基因名；不要用通用“日期 → 月日”猜测。
- 同时接受大小写变体和前导零形式（如 `01-Mar`）。
- 回补后再合并重复基因，所以 `1-Mar` 和 `MARCH1` 会按 mean 合并。
- 其他日期样式（如 `13-Mar`、`1-Apr`）不在固定表内：保留原值，并在 QC 的 `n_unresolved_date_like` 和 warning 中报告，需人工确认。
- 如果 gene 列被 Excel 读成真正的 Date/序列号（不是 `1-Mar` 文本），先重导出为文本；脚本不会猜测序列号对应哪个基因。

### 3.4 输出与检查

```text
processed/GSE50760/
  GSE50760_expression.RData
  GSE50760_qc_summary.tsv
  GSE50760_sample_qc.tsv
  GSE50760_metadata.tsv
  GSE50760_processing_log.txt
  GSE50760_sessionInfo.txt
```

RData 契约：

```r
load("processed/GSE50760/GSE50760_expression.RData")
# rt:           data.frame
#               rownames = gene symbols (unique, no whitespace, no "///")
#               colnames = sample names only (no gene_symbol column)
#               every column numeric, no NA/Inf
# metadata:     data.frame, one row per sample, sample_id first
# dataset_info: list(accession, source, platform, technology, species,
#                    input_value_type, output_value_type, output_class,
#                    gene_id_rule, date_gene_fix, collapse, n_genes, n_samples,
#                    transformations, negative_shift, log2_applied, knn,
#                    created_at, r_version, input_files, input_md5)
stopifnot(
  is.data.frame(rt),
  all(vapply(rt, is.numeric, logical(1))),
  !anyNA(rt),
  nrow(metadata) == ncol(rt),
  identical(colnames(rt), as.character(metadata$sample_id)),
  !anyDuplicated(rownames(rt)),
  !any(grepl("[[:space:]]", rownames(rt))),
  !any(grepl("///", rownames(rt), fixed = TRUE))
)
rt[1:5, 1:5]                     # 行名是 gene symbol，列名是样本
dataset_info$output_value_type   # 例如 "log2(TPM + 1)"
dataset_info$gene_id_rule        # 记录 /// 和空格处理规则
```

QC 至少报告：输入/输出维度、NA 数、删除的全 NA 行、负数最小值与 shift、log2 是否执行、输出范围、KNN 参数、映射率/折叠数，以及 `n_date_genes_fixed` / `n_unresolved_date_like`。

## 4. 必须暂停并追问（硬性 gate）

全表达 data.frame `*_expression.RData` 生成并通过上面的检查后，**必须停下来**，用 `ask_user_question` 提问，不能直接继续删样本或编临床表。一次问两个问题（工具支持多问题时），否则逐条问：

### Q1：是否删除 normal 组织？

> 全表达 data.frame 已生成：`<path>`。是否删除 normal / adjacent / 非肿瘤样本，只保留肿瘤？

选项建议：

- **先不删除，保留全部样本（推荐先保留）**：便于核对分组和配对；后续要删再跑。
- **删除 normal/adjacent，仅保留肿瘤**：需要 metadata 中有可用的组织列。
- **按我提供的样本列表删除**：用户给 `keep-file`。

执行前：

1. 先看 `metadata` 里可用的组织列：`tissue_type`、`tissue`、`source_name_ch1`、`characteristics_ch1` 等；打印 `table()`。
2. 先 dry-run，列出将删除的样本：

```bash
Rscript scripts/subset_expression.R \
  --rdata processed/GSE50760/GSE50760_expression.RData \
  --outdir processed/GSE50760/tumor_only \
  --column tissue \
  --drop-regex "normal|adjacent|non-tumor|healthy|benign" \
  --label tumor_only \
  --dry-run
```

3. 用户确认后再去掉 `--dry-run` 执行。不要根据标题自行猜测 normal；metadata 没有组织列时，要求用户给 `sample_id -> tissue_type` 映射。
4. 明确“转移灶算不算肿瘤”“是否只要 primary”“多样本患者保留哪个样本”；默认保留所有恶性样本，并在 summary 里写清楚。
5. 结果：`<accession>_tumor_only.RData` + `_removed_samples.tsv` + `_summary.tsv`。`dataset_info$subset` 记录规则和删除列表。

### Q2：是否需要整理预后/临床数据？

> 是否需要我额外整理预后或临床数据？如果已有临床文件，请给路径；如果没有，我会给你模板和字段格式。

选项建议：

- **需要，按标准模板整理**：发 `assets/clinical_template.csv` 和 `assets/prognosis_template.csv`，让用户填或提供原始表。
- **不需要，只要表达矩阵**：结束，交付表达 RData + QC。
- **已有临床文件，帮我校验并合并**：读取并校验，输出 `clinical.RData`、`clinical_clean.csv`、`sample_map.tsv`。

如果用户要临床数据：

1. 先确认 endpoint：OS、PFS、DFS、DSS 中的哪些；时间单位是月还是天；事件编码。
2. 用 `references/clinical-schema.md` 的字段解释，不要自创字段名。
3. 有原始文件就用：

```bash
Rscript scripts/validate_clinical.R \
  --clinical clinical/PROJECT_clinical.csv \
  --rdata processed/GSE50760/GSE50760_expression.RData \
  --outdir processed/GSE50760/clinical \
  --time-unit months
```

4. 校验失败（事件不是 0/1、时间为负、event=1 但无时间、patient_id 重复）必须修表，不要静默改。
5. 如果样本 ID 对不上，输出 `sample_map.tsv`；让用户决定是删、改 ID，还是保留未匹配样本。

## 5. 交付前检查清单

- [ ] `rt` 是 `data.frame`，genes × samples，所有列 numeric，无 NA/Inf。
- [ ] `rownames(rt)` 是唯一 gene symbols，无空格、无 `///`；重复已按 mean 合并。
- [ ] 日期型 gene symbols 已按固定表回补；`n_unresolved_date_like` 为 0，或已在回复中说明。
- [ ] `colnames(rt)` 只有样本名称，没有额外 `gene_symbol` 列。
- [ ] `metadata$sample_id` 与 `colnames(rt)` 完全一致且同序。
- [ ] `dataset_info$output_value_type` 明确（`TPM`、`log2(TPM + 1)`、`log2(FPKM + 1)`、`log2(intensity)` 等）。
- [ ] `dataset_info$transformations` 记录 KNN/shift/log2/subset；没有转换则 `"none"`。
- [ ] 原始文件、候选矩阵、fetch manifest、QC、processing log、sessionInfo 都在。
- [ ] 如果做了 normal 删除：removed-samples 表、规则、删除前后 n 都在。
- [ ] 如果做了临床：时间单位带后缀、event 为 0/1、缺失留空、样本映射表在。
- [ ] 最终回复里明确说明：数据来源、平台、value type、样本数/基因数、是否 log2、是否删除 normal、临床单位，以及任何未解决的 caveat。

## 6. 常见失败与处理

| 现象 | 处理 |
|---|---|
| GEO 403/429 | 降低频率（`--delay`）、稍后重试；不要循环猛抓。仍失败时让用户用浏览器下载 series matrix 后走 `--input` |
| series matrix 没有表 | 看 `per_sample/`、`suppl/`、`fetch_manifest.json`；必要时 `--suppl all`，但先确认不是 FASTQ/BAM |
| 只有 raw counts | 不直接当 TPM；找作者 TPM/recount3/ARCHS4，或先定量 |
| 多个 GPL 平台 | 分平台处理或加 batch；不要无脑合并 |
| 探针映射率 <50% | 检查 GPL 注释和探针 ID；用 `--annotation-probe-col`/`--annotation-gene-col` 显式指定 |
| 一个样本 >80% NA | KNN 的 `colmax` 会拒绝；明确删除该样本或提高 `--knn-colmax`，并在 QC 说明 |
| 某基因 >50% NA | `impute.knn` 用 row mean 兜底；在 QC 报告，不要声称全是 KNN |
| 全是 NA 的基因 | 删除并记录；KNN 无法补 |
| gene symbol 有空格或 `A /// B` | 自动删除空格并只保留 `A`；检查 `qc_summary.tsv` 的 `n_multi_symbol_rows` |
| gene symbol 重复 | 默认按 mean 合并；检查 `n_duplicate_genes_collapsed` 和 `dataset_info$collapse` |
| gene symbol 是 `1-Mar`/`1-Sep`/`1-Dec` | 按固定表回补为 MARCH/SEPT/DEC；检查 `n_date_genes_fixed` |
| 其他日期型 gene symbol | 不回补，保留并在 `n_unresolved_date_like` 警告；人工确认 |
| `impute` 未安装 | `BiocManager::install("impute")` |
| 负值 shift 后仍为负 | 只可能是浮点误差；脚本 clamp 到 0；否则检查 Inf/异常值 |
| 已经是 log2 又触发 log2 | 检查 `--value-type array_log2` 与输入范围；必要时人工确认后再处理 |
| 样本 ID 对不上 | 先输出 mapping，不要自动改样本名；让用户确认 |
| 临床 event=1 但 time 缺失 | 校验失败；回到原始表修，不能填 0 |

## 7. 数据库适配速查

| 来源 | 取什么 | 注意 |
|---|---|---|
| GEO | series matrix / per-sample processed files / `--annot` | 先看 `matrix_candidates.tsv`；RNA-seq 常是 FPKM/counts，不一定是 TPM |
| ArrayExpress/BioStudies | processed matrix + sample annotation | 以 BioStudies 文件为准；不要手抄 |
| TCGA/GDC | STAR counts 文件中的 `tpm` 列；clinical TSV/XML | `unstranded` counts ≠ TPM；临床字段按 schema |
| recount3 | 文档明确的 TPM/processed 层 | 确认单位后再 `--input` |
| ARCHS4 | 文档明确的 processed 层 | 确认 human/mouse 与基因 ID |
| 本地文件 | `--input` | 仍需确认方向、ID、value type |

更细的下载路径、GDC/ArrayExpress 说明见 `references/database-adapters.md`；临床字段定义见 `references/clinical-schema.md`。

## 8. 本 skill 的文件

```text
geo-expression-matrix/
  SKILL.md
  README.md
  scripts/
    geo_fetch.py
    prepare_expression.R
    subset_expression.R
    validate_clinical.R
  references/
    database-adapters.md
    clinical-schema.md
  assets/
    clinical_template.csv
    prognosis_template.csv
    samples_template.csv
    excel_date_gene_map.csv
```

最终回复使用用户的语言，先给结论和文件路径，再给关键 QC 数字和需要用户确认的事项；不要把中间命令日志整段贴给用户。
