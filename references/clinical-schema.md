# Clinical / prognosis schema

Use this schema when the user says “也帮我整理预后/临床数据” or supplies an existing clinical table. The validator is `scripts/validate_clinical.R`.

## Principles

1. **One row = one patient** by default. If a patient has multiple samples, add `sample_id` and run with `--allow-multiple-samples`; otherwise duplicate `patient_id` is an error.
2. **Do not invent missing values.** Blank means unknown. Never encode unknown as `0`, and never fill survival events from text.
3. **Time units are part of the column name.** Use `OS_time_months` or `OS_time_days`; never a bare `OS_time` in the final output. Input bare `OS_time` is accepted only with `--time-unit months|days`.
4. **Events are 0/1.** `1` = event (death, progression, recurrence); `0` = censored. Do not use `dead/alive`, `yes/no`, or `1/2`.
5. **Keep patient identifiers pseudonymous.** Do not put names, medical record numbers, or addresses in the table.
6. **Record the source and extraction date** in the final README, not as patient rows.

## Required identifiers

| Column | Required | Meaning |
|---|---|---|
| `patient_id` | yes | Stable pseudonymous patient identifier; unique unless multiple samples per patient |
| `sample_id` | recommended | Expression data.frame column name; required for multiple samples per patient and for sample-level matching |

## Standard clinical fields

| Column | Type | Notes |
|---|---|---|
| `age_years` | numeric | 0–120; `age` is also accepted and checked |
| `sex` | text | e.g. `female`, `male`, `unknown`; do not use numeric codes without a legend |
| `histology` | text | e.g. `adenocarcinoma` |
| `stage` | text | Include staging system in the README, e.g. `AJCC 8th`, `FIGO 2018` |
| `grade` | text | e.g. `G2` |
| `tissue_type` | text | `tumor`, `normal`, `adjacent_normal`, `metastasis`, `benign`, etc. |
| `tumor_purity` | numeric | 0–1; optional |
| `treatment` | text | Surgery/chemotherapy/radiation/immunotherapy; define combinations |
| `response` | text | RECIST/other criteria must be stated in README |
| `batch` | text | Sequencing/array batch or site; useful for confounding checks |
| `follow_up_months` | numeric | Optional overall follow-up |

Do not add redundant encodings such as `sex_male` unless the user's downstream model needs them. The validator preserves extra columns.

## Survival endpoints

| Endpoint | Definition | Event = 1 |
|---|---|---|
| `OS` | Overall survival from diagnosis/randomization to death | Death from any cause |
| `PFS` | Progression-free survival from treatment start to progression or death | Progression or death |
| `DFS` | Disease-free survival from resection/remission to recurrence or death | Recurrence or death |
| `DSS` | Disease-specific survival to death from the cancer | Cancer-specific death |

For each endpoint:

```text
<endpoint>_time_months   numeric, > 0 (0 allowed but flagged)
<endpoint>_event         0 or 1
```

or `_time_days`. The suffix must match the actual unit. If the source reports weeks, convert to months/days and document the conversion; do not leave a bare `_time`.

### Censoring and time origin

- State the time origin in the README: diagnosis date, surgery date, treatment start, etc.
- Censored patients need a non-missing follow-up time; `event=0, time=NA` is flagged.
- `event=1, time=NA` is an error.
- Zero follow-up is allowed by the validator but flagged; decide whether to exclude and document.
- Do not mix prevalent and incident cases without considering immortal time bias.

## Template: broad clinical table

`assets/clinical_template.csv`

```csv
patient_id,sample_id,age_years,sex,histology,stage,grade,tissue_type,tumor_purity,treatment,response,batch,OS_time_months,OS_event,PFS_time_months,PFS_event,DFS_time_months,DFS_event,DSS_time_months,DSS_event,notes
P001,S1,58,female,adenocarcinoma,IIA,G2,tumor,0.8,surgery,responder,batch1,36.5,1,30.0,1,30.0,1,36.5,1,example
P002,S2,63,male,adenocarcinoma,IIIB,G3,tumor,0.6,surgery+chemo,non-responder,batch1,18.0,0,12.0,1,NA,NA,NA,NA,example
```

## Template: survival-only table

`assets/prognosis_template.csv`

```csv
patient_id,sample_id,OS_time_months,OS_event,PFS_time_months,PFS_event,DFS_time_months,DFS_event,DSS_time_months,DSS_event,follow_up_months
P001,S1,36.5,1,30.0,1,30.0,1,36.5,1,36.5
P002,S2,18.0,0,12.0,1,NA,NA,NA,NA,18.0
```

## Template: sample mapping

`assets/samples_template.csv`

```csv
sample_id,patient_id,tissue_type,tumor_type,batch,group,include_in_analysis,notes
S1,P001,tumor,adenocarcinoma,batch1,T,TRUE,primary
S2,P002,tumor,adenocarcinoma,batch1,T,TRUE,primary
```

Use this when one patient contributes multiple samples (primary/metastasis/normal) or when the expression matrix's sample IDs need a separate clinical key.

## Validation behavior

`validate_clinical.R` checks:

- `patient_id` exists, non-empty, unique (unless `--allow-multiple-samples`).
- Survival columns are paired; time numeric and non-negative; events 0/1.
- `event=1` with missing time is an error; `event=0` with missing time is a warning.
- If `--rdata` is supplied, `sample_id` values are matched to `colnames(rt)`; unmatched samples are reported in `sample_map.tsv`.
- `age`/`age_years` outside 0–120 is a warning.
- Missing values are allowed and preserved as blank/`NA`.

Outputs:

```text
clinical/clinical.RData        clinical, clinical_info, optional sample_map
clinical/clinical_clean.csv    standardized columns + original extras
clinical/clinical_validation.tsv
clinical/sample_map.tsv        when --rdata is supplied
clinical/survival_analysis.csv when a standard endpoint is present
```

## Common pitfalls

- **Units:** months vs days is the most common silent error. Always inspect the source's unit and put it in the column name.
- **Event coding:** some tables use 1=censored, 2=dead. The validator rejects non-0/1 values; recode after checking the source legend.
- **Multiple samples per patient:** do not duplicate the same survival outcome across samples without accounting for within-patient correlation.
- **Stage/grade:** keep the staging system and version in the README; do not merge AJCC 7th and 8th without conversion.
- **Treatment:** “surgery+chemo” is not the same as “surgery” vs “chemo”; define categories before modeling.
- **Missing not at random:** survival time missing more often in sick patients can bias analyses; report missingness.
- **Immortal time:** do not start follow-up at a later treatment landmark unless the design explicitly accounts for it.
