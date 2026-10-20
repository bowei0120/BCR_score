# BCR_score

Disease prediction pipeline based on bulk RNA-seq immune repertoires (BCR light chains IGL/IGK).
Using the CDR3 amino acid sequences in sample clonotypes and a pre-trained Beta-Binomial
Bayesian model, it computes a sample's propensity scores (BCR score) for ulcerative colitis (UC)
and Crohn's disease (CD).

## Directory structure

```
BCR_score/
├── combine_result.py        # Merge MIXCR and TRUST4 clonotype results
├── predict.py               # Compute UC/CD disease BCR scores from a clonotype table
├── fisher_sets.json         # Precomputed per-chain universe/UC/CD CDR3aa reference sets
├── params_full.json         # Pre-trained Beta-Binomial parameters (a, b, N) per model
├── bcr_core.py              # Training library: cohort IO, Fisher tables, coverage, BB fitting
├── feature_selection.py     # Train: Fisher exact CDR3 selection -> fisher_sets.json
├── estimate_params.py       # Train: coverage + Beta-Binomial fit -> params.json / scales.tsv
├── cross_validate.py        # Train: leakage-free leave-one-patient-out CV -> AUC/ROC
├── generate_random_labels.py # Train: generate permuted-label cohort tables (step 1 of the permutation test)
├── resample_patients.py     # Train: one random specimen per patient, N iterations
└── example/
    ├── mixcr/               # MIXCR output example (SRR18509019)
    ├── trust4/              # TRUST4 output example (SRR18509019)
    └── SRR18509019/         # Pipeline output example (combined clonotypes + metrics + prediction)
```

## Dependencies

- Python 3
- pandas
- numpy
- scipy

## Input file generation (MIXCR / TRUST4)

The inputs of `combine_result.py` are the outputs of two immune-repertoire analysis tools run
on a sample's FASTQ reads. Paired-end example below (`${reads[0]}` / `${reads[1]}` are the
sample's paired-end FASTQ files, `${prefix}` is the sample prefix, `${task.cpus}` is the number
of threads; written as in a Nextflow pipeline):

### MIXCR (clonotype analysis)

```bash
mixcr analyze shotgun \
    -s hsa \
    --starting-material rna \
    --align "-OallowPartialAlignments=true" \
    ${reads[0]} \
    ${reads[1]} \
    mixcr.${prefix}
```

Parameter notes:

- `-s hsa`: species (human)
- `--starting-material rna`: starting material is RNA
- `--align "-OallowPartialAlignments=true"`: allow partial alignments

Output file `mixcr.${prefix}.clonotypes.ALL.txt`, which is the `-m` input of `combine_result.py`.

### TRUST4 (assembly and annotation)

```bash
run-trust4 \
    -t ${task.cpus} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -f /opt2/TRUST4/hg38_bcrtcr.fa \
    --ref /opt2/TRUST4/human_IMGT+C.fa \
    -o trust4.${prefix}
```

Parameter notes:

- `-t`: number of threads
- `-1` / `-2`: paired-end FASTQ files (use `-u <fastq>` for single-end data)
- `-f`: reference genome sequences (hg38 BCR/TCR loci)
- `--ref`: IMGT gene reference
- `-o`: output prefix

Output file `trust4.${prefix}_report.tsv`, which is the `-t` input of `combine_result.py`
(the file in the example directory is named `trsut4.SRR18509019_report.tsv`; the `trsut4`
prefix is the spelling produced by the pipeline).

## Usage

### 1. Merge MIXCR and TRUST4 results

```bash
python3 combine_result.py \
    -m <mixcr_clonotypes.txt> \
    -t <trust4_report.tsv> \
    -i <sample_id> \
    -o <output_dir>
```

Outputs:

- `<sample_id>.BCR_clonotypes.txt`: clonotype table merged by CDR3nt + CDR3aa
- `<sample_id>_metrics_manifest.yaml`: clonotype counts per gene type (IGH/IGK/IGL/TRA/TRB/TRD/TRG)

### 2. Compute disease BCR scores

```bash
python3 predict.py --table <samples.tsv> --out <result.tsv> \
    --fisher-json fisher_sets.json --params params_full.json
```

`--fisher-json` and `--params` are REQUIRED, with no built-in defaults; they must
point to the precomputed `fisher_sets.json` and `params_full.json`.

`--table` is a TSV with two columns:

- `sample`: sample name
- `clono_path`: path to the sample's clonotype table (i.e. the `<sample_id>.BCR_clonotypes.txt`
  produced by the merge step)

Note: the header line must NOT start with `#` (the script treats it as a comment line and skips it).

Output columns:

| Column | Meaning |
|---|---|
| `scale_UC` / `scale_CD` | Total UC / CD BCR scores (IGL + IGK) |
| `scale_IGL_UC` / `scale_IGL_CD` / `scale_IGK_UC` / `scale_IGK_CD` | Per-chain BCR scores |
| `n_IGL` / `n_IGK` | Number of sample CDR3aa found in the universe reference set |
| `k_IGL_UC` / `k_IGL_CD` / `k_IGK_UC` / `k_IGK_CD` | Number of sample CDR3aa found in the UC / CD reference sets |

## Method

For each sample and each chain (IGL / IGK), count the hits of its CDR3aa against the reference sets:

- `n`: size of the intersection with the universe set (all CDR3aa observed for that chain in the training cohort, counted when carried by any sample regardless of label)
- `k_UC` / `k_CD`: size of the intersection with the UC / CD-associated CDR3aa sets

**Important:** `n` and `k` are computed from **ALL unique CDR3aa of the sample** (chain-unfiltered), by intersecting that full sample CDR3aa set with each chain's reference sets (universe/UC/CD). This matches `predict.py` and the shared path in `bcr_core.py`.

Each model describes the hit probability with a Beta-Binomial distribution,
`k ~ Beta-Binomial(n, a, b)`, where `N` is the number of samples used to train the model.
The BCR score is the log Bayes factor:

```
log P(data | disease model) / P(data | healthy model)
```

To account for within-patient dependence (patients may contribute multiple specimens), training uses **one sample per patient** (via `resample_patients.py`) and evaluation uses **leave-one-patient-out cross-validation (LOPOCV)** with patient-wise grouping (via `cross_validate.py`).

## Training (reproducing the model)

The scripts below rebuild `fisher_sets.json` and `params_full.json` from a labeled
cohort. They all take a **cohort TSV** with four columns:

| Column | Meaning |
|---|---|
| `sample` | sample id (one clonotype file per row) |
| `patient` | patient id (a patient may contribute several specimens) |
| `label` | `HT`, `UC`, or `CD` |
| `clono_path` | path to the sample's clonotype table (the `<sample_id>.BCR_clonotypes.txt` from the merge step) |

Note: the header line must NOT start with `#` (the scripts skip it as a comment line).

### Handling within-patient non-independence

Since a single patient may contribute multiple biopsy specimens/samples, samples from the same patient are not independent. The pipeline explicitly accounts for this at different stages:

- **Downsampling to one sample per patient**: `resample_patients.py` randomly selects **exactly one** specimen per `patient` to form a "one-sample-per-patient" dataset. This is used to avoid over-representation of patients and to reduce within-patient correlation in feature selection and permutation/resampling analyses.
- **Patient-wise cross-validation**: `cross_validate.py` performs **Leave-One-Patient-Out Cross-Validation (LOPOCV)**. In each fold, **all samples belonging to the held-out patient are removed together** from feature selection and model fitting, and only the held-out patient's samples are scored with that fold's model.
- **Fisher counting**: `bcr_core.fisher_tables` counts carriers at the **sample level**, but is designed to be run on tables where patient-level dependence has already been addressed (e.g. the one-sample-per-patient tables produced by `resample_patients.py`).

In practice, the recommended workflow is: first generate one or more resampled cohorts with `resample_patients.py` (each with 1 sample/patient), then run `feature_selection.py`, `estimate_params.py`, or `cross_validate.py` on those resampled cohorts.

### 0. (Recommended) Resample to one specimen per patient

Before feature selection, parameter estimation, or cross-validation, it is strongly recommended to downsample to avoid within-patient leakage/duplication:

```bash
python3 resample_patients.py   --table cohort.tsv   --out train/resampled/   --seed 12345   --iter 100
```

This will generate `iter_001.tsv`, `iter_002.tsv`, ..., `iter_100.tsv` under `train/resampled/`, where each file contains **exactly one randomly selected sample per patient** (same `patient` set, possibly different sample choices across iterations). Use these resampled tables for training, permutation analyses, and LOPOCV.

### 1. Fisher feature selection

```bash
python3 feature_selection.py --table cohort.tsv --out-dir train/ \
    [--chains IGL IGK] \
    [--p 0.01] [--or-thresh 1.5] [--min-count 5]
```

For each chain, counts how many samples carrying each CDR3aa belong
to each label and runs a two-sided Fisher exact test per CDR3aa, with a
disease-oriented odds ratio (OR > 1 = enriched in disease).

- Counting unit is the sample: within-patient non-independence is handled
  upstream by `resample_patients.py`, which downsamples each patient to one
  specimen, so permutation testing, resampling, and normal analysis all share
  one sample-level pipeline.
- Outputs `IGL_fisher.tsv`, `IGK_fisher.tsv` (full tables) and `fisher_sets.json`
  in the exact schema of the released reference file (passable to
  `predict.py --fisher-json`).

### 2. Parameter estimation

```bash
python3 estimate_params.py --table cohort.tsv \
    --fisher-json train/fisher_sets.json --out-dir train/
```

- `coverage.tsv`: per sample and chain, `n` = |all sample CDR3aa ∩ universe set|,
  `k_UC` / `k_CD` = |all sample CDR3aa ∩ UC/CD sets| (sample side is
  chain-unfiltered, matching `predict.py`).
- `params.json`: Beta-Binomial parameters per chain × key in
  {HT_UC, UC, HT_CD, CD} — same schema as `params_full.json`, so usable via
  `predict.py --params`.
- `scales.tsv`: in-sample BCR scores for the training cohort.

### 3. Cross-validation (Leave-One-Patient-Out)

```bash
# After creating one-sample-per-patient cohorts
python3 cross_validate.py \
  --table train/resampled/iter_001.tsv \
  --out-dir train/cv/iter_001/ \
  [--chains IGL IGK] \
  [--p 0.01] [--or-thresh 1.5] [--min-count 5]
```

Implementation details (matching the manuscript):

1. **Patient grouping**: `cross_validate.py` groups samples by `patient`. For each distinct patient, it performs a fold in which **all samples of that patient are held out as the test set** (leave-one-patient-out).
2. **Leakage-free training per fold**: Within each fold, Fisher feature selection, reference set construction (`fisher_sets.json` equivalent), and Beta-Binomial parameter fitting are recomputed using **only the training patients**. No information from the held-out patient is used for feature selection or parameter estimation.
3. **Scoring**: The fold-specific model is applied to the held-out patient's samples to obtain out-of-sample BCR scores, which are then aggregated across all folds.
4. **With one-sample-per-patient tables**: If the input is produced by `resample_patients.py` (exactly one sample per patient), LOPOCV reduces to effectively "leave-one-patient = leave-one-sample" out, but the grouping is still enforced by `patient` to be robust to future cases with multiple specimens.

Outputs (per fold/summary) include out-of-sample scores and cross-validated AUCs (consistent with the manuscript's LOPOCV evaluation).
```

Leakage-free leave-one-sample-out CV: for every held-out sample, both the Fisher
feature selection and the Beta-Binomial parameters are refit on the remaining
samples ONLY, then the held-out sample is scored. One model is trained per
held-out sample, so on large cohorts this is the slowest step.
Outputs `loo_scores.tsv` (one out-of-sample score row per sample + `held_out`),
`auc.tsv` (UC vs HT, CD vs HT), `roc.tsv` (pooled ROC curve points).

### 4. Permutation test (empirical FDR, Supp. Fig. S1)

The permutation test's downstream analysis is identical to the normal analysis; the
ONLY difference is the initial label randomization. The repo therefore ships only the
label generator; once permuted tables exist, run the normal pipeline on them.

```bash
python3 generate_random_labels.py --table cohort.tsv --out-base-dir train/perm/ \
    [--n-files 200] [--seed 42]
```

Writes `random_CD/` and `random_UC/` under `--out-base-dir`, each containing
`--n-files` cohort tables (`CD_001.txt` ... / `UC_001.txt` ...). Labels are shuffled
only within the domain (rows with non-empty `clono_path` whose file exists on disk),
so the case:control ratio of domain rows is preserved exactly; non-domain rows keep
their original labels; the output keeps the same columns/order as the input table.
A fixed `--seed` makes the run reproducible.

For the null distribution, run the exact same `feature_selection.py` on every
permuted cohort table (compare significant CDR3 counts per chain / per group ->
empirical FDR, the analysis of Supp. Fig. S1) or `cross_validate.py` (compare
out-of-sample AUC).

### 5. Patient resampling

```bash
python3 resample_patients.py --table cohort.tsv --out-dir train/res/ \
    [--n-iter 1000] [--seed 0]
```

Patients with multiple specimens are downsampled to exactly ONE uniformly random
specimen per patient, `--n-iter` times (iteration i uses seed `--seed`+i, so runs
are deterministic). Writes `iter_XXX.tsv` files with the same cohort schema,
reusable as `--table` inputs for the scripts above — addressing the
within-patient non-independence. Resampling only regenerates the cohort table
(`iter_XXX.tsv`); the downstream analysis is exactly the normal pipeline: same
principle as the permutation test, only the input table changes.

Shared logic lives in `bcr_core.py`; all scripts need only pandas / numpy / scipy.

## Example: SRR18509019

```bash
# Merge (results in example/SRR18509019/)
python3 combine_result.py \
    -m example/mixcr/mixcr.SRR18509019.clonotypes.ALL.txt \
    -t example/trust4/trsut4.SRR18509019_report.tsv \
    -i SRR18509019 \
    -o example/SRR18509019

# Score
printf 'sample\tclono_path\nSRR18509019\t%s\n' \
    "$PWD/example/SRR18509019/SRR18509019.BCR_clonotypes.txt" > samples.tsv
python3 predict.py --table samples.tsv --out result.tsv \
    --fisher-json fisher_sets.json --params params_full.json
```
