# BCR_score

Disease prediction pipeline based on bulk RNA-seq immune repertoires (BCR light chains IGL/IGK).
Using the CDR3 amino acid sequences in sample clonotypes and a pre-trained Beta-Binomial
Bayesian model, it computes a sample's propensity scores (F scale) for ulcerative colitis (UC)
and Crohn's disease (CD).

## Directory structure

```
predict/
├── combine_result.py        # Merge MIXCR and TRUST4 clonotype results
├── predict.py               # Compute UC/CD disease F scales from a clonotype table
├── fisher_sets.json         # Precomputed per-chain cover/UC/CD CDR3aa reference sets
├── params_full.json         # Pre-trained Beta-Binomial parameters (a, b, N) per model
└── example/
    ├── mixcr/               # MIXCR output example (SRR18509019)
    ├── trust4/              # TRUST4 output example (SRR18509019)
    └── SRR18509019/         # Pipeline output example (combined clonotypes + metrics + prediction)
```

## Dependencies

- Python 3
- pandas
- scipy

## Input file generation (MIXCR / TRUST4)

The inputs of `combine_result.py` are the outputs of two immune-repertoire analysis tools run
on a sample's FASTQ reads. Paired-end example below (`${reads[0]}` / `${reads[1]}` are the
sample's paired-end FASTQ files, `${prefix}` is the sample prefix, `${task.cpus}` is the number
of threads; written as in a Nextflow pipeline):

### MIXCR (clonotype analysis)

```bash
mixcr analyze amplicon \
    -s hsa \
    --5-end v-primers \
    --3-end j-primers \
    --starting-material rna \
    --adapters no-adapters \
    --align "-OallowPartialAlignments=true" \
    ${reads[0]} \
    ${reads[1]} \
    mixcr.${prefix}
```

Parameter notes:

- `-s hsa`: species (human)
- `--5-end v-primers` / `--3-end j-primers`: amplicon mode with 5' V primers and 3' J primers
- `--starting-material rna`: starting material is RNA
- `--adapters no-adapters`: no adapter sequences
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

- `<sample_id>.TCR_clonotypes.txt`: clonotype table merged by CDR3nt + CDR3aa
- `<sample_id>_metrics_manifest.yaml`: clonotype counts per gene type (IGH/IGK/IGL/TRA/TRB/TRD/TRG)

### 2. Compute disease F scales

```bash
python3 predict.py --table <samples.tsv> --out <result.tsv>
```

Optional arguments: `--fisher-json fisher_sets.json`, `--params params_full.json`
(defaults to the files in this project; usually no need to specify).

`--table` is a TSV with two columns:

- `sample`: sample name
- `clono_path`: path to the sample's clonotype table (i.e. the `<sample_id>.TCR_clonotypes.txt`
  produced by the merge step)

Note: the header line must NOT start with `#` (the script treats it as a comment line and skips it).

Output columns:

| Column | Meaning |
|---|---|
| `scale_UC` / `scale_CD` | Total UC / CD F scales (IGL + IGK) |
| `scale_IGL_UC` / `scale_IGL_CD` / `scale_IGK_UC` / `scale_IGK_CD` | Per-chain F scales |
| `n_IGL` / `n_IGK` | Number of sample CDR3aa found in the cover reference set |
| `k_IGL_UC` / `k_IGL_CD` / `k_IGK_UC` / `k_IGK_CD` | Number of sample CDR3aa found in the UC / CD reference sets |

## Method

For each sample and each chain (IGL / IGK), count the hits of its CDR3aa against the reference sets:

- `n`: size of the intersection with the cover set (CDR3aa shared by healthy and disease samples in the training cohort)
- `k_UC` / `k_CD`: size of the intersection with the UC / CD-associated CDR3aa sets

Each model describes the hit probability with a Beta-Binomial distribution,
`k ~ Beta-Binomial(n, a, b)`, where `N` is the number of samples used to train the model.
The F scale is the log Bayes factor:

```
log P(data | disease model) / P(data | healthy model)
```

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
    "$PWD/example/SRR18509019/SRR18509019.TCR_clonotypes.txt" > samples.tsv
python3 predict.py --table samples.tsv --out result.tsv
```
