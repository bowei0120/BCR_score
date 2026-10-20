# BCR_score

基于 bulk RNA-seq 免疫组库（BCR 轻链 IGL/IGK）的疾病预测流水线。
根据样本克隆型中的 CDR3 氨基酸序列，利用预训练的 Beta-Binomial 贝叶斯模型，
计算样本对溃疡性结肠炎（UC）和克罗恩病（CD）的倾向分数（F scale）。

## 目录结构

```
BCR_score/
├── combine_result.py        # 合并 MIXCR 与 TRUST4 的克隆型结果
├── predict.py               # 基于克隆型表计算 UC/CD 疾病 F 分数
├── fisher_sets.json         # 每链 universe/UC/CD 的 CDR3aa 参考集合（预计算）
├── params_full.json         # 每链每个模型的 Beta-Binomial 参数 (a, b, N)
├── bcr_core.py              # 训练库：队列读入、Fisher 列联表、覆盖度、BB 拟合
├── feature_selection.py     # 训练：Fisher 精确检验 CDR3 特征选择 -> fisher_sets.json
├── estimate_params.py       # 训练：覆盖度 + Beta-Binomial 拟合 -> params.json / scales.tsv
├── cross_validate.py        # 训练：无泄漏留一患者交叉验证 -> AUC/ROC
├── generate_random_labels.py # 训练：生成置换标签队列表（置换检验的第一步）
├── resample_patients.py     # 训练：每患者随机抽取一个标本，重复 N 次
└── example/
    ├── mixcr/               # MIXCR 输出示例（SRR18509019）
    ├── trust4/              # TRUST4 输出示例（SRR18509019）
    └── SRR18509019/         # 流水线输出示例（合并克隆型表 + 指标 + 预测结果）
```

## 依赖

- Python 3
- pandas
- numpy
- scipy

## 输入文件生成（MIXCR / TRUST4）

`combine_result.py` 的输入是两个免疫组库分析工具对样本 FASTQ 的输出。
以双端数据为例（`${reads[0]}` / `${reads[1]}` 为该样本的双端 FASTQ，`${prefix}` 为样本前缀，
`${task.cpus}` 为线程数，以下为 nextflow 流程中的写法）：

### MIXCR（克隆型分析）

```bash
mixcr analyze shotgun \
    -s hsa \
    --starting-material rna \
    --align "-OallowPartialAlignments=true" \
    ${reads[0]} \
    ${reads[1]} \
    mixcr.${prefix}
```

参数说明：

- `-s hsa`：物种（human）
- `--starting-material rna`：起始材料为 RNA
- `--align "-OallowPartialAlignments=true"`：允许部分比对

输出文件 `mixcr.${prefix}.clonotypes.ALL.txt`，即 `combine_result.py` 的 `-m` 输入。

### TRUST4（组装与注释）

```bash
run-trust4 \
    -t ${task.cpus} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -f /opt2/TRUST4/hg38_bcrtcr.fa \
    --ref /opt2/TRUST4/human_IMGT+C.fa \
    -o trust4.${prefix}
```

参数说明：

- `-t`：线程数
- `-1` / `-2`：双端 FASTQ（单端数据用 `-u <fastq>`）
- `-f`：参考基因组序列（hg38 的 BCR/TCR 位点）
- `--ref`：IMGT 基因参考
- `-o`：输出前缀

输出文件 `trust4.${prefix}_report.tsv`，即 `combine_result.py` 的 `-t` 输入
（示例目录中的文件名为 `trsut4.SRR18509019_report.tsv`，前缀 `trsut4` 为流水线输出时的拼写）。

## 使用流程

### 1. 合并 MIXCR 与 TRUST4 结果

```bash
python3 combine_result.py \
    -m <mixcr_clonotypes.txt> \
    -t <trust4_report.tsv> \
    -i <sample_id> \
    -o <output_dir>
```

输出：

- `<sample_id>.BCR_clonotypes.txt`：以 CDR3nt+CDR3aa 为键合并后的克隆型表
- `<sample_id>_metrics_manifest.yaml`：各基因类型（IGH/IGK/IGL/TRA/TRB/TRD/TRG）的克隆型计数

### 2. 计算疾病 F 分数

```bash
python3 predict.py --table <samples.tsv> --out <result.tsv> \
    --fisher-json fisher_sets.json --params params_full.json
```

`--fisher-json` 和 `--params` 为必填参数，无内置默认值，必须分别指向预计算的
`fisher_sets.json` 和 `params_full.json`。

`--table` 为 TSV，需包含两列：

- `sample`：样本名
- `clono_path`：该样本的克隆型表路径（即合并步骤的 `<sample_id>.BCR_clonotypes.txt`）

注意：表头行不能以 `#` 开头（脚本会将其按注释行跳过）。

输出列说明：

| 列 | 含义 |
|---|---|
| `scale_UC` / `scale_CD` | UC / CD 总 F 分数（IGL + IGK 之和） |
| `scale_IGL_UC` / `scale_IGL_CD` / `scale_IGK_UC` / `scale_IGK_CD` | 各链的 F 分数 |
| `n_IGL` / `n_IGK` | 样本中落入 universe 参考集的 CDR3aa 数量 |
| `k_IGL_UC` / `k_IGL_CD` / `k_IGK_UC` / `k_IGK_CD` | 落入 UC / CD 参考集的 CDR3aa 数量 |

## 方法

对每个样本、每条链（IGL / IGK），统计其 CDR3aa 与参考集的命中数：

- `n`：与 universe 集（训练队列中该链观测到的全部 CDR3aa，任一标签的样本携带即计入）的交集数量
- `k_UC` / `k_CD`：与 UC / CD 关联 CDR3aa 集合的交集数量

每个模型用 Beta-Binomial 分布描述命中概率，`k ~ Beta-Binomial(n, a, b)`，
`N` 为训练该模型的样本数。注意 `n`/`k` 由样本的全部 CDR3aa 计算（不按链过滤），
与 `predict.py` 完全一致。F 分数为 log Bayes factor：

```
log P(数据 | 疾病模型) / P(数据 | 健康模型)
```

## 训练（模型复现）

以下脚本用于从带标签的队列重建 `fisher_sets.json` 和 `params_full.json`，
均接受一个四列的**队列 TSV**：

| 列 | 含义 |
|---|---|
| `sample` | 样本 id（每行一个克隆型文件） |
| `patient` | 患者 id（同一患者可有多个标本） |
| `label` | `HT`、`UC` 或 `CD` |
| `clono_path` | 该样本克隆型表路径（合并步骤产出的 `<sample_id>.BCR_clonotypes.txt`） |

注意：表头行不能以 `#` 开头（脚本会按注释行跳过）。

### 1. Fisher 特征选择

```bash
python3 feature_selection.py --table cohort.tsv --out-dir train/ \
    [--chains IGL IGK] \
    [--p 0.01] [--or-thresh 1.5] [--min-count 5]
```

对每条链，统计携带每个 CDR3aa 的样本在各标签下的数量，并对每个
CDR3aa 做双侧 Fisher 精确检验，OR 为疾病朝向（OR > 1 = 在疾病中富集）。

- 计数单位为样本（sample-level）——患者内非独立性由 `resample_patients.py`
  在上游解决（每患者重采样为一个标本），因此置换检验、重采样、正常分析
  全部使用同一套样本级分析流程。
- 输出 `IGL_fisher.tsv`、`IGK_fisher.tsv`（完整列联表）和 `fisher_sets.json`，
  schema 与发布的参考文件完全一致（可直接传给 `predict.py --fisher-json`）。

### 2. 参数估计

```bash
python3 estimate_params.py --table cohort.tsv \
    --fisher-json train/fisher_sets.json --out-dir train/
```

- `coverage.tsv`：每样本每链，`n` = |样本全部 CDR3aa ∩ universe 集|，
  `k_UC` / `k_CD` = |样本全部 CDR3aa ∩ UC/CD 集|（样本侧不按链过滤，
  与 `predict.py` 一致）。
- `params.json`：每链 × {HT_UC, UC, HT_CD, CD} 的 Beta-Binomial 参数，
  schema 与 `params_full.json` 相同，可通过 `predict.py --params` 使用。
- `scales.tsv`：训练队列的样本内 F 分数。

### 3. 交叉验证

```bash
python3 cross_validate.py --table cohort.tsv --out-dir train/cv/ \
    [--p 0.01] [--or-thresh 1.5] [--min-count 5]
```

无泄漏的留一样本交叉验证：对每个留出样本，Fisher 特征选择与 Beta-Binomial
参数都仅在其余样本上重新拟合，再给留出样本打分。每个留出样本训练一个模型
（完整留一样本），因此在大队列上是最慢的一步。
输出 `loo_scores.tsv`（每样本一行的 out-of-sample 分数 + `held_out`）、
`auc.tsv`（UC vs HT、CD vs HT）、`roc.tsv`（合并 ROC 曲线点）。

### 4. 置换检验（经验 FDR，补充图 S1）

置换检验的下游分析与正常分析完全一致，唯一区别是起始的标签随机化——因此
仓库只提供标签生成脚本；置换表生成后，直接对它们运行正常流程。

```bash
python3 generate_random_labels.py --table cohort.tsv --out-base-dir train/perm/ \
    [--n-files 200] [--seed 42]
```

在 `--out-base-dir` 下写出 `random_CD/` 和 `random_UC/`，各含 `--n-files` 个
队列表（`CD_001.txt` … / `UC_001.txt` …）。仅在 domain（`clono_path` 非空且
文件在磁盘上存在）内打乱标签，病例:对照比例严格保持不变；非 domain 行保持
原标签；输出与输入表列名/列序相同。固定 `--seed` 可复现。

对每个置换队列表运行与正常分析完全相同的 `feature_selection.py`（比较各链/
各组显著 CDR3 数 → 经验 FDR，即补充图 S1 的分析）或 `cross_validate.py`
（比较 out-of-sample AUC），即得零假设分布。

### 5. 患者重采样

```bash
python3 resample_patients.py --table cohort.tsv --out-dir train/res/ \
    [--n-iter 1000] [--seed 0]
```

将有多个标本的患者下采样为恰好一个均匀随机标本，重复 `--n-iter` 次
（第 i 次迭代使用种子 `--seed`+i，结果确定且各次迭代相互独立）。
输出与队列同 schema 的 `iter_XXX.tsv`，可作为上述脚本的 `--table` 输入——
解决患者内非独立性问题。重采样只重新生成队列表（`iter_XXX.tsv`），其后的
分析与正常分析完全一致——与置换检验同理，只有输入表变化。

公共逻辑位于 `bcr_core.py`；所有脚本仅依赖 pandas / numpy / scipy。

## 示例：SRR18509019

```bash
# 合并（结果见 example/SRR18509019/）
python3 combine_result.py \
    -m example/mixcr/mixcr.SRR18509019.clonotypes.ALL.txt \
    -t example/trust4/trsut4.SRR18509019_report.tsv \
    -i SRR18509019 \
    -o example/SRR18509019

# 打分
printf 'sample\tclono_path\nSRR18509019\t%s\n' \
    "$PWD/example/SRR18509019/SRR18509019.BCR_clonotypes.txt" > samples.tsv
python3 predict.py --table samples.tsv --out result.tsv \
    --fisher-json fisher_sets.json --params params_full.json
```
