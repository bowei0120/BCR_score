# BCR_score

基于 bulk RNA-seq 免疫组库（BCR 轻链 IGL/IGK）的疾病预测流水线。
根据样本克隆型中的 CDR3 氨基酸序列，利用预训练的 Beta-Binomial 贝叶斯模型，
计算样本对溃疡性结肠炎（UC）和克罗恩病（CD）的倾向分数（BCR score）。

## 目录结构

```
BCR_score/
├── combine_result.py        # 合并 MIXCR 与 TRUST4 的克隆型结果
├── predict.py               # 基于克隆型表计算 UC/CD 疾病 BCR score
├── fisher_sets.json         # 每链 universe/UC/CD 的 CDR3aa 参考集合（预计算）
├── params_full.json         # 每链每个模型的 Beta-Binomial 参数 (a, b, N)
├── bcr_core.py              # 训练库：队列读入、Fisher 列联表、覆盖度、BB 拟合
├── feature_selection.py     # 训练：Fisher 精确检验 CDR3 特征选择 -> fisher_sets.json
├── estimate_params.py       # 训练：覆盖度 + Beta-Binomial 拟合 -> params.json / scores.tsv
├── cross_validate.py        # 训练：无泄漏留一患者交叉验证 -> AUC/ROC
├── generate_random_labels.py # 训练：生成置换标签队列表（置换检验的第一步）
├── resample_patients.py     # 训练：每患者随机抽取一个标本，重复 N 次
├── ROC.R                    # 评估：ROC 分析（AUC/CI、Youden 最佳 cutoff），参数化脚本
├── BCR_score_result/        # 各队列的 BCR score 结果表（ROC.R 的输入）
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
- R（仅 `ROC.R` 需要）+ pROC 包（`install.packages("pROC")`）

## 输入文件生成（MIXCR / TRUST4）

`combine_result.py` 的输入是两个免疫组库分析工具对样本 FASTQ 的输出。
以双端数据为例（`${reads[0]}` / `${reads[1]}` 为该样本的双端 FASTQ，`${prefix}` 为样本前缀：

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

### 2. 计算疾病 BCR score

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
| `score_UC` / `score_CD` | UC / CD 总 BCR score（IGL + IGK 之和） |
| `score_IGL_UC` / `score_IGL_CD` / `score_IGK_UC` / `score_IGK_CD` | 各链的 BCR score |
| `n_IGL` / `n_IGK` | 样本中落入 universe 参考集的 CDR3aa 数量 |
| `k_IGL_UC` / `k_IGL_CD` / `k_IGK_UC` / `k_IGK_CD` | 落入 UC / CD 参考集的 CDR3aa 数量 |

### 3. ROC 分析

`BCR_score_result/` 存放各队列的打分结果表（`score_<队列>.csv`，列为
`label`、`score_UC`、`score_CD`，行名为样本 id），由 `ROC.R` 对其做 ROC 分析。
`ROC.R` 为参数化脚本，全部通过命令行传参：

```bash
Rscript ROC.R --table BCR_score_result/score_PRJNA797175.csv \
    --case UC --control Control --score score_UC \
    [--label label] [--cutoff 3.258] [--out roc.pdf]
```

| 参数 | 说明 |
|---|---|
| `--table` | **必填**：分数表路径（CSV 按 `row.names=1` 读取；TSV 需表头 + 行名列） |
| `--case` | **必填**：标签列中的病例标签（如 `UC`、`CD`） |
| `--control` | **必填**：标签列中的对照标签（如 `Control`、`HT`） |
| `--score` | **必填**：使用的分数列名（如 `score_UC`、`score_CD`） |
| `--label` | 标签列名，默认 `label` |
| `--cutoff` | 可选：固定阈值，输出该点的敏感度/特异度 |
| `--out` | 可选：输出图片路径，按扩展名选设备（`.pdf`/`.png`/`.jpg`/`.svg`）；省略则交互式绘图 |

输出到标准输出：AUC（95% CI）、Youden 最佳 cutoff 及其 sens/spec，
若指定 `--cutoff` 则追加该阈值下的 sens/spec。

示例：

```bash
# UC vs Control，固定阈值 3.258
Rscript ROC.R --table BCR_score_result/score_PRJNA797175.csv \
    --case UC --control Control --score score_UC --cutoff 3.258

# CD vs Control，并输出 PDF 图
Rscript ROC.R --table BCR_score_result/score_PRJEB24645.csv \
    --case CD --control Control --score score_CD --out roc_CD.pdf
```

## 方法

对每个样本、每条链（IGL / IGK），统计其 CDR3aa 与参考集的命中数：

- `n`：与 universe 集（训练队列中该链观测到的全部 CDR3aa，任一标签的样本携带即计入）的交集数量
- `k_UC` / `k_CD`：与 UC / CD 关联 CDR3aa 集合的交集数量

**重要说明：** `n` 和 `k` 是从**样本的全部唯一 CDR3aa**（不按链过滤）计算的，即将样本全部的 CDR3aa 集合分别与各链的参考集（universe/UC/CD）做交集来获得。这一做法与 `predict.py` 以及 `bcr_core.py` 的实现保持一致。

每个模型用 Beta-Binomial 分布描述命中概率，`k ~ Beta-Binomial(n, a, b)`，
`N` 为训练该模型的样本数。BCR score 为 log Bayes factor：

```
log P(数据 | 疾病模型) / P(数据 | 健康模型)
```

为了解决同一患者可能贡献多个标本（样本间非独立）的问题，本流程在训练时采用**每患者只保留一个样本**（使用 `resample_patients.py`），在评估时采用**按患者分组的留一患者交叉验证（Leave-One-Patient-Out CV，LOPOCV）**（使用 `cross_validate.py`）。

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

### 同一患者多个样本的处理方式

由于同一患者可能贡献多个活检标本/样本，这些样本之间并非独立。为此，训练流程在不同阶段做了相应处理：

- **每患者只保留一个样本（下采样）**：`resample_patients.py` 会为每个 `patient` 随机选取**恰好 1 个**样本，构建“每患者一个样本”（one-sample-per-patient）的数据集，用于降低同一患者样本之间的相关性，避免在特征选择、置换检验和重采样分析中出现患者过度代表的问题。
- **按患者分组的留一患者交叉验证（LOPOCV）**：`cross_validate.py` 执行**留一患者交叉验证（Leave-One-Patient-Out Cross-Validation, LOPOCV）**。在每一折中，**被留出的整个患者的所有样本一起从训练集中剔除**，仅用于评分，特征选择和 Beta-Binomial 参数拟合均仅基于其余患者的样本，不使用留出患者的任何信息。
- **Fisher 计数**：`bcr_core.fisher_tables` 在样本层面（sample-level）统计携带者数量，但推荐在**每患者一个样本**的数据集上运行，以更好地控制患者内相关性。

推荐的训练流程：先用 `resample_patients.py` 生成若干组“每患者一个样本”的队列表，再在这些重采样后的表上运行 `feature_selection.py`、`estimate_params.py` 或 `cross_validate.py`。

### 0. （推荐）每患者随机抽取一个样本

在进行特征选择、参数估计或交叉验证前，强烈建议先进行下采样以避免同患者样本造成的泄漏或重复偏倚：

```bash
python3 resample_patients.py --table cohort.tsv --out-dir train/resampled/ \
    [--n-iter 100] [--seed 12345]
```

该命令会在 `train/resampled/` 下生成 `iter_000.tsv`、`iter_001.tsv`、...、`iter_099.tsv`。每个文件都包含**每个患者恰好一个随机选取的样本**（患者集合相同，不同迭代的样本选择可能不同）。建议使用这些重采样后的表进行训练、置换检验和 LOPOCV。

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
- `scores.tsv`：训练队列的样本内 BCR score。

### 3. 交叉验证（留一患者交叉验证，LOPOCV）

```bash
# 基于“每患者一个样本”的重采样队列表进行交叉验证
python3 cross_validate.py \
  --table train/resampled/iter_001.tsv \
  --out-dir train/cv/iter_001/ \
  [--chains IGL IGK] \
  [--p 0.01] [--or-thresh 1.5] [--min-count 5]
```

实现要点（与论文方法一致）：

1. **按患者分组**：`cross_validate.py` 按 `patient` 字段分组。对每个患者执行一折，**将该患者的全部样本一起作为测试集**（留一患者）。
2. **各折严格无泄漏**：每一折内的 Fisher 特征选择、参考集构建（等价于 `fisher_sets.json`）以及 Beta-Binomial 参数拟合均**仅使用训练集患者**的数据完成，完全不使用留出患者的信息。
3. **独立评分**：用该折训练得到的模型对留出患者的样本进行评分，最后将所有折的结果汇总，计算交叉验证 AUC。
4. **与“每患者一个样本”配合**：如果输入已经是 `resample_patients.py` 产生的“每患者一个样本”表，则 LOPOCV 实际上等效于“留一患者=留一样本”，但脚本仍按 `patient` 分组以保持逻辑通用性。

输出包括各折的 out-of-sample 评分结果及交叉验证汇总的 AUC 等指标。

在大队列上这是最慢的一步（每个留出患者训练一个模型）。
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

此外，当预测结果中包含贡献了多个标本的患者时，样本级 ROC 会因同一患者被重复
计入而偏乐观。此时可通过患者重采样进一步评估 **patient-level 的 ROC**：对每个
迭代（每位患者恰好一个标本）独立计算 ROC/AUC，再汇总各次迭代的 AUC 分布，
得到不依赖于具体标本选择的患者级判别性能估计。

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
