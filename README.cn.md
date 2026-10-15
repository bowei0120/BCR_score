# BCR_score

基于 bulk RNA-seq 免疫组库（BCR 轻链 IGL/IGK）的疾病预测流水线。
根据样本克隆型中的 CDR3 氨基酸序列，利用预训练的 Beta-Binomial 贝叶斯模型，
计算样本对溃疡性结肠炎（UC）和克罗恩病（CD）的倾向分数（F scale）。

## 目录结构

```
predict/
├── combine_result.py        # 合并 MIXCR 与 TRUST4 的克隆型结果
├── predict.py               # 基于克隆型表计算 UC/CD 疾病 F 分数
├── fisher_sets.json         # 每链 cover/UC/CD 的 CDR3aa 参考集合（预计算）
├── params_full.json         # 每链每个模型的 Beta-Binomial 参数 (a, b, N)
└── example/
    ├── mixcr/               # MIXCR 输出示例（SRR18509019）
    ├── trust4/              # TRUST4 输出示例（SRR18509019）
    └── SRR18509019/         # 流水线输出示例（合并克隆型表 + 指标 + 预测结果）
```

## 依赖

- Python 3
- pandas
- scipy

## 输入文件生成（MIXCR / TRUST4）

`combine_result.py` 的输入是两个免疫组库分析工具对样本 FASTQ 的输出。
以双端数据为例（`${reads[0]}` / `${reads[1]}` 为该样本的双端 FASTQ，`${prefix}` 为样本前缀，
`${task.cpus}` 为线程数，以下为 nextflow 流程中的写法）：

### MIXCR（克隆型分析）

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

参数说明：

- `-s hsa`：物种（human）
- `--5-end v-primers` / `--3-end j-primers`：5' 端 V 引物、3' 端 J 引物的扩增子（amplicon）模式
- `--starting-material rna`：起始材料为 RNA
- `--adapters no-adapters`：无接头序列
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

- `<sample_id>.TCR_clonotypes.txt`：以 CDR3nt+CDR3aa 为键合并后的克隆型表
- `<sample_id>_metrics_manifest.yaml`：各基因类型（IGH/IGK/IGL/TRA/TRB/TRD/TRG）的克隆型计数

### 2. 计算疾病 F 分数

```bash
python3 predict.py --table <samples.tsv> --out <result.tsv>
```

可选参数：`--fisher-json fisher_sets.json`、`--params params_full.json`
（默认即本项目内文件，一般无需指定）。

`--table` 为 TSV，需包含两列：

- `sample`：样本名
- `clono_path`：该样本的克隆型表路径（即合并步骤的 `<sample_id>.TCR_clonotypes.txt`）

注意：表头行不能以 `#` 开头（脚本会将其按注释行跳过）。

输出列说明：

| 列 | 含义 |
|---|---|
| `scale_UC` / `scale_CD` | UC / CD 总 F 分数（IGL + IGK 之和） |
| `scale_IGL_UC` / `scale_IGL_CD` / `scale_IGK_UC` / `scale_IGK_CD` | 各链的 F 分数 |
| `n_IGL` / `n_IGK` | 样本中落入 cover 参考集的 CDR3aa 数量 |
| `k_IGL_UC` / `k_IGL_CD` / `k_IGK_UC` / `k_IGK_CD` | 落入 UC / CD 参考集的 CDR3aa 数量 |

## 方法

对每个样本、每条链（IGL / IGK），统计其 CDR3aa 与参考集的命中数：

- `n`：与 cover 集（训练集中健康与疾病样本共有的 CDR3aa）的交集数量
- `k_UC` / `k_CD`：与 UC / CD 关联 CDR3aa 集合的交集数量

每个模型用 Beta-Binomial 分布描述命中概率，`k ~ Beta-Binomial(n, a, b)`，
`N` 为训练该模型的样本数。F 分数为 log Bayes factor：

```
log P(数据 | 疾病模型) / P(数据 | 健康模型)
```

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
    "$PWD/example/SRR18509019/SRR18509019.TCR_clonotypes.txt" > samples.tsv
python3 predict.py --table samples.tsv --out result.tsv
```
