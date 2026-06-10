# P8/P10 参数选择与跨数据一致性分析计划

日期：2026-05-23

这份文档记录当前重新整理后的 P8/P10 分析目标。核心思想是把：

- 最后要选择的参数；
- 需要分层查看、但不是最终选择目标的变量；
- 需要保留/排除的数据；
- xcorr top hit 的 biological interpretation；
- activation map / ROI consistency；

分开处理，避免再把所有东西混在一个不可解释的 top5 平均分数里。

## 当前主要问题

之前的 P8/P10 cross-session summary 有几个问题：

1. 过度依赖 `mean(top5 peak_abs_corr)`，但我们真正关心的不是 topN 平均值，而是 top hit 的身份、density 类型、selectivity label、timescale、lag、ROI/activation 是否跨数据一致。
2. 参数选择变量和分层查看变量混在一起，导致图很多但不回答明确问题。
3. P8/P10 没有足够清楚地区分：
   - BLP observable type；
   - BOLD observable type；
   - BLP density class；
   - BOLD efun type；
   - raw vs standardized training/provenance；
   - P10 中 BOLD-side dimred method/k 和 BLP-side dimred method/k。
4. 目前 P10 还有 partial broad run 和 clean first-pass run 混用同一 tag 的历史问题，所以 P10 结果需要更严格 provenance filter 或重新用 clean tag 跑。

## 需要选择的参数

这些是最终希望通过分析帮助选择的参数。

| 参数 | 当前选项 | 说明 |
|---|---|---|
| BLP observable type | `abs_projected_vlambda`, `complex_split_projected_vlambda` | 第一层主要比较 abs vs complex-split/csplit。 |
| BLP training/input normalization branch | `raw`, `standardized` | 必须作为主线 provenance 分层变量。不能把 raw 和 standardized 混在同一结论里。 |
| BOLD observable type | `HP_svd100`, `global_svd100`, `gsvd100_ds`, `roi_mean` | 当前主线 BOLD observable set。 |
| BOLD training/input normalization branch | `raw`, `standardized` | P7 找 best checkpoint、P8/P9/P10 provenance 都必须区分 raw vs standardized。 |
| BLP dimred method | `SVD`, `NMF`, `MDS`, `UMAP` | 只对 dimred BLP efun density 生效。 |
| BLP dimred component number | `k03:k08` | 不要混在一起先平均；应作为横轴或 facet。 |
| P10 BOLD-side dimred method | `SVD`, `NMF`, `MDS`, `UMAP` | P10 额外参数，来自 P9。不能和 BLP dimred method 混淆。 |
| P10 BOLD-side dimred component number | `k03:k08` | P10 额外参数，来自 P9。 |

## 需要分类查看的变量

这些变量不是最终要单独“挑选”的参数，但必须作为统计分层或解释维度。

| 分类变量 | 选项/含义 |
|---|---|
| pipeline | `P8`, `P10` |
| BLP density class | `event_density`, `raw_BLP_efun_density`, `dimred_BLP_efun_density` |
| BLP density condition | `abs`, `csplit` |
| BOLD efun type | `efun`, `deconv_efun` |
| activity/density version | old/samplewise/maxabs vs current `rmsenv_adaptive` |
| topN | `1`, `3`, `5`, `10`, `20`, `50` or larger sensitivity range |
| lag | peak lag in seconds; should be summarized separately from magnitude |
| sign | signed correlation is useful for QC, but primary similarity score uses `peak_abs_corr` |
| dimred BLP component label | theta/ripple/gamma/SWR-compatible/mixed/nonselective labels from P5 selectivity |
| raw BLP eigenfunction metadata | raw efun index, eigenvalue-derived timescale, frequency, decay rate |
| activation/ROI output | BOLD efun/deconv efun activation map and ROI vector for selected top hits |

## 需要保留的数据

当前主线目标是 9 个数据，以后可以扩展。

```text
e10aw1
e10bv1
e10fV1
e10gb1
e10gh1
e10gw1
f12m01
k13m17
k13m23
```

不应进入当前主线：

```text
f12m01old
_legacy_quarantine
legacy_pipeline4plus_20260511
summary/cache-only folders
```

每个 dataset 应该有状态：

```text
ready
missing_upstream
missing_p5_density
missing_p7_best_checkpoint
missing_p8
missing_p9
missing_p10
invalid_or_stale
legacy_only
```

## 主要分析单位

后续所有统计应该先生成一个干净的 long table，而不是一开始就画 summary heatmap。

建议核心行单位：

```text
pipeline
dataset
BLP_training_branch
BOLD_training_branch
BOLD_observable
BOLD_efun_type
BLP_observable
BLP_density_class
BLP_density_condition
BLP_method
BLP_k
BLP_component_or_raw_mode
P10_BOLD_dimred_method
P10_BOLD_dimred_k
peak_abs_corr
peak_corr
lag_sec
rank_within_context
topN_membership
selectivity_label
raw_efun_timescale_sec
raw_efun_frequency_hz
activation_map_file
roi_summary_file
roi_vector_file
source_csv
source_mat
provenance_status
```

这里的 `rank_within_context` 很重要。context 至少应该包括：

```text
dataset
pipeline
BLP_training_branch
BOLD_training_branch
BOLD_observable
BOLD_efun_type
```

P10 的 context 还要包括：

```text
P10_BOLD_dimred_method
P10_BOLD_dimred_k
```

## 结论 1：raw/dimred BLP efun density 是否比 event density 更强

问题：

```text
raw BLP efun density / dimred BLP efun density
是否比 event density 更能预测或对应 BOLD efun/deconv efun？
```

不要只看 `mean(top5)`。建议看以下统计：

1. Winner fraction：
   - 每个 context 里，`event_density`、`raw_BLP_efun_density`、`dimred_BLP_efun_density` 谁的 rank 最高。
   - 统计每类 density 在 9 个 dataset 中赢了几次。

2. Best-hit delta：

```text
best_raw_abs_corr - best_event_abs_corr
best_dimred_abs_corr - best_event_abs_corr
best_dimred_abs_corr - best_raw_abs_corr
```

3. TopN membership curve：

```text
N = 1, 3, 5, 10, 20, 50
```

看每类 density 在 topN 中出现比例，而不是把 topN 的相关系数平均掉。

4. Candidate-count correction：
   - raw BLP efun density 和 dimred density 的 candidate 数量比 event density 多。
   - 只比较 maximum 会偏向 candidate 更多的类型。
   - 因此必须同时报告：
     - 每类 density 的 candidate count；
     - topN enrichment；
     - per-dataset consistency；
     - permutation 或 matched-count sensitivity（后续可选）。

推荐图：

- density class winner fraction bar plot；
- topN membership curve；
- paired delta violin/box plot；
- dataset x density class best score heatmap；
- raw vs dimred vs event scatter plot。

## 结论 2：高相似度 dimred BLP efun density 的 selectivity label

问题：

```text
和 BOLD deconv efun 相似度高的 dimred BLP efun density，
它们对应 theta/ripple/gamma/SWR/mixed/nonselective 中哪类 subprocess？
```

筛选：

```text
BOLD_efun_type = deconv_efun
BLP_density_class = dimred_BLP_efun_density
```

需要 join P5 dimred component selectivity label。当前 label 应该以后基于：

```text
abs + adaptive RMS envelope activity
```

而不是旧的 signed 或 maxabs prototype。

推荐统计：

- topN 内 label composition；
- 每个 dataset 的 winner label；
- 每个 `BLP observable x method x k x BOLD observable` 的 label enrichment；
- label 是否跨 dataset 一致；
- label 是否依赖 BOLD efun vs deconv efun。

推荐图：

- `method x k` label stacked bar；
- dataset x label heatmap；
- topN label composition curve；
- label winner table；
- abs vs csplit label comparison。

## 结论 3：raw BLP efun density 的 timescale distribution

问题：

```text
和 BOLD efun 相似度高的 raw BLP efun density
和 BOLD deconv efun 相似度高的 raw BLP efun density
是否来自完全不同的 raw efun timescale distribution？
```

筛选：

```text
BLP_density_class = raw_BLP_efun_density
```

分组：

```text
BOLD_efun_type = efun
BOLD_efun_type = deconv_efun
```

需要 metadata：

```text
raw_efun_index
eigenvalue_real_part
eigenvalue_imag_part
timescale_sec
frequency_hz
decay_rate
```

推荐统计：

- raw efun index histogram；
- log-timescale distribution；
- median/quantile timescale per dataset；
- deconv vs efun paired timescale delta；
- topN sensitivity；
- fast/slow bin enrichment。

推荐图：

- raw efun index distribution；
- log-timescale violin/ridge plot；
- ECDF of timescale；
- dataset x median timescale heatmap；
- efun vs deconv efun paired scatter。

## 结论 4：P10 是否和 P8 有类似结论

P10 不能简单和 P8 混在一起，因为 P10 多了 BOLD-side dimred 参数：

```text
P9/P10 BOLD dimred method
P9/P10 BOLD dimred k
```

P10 应该 mirror P8 的问题：

1. density class competition：
   - event vs raw BLP efun vs dimred BLP efun；
2. dimred BLP label enrichment：
   - top P10 hits 里的 BLP dimred component 是 theta/ripple/mixed 还是 nonselective；
3. raw BLP efun timescale distribution：
   - P10 是否也偏向同一类 fast/slow raw modes；
4. P8 vs P10 conclusion agreement：
   - winner density class 是否一致；
   - dimred label distribution 是否一致；
   - raw timescale distribution 是否一致；
   - abs/csplit、method、k 的 preference 是否一致。

推荐图：

- P8 vs P10 density class winner agreement heatmap；
- P8 vs P10 topN membership curves side by side；
- P8 vs P10 dimred label composition comparison；
- P8 vs P10 raw timescale ECDF overlay；
- P10 BOLD-side method/k x BLP-side method/k score/label matrix。

## 结论 5：高相似度 BOLD efun/deconv efun 的 activation map 和 ROI summary 是否跨数据一致

这个应该是第二阶段分析，不应该默认对所有 top hits 大量画 activation map。

第一阶段先用前四个问题选出少数 candidate groups，例如：

```text
BLP_training_branch = standardized
BOLD_training_branch = standardized
BLP_observable = complex_split_projected_vlambda
BOLD_observable = roi_mean
BOLD_efun_type = deconv_efun
BLP_density_class = dimred_BLP_efun_density
BLP_method = NMF
BLP_k = k06
selectivity_label = theta/ripple/SWR-compatible
```

第二阶段再看这些 group 的 BOLD activation/ROI：

- 每个 dataset 取 top hit 或 top few hits；
- 找对应 BOLD efun/deconv efun activation map；
- 找对应 ROI vector；
- 做 cross-dataset ROI correlation；
- 做 top ROI overlap；
- 做 ROI profile clustering；
- activation map 作为浏览和 sanity check，不作为第一阶段筛选依据。

推荐图：

- dataset x dataset ROI correlation heatmap；
- top ROI overlap matrix；
- ROI profile line/bar overlay；
- selected activation map contact sheet；
- per-candidate-group ROI consistency score table。

## topN 策略

topN 不应该主要作为平均值。

推荐把 topN 当作 sensitivity axis：

```text
N = 1, 3, 5, 10, 20, 50
```

每个 N 计算：

- density class winner/membership；
- dimred selectivity label composition；
- raw efun timescale distribution；
- P8/P10 agreement；
- activation/ROI consistency候选是否稳定。

可保留的分数：

```text
best_peak_abs_corr
median_peak_abs_corr_in_topN
max_peak_abs_corr_in_topN
rank_of_best_event_density
rank_of_best_raw_efun_density
rank_of_best_dimred_efun_density
topN_membership_count
topN_membership_fraction
```

不建议作为主结论：

```text
mean(top5 peak_abs_corr)
```

原因：

- 它会丢掉 top hit identity；
- 它会混合不同 density class；
- 它不告诉我们 selectivity label；
- 它不告诉我们 raw efun timescale；
- 它对 candidate 数量不敏感；
- 它不容易映射到 activation/ROI consistency。

## 推荐 P11 新增输出

P11 不应该只平铺图，而应该生成当前分析需要的 canonical long tables 和少量 summary figures。

建议新增：

```text
results\pipeline11_current_mainline_check\P8_P10_xcorr_hit_long_table.csv
results\pipeline11_current_mainline_check\P8_P10_density_class_competition.csv
results\pipeline11_current_mainline_check\P8_P10_topN_sensitivity.csv
results\pipeline11_current_mainline_check\P8_P10_dimred_selectivity_label_hits.csv
results\pipeline11_current_mainline_check\P8_P10_raw_efun_timescale_hits.csv
results\pipeline11_current_mainline_check\P8_P10_roi_consistency_candidates.csv
```

建议新增 summary figure folders：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\density_class_competition\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\topN_sensitivity\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\dimred_selectivity_labels\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\raw_efun_timescales\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\roi_consistency_candidates\
```

## 推荐图组设计

重要原则：

- 只使用当前 available/current/non-stale data。
- 数据不齐时，不强行补空条件；每张图必须标明 `n_available dataset` 和 dataset names。
- raw 和 standardized 必须分开画，不能混成一个结论。
- `efun` 和 `deconv_efun` 必须分开画，不能混成一个结论。
- P8 和 P10 各自都要有完整图组；P8/P10 agreement 只能作为最后的对照图，不能替代 P10 自己的分析。
- topN 不固定为 5，而是作为 sensitivity axis：

```text
topN = 1, 3, 5, 10, 20, 50
```

推荐 summary figure 根目录：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\
```

推荐子目录：

```text
00_availability\
01_p8_density_class_competition\
02_p10_density_class_competition\
03_p8_blp_method_k\
04_p10_blp_method_k\
05_p10_bold_method_k\
06_p8_dimred_selectivity_labels\
07_p10_dimred_selectivity_labels\
08_p8_raw_efun_timescales\
09_p10_raw_efun_timescales\
10_p8_roi_activation_consistency\
11_p10_roi_activation_consistency\
12_p8_vs_p10_comparison\
```

### 00. Data Availability

每组图之前都应该先有可用性图，避免把数据缺失误解成 biological result。

图：

```text
condition availability heatmap
```

设计：

- x 轴：dataset；
- y 轴：condition；
- color：`available`, `missing`, `stale`, `invalid`, `legacy_only`；
- condition 至少包括：

```text
pipeline
BLP raw/standardized branch
BOLD raw/standardized branch
BLP observable type
BOLD observable type
BOLD efun type
P8/P10 status
```

P10 的 availability 还要包括：

```text
P9/BOLD dimred method
P9/BOLD dimred k
```

每张后续分析图都要在标题或 subtitle 写：

```text
n_available = ...
datasets = ...
```

### 01. P8 Density Class Competition

问题：

```text
P8 中 raw BLP efun density / dimred BLP efun density / event density 谁更强？
```

图 A：

```text
P8 topN density-class membership curve
```

设计：

- x 轴：topN = `1, 3, 5, 10, 20, 50`；
- y 轴：topN 内属于某个 density class 的比例；
- color：`event_density`, `raw_BLP_efun_density`, `dimred_BLP_efun_density`；
- facet：
  - BLP raw/standardized；
  - BOLD raw/standardized；
  - BLP observable type；
  - BOLD observable type；
  - BOLD efun type；
- 每个 facet 标 `n_available`。

图 B：

```text
P8 density-class winner fraction
```

设计：

- 每个 available context 内，top1 winner 属于哪类 density；
- y 轴：available dataset fraction；
- color：event/raw/dimred；
- 每个 bar 标 `count / n_available`。

### 02. P10 Density Class Competition

P10 必须有和 P8 对称的图，不能只画 P8/P10 agreement。

图 A：

```text
P10 topN density-class membership curve
```

设计：

- x 轴：topN；
- y 轴：topN 内 density class 比例；
- color：event/raw/dimred；
- facet 同 P8，并额外包括：
  - P9/BOLD dimred method；
  - P9/BOLD dimred k。

图 B：

```text
P10 density-class winner fraction
```

设计：

- 同 P8；
- 每张图或每个 facet 固定 P9/BOLD method/k；
- 不把所有 P9 method/k 混在一起。

### 03. P8 BLP Dimred Method/k Parameter Heatmap

问题：

```text
P8 里哪个 BLP dimred method/k 更稳定进入 top hits 或成为 winner？
```

图：

```text
P8 BLP dimred method-k robustness heatmap
```

设计：

- x 轴：`k03:k08`；
- y 轴：`SVD`, `NMF`, `MDS`, `UMAP`；
- color：在 available datasets 中进入 topN 或 top1 winner 的比例；
- 每格标：`count / n_available`；
- facet：
  - BLP observable type；
  - BLP raw/standardized；
  - BOLD observable type；
  - BOLD raw/standardized；
  - BOLD efun type。

### 04. P10 BLP Dimred Method/k Parameter Heatmap

问题：

```text
P10 里哪个 BLP dimred density method/k 更稳定？
```

图：

```text
P10 BLP dimred method-k robustness heatmap
```

设计：

- x 轴：BLP `k03:k08`；
- y 轴：BLP `SVD/NMF/MDS/UMAP`；
- color：topN 或 winner fraction；
- facet：
  - BLP observable type；
  - BLP raw/standardized；
  - BOLD observable type；
  - BOLD raw/standardized；
  - BOLD efun type；
  - P9/BOLD method/k。

### 05. P10 BOLD-side Dimred Method/k Parameter Heatmap

问题：

```text
P10 里哪个 P9/BOLD-side dimred method/k 更稳定？
```

图 A：

```text
P10 BOLD-side dimred method-k robustness heatmap
```

设计：

- x 轴：P9/BOLD `k03:k08`；
- y 轴：P9/BOLD `SVD/NMF/MDS/UMAP`；
- color：topN 或 winner fraction；
- facet：
  - BLP density class；
  - BLP observable type；
  - BOLD observable type；
  - BOLD efun type；
  - raw/standardized branches。

图 B，可选但非常有用：

```text
P10 BLP-method-k x BOLD-method-k matrix
```

设计：

- x 轴：BLP method_k；
- y 轴：P9/BOLD method_k；
- color：available dataset winner/topN fraction；
- 用于判断 P10 里 BLP dimred 和 BOLD dimred 是否有稳定组合。

### 06. P8 Dimred Selectivity Label

问题：

```text
P8 里和 BOLD efun/deconv efun 高相似的 dimred BLP component 是 theta/ripple/什么？
```

图 A：

```text
P8 dimred selectivity label composition
```

设计：

- x 轴：BLP method_k；
- y 轴：topN hit proportion；
- color：
  - theta；
  - ripple；
  - theta+ripple；
  - SWR-compatible；
  - gamma；
  - mixed；
  - nonselective；
- facet：
  - BOLD observable；
  - BOLD efun type；
  - BLP observable；
  - raw/standardized branches。

图 B：

```text
P8 dataset-by-label tile
```

设计：

- x 轴：dataset；
- y 轴：BLP method_k；
- color：top hit selectivity label；
- 每张图固定一个 high-level condition。

### 07. P10 Dimred Selectivity Label

P10 也要画自己的 selectivity label 图。

图 A：

```text
P10 dimred selectivity label composition
```

设计：

- x 轴：BLP method_k；
- y 轴：topN hit proportion；
- color：selectivity label；
- facet：
  - BOLD observable；
  - BOLD efun type；
  - BLP observable；
  - raw/standardized branches；
  - P9/BOLD method_k。

图 B：

```text
P10 dataset-by-label tile
```

设计：

- x 轴：dataset；
- y 轴：BLP method_k 或 P9/BOLD method_k；
- color：top hit selectivity label；
- P9/BOLD method/k 必须写在标题或 facet 中。

### 08. P8 Raw Efun Timescale

问题：

```text
P8 中 raw BLP efun density top hits 的 timescale 是否随 BOLD efun/deconv efun 改变？
```

图 A：

```text
P8 raw efun timescale ECDF/ridge
```

设计：

- x 轴：`log10(timescale_sec)`；
- color：`efun` vs `deconv_efun`；
- facet：
  - BOLD observable；
  - BLP observable；
  - raw/standardized branches；
  - topN。

图 B：

```text
P8 raw efun index histogram
```

设计：

- x 轴：raw efun index；
- y 轴：top hit count；
- 分开 efun/deconv efun；
- 标明 topN 和 n_available。

### 09. P10 Raw Efun Timescale

P10 也要画自己的 raw efun timescale 图。

图 A：

```text
P10 raw efun timescale ECDF/ridge
```

设计：

- x 轴：`log10(timescale_sec)`；
- color：`efun` vs `deconv_efun`；
- facet：
  - BOLD observable；
  - BLP observable；
  - raw/standardized branches；
  - P9/BOLD method_k；
  - topN。

图 B：

```text
P10 raw efun index histogram
```

设计：

- x 轴：raw efun index；
- y 轴：top hit count；
- 分开 efun/deconv efun；
- 标明 P9/BOLD method/k、topN、n_available。

### 10. P8 ROI / Activation Consistency

这一步只对 P8 筛出的少数 candidate groups 画，不对所有 top hits 画。

图：

```text
P8 ROI profile correlation heatmap
P8 ROI profile overlay
P8 activation map contact sheet
```

每张图必须标明：

```text
candidate group
n_available
dataset list
topN rule
selected hit rule
```

### 11. P10 ROI / Activation Consistency

P10 也要有自己的 ROI/activation consistency 图。

图：

```text
P10 ROI profile correlation heatmap
P10 ROI profile overlay
P10 activation map contact sheet
```

每张图额外标明：

```text
P9/BOLD method_k
BLP density source
BLP method_k if dimred
```

### 12. P8 vs P10 Comparison

最后才画 P8/P10 对照图。只在 P8 和 P10 都有同一 condition 的时候画。

图：

```text
P8 vs P10 density-class comparison
P8 vs P10 selectivity-label comparison
P8 vs P10 raw-timescale comparison
P8 vs P10 ROI-consistency comparison
```

设计：

- x 轴：P8 conclusion；
- y 轴：P10 conclusion；
- color：dataset count 或 agreement fraction；
- 没有 paired data 的 condition 不画，或者明确标为 missing。

## 需要在代码里注意的 provenance 规则

1. raw 和 standardized 必须分开。
   - P4 best checkpoint resolution 要按 `dataset x observable x normalization branch` 找。
   - P7/P8/P9/P10 都要继承这个 branch。

2. P8/P10 必须继承 current P7/P9 provenance。
   - 不能只看 folder 是否存在。
   - 必须确认对应 current BOLD_POST / P9 result。

3. P5 density version 必须写清楚。
   - old/maxabs/samplewise 与 `rmsenv_adaptive` 不能混。
   - 当前新分析优先使用 `rmsenv_adaptive`。

4. BLP dimred result 与 dimred density 要分清。
   - `pipeline5_eigenfunction_reduction` 保存 component 本体。
   - `pipeline5_dimred_thresholded_density` 是 derived density。

5. P10 必须避免同一个 save tag 混入 partial broad run。
   - 后续最好换 clean tag 或在 summary 脚本里严格 filter：

```text
P9 feature
P9 method/k
BLP density condition
BLP density method/k
activity version
dataset provenance
```

## 当前建议优先级

第一步：生成干净 long table。

第二步：做 density class competition 和 topN sensitivity。

第三步：join P5 dimred selectivity label。

第四步：join raw efun timescale metadata。

第五步：比较 P8 vs P10。

第六步：只对少数 candidate groups 做 ROI/activation consistency。

这样分析才能服务最终参数选择：

```text
BLP observable type
BOLD observable type
BLP dimred method
BLP dimred component number
raw vs standardized branch
```

## 2026-05-23 当前图生成记录

已新增图优先脚本：

```text
scripts/plot_p8_p10_parameter_selection_figures.py
```

这版脚本不再使用旧的 top5 `*_top.csv` 作为唯一来源，而是读取完整
`*_peaks.csv` 文件，重新构建：

```text
topN = 1, 3, 5, 10, 20, 50
```

因此 `topN > 5` 的图不需要重新计算 xcorr，只要对应的 `peaks.csv`
存在即可。

当前输出：

```text
results\p8_p10_parameter_selection_current_rmsenv_adaptive\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection\
```

当前 available 数据：

```text
e10fV1
e10gb1
e10gh1
e10gw1
f12m01
```

当前缺少 P8/P10 peaks 的数据：

```text
e10aw1
e10bv1
k13m17
k13m23
```

当前已生成图组：

```text
00_availability
01_p8_density_class_competition
02_p10_density_class_competition
03_p8_blp_method_k
04_p10_blp_method_k
05_p10_bold_method_k
06_p8_dimred_selectivity_labels
07_p10_dimred_selectivity_labels
08_p8_raw_efun_timescales
09_p10_raw_efun_timescales
10_p8_roi_activation_consistency
11_p10_roi_activation_consistency
12_p8_vs_p10_comparison
```

当前重要限制：

- raw/standardized branch 还没有从现有 P8/P10 xcorr filenames 中稳定解析
  出来；当前图代表 current available branch，不能混称为 raw 或
  standardized。
- raw efun index histogram 已经可以画，因为 peaks CSV 里有
  `density_index`。
- raw efun timescale ECDF 还不能画真实数值，因为 peaks CSV 里没有
  `timescale_sec`；需要后续从 P5 raw density `mode_metadata` join/export 到
  P8/P10 hit rows。
- ROI/activation consistency 目前只生成 candidate-selection-needed 占位图；
  真正的 ROI/map consistency 应该在从前面图组筛出少数 candidate groups
  后再跑。

### 2026-05-23 图 v1 需要修改的问题

当前 `p8_p10_parameter_selection` 第一版图已经能从 full `peaks.csv`
重建 `topN > 5`，但是有几类图的统计定义还需要修改。

#### 03/04 BLP method-k robustness

当前 `03_p8_blp_method_k` 和 `04_p10_blp_method_k` 画的是：

```text
top10 membership fraction
```

也就是某个 BLP dimred `method/k` 是否出现在 top10 hits 中：

```text
出现的 context 数 / available context 总数
```

它不是 `mean(top10 peak_abs_corr)`。

这个定义比 topN 平均值更接近当前目标，但仍然不够，因为它把 rank1
和 rank10 当成一样。下一版应改成三个互补图：

```text
top1 winner fraction
topN membership curve: N = 1, 3, 5, 10, 20, 50
median best rank heatmap
```

其中：

- `top1 winner fraction` 回答哪个 method/k 经常赢；
- `topN membership curve` 回答哪个 method/k 随 topN 放宽仍然稳定出现；
- `median best rank heatmap` 回答它通常排在多靠前，而不只是是否进入 top10。

#### 05 P10 BOLD-side method-k robustness

当前 `05_p10_bold_method_k` 的数值不合理。它现在更接近：

```text
某个 P9/BOLD method-k 占全部 P10 context 的比例
```

所以会出现类似 `0.11` 或 `0.01` 这种很小的数值。这不能解释为
performance、stability 或 biological preference。

下一版应废弃这张 v1 图，改成：

```text
P10 BOLD-side method-k top1 winner fraction
P10 BOLD-side method-k topN membership curve
P10 BOLD-side method-k median best rank
```

并且 P10 必须分层或固定：

```text
BLP density class
BLP observable type
BOLD observable type
BOLD efun type
raw/standardized branch
```

不能把所有 P9 method/k、BLP density class 和 BOLD efun/deconv efun 混在
一起算一个比例。

#### 06/07 Dimred selectivity label composition

当前 `06_p8_dimred_selectivity_labels` 和
`07_p10_dimred_selectivity_labels` 中的数值是比例：

```text
在 top10 dimred BLP density hits 中，
某个 BLP method/k 被标成某个 selectivity label 的比例
```

这些数值不是 xcorr 强度。

当前 label 含义：

```text
partial_or_inactive
```

表示 component 对某些 event family 有 activity，但没有满足干净的
theta/ripple selectivity 标准。

当前：

```text
unlabeled
```

混合了两种情况：

1. P5 label table 没有成功 join 到该 component；
2. P5 label 认为该 component 没有明确 event-family activity。

这会混淆“数据缺失”和“生物学上 nonselective”。下一版必须拆开：

```text
label_missing
no_event_activity
partial_or_inactive
```

建议最终 label set 至少包括：

```text
theta_selective_similar
theta_selective_unequal
ripple_selective_similar
ripple_selective_unequal
theta_ripple_joint
mixed_theta_ripple
mixed_theta_gamma
mixed_ripple_gamma
gamma_selective
pan_event
partial_or_inactive
no_event_activity
label_missing
nonselective
```

图标题和 caption 必须写清楚：

```text
value = fraction of topN dimred BLP density hits
```

#### 08/09 Raw efun timescale distribution

当前 `08_p8_raw_efun_timescales` 和 `09_p10_raw_efun_timescales`
没有真正回答原始问题。

原始问题是：

```text
和 BOLD efun 相似度高的 raw BLP efun density
vs
和 BOLD deconv efun 相似度高的 raw BLP efun density

是否来自完全不同的 raw efun timescale distribution？
```

当前 v1 只画了 pooled raw efun index histogram。它没有：

- 分开 `efun` 和 `deconv_efun`；
- 画 timescale distribution；
- 对比两个 distribution。

原因是当前 P8/P10 `peaks.csv` 里有：

```text
density_index
```

但没有：

```text
timescale_sec
frequency_hz
lambda_continuous_real
lambda_continuous_imag
```

下一版必须先把 P5 raw density `mode_metadata` join/export 到 P8/P10 hit
rows，然后画：

```text
raw efun timescale ECDF: efun vs deconv_efun
raw efun timescale violin/ridge: efun vs deconv_efun
raw efun index histogram: efun and deconv_efun separated
dataset-level median timescale paired plot
```

这些图应分别对 P8 和 P10 生成，而不是只画一个 pooled histogram。

#### 当前 v1 图的使用建议

可以暂时参考：

```text
01/02 density class topN curves
03/04 top10 membership idea
06/07 label composition rough view
```

需要重画或废弃：

```text
05 P10 BOLD-side method-k robustness
08/09 raw efun timescale plots
```

需要修正后再解释：

```text
06/07 label composition: split label_missing vs no_event_activity
03/04 method-k robustness: add winner fraction and median rank
```

### 2026-05-23 v2 实现记录

本轮已经把上面这些修改落到脚本和新输出目录里。旧版图没有覆盖，v2 单独保存。

代码入口：

```text
scripts/plot_p8_p10_parameter_selection_figures.py
scripts/export_p5_raw_density_mode_metadata_for_p8_p10.m
```

v2 数值结果：

```text
results/p8_p10_parameter_selection_current_rmsenv_adaptive_v2/
```

v2 图：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection_v2\
```

本轮生成状态：

```text
P8 peaks rows: 184000
P10 peaks rows: 444213
P5 raw density mode_metadata rows: 2512
raw density hit rows with finite positive timescale: 32182
TopN values: 1, 3, 5, 10, 20, 50
main heatmap TopN: 10
```

v2 图的定义如下。

#### 03/04 BLP method-k v2

现在不再只画 “top10 里出现过的比例”。每个 P8/P10 的 BLP method-k 都画三类图：

```text
topN_membership:
  topN = 1,3,5,10,20,50 时，该 method-k 是否进入 topN xcorr hits

winner_fraction:
  在 top10 xcorr hits 里，dimred BLP density 中 rank 最靠前的 method-k
  是该 method-k 的 context 比例

median_rank:
  该 method-k 在进入 top10 时的 median best rank
```

这些图回答的问题分别是：

```text
这个 method-k 是否经常出现？
这个 method-k 是否经常赢？
这个 method-k 赢不了时通常排多靠前？
```

输出目录：

```text
03_p8_blp_method_k_v2
04_p10_blp_method_k_v2
```

#### 05 P10 BOLD-side P9 method-k v2

旧版 05 的定义已经废弃。v2 按 P10 upstream P9 method-k 重新统计：

```text
parent context = dataset × BOLD observable × BOLD efun family × condition_scope
score per P9 method-k = that P9 method-k run 中 topN xcorr hits 的最大 |r|
winner_fraction = 在 parent context 内成为最高 score 的 P9 method-k 比例
median_rank = 在 parent context 内的 P9 method-k 排名中位数
```

这里的 `winner_fraction` 本来就会比较小，因为候选 P9 method-k 很多；它不能和
BLP density class 的 winner fraction 直接比较。

输出目录：

```text
05_p10_bold_method_k_v2
```

2026-05-23 correction: the first v2 `topN_presence` plot was still misleading.
It measured whether a P9 method-k run existed for a parent context after taking
topN density hits inside that run. Because every available run has at least one
hit, changing density topN did not change this value, so the plot looked
artificially regular. Treat `p10_bold_method_k_topN_presence__*.png` as an
availability diagnostic only, not as a scientific result.

The corrected v2.1 output is:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection_v2_1\05_p10_bold_method_k_v2_1\
```

v2.1 definitions:

```text
density_top10:
  use the top10 density hits inside each P10 run to score that P9 method-k run

method_topK_fraction:
  after ranking P9 method-k candidates inside each parent context, fraction of
  contexts where the method-k is in method-rank topK

winner_fraction:
  fraction of contexts where the method-k is rank 1 among P9 method-k candidates

median_method_rank:
  median rank among available P9 method-k candidates. This can be larger than
  10 because it is not the rank within top10 density hits; it is the rank among
  up to 24 current P9 method-k candidates.

median_score:
  median max |r| from density_top10 hits for that P9 method-k
```

#### 06/07 Dimred selectivity label v2

现在把旧的 `unlabeled` 拆成：

```text
label_missing:
  P8/P10 hit row 没有成功 join 到 P5 dimred component label table

no_event_activity:
  成功 join 到 P5 label table，但 P5 selectivity 判断为没有明确 event-family activity
```

因此这些图现在可以区分 “数据/metadata 没接上” 和 “component 本身不 selective”。

输出目录：

```text
06_p8_dimred_selectivity_labels_v2
07_p10_dimred_selectivity_labels_v2
```

#### 08/09 Raw efun timescale v2

现在先从 P5 raw density MAT 中导出 `mode_metadata`，再 join 回 P8/P10 hit rows。

新增图：

```text
raw efun index histogram: efun vs deconv_efun 分开显示
raw efun timescale ECDF: efun vs deconv_efun 分开显示
```

这直接回答：

```text
和 BOLD efun 相似度高的 raw BLP efun density
vs
和 BOLD deconv efun 相似度高的 raw BLP efun density
是否来自不同 raw efun timescale distribution？
```

注意：本轮 raw metadata join 是成功的；少数 hit 没有有限正 timescale，
是 eigenvalue real part 导致 timescale 本身不可定义，不是 path join 失败。

输出目录：

```text
08_p8_raw_efun_timescales_v2
09_p10_raw_efun_timescales_v2
```

当前 v2/v2.1 的初步解读：

```text
和 BOLD efun 相似度高的 raw BLP efun density
vs
和 BOLD deconv_efun 相似度高的 raw BLP efun density

当前趋势是 deconv_efun 相关 raw efun 稍微更偏长 timescale，P10 里更明显。
```

但是这个结论不能只靠 `top10`。下一步必须做 topN sensitivity：

```text
topN = 1,3,5,10,20,50
每个 topN 分别画 raw efun timescale ECDF: efun vs deconv_efun
再画 median log10(timescale_sec) 随 topN 变化的 summary line plot
输出 raw_timescale_topN_summary.csv，记录 n / median / q25 / q75 / mean
```

解释规则：

```text
如果 deconv_efun 在 top1/top3/top5/top10/top20/top50 都稳定右移，
才可以说 deconv_efun 相关 raw BLP efun 偏长 timescale。

如果只在 top10 或 top20 出现，而 top1/top3 不出现，
则只能说这是中等 topN pooled hit 的趋势，不能作为强结论。

后续还需要 dataset-stratified 版本，避免被某一个 dataset 带偏。
```

Implementation status, 2026-05-23:

```text
Completed in scripts/plot_p8_p10_parameter_selection_figures.py.
Generated under:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_parameter_selection_v2_1\

CSV:
results\p8_p10_parameter_selection_current_rmsenv_adaptive_v2_1\raw_efun_timescale_topN_summary.csv

Figures:
08_p8_raw_efun_timescales_v2\p8_raw_efun_timescale_median_by_topN.png
08_p8_raw_efun_timescales_v2\topN_sensitivity\p8_raw_efun_timescale_ecdf_by_family_top*.png
09_p10_raw_efun_timescales_v2\p10_raw_efun_timescale_median_by_topN.png
09_p10_raw_efun_timescales_v2\topN_sensitivity\p10_raw_efun_timescale_ecdf_by_family_top*.png
```

First read after adding topN sensitivity:

```text
P8:
  top3/top5/top10 show deconv_efun raw hits shifted to longer timescale,
  but the separation weakens by top20/top50.

P10:
  top1 does not support deconv_efun being longer.
  top3/top5/top10/top20/top50 increasingly support deconv_efun being longer.

Therefore the current wording should be:
  "deconv_efun-related raw BLP efun hits tend to include longer-timescale modes
   once more than the single best hit is considered, especially in P10."

It should not be stated as:
  "the single strongest deconv_efun hit is always longer-timescale."
```

Interpretation note:

```text
Longer timescale here means slower Koopman decay mode.

timescale_sec is interpreted from the continuous-time eigenvalue real part:
timescale_sec = -1 / real(lambda_continuous), when real(lambda_continuous) < 0.

Therefore, larger timescale_sec means the mode decays more slowly.
This is not the same as oscillation frequency; oscillatory frequency is tied
to the eigenvalue imaginary part.
```

Current biological/analysis interpretation:

```text
The topN sensitivity result is consistent with the prior observation that
BOLD deconv_efun coupling tends to involve slower raw BLP Koopman modes than
BOLD efun coupling, especially in P10 and when considering topN hits beyond
the single best hit.

This should be treated as an interpretation hypothesis to preserve and test
with dataset-stratified plots, not yet as a final cross-session conclusion.
```
