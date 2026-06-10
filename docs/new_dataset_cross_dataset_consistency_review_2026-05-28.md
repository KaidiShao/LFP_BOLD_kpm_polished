# 新数据 cross-dataset consistency 分析回顾

日期：2026-05-28

目的：为后续新开一份“一个新 standardized complex-split dataset 应该怎么分析，才能判断是否复现 E10gb1 hypothesis”的正式设计文档做准备。

这份文档先回顾已有材料和已有分析，不直接给最终规则。重点是把已经发生过的口径变化写清楚，尤其是 `ripple-gamma` label 的定义问题。

## 1. 当前要验证的 hypothesis

当前工作性假设来自 E10gb1 standardized complex-split 分析：

```text
theta-like slow/state BLP subprocess
    -> BOLD efun / slow BOLD state variable

ripple-gamma-like fast/intrinsic-trigger BLP subprocess
    -> BOLD deconv_efun / residual perturbation
```

这里的核心不是简单找最大 xcorr，而是判断：

1. BLP eigenfunction-derived density 是否比 event density 更能解释 BOLD efun/deconv_efun。
2. P5 dimred BLP components 是否能在同一个 method/k 下分出 theta-like 和 ripple-gamma-like 两个 subprocess。
3. P8/P10 的 top xcorr hits 是否分别落在符合解释的 BLP subprocess 上。
4. raw BLP efun density 的 top hits 是否在 BOLD efun vs deconv_efun 之间呈现 slow/fast timescale 分离。
5. 这些 pattern 在新 dataset 里是否和 E10gb1 一致，而不是数值完全一样。

## 2. 已回顾的核心文档

### 2.1 `docs/p8_p10_parameter_selection_analysis_plan_2026-05-23.md`

这份文档把 P8/P10 分析重新整理成参数选择问题。

保留下来的关键原则：

- 不再依赖 `mean(top5 peak_abs_corr)` 作为主结论。
- `topN = 1, 3, 5, 10, 20, 50` 应该作为 sensitivity axis。
- `raw` 和 `standardized` 必须分开。
- `efun` 和 `deconv_efun` 必须分开。
- P8 和 P10 要各自完整分析，最后再比较 agreement。
- P10 有额外的 BOLD-side dimred method/k，不能和 BLP-side dimred method/k 混淆。
- ROI/activation 是第二阶段，只对少数 candidate groups 做，不应该在第一阶段大量默认画图。

这份文档里已经定义了需要保留的核心 long-table 字段，例如：

```text
dataset
pipeline
BOLD_observable
BOLD_efun_type
BLP_density_class
BLP_method
BLP_k
peak_abs_corr
lag_sec
rank_within_context
selectivity_label
raw_efun_timescale_sec
provenance_status
```

### 2.2 `docs/eigenfunction_activity_magnitude_policy_2026-05-22.md`

这份文档定义了当前 P5/P8/P10 应该使用的 activity 口径。

关键结论：

- 不应该直接用 signed eigenfunction 做 threshold，因为 eigenfunction 有 sign/polarity ambiguity。
- 当前主线应该用 magnitude/envelope：

```text
activity(t, mode) = sqrt(movmean(abs(phi(t, mode)).^2, adaptive_window))
```

- adaptive window 应该根据 Koopman eigenvalue 的 decay timescale 决定。
- raw efun slow/fast 解释现在优先用 discrete Koopman eigenvalue：

```text
tau = -dt / log(abs(lambda_discrete))
```

旧的 continuous/bilinear estimate 只能作为 legacy/fallback provenance，不作为主解释依据。

### 2.3 `docs/e10gb1_p5_standardize_condition_comparison_archive_2026-05-24.md`

这份文档记录 E10gb1 raw vs standardized P5 probe。

关键结论：

- strict ripple-only selectivity 在 standardization 后并没有变好。
- strict ripple-only 不是当前最合适的 biological target。
- 当前更合理的 subprocess 定义是：

```text
theta subprocess
ripple-gamma subprocess
```

- gamma 和 ripple 同时 active 是可以接受的。
- sharp-wave-ripple 不能当成 pure ripple-only，因为它可能含 theta/gamma 结构。
- pure theta activity 是判断 ripple-gamma component 时最需要警惕的 confound。
- E10gb1 里 standardized complex split 是目前最有用的分支。

E10gb1 standardized csplit 的两 subprocess candidates 很多：

```text
SVD:  k03-k08
NMF:  k03-k08
MDS:  k04-k08
UMAP: k05-k08
```

紧凑例子：

```text
standardized csplit, SVD k03:
    theta component = 1
    ripple-gamma no-pure-theta component = 2

standardized csplit, NMF k03:
    theta component = 2
    ripple-gamma no-pure-theta component = 1
```

### 2.4 `docs/e10gb1_bold_blp_slow_variable_intrinsic_trigger_model_2026-05-25_zh.md`

这份文档记录了 E10gb1 hypothesis 的核心叙事。

主要实验回顾：

1. 四个 BLP condition 对比：

```text
abs_projected_vlambda
complex_split_projected_vlambda
abs_projected_vlambda_standardize
complex_split_projected_vlambda_standardize
```

当前最有解释力的是：

```text
complex_split_projected_vlambda_standardize
```

2. raw BLP efun density 的 top xcorr 时间序列：

- `BOLD efun x raw BLP density` 常常对应 slow state / session-wise drift。
- `BOLD deconv_efun x raw BLP density` 更容易对应 fast activity / burst-like activity。

3. dimred BLP density 的 top xcorr label：

- E10gb1 P8/P10 的 dimred top hits 包含 `theta_strict` 和 `ripple_gamma_no_pure_theta`。
- dimred top5 time-series 也支持 slow variable / intrinsic trigger 的解释。

### 2.5 `docs/pipeline12_standardized_csplit_consistency_design_2026-05-25.md`

这份文档是 P12 的第一版设计。

P12 的第一版目标：

```text
只看 standardized complex split
只看 RMS-envelope adaptive density
检查新 dataset 是否复现 E10gb1 的 slow/state vs trigger pattern
```

必要输入：

```text
P4 full EDMD outputs with efuns
P5 raw thresholded density for standardized csplit
P5 dimred thresholded density for standardized csplit, SVD/NMF/MDS/UMAP k03:k08
P5 dimred component process labels
P8 xcorr peak tables for standardized csplit
P10 xcorr peak tables for standardized csplit
raw efun mode metadata with preferred discrete-log timescale
```

P12 的已有判据：

1. standardized csplit P5 能找到 theta component 和 ripple-gamma component。
2. dimred density 在 P8/P10 中相对 event density 有竞争力。
3. 高 xcorr dimred density 的 label 随 BOLD-side target 分化。
4. raw density 的 topN timescale 在 efun vs deconv_efun 间有分化。
5. 以上 pattern 按 BOLD observable 分开后仍可解释。

## 3. 已回顾的核心结果目录

### 3.1 E10gb1 P5 two-subprocess candidates

结果：

```text
results\e10gb1_p5_two_subprocess_candidates\summary.md
results\e10gb1_p5_two_subprocess_candidates\method_k_two_subprocess_candidates.csv
results\e10gb1_p5_two_subprocess_candidates\component_subprocess_flags.csv
```

结论：

- E10gb1 standardized csplit 很多 method/k 同时有 theta component 和 ripple-gamma no-pure-theta component。
- 这是当前 hypothesis 的主要 anchor dataset。

### 3.2 E10gb1 standardized csplit P8/P10 probe

结果：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\
```

summary 中的关键点：

- P8 winner source: `dimred=8/8`。
- P10 winner source: `dimred=56`, `raw=32`, `event=8`。
- P8 efun strongest dimred hits 多为 `theta_strict`。
- P10 中也有 `ripple_gamma_no_pure_theta` top hits。

重要解释：

- E10gb1-only probe 支持 hypothesis，但不是 cross-session conclusion。

### 3.3 P12 standardized csplit 7-dataset results

结果：

```text
results\pipeline12_standardized_csplit_consistency\
```

当前 7 个 dataset：

```text
e10gb1
e10fV1
e10gh1
e10gw1
f12m01
k13m17
k13m23
```

已完成：

- P5 standardized csplit outputs complete。
- P8/P10 standardized csplit outputs complete for available BOLD observables。
- P12 `missing_requirements.csv` 均为空。

当前可用 summary figures：

```text
results\pipeline12_standardized_csplit_consistency\topn_comparison_7datasets\
results\pipeline12_standardized_csplit_consistency\matched_topn_panels_7datasets\
E:\DataPons_processed\summary_figures\pipeline12_standardized_csplit_consistency_7datasets\
```

## 4. 重要修正：P12 当前 label 图会低估 ripple-gamma subprocess

这是这次回顾最重要的一点。

当前 P12 的 dimred label 图主要使用：

```text
primary_process_label
```

这个 label 偏严格，强调 clean selectivity。它适合回答：

```text
top xcorr hits 里有多少是 clean theta-selective 或 clean ripple/gamma-selective？
```

但它不适合回答：

```text
top xcorr hits 里有多少包含 ripple-gamma subprocess activity？
```

E10gb1 的例子很清楚。当前 standardized csplit P5 label 表中：

```text
total dimred components = 132
ripple_active = 79 / 132 = 59.8%
gamma_active = 66 / 132 = 50.0%
ripple_gamma_joint = 52 / 132 = 39.4%
ripple_gamma_no_theta_active = 19 / 132 = 14.4%
theta_selective = 24 / 132 = 18.2%
```

但是这些 `ripple_gamma_joint` components 在 `primary_process_label` 中被归为：

```text
pan_event = 33
partial_or_inactive = 19
```

因此，P12 当前三分类图里几乎看不到 E10gb1 的 ripple-gamma activity，并不代表 E10gb1 没有 ripple-gamma components。

后续正式文档必须区分两层 label：

### 4.1 selectivity label

用于保守判断 clean component：

```text
theta_selective
ripple_gamma_selective / ripple_gamma_no_pure_theta
other
```

### 4.2 activity/subprocess label

用于判断是否存在 subprocess activity：

```text
theta_active
theta_selective
ripple_active
gamma_active
ripple_gamma_joint
ripple_gamma_no_theta_active
ripple_gamma_no_theta_selective
pan_event
partial_or_inactive
```

正式的新 dataset consistency 分析不能只看 selectivity label。必须同时看 subprocess activity label。

## 5. 当前可以保留的结论

### 5.1 比较稳的结论

BLP eigenfunction-derived density 比 event density 更有竞争力。

在 7-dataset P12 topN 中，`event_density` 的 topN membership 很低；`raw_efun_density` 和 `dimred_efun_density` 更常进入 top hits。

这个结论是目前最稳的 cross-dataset pattern。

### 5.2 中等可信、需要分 dataset/observable 验证的结论

E10gb1 支持：

```text
BOLD efun
    更像 slow/state variable

BOLD deconv_efun
    更像 intrinsic perturbation / fast trigger
```

但 pooled 7-dataset P12 不应该直接用来宣布这个结论完成。必须按：

```text
dataset
BOLD observable
pipeline P8/P10
topN
raw vs dimred
```

分开看。

### 5.3 目前不能直接宣布的结论

不能只根据当前 P12 `primary_process_label` 图说：

```text
7 datasets 没有 ripple-gamma components
```

这个说法不成立。当前 P12 图低估了 ripple-gamma activity，因为它没有画 activity/joint label。

更准确的说法是：

```text
clean ripple-gamma-selective top hits 在当前 strict primary label 下不强；
但 ripple-gamma-active/joint components 需要用 activity-label 层重新统计。
```

## 6. 正式新文档应该继承的分析顺序

后续正式设计文档建议按以下顺序组织。

### Step 0. Availability / provenance audit

先确认：

- P4 full EDMD chunks 是否存在。
- P5 standardized csplit reduction/density/peak stats 是否完整。
- P5 label 表是否使用 RMS-envelope adaptive。
- P8/P10 是否绑定 current-best P7/P9，不混 stale runs。
- P10 是否有 current-best P9 `efun_real` 和 `deconv_real`。

### Step 1. P5 subprocess existence

先问新 dataset 自己是否有候选 subprocess：

```text
同一个 method/k 下是否同时存在：
1. theta-like component
2. ripple-gamma-like component
```

这里必须同时报告：

- strict selectivity label；
- relaxed activity/joint label。

### Step 2. Density source competition

比较：

```text
event_density
raw_BLP_efun_density
dimred_BLP_efun_density
```

不要只看 mean top5。要看：

- winner source；
- topN membership；
- best-hit delta；
- candidate-count caveat。

### Step 3. BOLD efun vs deconv_efun split

分别看：

```text
P8 efun
P8 deconv_efun
P10 efun
P10 deconv_efun
```

重点：

- raw efun index/timescale distribution；
- topN sensitivity；
- top-hit time series；
- lag direction。

### Step 4. Dimred top-hit subprocess label

对 top dimred BLP density hits，join P5 labels。

必须同时画：

1. strict primary selectivity composition；
2. relaxed subprocess activity composition。

### Step 5. BOLD observable-specific interpretation

不能把 BOLD observable 混起来。

至少分开：

```text
HP_svd100
global_svd100 / gsvd100_ds
roi_mean
```

解释时需要允许：

- global/state-like BOLD efun 更偏 slow/theta；
- HP/local BOLD efun 可能直接对 local fast/ripple-gamma；
- deconv_efun 更偏 residual/trigger，但仍需 dataset 内验证。

### Step 6. Cross-dataset consistency score

正式文档需要定义新 dataset 和 E10gb1 的一致性不是数值一样，而是 pattern 一样。

候选评分维度：

- P5 是否有 two-subprocess candidates。
- P8/P10 是否都显示 raw/dimred density 优于 event density。
- P8/P10 是否都有 efun vs deconv_efun 分化。
- dimred top hits 是否有 theta/ripple-gamma subprocess label。
- raw top hits 是否显示 slow/fast split。
- 这些 pattern 是否在多个 BOLD observable 上可解释。

### Step 7. Second-stage ROI/activation review

只对少数通过 Step 1-6 的 candidate groups 看：

- BOLD activation map；
- ROI summary；
- cross-dataset ROI profile correlation；
- top ROI overlap。

这一步不应该替代前面的 numeric/process consistency。

## 7. 下一份正式文档需要修正的内容

1. P12 设计文档里的 `ripple-gamma` label 部分需要更新，不能只依赖 `primary_process_label`。
2. P12 图需要新增 activity/joint label panels。
3. E10gb1 anchor pattern 应该明确拆成：
   - P5 subprocess existence；
   - P8/P10 density competition；
   - raw timescale split；
   - dimred top-hit label；
   - top-hit time-series sanity check。
4. 新 dataset 不能只看 pooled 7-dataset 图；必须先做 per-dataset report，再做 cross-dataset summary。
5. cross-dataset consistency 应该是 “pattern matching”，不是 top xcorr 数值完全相等。

