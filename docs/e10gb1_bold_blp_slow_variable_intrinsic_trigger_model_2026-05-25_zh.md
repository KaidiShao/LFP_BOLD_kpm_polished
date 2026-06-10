# E10gb1 BOLD-BLP slow variable / intrinsic trigger 模型记录

日期：2026-05-25

状态：这是基于 E10gb1 条件比较得到的当前工作性解释模型，还不是
cross-session 最终结论。

## 1. 核心问题

现在真正关心的问题不是单纯问：

```text
哪个 BLP density 和 BOLD efun 的 xcorr 最大？
```

而是问：

```text
BLP Koopman eigenfunction / dimred eigenfunction 里面，
是不是存在可以解释 BOLD 不同动力学成分的 subprocess？
```

目前最重要的区分是：

- `BOLD efun x BLP density`
- `BOLD deconv efun x BLP density`

现在的解释是：

- `BOLD efun` 更像 BOLD 系统里的慢变量、latent state、state variable。
- `BOLD deconv efun` 更像从 BOLD 慢动力学里分离出来的 residual /
  intrinsic perturbation。

因此如果 BLP density 能分别和这两类 BOLD-side target 对上，就可能得到一个
很漂亮的 subprocess 模型：

```text
theta-like slow BLP subprocess
    -> BOLD slow state variable

ripple-gamma-like fast BLP subprocess
    -> BOLD residual / intrinsic trigger
```

这个解释不是只从 xcorr 数字来的，而是因为我们画了 top xcorr pair 的实际
时间序列之后，slow/state-like 和 fast/activity-like 的差别变得很清楚。

## 2. 比较过的四个 BLP condition

E10gb1 目前比较的是这四个 BLP 条件：

```text
abs_projected_vlambda
complex_split_projected_vlambda
abs_projected_vlambda_standardize
complex_split_projected_vlambda_standardize
```

也就是：

```text
BLP observable type: abs vs complex split
BLP normalization: raw vs standardized
```

可以写成一个 2 x 2：

| observable | normalization | condition |
| --- | --- | --- |
| abs | raw | `abs_projected_vlambda` |
| complex split | raw | `complex_split_projected_vlambda` |
| abs | standardized | `abs_projected_vlambda_standardize` |
| complex split | standardized | `complex_split_projected_vlambda_standardize` |

目前最有科学解释力、也最值得继续追的 condition 是：

```text
complex_split_projected_vlambda_standardize
```

原因是它更容易在低维 dimred component 里同时看到：

- theta-like component
- ripple-gamma-like component

也就是说，它更有可能把两个 subprocess 分开。

## 3. 当前 activity / density 口径

当前 P5/P8/P10 相关 density 使用的是 RMS-envelope adaptive activity
口径，不再直接用原始 eigenfunction 正负值做 threshold。

原因是 eigenfunction 有 sign / polarity ambiguity。如果直接 threshold
原始 eigenfunction，sign 翻转会改变 active/inactive 的方向。

当前工作规则是：

```text
activity(t, mode) = sqrt(movmean(abs(phi(t, mode)).^2, window_samples))
```

window 根据 eigenvalue 对应的 timescale 自适应决定。

重要修正：

raw efun 的 slow/fast timescale 现在优先用 discrete Koopman eigenvalue：

```text
tau = -dt / log(abs(lambda_discrete))
```

之前一些早期图里使用过 bilinear continuous estimate，这个现在只作为
legacy / fallback provenance，不作为主要 slow/fast 解释依据。

相关文档：

```text
docs\eigenfunction_activity_magnitude_policy_2026-05-22.md
docs\rmsenv_adaptive_p5_p8_p10_run_log_20260523.md
```

## 4. 实验一：P5 dimred selectivity 比较

目的：

看每个 BLP condition 下，P5 dimred eigenfunction component 是否能分出
theta subprocess 和 ripple-gamma subprocess。

比较范围：

```text
condition: abs / csplit x raw / standardized
method: SVD / NMF / MDS / UMAP
k: 03:08
```

主要输出：

```text
results\e10gb1_p5_two_subprocess_candidates\summary.md
results\e10gb1_p5_two_subprocess_candidates\method_k_two_subprocess_candidates.csv
results\e10gb1_p5_two_subprocess_candidates\component_subprocess_flags.csv
```

主要脚本：

```text
scripts\export_e10gb1_standardized_edmd_outputs.py
scripts\script_run_e10gb1_standardized_p5_ripple_probe.m
scripts\compare_e10gb1_p5_standardize_ripple_selectivity.py
scripts\compare_e10gb1_p5_relaxed_ripple_gamma.py
```

### 4.1 selectivity 定义的变化

一开始看过严格 ripple-only selectivity，但这个定义太窄。

现在更合理的定义是两个 subprocess：

```text
theta subprocess
ripple-gamma subprocess
```

新的解释规则：

- ripple 和 gamma 同时 active 是可以接受的。
- sharp-wave-ripple 不能当作纯 ripple-only，因为它可能包含 theta/gamma。
- ripple-gamma component 最重要的是不要在 pure theta 上 active。
- gamma 本身不是必须单独解释的目标，但 gamma 和 ripple 一起出现是有意义的。

### 4.2 P5 主要发现

在 E10gb1 里，`complex_split_projected_vlambda_standardize` 最有希望。

它在很多 method/k 下都能同时找到：

- 一个 theta-like component
- 一个 ripple-gamma no-pure-theta component

比较紧凑的例子：

```text
standardized csplit, SVD k03:
    theta component = 1
    ripple-gamma no-pure-theta component = 2

standardized csplit, NMF k03:
    theta component = 2
    ripple-gamma no-pure-theta component = 1
```

更广泛地看，standardized csplit 下有两 subprocess candidate 的范围包括：

```text
SVD:  k03-k08
NMF:  k03-k08
MDS:  k04-k08
UMAP: k05-k08
```

解释：

standardization 似乎增强了 fast / ripple-gamma-like component 的可见性，
而 complex split 又保留了足够的相位/复数结构信息，所以两者结合比较有用。

## 5. 实验二：raw BLP efun density 的 P8/P10 xcorr 分布

目的：

比较 raw BLP efun density 和 BOLD efun / BOLD deconv efun 的 xcorr。

重点不是只看 top5 平均值，而是看：

- top1/top3/top5/top10/top20/top50 分布
- standardized vs nonstandard
- efun vs deconv_efun
- raw efun index / timescale distribution

主要输出：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_xcorr_distribution_std_vs_nonstandard_long.csv
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_xcorr_distribution_std_vs_nonstandard_summary.csv
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_xcorr_topN_compact_interpretation.csv
```

主要图：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_ecdf_std_vs_nonstandard_top50.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_topN_median_std_vs_nonstandard.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_topN_delta_std_minus_nonstandard.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_threshold_fraction_by_topN.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_exact_rank_profile_std_vs_nonstandard.png
```

### 5.1 raw density xcorr 的现象

目前 E10gb1 csplit nonstandard vs standardized 的结果大概是：

- P8 efun：standardized 比 nonstandard 弱，尤其 top1/top3/top5。
- P8 deconv_efun：standardized 和 nonstandard 接近。
- P10 efun：standardized 更强。
- P10 deconv_efun：standardized 更弱。

这说明 xcorr 大小本身不能单独解释。

更重要的是：

```text
efun 和 deconv_efun 可能对应不同机制。
```

## 6. 实验三：raw top5 xcorr pair 的时间序列

这是这次最关键的可视化之一。

目的：

直接看 top5 xcorr pair 的时间序列，确认高 xcorr 到底来自什么形态。

主脚本：

```text
scripts\plot_e10gb1_raw_csplit_top5_xcorr_time_series.m
```

输出表：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_top5_xcorr_timeseries_rows.csv
```

输出图目录：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_top5_xcorr_timeseries\
```

图包括：

```text
raw_csplit_top5_xcorr_timeseries__nonstandard__P8__efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P8__efun.png
raw_csplit_top5_xcorr_timeseries__nonstandard__P8__deconv_efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P8__deconv_efun.png
raw_csplit_top5_xcorr_timeseries__nonstandard__P10__efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P10__efun.png
raw_csplit_top5_xcorr_timeseries__nonstandard__P10__deconv_efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P10__deconv_efun.png
```

图的读法：

- 蓝线：raw BLP efun density
- 红线：对应的 BOLD-side target，按 peak lag shift 后画出
- 两条线都 z-score
- 每张图是一个 condition x pipeline x feature family 的 top5 xcorr pair

### 6.1 raw top5 时间序列结论

`BOLD efun x raw BLP density`：

- 常常对应 slow state / session-wise drift。
- 图上更像一个慢变量，或者 session 级别的状态变化。
- 这个 slow component 和 theta/state-like 解释比较吻合。

`BOLD deconv efun x raw BLP density`：

- 更容易对应 fast activity / burst-like activity。
- 这更像 residual 里的 intrinsic perturbation。
- 这个 fast component 和 ripple-gamma / trigger-like 解释比较吻合。

这个图让当前模型从“top xcorr 数字”变成了“可以被时间序列支持的机制解释”。

## 7. 实验四：dimred BLP efun density 的 label 检查

目的：

确认 raw efun 里看到的 slow/fast 区分，在 dimred efun density 里是否也存在。

主要表：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\top_dimred_rows.csv
results\e10gb1_standardized_csplit_p8_p10_probe\best_dimred_method_k_rows.csv
results\e10gb1_standardized_csplit_p8_p10_probe\all_top_rows.csv
```

关键图：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\p8_p10_dimred_two_subprocess_label_counts_by_family_top50.png
```

### 7.1 dimred label 结果

dimred top hits 也显示了 feature-family 分离：

- `efun` top hits 里有明显的 `theta_strict`，但也会混入
  `ripple_gamma_no_pure_theta` 和 `other`。
- `deconv_efun` top hits 很集中地落在
  `ripple_gamma_no_pure_theta`。

这说明：

```text
raw efun 的 slow/fast 现象不是 raw mode 偶然结果；
dimred 之后仍然保留了这个 subprocess 区分。
```

## 8. 实验五：dimred top5 xcorr pair 的时间序列

目的：

直接画 dimred BLP efun density 的 top5 xcorr pair，检查它们是否支持同样的
slow variable / intrinsic trigger 模型。

主脚本：

```text
scripts\plot_e10gb1_dimred_csplit_top5_xcorr_time_series.m
```

输出表：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\dimred_csplit_top5_xcorr_timeseries_rows.csv
```

输出图目录：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\dimred_csplit_top5_xcorr_timeseries\
```

图包括：

```text
dimred_csplit_top5_xcorr_timeseries__P8__efun.png
dimred_csplit_top5_xcorr_timeseries__P8__deconv_efun.png
dimred_csplit_top5_xcorr_timeseries__P10__efun.png
dimred_csplit_top5_xcorr_timeseries__P10__deconv_efun.png
```

### 8.1 dimred top5 结果

P8 efun top5：

- 全部是 `mds k08 component 3`
- label = `theta_strict`
- 时间序列视觉上是 slow/state-like

P10 efun top5：

- 前三个是 `mds k08 component 3`
- label = `theta_strict`
- 后两个是 `nmf k03 component 1`
- label = `ripple_gamma_no_pure_theta`
- 解释：efun 侧以 theta/state 为主，但可以混入 ripple-gamma 成分。

P8 deconv_efun top5：

- 主要是 `umap k03/k04 component 2` 和 `mds k03 component 3`
- label = `ripple_gamma_no_pure_theta`
- 这边非常干净地偏向 ripple-gamma。

P10 deconv_efun top5：

- 主要是 `umap k03/k04 component 2`
- label = `ripple_gamma_no_pure_theta`
- 同样支持 intrinsic trigger / fast subprocess 解释。

### 8.2 dimred 结论

dimred BLP efun density 明确支持：

```text
BOLD efun        -> theta/state-like slow subprocess
BOLD deconv efun -> ripple-gamma-like fast subprocess
```

其中 deconv 侧比 efun 侧更干净。

## 9. 当前最漂亮的解释模型

当前模型可以写成：

```text
theta-like slow BLP subprocess
    -> tracks / modulates BOLD latent state variable
    -> captured by BOLD efun x BLP density

ripple-gamma-like fast BLP subprocess
    -> aligns with BOLD residual / intrinsic perturbation
    -> captured by BOLD deconv efun x BLP density
```

换句话说：

```text
BOLD efun:
    slow variable
    state-like
    theta-related

BOLD deconv efun:
    residual / intrinsic perturbation
    trigger-like
    ripple-gamma-related
```

这比单纯说 “BOLD efun 对应 slow，BOLD deconv 对应 fast” 更好。

因为这里可以进一步解释为：

```text
一个 subprocess 是 slow variable；
另一个 subprocess 是 intrinsic trigger。
```

## 10. 为什么不能把 efun 和 deconv_efun 混在一个总排名里

如果把 `efun` 和 `deconv_efun` 混在一起做一个总的 xcorr ranking，
会把两个不同机制混掉：

```text
slow state-variable coupling
fast residual-trigger coupling
```

所以 P8/P10 后续应该首先分层：

```text
efun        -> slow variable / theta-state analysis
deconv_efun -> intrinsic trigger / ripple-gamma analysis
```

在这个分层之后，再比较：

- BLP observable：abs vs complex split
- BLP normalization：raw vs standardized
- dimred method：SVD / NMF / MDS / UMAP
- component number：k03-k08
- BOLD observable type
- raw density vs dimred density

## 11. 当前还不能过度解释的地方

1. 目前最强证据来自 E10gb1。

2. P5 selectivity 已经比较了四个 BLP condition，但 P8/P10 的最细时间序列
   检查目前主要集中在 csplit nonstandard vs csplit standardized。

3. top5 里有一些重复，比如同一个 BLP component 对多个 BOLD mode 或多个
   method/k setting 都出现。这说明结果可能稳定，但后续也需要 de-duplicate
   by density component。

4. 高 xcorr 不等于因果。当前模型是 phenomenological / mechanistic
   interpretation。

5. standardized csplit 很有希望，但还需要跨 dataset 验证后才能作为主线选择。

## 12. 下一步建议

### 12.1 cross-session 验证

对每个 dataset 检查：

```text
efun top hits 是否富集 theta/state label
deconv_efun top hits 是否富集 ripple-gamma label
```

这是最关键的下一步。

### 12.2 P11 需要加入的 summary

P11 后续应该加入一个专门的 slow variable / intrinsic trigger section：

- label composition by BOLD-side target
- topN xcorr strength by BOLD-side target
- raw efun timescale distribution by BOLD-side target
- dimred component label distribution by BOLD-side target
- top5 或 de-duplicated top component 时间序列快照

### 12.3 四个 condition 的完整 P8/P10 对比

后续应该对四个 condition 做同等深度的 P8/P10 检查：

```text
abs raw
csplit raw
abs standardized
csplit standardized
```

但是输出必须按 BOLD-side target 分开：

```text
efun
deconv_efun
```

不要混在一个总图里。

### 12.4 参数选择应该怎么看

这个模型可以帮助选择参数：

- 如果目标是 slow state / theta variable：
  - 看 `BOLD efun x BLP dimred density`
  - 关注 theta label、slow timescale、cross-session consistency

- 如果目标是 intrinsic trigger / ripple-gamma perturbation：
  - 看 `BOLD deconv efun x BLP dimred density`
  - 关注 ripple-gamma label、fast timescale、burst-like time series

理想参数应该能在同一个 condition 里同时得到：

```text
一个稳定 theta/state component
一个稳定 ripple-gamma/trigger component
```

目前 E10gb1 里最值得继续追的是：

```text
complex_split_projected_vlambda_standardize
```

尤其是低 k 的 SVD/NMF，以及在 top xcorr 中表现很强的 MDS/UMAP components。

## 13. 实验六：四种 BOLD observable 分开检查

### 13.1 为什么要单独看 BOLD observable

前面的模型把 `BOLD efun` 和 `BOLD deconv_efun` 分开以后已经清楚很多，
但是后续又发现一个重要问题：

```text
BOLD observable type 本身也不能混在一起解释。
```

尤其是 `HP_svd100` 和其他几个 observable 不太一样。`HP_svd100` 更偏 local /
high-pass BOLD 表征；而 `global_svd100`、`gsvd100_ds`、`roi_mean` 更偏 global /
state-like 或 ROI-summary 表征。因此 “BOLD efun 是 slow variable” 这个说法不能
不分 observable 地套用。

另外，当前 P8/P10 结果里的 `pv` 不是新的生物学 observable，它只是
`projected_vlambda` 的 run-tag 缩写。当前分析里实际对应的 BOLD observable 是：

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

### 13.2 新增脚本和输出

主脚本：

```text
scripts\plot_e10gb1_observable_separated_subprocess_probe.m
```

输入表：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\all_top_rows.csv
```

输出表：

```text
results\e10gb1_standardized_csplit_p8_p10_probe\observable_separated_dimred_label_summary_top50.csv
results\e10gb1_standardized_csplit_p8_p10_probe\observable_separated_raw_timescale_summary_top50.csv
results\e10gb1_standardized_csplit_p8_p10_probe\observable_separated_dimred_top5_rows.csv
```

输出图：

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\observable_separated_subprocess_probe\observable_separated_dimred_label_fraction_top50.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\observable_separated_subprocess_probe\observable_separated_dimred_xcorr_strength_top50.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\observable_separated_subprocess_probe\observable_separated_raw_timescale_top50.png
```

这里主要看 top50 的 label composition、xcorr strength 和 raw efun timescale。
top1/top5 仍然有解释价值，但 top50 更适合看一个 observable 的整体倾向。

### 13.3 当前主要发现

#### A. `BOLD deconv_efun -> ripple-gamma / intrinsic trigger` 大体成立

这个结论在很多 observable 里仍然成立，尤其是：

```text
P8 deconv_efun:
    gsvd100_ds
    HP_svd100
    roi_mean

P10 deconv_efun:
    HP_svd100
```

这些条件下，高 xcorr 的 dimred BLP density 大多是
`ripple_gamma_no_pure_theta` label。这个支持之前的解释：

```text
BOLD deconv_efun 更像 residual / intrinsic perturbation；
它更容易对上 ripple-gamma-like fast subprocess。
```

但有一个重要例外：

```text
P8 global_svd100 deconv_efun
```

它在当前 top50 里反而 theta / other 成分更多，所以不能把 deconv 侧简单写成
“所有 observable 都是 ripple-gamma”。

#### B. `BOLD efun -> slow theta/state variable` 只对部分 observable 成立

这个结论最清楚地出现在：

```text
global_svd100 efun
```

它更符合 slow / theta / state-like variable 的解释。

但是：

```text
HP_svd100 efun
```

并不像 slow state variable。它更偏 local / fast / ripple-gamma-like。
所以 HP 不能和 global observable 混着看。

`roi_mean` 和 `gsvd100_ds` 处于中间状态：

- P8 里有时 top1 仍然是 theta/state-like。
- 但在 top50 或 P10 里，ripple-gamma 成分可以明显增加。

因此现在更准确的说法是：

```text
global_svd100 / 部分 roi_mean 的 efun:
    更适合解释 slow theta/state variable

HP_svd100 的 efun:
    更适合解释 local / fast / ripple-gamma-like mode

deconv_efun:
    多数情况下更适合解释 intrinsic trigger，
    但仍要按 observable 分开检查。
```

#### C. raw efun timescale 也要按 observable 分开解释

raw efun timescale 的 top50 summary 显示，median top50 timescale 不等于 top1
timescale。某些 `global_svd100 efun` top1 hit 可以是很慢的 mode，但 top50 median
会因为混入大量 fast / repeated rows 而变短。

所以后续不能只看一个总的 median：

```text
需要同时看：
1. top1 / top5 的代表性 slow mode
2. top10 / top20 / top50 的整体分布
3. 按 BOLD observable 分层的分布
```

### 13.4 对工作模型的修正

之前的简化模型是：

```text
BOLD efun        -> slow theta/state variable
BOLD deconv_efun -> fast ripple-gamma/intrinsic trigger
```

现在更准确的模型是：

```text
global/state-like BOLD observable 的 efun:
    更可能对应 slow theta/state variable

local/HP BOLD observable 的 efun:
    不一定是 slow variable；
    它可以直接对应 local fast / ripple-gamma-like activity

BOLD deconv_efun:
    更常对应 residual / intrinsic perturbation / trigger，
    但也必须按 observable 单独验证
```

因此后续 P8/P10 不能只做一个总图。至少要分三层：

```text
pipeline: P8 vs P10
BOLD feature family: efun vs deconv_efun
BOLD observable: global_svd100 / gsvd100_ds / HP_svd100 / roi_mean
```

在这个分层之后，才比较：

```text
BLP observable: abs vs csplit
BLP normalization: raw vs standardized
dimred method: SVD / NMF / MDS / UMAP
component number: k03-k08
```

### 13.5 目前最重要的解释结论

现在可以保留的结论不是一句简单的 “efun 慢、deconv 快”，而是：

```text
不同 BOLD observable 把 BOLD dynamics 投影到了不同层面。

global/state-like BOLD efun 更容易看到 slow theta/state subprocess。

HP/local BOLD efun 和很多 deconv_efun 更容易看到 ripple-gamma / intrinsic-trigger subprocess。

complex_split_projected_vlambda_standardize 仍然是 E10gb1 里最值得继续追的 BLP condition，
因为它在 P5 里能同时给出 theta component 和 ripple-gamma component，
并且在 P8/P10 里可以把这两类 BOLD-side target 分开解释。
```

这个版本比之前更严谨：它把 BOLD observable type 本身也纳入了解释，而不是默认所有
BOLD efun 都代表同一种 slow state。
