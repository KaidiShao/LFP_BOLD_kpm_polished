# 当前 SOP 图设计候选清单

日期：2026-06-10

本文件不是执行 SOP，而是图设计审查清单。目的：

```text
把当前 big picture、Koopman event narrative 和最小 SOP
翻译成“到底要画哪些图来支撑每一层证据”。
```

相关背景文档：

```text
notes/koopman_lfp_bold_framework.md
notes/koopman_event_narrative.md
docs/current_minimal_cross_dataset_consistency_sop_2026-06-10_zh.md
docs/roi_real_mean_mainline_rationale_2026-06-10_zh.md
```

当前核心解释：

```text
LFP/BLP 侧：
    standardized complex-split BLP dimred efun density
    应该能拆出 theta-like slow subprocess
    和 ripple/gamma-no-theta fast subprocess。

BOLD 侧：
    roi_mean BOLD efun 更像 global slow state variable / x^(g) readout。
    roi_mean BOLD deconv_efun 更像 BOLD modal innovation / eta^B readout。

LFP-BOLD coupling：
    如果 RG-no-theta subprocess density 更强地耦合 deconv_efun，
    则支持 fast hippocampal subprocess 作为
    intrinsic perturbation to global BOLD dynamics 的解释。
```

下面每一节都列：

```text
问题
建议图
已有图 / 数据
建议状态
```

---

## 1. P5 subprocess gate：先证明 LFP 侧能拆出两个 subprocess

### 正式代码入口

P5 subprocess gate 的正式最小代码已经单独整理到：

```text
scripts/pipeline5_sop/
```

一键入口：

```text
scripts/pipeline5_sop/run_p5_sop.py
```

默认口径：

```text
condition = complex_split_projected_vlambda_standardize
transform = adaptive_envelope
k = 03:16
strict labels = P2 theta/gamma/ripple band-event selectivity
```

默认只生成 P5 gate 主图、P5-pair-pass method-k diagnostic plate，以及可选的
theta/RG density correlation SOP 图；不默认生成 abs sensitivity、top30 window QC
或 consensus trajectory showcase。

### 问题

```text
不先看 BOLD，只看 P5：
standardized complex-split BLP dimred space
是否能稳定找到 theta-like slow subprocess
和 ripple/gamma-no-theta fast subprocess？
```

现在不需要重新纠结传统 event family。P5 band-event strict label 已经用 P2
theta/gamma/ripple windows 和 non-event baseline 做了 event-vs-baseline
effect z，因此能回答：

```text
这个 dimred component 是否对 theta / gamma / ripple band events
有相对 baseline 的显著或强响应？
```

### 建议图

主图候选：

```text
1. strict theta-vs-ripple/gamma effect scatter
2. strict label grid by method-k
3. strict label composition by method-k
4. strict two-subprocess candidate map
```

当前 SOP 默认只画：

```text
transform = adaptive_envelope
```

`abs` 不进默认 SOP，只在需要做 sensitivity / debugging 时手动补画。

图类型和文件对应关系：

```text
1. strict theta-vs-ripple/gamma effect scatter
   文件名模式：
       04_<transform>_strict_theta_vs_ripple_gamma_effect_scatter.png
   当前例子：
       04_adaptive_envelope_strict_theta_vs_ripple_gamma_effect_scatter.png
   作用：
       每个点是一个 dimred component。
       横轴 theta event effect z，
       纵轴 max(gamma, ripple) event effect z。
       用来看 theta-like 和 RG-like component 是否在二维 effect space 中分开。

2. strict label grid by method-k
   文件名模式：
       02_<transform>_strict_label_grid_by_method_k.png
   作用：
       每个 method-k 下面逐个 component 显示 strict label：
       theta_selective / gamma_selective / ripple_selective /
       ripple_gamma_no_theta / mixed_or_partial / inactive。
       用来看每个 method-k 的 component-level label 结构。

3. strict label composition by method-k
   文件名模式：
       03_<transform>_strict_label_composition_by_method_k.png
   作用：
       每个 method-k 的 strict label 组成比例。
       用来看某个 method-k 是 clean subprocess 多，
       还是 mixed / inactive 多。

4. strict two-subprocess candidate map
   文件名模式：
       05_<transform>_strict_two_subprocess_candidate_map.png
   当前例子：
       05_adaptive_envelope_strict_two_subprocess_candidate_map.png
   作用：
       每个格子是一个 method-k。
       标记是否同时有 theta_selective 和 ripple_gamma_no_theta component，
       并写出代表 component index，例如 T:C1, RG:C3。
       这是 P5 subprocess gate 的最直接 pass map。
```

解释图候选：

```text
5. selected method-k diagnostic plate
```

解释图和文件对应关系：

```text
5. selected method-k diagnostic plate
   文件名模式：
       <dataset>_<transform>_<method>_k<kk>_thetaC<cc>_rgC<cc>_*diagnostic_plate.png
   作用：
       左边画 selected theta/RG component 的 peri-event mean +/- ribbon；
       右边画 top windows，并叠加 P2 theta/gamma/ripple event rugs。
       用于 P5-pair-pass method-k 的机制解释和 QC。
       当前默认只画 adaptive_envelope，并且只画同时有
       theta_selective 与 ripple_gamma_no_theta / RG-like candidate 的 method-k。
```

showcase only，不进默认 SOP：

```text
6. consensus state trajectory in two-subprocess coordinates
   x-axis = selected theta-like dimred component
   y-axis = selected RG-no-theta dimred component
   z-axis = residual / leftover variability / optional third component
```

showcase 图对应关系：

```text
6. consensus state trajectory in two-subprocess coordinates
   文件名模式：
       待标准化；不作为当前默认 SOP 输出。
   作用：
       用 selected theta-like component 和 selected RG-no-theta component
       作为前两个坐标轴来画 consensus-state trajectory。
       只用于最后挑 1--2 个漂亮 dataset 做机制展示。
```

这个 trajectory 图适合最后挑 1--2 个漂亮 dataset 做机制展示，不适合默认每个
dataset 都画。

### 建议状态

```text
默认 SOP 必画：
    1, 2, 3, 4
    transform = adaptive_envelope only

解释 / supplement：
    5
    只对 P5-pair-pass method-k 画 adaptive_envelope diagnostic plate

showcase only：
    6 consensus trajectory
```

---

## 2. P5 theta/RG density correlation：证明不是同一个 subprocess

```text
即使 P5 label 能找到 theta-like 和 RG-no-theta component，
这两个 density 会不会其实是同一个 time series 的不同标签？
```

默认 SOP 只保留一张主图：

```text
dataset x method-k canonical theta/RG density correlation heatmap
```

图的每个格子是一个 dataset × method-k 的 canonical theta/RG density correlation。
优先看 session-demeaned correlation。接近 0 表示 theta-like 和 RG-no-theta
density 在时间上相对独立；绝对值很高则说明这两个 labeled subprocess 在该
dataset/method-k 里可能没有真正分开。

正式代码入口：

```text
scripts/pipeline5_sop/run_p5_sop.py --refresh-density-correlation
scripts/pipeline5_sop/plot_theta_rg_density_correlation.py
```

---

## 3. Density-source competition：证明 subprocess 比 simple event 更像 latent mechanism

### 问题

```text
P8/P10 中，真正解释 BOLD efun/deconv_efun 的 density source 是什么？

simple event_density 是否已经足够？
raw BLP efun density 是否更好？
dimred BLP efun density 是否最好？
```

这是 big picture 里的关键检验：

```text
frequency-defined event density
vs raw Koopman eigenfunction activity density
vs low-dimensional Koopman subprocess density
```

如果：

```text
dimred density > raw density >> event density
```

则支持：

```text
P5 dimred subprocess representation
比传统 event count 更接近 BOLD-scale latent action variable。
```

默认 SOP 只看 `roi_mean` observable，并且只保留这四张图：

```text
01_density_source_competition_by_observable/roi_mean/
    p8_density_class_fraction__roi_mean__efun.png
    p8_density_class_fraction__roi_mean__deconv_efun.png
    p10_density_class_fraction__roi_mean__efun.png
    p10_density_class_fraction__roi_mean__deconv_efun.png
```

读图逻辑：

```text
event_density 低：
    simple event count 不足以解释 BOLD efun/deconv_efun。

raw / dimred BLP efun density 高：
    Koopman eigenfunction activity 比传统 event density 更像 latent mechanism。

dimred density 高于 raw density：
    P5 dimred subprocess representation 有额外价值。

deconv_efun 比 efun 更偏 dimred/RG-like source：
    支持 fast intrinsic perturbation / trigger readout。
```

---

## 4. roi_mean BOLD deconv_efun 与 dimred RG-like density 的关系

### 问题

```text
roi_mean BOLD 代表 global pattern。
efun 是 latent state readout。
deconv_efun 是 BOLD modal innovation / intrinsic perturbation readout。

当前主线更关心 event-triggered activities：
哪一类 LFP/BLP dimred subprocess density
更像 deconv_efun？
```

如果 deconv_efun 的 top xcorr hits 更偏 `ripple_gamma_no_theta`，
则支持：

```text
RG-like hippocampal subprocess
与 event-triggered global BOLD modal innovation 相关。
```

这可以连接到：

```text
ripple-triggered system reorganization
memory consolidation / replay-related perturbation
```

但不能写成 causal proof。

### 建议图

主图候选：

```text
1. roi_mean deconv-vs-efun RG enrichment by topN
2. roi_mean efun/deconv_efun strict-label fraction by topN
3. roi_mean method-k hit map for ripple_gamma_no_theta
4. roi_mean strict effect scatter for top hits
```

### 已有图 / 数据

CSV：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    p8_p10_strict_band_coupling/
        deconv_vs_efun_rg_enrichment.csv
        deconv_vs_efun_rg_enrichment_by_observable.csv
        deconv_vs_efun_rg_enrichment_matched_units.csv
        dimred_strict_label_topN_fraction.csv
        dimred_strict_label_topN_fraction_by_observable.csv
```

现有图：

```text
E:/DataPons_processed/summary_figures/pipeline11_current_analysis_summary/
standardized_csplit_k03_k16_all_current_20260607/p8_p10_strict_band_coupling/
    00_deconv_vs_efun_rg_enrichment_by_observable/roi_mean/
        p8_deconv_vs_efun_rg_enrichment__roi_mean.png
        p10_deconv_vs_efun_rg_enrichment__roi_mean.png

    02_dimred_strict_label_composition_by_observable/roi_mean/
        p8_dimred_strict_label_fraction__roi_mean__efun.png
        p8_dimred_strict_label_fraction__roi_mean__deconv_efun.png
        p10_dimred_strict_label_fraction__roi_mean__efun.png
        p10_dimred_strict_label_fraction__roi_mean__deconv_efun.png

    05_method_k_strict_label_hit_maps_by_observable/roi_mean/
        p8_ripple_gamma_no_theta_method_k_hits_top10__roi_mean__deconv_efun.png
        p8_ripple_gamma_no_theta_method_k_hits_top10__roi_mean__efun.png
        p10_ripple_gamma_no_theta_method_k_hits_top10__roi_mean__deconv_efun.png
        p10_ripple_gamma_no_theta_method_k_hits_top10__roi_mean__efun.png
```

### 建议状态

```text
默认 SOP 必画。

重点看 P8；P10 作为 robustness / extension。
```

---

## 5. BOLD efun 是否对应 global slow variable

### 问题

当前推测：

```text
BOLD efun      -> global slow state variable
BOLD deconv    -> fast event-triggered modal innovation
```

理想证据：

```text
和 BOLD efun xcorr 高的 dimred efun density 更偏 theta-like；
和 BOLD deconv_efun xcorr 高的 dimred density 更偏 RG-like。
```

但如果 dimred label 不能强烈支持 efun-theta，那么还有第二条证据：

```text
raw BLP efun density 与 BOLD efun/deconv_efun 的 top hits
是否来自不同 BLP timescale distribution？

P8 selected BOLD modes：
    efun selected modes 是否更慢、更靠前、更接近单位圆？
    deconv selected modes 是否更快、更靠后？
```

### 建议图

主图候选：

```text
1. P8 selected BOLD mode magnitude tau by feature
2. P8 selected BOLD mode magnitude tau by target
3. P8 selected sorted mode index by feature
4. P8 selected mode lambda plane
5. topN sensitivity curve for selected-mode tau
```

补充图：

```text
6. raw BLP efun density top hit timescale ECDF
7. raw BLP efun density topN median timescale by topN
8. raw BLP efun density whole-trace examples
```

### 已有图 / 数据

P8 selected BOLD mode timescale，P5-pair-pass strict subset：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_p8_selected_bold_mode_timescales_p5_pair_pass/
        p8_selected_bold_mode_timescales_long.csv
        p8_selected_bold_mode_timescale_summary.csv
        README_p8_selected_bold_mode_timescales.md
```

现有图：

```text
E:/DataPons_processed/summary_figures/pipeline11_current_analysis_summary/
standardized_csplit_k03_k16_all_current_20260607/
roi_mean_p8_selected_bold_mode_timescales_p5_pair_pass/
    01_p8_selected_mode_magnitude_tau_by_feature__top3.png
    01_p8_selected_mode_magnitude_tau_by_feature__top5.png
    01_p8_selected_mode_magnitude_tau_by_feature__top10.png
    01_p8_selected_mode_magnitude_tau_by_feature__top20.png
    02_p8_selected_mode_magnitude_tau_by_target__top10.png
    03_p8_selected_mode_sorted_index_by_feature__top10.png
    04_p8_selected_mode_lambda_plane__top10.png
    05_p8_selected_mode_tau_by_topn_feature.png
```

Raw BLP efun density timescale probe：

```text
results/bold_efun_raw_density_slow_state_probe_current_real_fixed_zoom_signed_aligned/figures/
    01_raw_topN_timescale_ecdf_by_pipeline_family.png
    02_raw_topN_median_timescale_by_topN.png
    03_raw_topN_slow_ge10s_fraction_by_topN.png
    04_raw_efun_index_hist_raw_top10.png
    05_p8_bold_vs_raw_density_whole_trace_examples_raw_top10__wide.png
```

### 建议状态

```text
默认 SOP 必画：
    selected BOLD mode tau/index/lambda/topN figures

补充：
    raw BLP efun density timescale figures

whole-trace examples：
    showcase only，不作为主证据
```

---

## 6. ROI real_mean：BOLD spatial footprint 是否区分 efun/deconv 或 theta/RG

### 问题

不能只看时间上的相似性，还要看 selected BOLD modes 的 spatial footprint。

当前 ROI 主线口径：

```text
ROI profile mode = real_mean
```

核心检验：

```text
4 target combinations:
    efun x theta
    efun x RG
    deconv_efun x theta
    deconv_efun x RG

是否 efun_theta 和 deconv_RG 最不一样？
```

如果成立，可以支持：

```text
slow global state readout
和 fast event-triggered innovation readout
不是同一个 spatial footprint。
```

如果不成立，也仍然可能说明：

```text
ROI 层面主要支持 efun vs deconv branch difference，
而不是 theta/RG spatial separation。
```

### 建议图

主图候选：

```text
1. real_mean 4-target pairwise ROI profile correlation matrix
2. real_mean 4-target pairwise top ROI Jaccard matrix
3. P8 selected BOLD mode exact overlap matrix
4. P8 selected BOLD mode adjacent/nearby overlap matrix
5. all-ROI-by-dataset signed profile, anatomical ROI order, wide layout
6. top positive / top negative ROI labels by target, derived from real_mean
```

### 已有图 / 数据

CSV：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_signed_target_pairwise_tests/
        pairwise_roi_contrasts.csv
        pairwise_selected_unit_overlap.csv
        target_profile_availability.csv
```

现有图，注意这些是旧版多 ROI mode 输出；现在主线只看 `real_mean`：

```text
E:/DataPons_processed/summary_figures/pipeline11_current_analysis_summary/
standardized_csplit_k03_k16_all_current_20260607/
roi_mean_signed_target_pairwise_tests/
    01_pairwise_roi_profile_corr__top3__P8__real_mean.png
    01_pairwise_roi_profile_corr__top5__P8__real_mean.png
    01_pairwise_roi_profile_corr__top10__P8__real_mean.png
    01_pairwise_roi_profile_corr__top20__P8__real_mean.png
    02_pairwise_top_roi_jaccard__top10__P8__real_mean.png
    03_pairwise_p8_selected_mode_exact_jaccard__top10.png
    04_pairwise_p8_selected_mode_adjacent_fraction__top10.png
```

### 建议状态

```text
默认 SOP 应该保留 real_mean pairwise matrix 和 selected-mode overlap。

但需要重新生成 real_mean-only 版本，
不要再把 mean_abs / positive_real / negative_real 作为并列 condition。

all-ROI-by-dataset signed profile 很可能很有用，
但应该作为 interpretive figure，不一定进第一版最小 SOP。
```

---

## 7. P10 robustness：BOLD dimred component number / dictionary size / topN

### 问题

P10 的解释比 P8 更复杂，因为 P10 经过 BOLD-side P9 reduction。

需要关注：

```text
1. P10 是否重复 P8 的方向：
   deconv_efun 更偏 RG-like dimred density？

2. P10 是否受 BOLD dimred component number 影响？

3. BOLD ResKoopNet dictionary size / checkpoint / P9 reduction setting
   是否影响结果？

4. 结论是否只在 top10 成立，还是 top3/top5/top20 也成立？
```

### 建议图

主图候选：

```text
1. P10 deconv-vs-efun RG enrichment by topN
2. P10 density-source competition by topN
3. P10 dimred strict-label fraction by topN
```

future / parameter sweep：

```text
4. P10 result vs BOLD-side P9 component number
5. P10 result vs BOLD ResKoopNet dictionary size
6. P10 result vs P7 checkpoint / best checkpoint version
```

### 已有图 / 数据

已有 P10 topN / label / density-source 图：

```text
E:/DataPons_processed/summary_figures/pipeline11_current_analysis_summary/
standardized_csplit_k03_k16_all_current_20260607/p8_p10_strict_band_coupling/
    00_deconv_vs_efun_rg_enrichment_by_observable/roi_mean/
        p10_deconv_vs_efun_rg_enrichment__roi_mean.png

    01_density_source_competition_by_observable/roi_mean/
        p10_density_class_fraction__roi_mean__deconv_efun.png
        p10_density_class_fraction__roi_mean__efun.png

    02_dimred_strict_label_composition_by_observable/roi_mean/
        p10_dimred_strict_label_fraction__roi_mean__deconv_efun.png
        p10_dimred_strict_label_fraction__roi_mean__efun.png
```

### 建议状态

```text
默认 SOP：
    P10 作为 robustness，不阻塞 P8 主结论。
    必须显示 topN sensitivity。

不应现在写死：
    BOLD P9 component number / dictionary size 已经充分验证。

需要后续补：
    P10 vs BOLD dimred component number / dictionary size 的 systematic sweep。
```

---

## 8. 我建议的最小默认图组

如果要极简、可跨数据重复，我建议默认 SOP 只保留这 8 类：

```text
LFP/P5:
1. P5 strict two-subprocess candidate map
2. P5 theta/RG density correlation heatmap

LFP density source:
3. P8/P10 event vs raw vs dimred density-source competition

BOLD coupling:
4. P8 roi_mean deconv-vs-efun RG enrichment by topN
5. P8 roi_mean strict-label composition by topN

BOLD slow vs fast evidence:
6. P8 selected BOLD mode timescale/index/lambda, topN sensitivity

Spatial sanity:
7. real_mean ROI 4-target pairwise matrix
8. P8 selected BOLD mode overlap matrix
```

Showcase / paper-figure candidates：

```text
1. two-subprocess consensus trajectory
2. selected method-k diagnostic plate
3. raw/dimred density vs BOLD efun/deconv whole-trace examples
4. all-ROI-by-dataset real_mean signed profile
```

暂时不建议默认画：

```text
1. every method-k paired top-window sheet
2. all observable panels beyond roi_mean
3. all ROI modes mean_abs / real_mean / positive / negative as parallel conditions
4. P10 component-index overlap as direct interpretability evidence
```

---

## 9. 当前需要你判断的图

我建议你先判断以下几类是否进入正式 SOP：

```text
1. P5 density correlation：
   只要 heatmap，还是也要 distribution？

2. density-source competition：
   roi_mean 四张 fraction 图是否够清楚？

3. P8 coupling：
   enrichment by topN 和 strict-label composition 是否都保留？

4. ROI：
   real_mean pairwise matrix 是否足够？
   是否还要 all-ROI-by-dataset wide signed profile？

5. showcase：
   two-subprocess consensus trajectory 是否只挑 e10gb1/f12m05 这种漂亮数据？
```
