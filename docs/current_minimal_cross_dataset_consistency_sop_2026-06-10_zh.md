# 当前 cross-dataset consistency 最小 SOP

日期：2026-06-10

本文件是当前执行版 SOP。它替代 2026-06-04 版作为后续新数据检查入口。
旧版和更长的历史记录仍保留为 provenance：

```text
docs/current_minimal_cross_dataset_consistency_sop_2026-06-04_zh.md
docs/current_minimal_cross_dataset_consistency_rationale_2026-06-04_zh.md
```

## 1. 当前 working hypothesis

现在的主叙事已经收敛到一条最小、可跨数据验证的链：

```text
standardized complex-split BLP dimred density
可以拆出 theta-like slow/state subprocess
和 ripple/gamma-no-theta fast subprocess。

roi_mean BOLD efun 更像 slow/global state readout。
roi_mean BOLD deconv_efun 更像 fast intrinsic perturbation readout。

P8/P10 中，deconv_efun 比 efun 更倾向于和
ripple_gamma_no_theta / RG-like fast subprocess coupling。

P8 selected BOLD mode eigenvalue/timescale 进一步显示：
efun selected modes 更慢、更靠前、更接近单位圆；
deconv_efun selected modes 更快、更靠后。
```

当前最稳的生物学表述是：

```text
BLP dimred density 中存在 theta-like slow state 和 RG-like fast subprocess；
BOLD efun 更像 slow/state readout；
BOLD deconv_efun 更像 fast intrinsic perturbation readout；
deconv branch preferentially couples to RG-no-theta subprocess。
```

可以作为 working hypothesis 联系到：

```text
ripple/gamma-triggered intrinsic perturbation
或 replay / memory-consolidation 相关 network reorganization
```

但现阶段不要写成因果结论。

## 2. 当前主线口径

第一阶段固定口径，不把所有历史版本、observable、density source 混在一起：

```text
BLP observable:
    complex_split_projected_vlambda_standardize

BLP activity / density:
    adaptive RMS-envelope 为主
    abs 只作为 sensitivity / QC

P5 label:
    只用 P2 theta/gamma/ripple band-event strict label
    暂时不使用 dominant label 作为主结论

BOLD observable:
    roi_mean

BOLD feature family:
    efun
    deconv_efun

P8/P10 density source:
    dimred_efun_density, csplit, standardized RMS-envelope adaptive
    raw_csplit_q070 作为辅助比较
    event_density 作为 baseline / negative control
```

核心 strict labels：

```text
theta_selective
ripple_gamma_no_theta
```

RG-like fast label 第一阶段以 `ripple_gamma_no_theta` 为主。
`gamma_selective` / `ripple_selective` 可以记录，但不要和主结论混在一起。

## 3. Dataset-level validation gates

每个新 dataset 不再直接混进总图，而是按同一套 gate 逐项判断。

最终每个 dataset 给出一个状态：

```text
P5 pair pass: yes / no / weak / missing
deconv RG enrichment: yes / weak / no / missing
efun slower than deconv: yes / weak / no / missing
ROI sanity: compatible / unclear / contradictory / missing
overall: supports / partial / fails / missing
```

### Gate 0: 数据和主线结果完整性

最低要求：

```text
P2 theta/gamma/ripple event detection exists
P5 standardized complex-split reduction exists
P5 adaptive RMS-envelope dimred density exists
P8 roi_mean xcorr exists
P7 roi_mean BOLD_POST / ROI profile exists
```

P10 是增强证据，不应该阻塞 P8 主结论。
P7 ROI profile 缺失时，先跳过 ROI sanity，不影响 P5/P8/timescale 主判断。

需要标注但不自动排除的问题：

```text
raw signal saturation / clipping
P2 event detection abnormality
P5 strict theta component missing
P8/P10 only event_density hits
P7 roi_mean BOLD_POST stale or not best checkpoint
```

### Gate 1: P5 subprocess gate

问题：

```text
这个 dataset 的 BLP dimred space 中，是否至少有一个 method-k
可以同时找到：

1. theta_selective component
2. ripple_gamma_no_theta component
```

严格 pass 单位：

```text
dataset x method_k
```

而不是只说整个 dataset pass。后续 P8 主分析应该优先只保留
P5-pair-pass 的 `dataset x method_k`。

主图：

```text
00_strict_label_counts_abs_vs_adaptive_envelope.png
04_<transform>_strict_theta_vs_ripple_gamma_effect_scatter.png
02_<transform>_strict_label_grid_by_method_k.png
03_<transform>_strict_label_composition_by_method_k.png
05_<transform>_strict_two_subprocess_candidate_map.png
```

解释图：

```text
selected component peri-event mean by band
method-k diagnostic plate:
    left: theta/RG selected components' event-aligned mean +/- ribbon
    right: top windows with theta-like and RG-like traces plus P2 event rugs
```

QC / appendix：

```text
01_<transform>_component_p2_band_event_effect_heatmap.png
P2 theta/gamma/ripple top30 event-window figures for suspicious datasets
```

必须保留的数值表：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    p8_p10_strict_band_coupling/
        p5_strict_band_label_background.csv
```

当前 P5-pair-pass 结果示例：

```text
通过 P5 pair-pass 的 dataset 包括：
e10gb1, e10fV1, e10gh1, e10gw1,
f12m01, f12m02, f12m03, f12m05,
k13m18, k13m23

若某些 K13 dataset 没有稳定 theta strict component，
不要继续把它们作为 theta/RG two-subprocess 支持证据。
```

### Gate 1b: P5 theta/RG density independence

问题：

```text
theta_selective density 和 ripple_gamma_no_theta density
是不是只是同一个 trace 的不同标签？
```

主图 / 表：

```text
canonical theta/RG density correlation
all theta/RG pair density correlation
session-demeaned density correlation
```

脚本：

```text
scripts/analyze_p5_theta_rg_density_correlation.py
```

结果目录：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    p5_theta_rg_density_correlation/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        p5_theta_rg_density_correlation\
```

解释：

```text
如果 theta/RG density 大多低相关，说明 P5 层面的 subprocess split 是实质性的。
如果某个 dataset theta/RG density 高相关，它只能作为 weak / ambiguous。
```

### Gate 2: P8/P10 coupling preference

问题：

```text
在 P5-pair-pass method-k 中，
roi_mean deconv_efun 的 top xcorr hits
是否比 roi_mean efun 更偏向 ripple_gamma_no_theta？
```

主图：

```text
roi_mean x efun/deconv_efun x topN strict-label composition
roi_mean deconv-vs-efun RG enrichment by topN
roi_mean x efun/deconv_efun density-source competition
```

主指标：

```text
RG fraction:
    topN hits 中 ripple_gamma_no_theta 的比例

delta:
    deconv RG fraction - efun RG fraction

matched-unit win rate:
    matched context 中 deconv RG fraction > efun RG fraction 的比例

density source competition:
    dimred density vs raw density vs event density
```

当前解释规则：

```text
pooled delta > 0:
    支持 aggregate deconv RG enrichment

pooled delta > 0 且 matched-unit win rate > 0.5:
    才能写成更稳定的 context-level preference

event_density dominant:
    不支持当前 dimred BLP subprocess 叙事

dimred density > raw density:
    支持 P5 subprocess representation 有额外价值
```

当前已观察到的方向：

```text
P8 deconv:
    dimred density > raw density >> event density

P10 deconv:
    dimred density > raw density >> event density

这支持“dimred BLP subprocess density 比 simple event density 更能解释
BOLD deconv_efun coupling”。
```

### Gate 3: P8 selected BOLD mode eigenvalue / timescale

这是当前最重要的新证据之一。

问题：

```text
P8 top xcorr selected BOLD modes 是否显示：
efun selected modes 更慢 / 更靠前 / 更接近单位圆；
deconv_efun selected modes 更快 / 更靠后？
```

只对 P8 直接解释 BOLD mode index：

```text
P8 bold_mode_index = P7 sorted BOLD Koopman mode
P10 bold_component_index = P9 reduced component
    不能直接当同一个 BOLD mode index
    需要另做 weighted source-mode timescale
```

timescale 定义：

```text
lambda = P7 BOLD_POST EDMD discrete Koopman eigenvalue
dt     = BOLD_POST/dt 或 EDMD_outputs/dt

magnitude_tau_sec = dt / abs(log(abs(lambda)))

decay_tau_sec:
    只对 abs(lambda) < 1 有纯 decay 解释

growth_tau_sec:
    只对 abs(lambda) > 1 有 growth 解释
```

由于部分 BOLD eigenvalue 略大于 1，主文中应写：

```text
near-unit / magnitude timescale
```

不要把所有 mode 都解释成 decay time constant。

主脚本：

```text
scripts/analyze_roi_mean_p8_selected_bold_mode_timescales.py
```

broad 版输出：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_p8_selected_bold_mode_timescales/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        roi_mean_p8_selected_bold_mode_timescales\
```

P5-pair-pass strict 版输出：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_p8_selected_bold_mode_timescales_p5_pair_pass/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        roi_mean_p8_selected_bold_mode_timescales_p5_pair_pass\
```

主图：

```text
01_p8_selected_mode_magnitude_tau_by_feature__top<N>.png
02_p8_selected_mode_magnitude_tau_by_target__top<N>.png
03_p8_selected_mode_sorted_index_by_feature__top<N>.png
04_p8_selected_mode_lambda_plane__top<N>.png
05_p8_selected_mode_tau_by_topn_feature.png
```

必须同时报告两种统计：

```text
hit-weighted:
    每个 top xcorr hit row 都计入
    回答“top hits 整体偏向什么 timescale”

unique selected BOLD mode:
    对 dataset x run_tag x bold_feature_family x bold_mode_index 去重
    target-level 对 dataset x run_tag x target x bold_mode_index 去重
    回答“真正被选中的 BOLD modes 本身偏向什么 timescale”
```

当前 strict P5-pair-pass 结果：

```text
unique mode, top5:
    efun tau        ~= 11.38 s
    deconv_efun tau ~= 4.22 s

target-level, top5 unique:
    efun_theta   ~= 64.43 s
    efun_RG      ~= 7.62 s
    deconv_theta ~= 5.00 s
    deconv_RG    ~= 3.88 s
```

当前解释：

```text
efun selected BOLD modes 比 deconv_efun selected modes 更慢；
efun_theta 是最慢的 slow-state branch；
deconv_RG 是最快的 fast perturbation branch。

这条证据比 ROI footprint 更适合作为
slow/state readout vs fast perturbation readout 的主证据。
```

### Gate 4: ROI spatial sanity

ROI 现在只作为 spatial sanity / falsification，不再作为主证据。

当前关键结论：

```text
mean_abs ROI profile 往往把差异抹平。
signed real / positive / negative ROI profile 能看到更多差异。

但是 ROI 差异主要反映 efun vs deconv_efun branch，
不稳定支持 theta vs RG 本身有完全独立 spatial footprint。
```

主脚本：

```text
scripts/export_p7_roi_mean_profiles_direct_from_bold_post.py
scripts/analyze_roi_mean_signed_target_modes.py
scripts/analyze_roi_mean_signed_target_pairwise.py
```

ROI value modes：

```text
mean_abs
real_mean
positive_real
negative_real
```

主图：

```text
efun_theta vs deconv_RG signed ROI profile corr by topN
efun_theta vs deconv_RG top ROI set Jaccard by topN
4-target pairwise ROI profile corr matrix:
    efun_theta
    efun_RG
    deconv_theta
    deconv_RG
4-target pairwise top ROI Jaccard matrix
P8 selected BOLD mode exact/adjacent overlap matrix
top positive/negative ROI labels by target
```

输出：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_signed_target_mode_tests/
    roi_mean_signed_target_pairwise_tests/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        roi_mean_signed_target_mode_tests\
        roi_mean_signed_target_pairwise_tests\
```

当前解释：

```text
同一个 BOLD feature family 内：
    efun_theta 和 efun_RG ROI profile 很像
    deconv_theta 和 deconv_RG ROI profile 很像

跨 BOLD feature family：
    efun vs deconv_efun signed ROI profile 更不同

因此 ROI 支持的是 efun/deconv branch difference，
不是 theta/RG spatial separation。
```

## 4. 需要画的正式图清单

### A. 每个 dataset 的 P5 gate 图

必画：

```text
strict theta-vs-ripple/gamma effect scatter
strict label grid by method-k
strict label composition by method-k
strict two-subprocess candidate map
```

可选：

```text
abs vs adaptive_envelope label counts
selected method-k diagnostic plate
P2 top30 event windows for suspicious data
```

### B. 每个 dataset 的 P8 coupling 图

必画：

```text
roi_mean efun/deconv_efun strict-label composition by topN
roi_mean deconv-vs-efun RG enrichment by topN
density-source competition by feature family
```

可选：

```text
whole-trace examples:
    raw density / dimred density vs BOLD efun/deconv_efun
    unsmoothed trace
    fixed x-axis
    session borders and masked regions marked
```

### C. 每个 dataset 的 P8 BOLD mode timescale 图

必画：

```text
selected mode magnitude tau by feature
selected mode magnitude tau by target
selected sorted mode index by feature
lambda plane scatter
```

跨 dataset 汇总：

```text
dataset x feature median selected-mode tau
dataset x target median selected-mode tau
topN sensitivity curve:
    efun tau vs deconv tau
    hit-weighted and unique-mode both shown
```

### D. ROI sanity 图

必画但只作为 sanity：

```text
signed ROI profile corr:
    efun_theta vs deconv_RG
    top3/top5/top10/top20

4-target pairwise ROI matrix:
    efun_theta
    efun_RG
    deconv_theta
    deconv_RG

P8 selected BOLD mode overlap matrix
```

可选：

```text
all ROI by dataset:
    keep anatomical ROI order
    per-dataset color scale

top positive / negative ROI labels by target
```

## 5. 历史叙事的当前状态

### 5.1 Event density alone

状态：降级为 baseline / negative control。

原因：

```text
P8/P10 top hits 中 event_density 很少成为主解释；
dimred density 和 raw density 更强。
```

### 5.2 Pure ripple component

状态：被 `ripple_gamma_no_theta` fast subprocess 取代。

原因：

```text
ripple 和 gamma 经常共同出现。
要求 pure ripple 太严格，且不符合当前 P5 label 结构。
```

### 5.3 theta-slow -> BOLD efun 的简单一对一叙事

状态：保留但弱化。

原因：

```text
BOLD efun 确实更慢；
efun_theta 是最慢 target；
但 efun top hits 不总是 clean theta-selective。

所以可以写：
    efun branch is slow/state-like

不要写：
    efun branch is purely theta-specific
```

### 5.4 theta/RG 有不同 ROI network

状态：不是主证据。

原因：

```text
同一个 BOLD feature family 内，theta/RG ROI profile 很像；
ROI 差异主要来自 efun vs deconv_efun feature family。
```

当前可以写：

```text
ROI signed profile is compatible with efun/deconv branch difference,
but does not by itself prove theta/RG spatial separation.
```

### 5.5 Zoom-window trace correlation

状态：只作为 mechanism illustration。

原因：

```text
局部窗口可能挑到过于好看的片段；
主证据必须来自 whole trace / topN distribution / cross-dataset repetition。
```

## 6. 新 dataset 的最小执行顺序

### Step 1: 主线结果检查

```text
P2 event detection
P5 standardized complex-split k03:k16 reduction/density
P7 roi_mean BOLD_POST
P8 roi_mean xcorr
P10 roi_mean xcorr, optional
```

### Step 2: 生成 P5 band selectivity

```text
python scripts/run_pipeline5_band_selectivity.py ...
```

默认不要批量生成 heavy paired top-window sheets。
只有需要解释某个 method-k 时再按需画 diagnostic plate。

### Step 3: 生成 minimal coupling summary

```text
/home/kdshao/anaconda3/bin/python scripts/run_minimal_new_dataset_consistency.py ...
```

### Step 4: 生成 signed ROI / pairwise control

```text
/home/kdshao/anaconda3/bin/python scripts/analyze_roi_mean_signed_target_modes.py \
    --top-ns 3 5 10 20

/home/kdshao/anaconda3/bin/python scripts/analyze_roi_mean_signed_target_pairwise.py \
    --top-ns 3 5 10 20
```

### Step 5: 生成 P8 selected BOLD mode timescale

Broad 版：

```text
/home/kdshao/anaconda3/bin/python scripts/analyze_roi_mean_p8_selected_bold_mode_timescales.py \
    --top-ns 3 5 10 20
```

P5-pair-pass strict 版：

```text
/home/kdshao/anaconda3/bin/python scripts/analyze_roi_mean_p8_selected_bold_mode_timescales.py \
    --top-ns 3 5 10 20 \
    --require-p5-pair-pass \
    --output-dir results/standardized_csplit_k03_k16_all_current_20260607/roi_mean_p8_selected_bold_mode_timescales_p5_pair_pass \
    --figure-dir /mnt/e/DataPons_processed/summary_figures/pipeline11_current_analysis_summary/standardized_csplit_k03_k16_all_current_20260607/roi_mean_p8_selected_bold_mode_timescales_p5_pair_pass
```

## 7. 主结论写作模板

支持时：

```text
Across datasets with a valid standardized complex-split P5 theta/RG pair,
roi_mean BOLD deconv_efun preferentially couples to RG-no-theta dimred BLP
density.  P8 selected-mode eigenvalue analysis further shows that efun-coupled
BOLD modes are slower and closer to the unit circle, whereas deconv-coupled
BOLD modes are faster.  This supports a working model in which BLP theta-like
activity tracks slow/state-like BOLD readout, while RG-no-theta activity couples
to a faster intrinsic-perturbation BOLD branch.
```

中文：

```text
在通过 standardized complex-split P5 theta/RG pair gate 的数据中，
roi_mean BOLD deconv_efun 比 efun 更倾向于耦合 RG-no-theta dimred BLP density。
同时，P8 selected BOLD mode eigenvalue/timescale 显示 efun 选中的 BOLD modes
更慢、更接近单位圆，而 deconv_efun 选中的 modes 更快。
这支持一个工作模型：theta-like BLP activity 更接近 slow/state-like BOLD readout，
而 RG-no-theta activity 更接近 fast intrinsic-perturbation BOLD branch。
```

必须保留的限制：

```text
This is cross-dataset consistency evidence, not causal proof.
ROI profiles are sanity checks and do not currently prove theta/RG spatial separation.
P10 selected components need separate weighted-source-mode timescale analysis.
```
