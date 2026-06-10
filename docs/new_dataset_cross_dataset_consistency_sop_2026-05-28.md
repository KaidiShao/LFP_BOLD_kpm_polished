# 新 dataset P5 band selectivity consistency SOP

日期：2026-05-29

目的：定义一个新 dataset 有 standardized complex-split P4/P5 结果之后，如何先在 P5 层面检查是否复现 E10gb1 的 two-subprocess 叙事：

```text
theta-like slow/state BLP subprocess
    -> later expected to relate more to BOLD efun / slow BOLD state variable

ripple-gamma-like fast/intrinsic-trigger BLP subprocess
    -> later expected to relate more to BOLD deconv_efun / residual perturbation
```

本 SOP 当前只保留 **P2 band event -> P5 strict band selectivity** 这一条主线，以及已经生成的结果总结。旧版 component-label 设计和旧 P8/P10 总图设计暂时全部移除。

## 0. 2026-06-01 当前主计划：少图、分 gate、先 P5 后 P8/P10

当前 SOP 的主目标不是继续生成更多 summary 图，而是把 cross-dataset consistency 压成少数几个可以判定的问题。

核心叙事：

```text
standardized complex-split BLP 里是否能稳定拆出两类 subprocess：

1. theta_selective component
   -> slow / session-state-like BLP variable

2. ripple-gamma selective component
   -> fast intrinsic-trigger / perturbation-like BLP variable

然后再问：
BOLD efun 和 BOLD deconv_efun 是否分别更像这两类 BLP subprocess。
```

正式阅读顺序必须是：

```text
Gate 0: dataset QC
Gate 1: P5 subprocess gate
Gate 2: P8/P10 coupling gate
Gate 3: cross-dataset consistency
```

### Gate 0: dataset QC

目的：先判断 dataset 是否可以进入主分析。

主图只看：

```text
P2 theta/gamma/ripple top-window figures
P5 top-window diagnostic plate
```

输出结论只分三类：

```text
usable
suspicious_but_usable
exclude_from_cross_dataset_conclusion
```

这一步主要用于识别 event detection 或原始信号异常，例如 K13m17 这种疑似 theta event saturation / clipping 的情况。

### Gate 1: P5 subprocess gate

核心问题：

```text
这个 dataset 的 dimred BLP efun space 里，是否能找到：

1. theta_selective component
2. ripple_gamma_no_theta / gamma_selective / ripple_selective component
```

主图只保留：

```text
strict selectivity plane
strict two-subprocess candidate map
每个候选 method-k 的 diagnostic plate
```

判定规则：

```text
同一个 method-k 同时有 theta component 和 RG component：
    P5 subprocess gate pass

多个 method-k pass：
    subprocess structure 稳定

只有一个 method-k 勉强 pass：
    结果可疑，需要 trace/window 解释图支持

P5 不 pass：
    暂停解释该 dataset 的 P8/P10 coupling
```

当前主线暂时不使用 dominant label。P5 主判断只使用 strict band-selectivity label。

### Gate 2: P8/P10 coupling gate

P5 pass 以后才看 P8/P10。

核心问题：

```text
BOLD efun 和 BOLD deconv_efun 是否偏向不同 BLP subprocess？
```

必须按 BOLD observable 分开看：

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

每个 observable 内再分：

```text
efun
deconv_efun
```

每一栏看：

```text
raw BLP density top hits:
    是 slow timescale 还是 fast timescale？

dimred BLP density top hits:
    strict label 是 theta group 还是 RG group？

whole trace:
    是否支持 slow-state / fast-trigger 的解释？
```

topN sensitivity 要保留，但不能把 `topN mean` 作为主证据。

默认 topN：

```text
top1, top3, top5, top10, top20
```

主证据优先级：

```text
1. top hits 的 strict-label composition
2. top hits 的 raw-efun timescale distribution
3. top hits 是否在不同 dataset / observable 中重复出现
4. representative whole-trace examples
```

zoom trace 只作为形状解释图，不能单独作为全局强度证据。

### Gate 3: cross-dataset consistency

cross-dataset 主结论只回答五个问题：

```text
1. 哪些 dataset 通过 P5 subprocess gate？
2. 哪些 method-k 在多个 dataset 都能找到 theta + RG pair？
3. BOLD efun top hits 是否更偏 theta / slow？
4. BOLD deconv_efun top hits 是否更偏 RG / fast？
5. 这个结论在哪些 BOLD observable 里成立，尤其 HP_svd100 是否不同于 global observables？
```

最终主图只保留少数几张：

```text
dataset x method-k 的 P5 pass map
observable x efun/deconv 的 subprocess preference map
efun vs deconv raw efun timescale distribution by topN
representative whole-trace examples, split by dataset and BOLD observable
```

其它图全部降级为 QC / appendix / debugging，不作为主结论入口。

### 2026-06-01 当前读表结论

当前这版结论来自已有的 7-dataset standardized complex-split 结果，不额外补算。解读时先执行 Gate 1，排除 P5 gate 不通过的数据，再看 P8/P10。

Gate 1 结论：

```text
e10gb1 : P5 pass, adaptive_envelope 24/24 method-k 有 theta + RG pair
e10fV1 : P5 pass, adaptive_envelope 18/24
e10gh1 : P5 pass, adaptive_envelope 24/24
e10gw1 : P5 pass, adaptive_envelope 18/24
f12m01 : P5 pass, adaptive_envelope 16/24
k13m23 : P5 pass but weaker, adaptive_envelope 11/24
k13m17 : P5 fail, adaptive_envelope 0/24; 当前主结论中排除
```

排除 `k13m17` 后，跨 dataset 最稳定的 P5 method-k 是：

```text
6/6 pass:
    mds_k04, mds_k05, mds_k06, mds_k07, mds_k08
    nmf_k04
    umap_k04, umap_k05, umap_k06, umap_k08
```

因此 P5 层面的主判断是：

```text
除 k13m17 外，多数 dataset 都能在 standardized complex-split dimred BLP efun space
里拆出 theta-like 和 ripple-gamma-like subprocess。
P5 two-subprocess 叙事成立。
```

Gate 2 结论：

```text
P8/P10 top hits 的主 density source 是 Koopman efun density，
不是 P2 event density。

event_density 通常只占 top hits 的 0-8%。
主要竞争是 raw BLP efun density 和 dimred BLP efun density。
```

当前不能说 dimred density 全面压倒 raw density。更准确的结论是：

```text
raw BLP efun density 和 dimred BLP efun density 是互补证据。
raw efun density 与 BOLD efun 的关系尤其重要，因为它更容易暴露
session-wise slow state change。
```

dimred strict-label coupling 的当前结论：

```text
BOLD deconv_efun:
    稳定偏 ripple-gamma / fast-trigger group。
    排除 k13m17 后，top10 dimred hits 里 RG group 通常约 0.63-0.94。

BOLD efun:
    不应简单解释为 theta_selective。
    top10 dimred hits 的 theta fraction 通常较低，更多是 mixed / RG / inactive。
    因此 BOLD efun 的 slow-state 证据不能只靠 P5 strict theta label。
```

P5-stable-method-k-only sensitivity check：

```text
目的：
    检查 deconv_efun -> RG 的结论是不是只来自全 method-k 搜索空间太大。

stable method-k set:
    mds_k04, mds_k05, mds_k06, mds_k07, mds_k08
    nmf_k04
    umap_k04, umap_k05, umap_k06, umap_k08

excluded dataset:
    k13m17
```

当前实现：

```text
scripts\summarize_p8_p10_stable_method_k_coupling.py

outputs:
results\pipeline8_10_strict_band_coupling_p5_stable_method_k_20260601\
```

top10, competed-topN 读法下，也就是 raw/event/stable-dimred 一起竞争后再看 dimred label：

```text
P8 efun:
    all method-k RG fraction       = 0.34
    P5-stable-only RG fraction     = 0.45

P8 deconv_efun:
    all method-k RG fraction       = 0.68
    P5-stable-only RG fraction     = 0.71

P10 efun:
    all method-k RG fraction       = 0.50
    P5-stable-only RG fraction     = 0.51

P10 deconv_efun:
    all method-k RG fraction       = 0.80
    P5-stable-only RG fraction     = 0.80
```

dimred-only top10 读法下，也就是只在 dimred density 内部排序：

```text
P8 efun:
    all method-k RG fraction       = 0.55
    P5-stable-only RG fraction     = 0.53

P8 deconv_efun:
    all method-k RG fraction       = 0.78
    P5-stable-only RG fraction     = 0.81

P10 efun:
    all method-k RG fraction       = 0.62
    P5-stable-only RG fraction     = 0.58

P10 deconv_efun:
    all method-k RG fraction       = 0.84
    P5-stable-only RG fraction     = 0.84
```

因此 stable-only 敏感性检查支持：

```text
deconv_efun -> RG / fast-trigger 的结论不是全搜索空间偶然造成的；
即使只保留跨 dataset P5 gate 最稳定的 method-k，这个结论仍然成立。

efun 仍然不是 clean theta-selective。
efun 的 slow-state 解释仍应主要依赖 raw efun timescale + whole-trace evidence。
```

raw efun timescale 的当前结论：

```text
BOLD deconv_efun top raw-density hits:
    几乎稳定来自 fast raw modes。
    P8 top10 median tau 多在约 0.001-0.003 sec，slow fraction 接近 0。

BOLD efun top raw-density hits:
    更容易包含少数 very slow / session-wise state modes。
    这些 slow modes 会显著拉高 mean tau，并且在 whole-trace 图里表现为
    session-wise baseline/state shift。
```

2026-06-01 已补正式版本：

```text
scripts\analyze_bold_efun_top_raw_density_hits.py

新增口径：
    --exclude-datasets k13m17
    real BOLD feature only
    xcorr-sign aligned BOLD trace
    whole-trace robust z normalization
    masked samples shaded gray
    session borders shown as dashed lines

P8 output:
results\bold_efun_raw_density_slow_state_probe_p5gate_pass_only_p8_20260601\

P10 output:
results\bold_efun_raw_density_slow_state_probe_p5gate_pass_only_p10_20260601\
```

排除 `k13m17` 后，raw top10 timescale 结果：

```text
P8 raw_top10:
    efun        weighted mean tau = 10.842 sec, slow >=10s fraction = 0.09
    deconv_efun weighted mean tau = 0.008 sec,  slow >=10s fraction = 0.00

P10 raw_top10:
    efun        weighted mean tau = 13.753 sec, slow >=10s fraction = 0.07
    deconv_efun weighted mean tau = 0.025 sec,  slow >=10s fraction = 0.00
```

按 observable 看，BOLD efun 的 slow tail 主要出现在：

```text
P8:
    global_svd100 slow>=10s fraction = 0.16
    HP_svd100     slow>=10s fraction = 0.10
    roi_mean      slow>=10s fraction = 0.11
    gsvd100_ds    slow>=10s fraction = 0.00

P10:
    global_svd100 slow>=10s fraction = 0.14
    HP_svd100     slow>=10s fraction = 0.09
    roi_mean      slow>=10s fraction = 0.03
    gsvd100_ds    slow>=10s fraction = 0.01
```

这说明：

```text
BOLD efun 与 raw BLP efun density 的关系确实更容易带出 slow/session-wise state modes，
但 slow signal 是 tail effect，不是 median effect。

BOLD deconv_efun 几乎没有 slow tail，稳定对应 fast raw modes。
```

同时已补 selected-method dimred trace：

```text
scripts\plot_bold_density_trace_by_dataset_observable.py

新增口径：
    --method-k-scope p5_stable_method_k

stable method-k set:
    mds_k04, mds_k05, mds_k06, mds_k07, mds_k08
    nmf_k04
    umap_k04, umap_k05, umap_k06, umap_k08

output:
results\p8_p10_dimred_stable_method_k_trace_by_observable_20260601\
```

这些图对每个 P5-pass dataset、P8/P10、每个 BOLD observable、efun/deconv，画：

```text
dimred_efun_density full trace
dimred_efun_density fixed zoom

selection:
    top3 unique BLP dimred density per observable x efun/deconv
    only P5-stable method-k
```

`k13m23` 当前只生成 HP_svd100 和 roi_mean 两个 observable 的 trace，因为 global_svd100 / gsvd100_ds 的 P8/P10 结果缺失。

因此当前推荐叙事应写成：

```text
BLP 侧能稳定拆出 theta slow-state subprocess 和 ripple-gamma fast-trigger subprocess。
BOLD deconv_efun 稳定耦合 ripple-gamma / fast-trigger。
BOLD efun 更像 session-wise slow-state coupling，
但这个 slow-state 不一定等价于 P5 strict theta label；
需要用 raw efun timescale + whole-trace examples 单独证明。
```

这也解释了为什么后续主图必须同时保留：

```text
1. dimred strict-label composition
2. raw efun timescale distribution
3. representative whole-trace examples
```

其中 raw efun density 与 BOLD efun 的关系需要作为主分析线，而不是附属 QC。

## 1. 固定主线口径

第一阶段只看一个 BLP branch：

```text
complex_split_projected_vlambda_standardize
```

默认 P5 dimred grid：

```text
method = SVD, NMF, MDS, UMAP
k      = 03:08
```

默认 activity transforms：

```text
abs
adaptive_envelope
```

当前 P5/P8/P10 density 口径仍然是：

```text
complex_split_projected_vlambda_standardize_rmsenv_adaptive
threshold ratio q070
activity transform = RMS envelope of abs eigenfunction/component
activity window = eigenvalue-adaptive
```

## 2. Availability / provenance audit

任何新 dataset 进入解释前，先检查：

1. P4 standardized csplit full EDMD outputs 是否存在，不只看 summary MAT。
2. P5 standardized csplit reduction 是否有 24/24 method-k。
3. P2 band event detection 文件是否存在。
4. P5 band-selectivity label CSV 是否由当前 P2/P5 结果生成。
5. P5 dimred density 是否覆盖 24/24 method-k。
6. 后续若看 P8/P10，必须确认其 density source 来自当前 P5，而不是 legacy/cache。

缺失或旧版本不能用历史图补上，必须显式标记：

```text
missing_p2_band_events
missing_p4_full_outputs
missing_p5_reduction
missing_p5_density
missing_p5_band_selectivity
stale_or_legacy_only
```

## 3. P5 band selectivity 定义

关键修正：P5 subprocess label 的来源必须是 P2 band event detection，而不是 overlapping event-family label。P2 直接给出三个基础 band：

```text
theta  = 2-15 Hz
gamma  = 30-90 Hz
ripple = 90-190 Hz
```

输入：

```text
E:\DataPons_processed\<dataset>\pipeline2_event_detection\<dataset>_bandpass_events_3bands.mat
E:\DataPons_processed\<dataset>\pipeline5_eigenfunction_reduction\complex_split_projected_vlambda_standardize\<method>_k<kk>\mat\*.mat
```

对每个 dataset、method、k、component：

1. 从 P2 `DetectResults` 生成 theta/gamma/ripple 三个 full-length event mask。
2. 每个 band 内合并所有 channel 的 event windows。
3. baseline 定义为不在 theta/gamma/ripple 任一 event window 内的 samples。
4. 读取 P5 `result/core/temporal_components_time_by_comp`。
5. 同时计算两种 activity：

```text
abs:
    A_c(t) = abs(component_c(t))

adaptive_envelope:
    A_c(t) = sqrt(movmean(abs(component_c(t))^2, adaptive_window_c))
```

adaptive envelope 必须 session-aware，不能跨 session boundary 平滑。

dimred component 的 adaptive window 用 source-mode timescale 决定：

```text
source tau_j     = preferred discrete-log eigenvalue timescale
component tau_c  = weighted_median(source tau_j, abs(component loading_jc))
window_c         = clamp(0.5 * component tau_c, 0.03 sec, 1.0 sec)
```

每个 band 的响应分数：

```text
effect_z_band =
    (mean(A_c inside band event windows) - mean(A_c inside non-event baseline))
    / sd(A_c inside non-event baseline)
```

当前正式分析只使用 strict label：

```text
theta_selective
gamma_selective
ripple_selective
ripple_gamma_no_theta
mixed_or_partial
inactive
```

当前 strict threshold：

```text
active effect threshold = 0.50 baseline SD
off-target fraction     = 0.60 of target response
target margin           = 0.50 baseline SD
```

P5 band selectivity 的主判断只看 `strict_label`。

## 4. P5 band selectivity 正式主图

每个 dataset、condition、activity transform 都应该按同一套图生成。默认 transform 为：

```text
abs
adaptive_envelope
```

主图包括：

```text
00_strict_label_counts_abs_vs_adaptive_envelope.png
04_<transform>_strict_theta_vs_ripple_gamma_effect_scatter.png
02_<transform>_strict_label_grid_by_method_k.png
03_<transform>_strict_label_composition_by_method_k.png
05_<transform>_strict_two_subprocess_candidate_map.png
```

这些图的用途：

1. `strict_label_counts`：比较 `abs` 和 `adaptive_envelope` 整体上哪个更容易产生 clean theta / ripple-gamma component。
2. `strict_theta_vs_ripple_gamma_effect_scatter`：主二维分离图。横轴是 theta event effect z，纵轴是 `max(gamma, ripple)` event effect z。用来看 theta 和 ripple-gamma component 是否分开。
3. `strict_label_grid_by_method_k`：每个格子是一个 component，标签是 strict label。用来看每个 method-k 里有没有 clean component。
4. `strict_label_composition_by_method_k`：stacked bar。用来看某个 method-k 是不是大部分都是 mixed，还是能稳定给出 clean theta 和 clean RG。
5. `strict_two_subprocess_candidate_map`：method x k 矩阵，标出是否同时有 theta candidate 和 ripple-gamma candidate，并写出 component index。这是 P5 参数选择的主图之一。

`strict_two_subprocess_candidate_map` 的通过标准是同一个 method-k 里同时有：

```text
at least one theta_selective component
at least one ripple_gamma_no_theta / gamma_selective / ripple_selective component
```

## 5. P5 band selectivity 解释图

除了 summary 主图，每个 method-k 还应该能生成用于解释的 component-pair 图：

```text
06_<transform>_selected_component_perievent_mean_by_band.png
<dataset>_<transform>_<method>_k<kk>_thetaC<cc>_rgC<cc>_<ribbon>_diagnostic_plate.png
<dataset>_<transform>_<method>_k<kk>_thetaC<cc>_rgC<cc>_paired_top_windows.png
```

`selected_component_perievent_mean_by_band`：

```text
对选出来的代表 theta component 和 RG component，
分别画 theta/gamma/ripple P2 event peak 对齐的 mean response。
```

正式版支持 variability ribbon：

```text
SEM = 默认主图，适合看均值结构是否稳定
STD = QC 图，适合看 event-to-event variability
IQR = robust QC 图，适合处理 outlier-heavy dataset
```

`diagnostic_plate` 是当前正式推荐的单 method-k 解释图。它把两类信息放在同一张大图里：

```text
左边：
    theta-like component 的 peri-event mean +/- ribbon
    ripple-gamma-like component 的 peri-event mean +/- ribbon
    三条线分别对齐 P2 theta/gamma/ripple event peak

右边：
    前 15 个 theta-component top windows
    后 15 个 ripple-gamma-component top windows
    每个小图同时画 theta-like 和 RG-like trace
    叠加 P2 theta/gamma/ripple event rugs
```

这张图同时回答两个问题：

```text
平均意义上：这个 component 对哪个 band event 响应？
真实窗口里：theta 和 RG 两条 component activity 是否真的分开？
```

`paired_top_windows` 是不带左侧 peri-event summary 的 component-pair window sheet。它应该对每个有合法 theta/RG pair 的 method-k 都生成，并且 `abs` 与 `adaptive_envelope` 分开生成。

## 6. P5 band selectivity QC 图

完整的 component x band response heatmap 只作为 QC/appendix：

```text
01_<transform>_component_p2_band_event_effect_heatmap.png
```

它可以帮助检查 broad failure 或异常 dataset，但不作为主参数选择图，因为所有 method-k/component 行挤在一起后不够直观。

对可疑数据，例如 K13m17，还需要生成 P2 top30 event-window figures：

```text
P2 theta top30 event windows
P2 gamma top30 event windows
P2 ripple top30 event windows
```

这一步用于判断 P2 event detection 或原始 BLP signal 是否异常，不直接替代 P5 component selectivity 结论。

## 7. 一键 pipeline

当前对应实现脚本：

```text
scripts\run_pipeline5_band_selectivity.py
scripts\analyze_e10gb1_p5_dimred_band_event_response_v2.py
scripts\plot_p5_event_diagnostic_plate.py
scripts\plot_p5_paired_subprocess_traces.py
scripts\plot_p2_p5_top30_windows.py
```

正式一键入口：

```powershell
wsl.exe -d Ubuntu-22.04 --cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished `
  /home/kdshao/anaconda3/bin/python scripts/run_pipeline5_band_selectivity.py `
  --dataset <dataset>
```

先检查数据是否齐：

```powershell
wsl.exe -d Ubuntu-22.04 --cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished `
  /home/kdshao/anaconda3/bin/python scripts/run_pipeline5_band_selectivity.py `
  --dataset <dataset> --dry-run
```

默认会执行：

```text
1. P2-band event response label table
2. P5 band-selectivity summary figures
3. SEM diagnostic plates for all method-k pairs that pass theta/RG candidate selection
4. paired top-window sheets for all valid method-k pairs
5. P2/P5 top30 window QC if the required P2 top-window table exists
```

常用选项：

```text
--ribbons sem std iqr
--summary-only
--skip-diagnostic-plates
--skip-paired-top-windows
--skip-p2-p5-top30
--dry-run
```

## 8. P5 subprocess 3D trajectory

在 P5 band selectivity 成立以后，可以进一步画 3D trajectory，但不再使用旧的 `first three dimred components`。正式 trajectory 应该围绕当前 two-subprocess 叙事定义：

```text
x = selected theta-selective component activity
y = selected ripple-gamma-no-theta/gamma/ripple component activity
z = residual PC1 after regressing x/y from the remaining components
```

第三维不能随便再挑一个看起来像 fast 的 component。原因是：

```text
1. arbitrary component choice 很难跨 dataset / method-k 比较。
2. NMF/MDS/UMAP/SVD component 不一定能用同一种 component index 解释。
3. 如果第三维是人为挑选的 secondary fast component，容易把结果解释成选择规则本身。
```

因此默认第三维必须是一个规则化 residual axis：

```text
input:
    A(t, c) = same method-k 的 dimred component activity
    默认 activity transform = adaptive_envelope

selected axes:
    x = best theta_selective component
    y = best ripple_gamma_no_theta/gamma_selective/ripple_selective component

residual construction:
    A_other = all components excluding x/y
    A_other_resid = A_other - projection(A_other onto [x, y])
    z = PC1(A_other_resid)
```

如果 `k03` 只剩一个 component，`z` 就是该剩余 component 在回归掉 `x/y` 后的 residual trace。这个仍然是合法定义。

标题和 manifest 必须记录：

```text
dataset
condition
activity_transform
method
k
theta_component
ripple_gamma_component
z_axis = residual_pc1_after_theta_rg_regression
z_explained_variance_ratio
correlation(z, x)
correlation(z, y)
```

解释原则：

```text
如果 z_explained_variance_ratio 很低：
    theta + RG 两维已经解释了主要可视化结构，3D 图主要是二维结构的辅助投影。

如果 z_explained_variance_ratio 较高：
    第三维代表 theta/RG 之外的 leftover variability。
    后续再检查它是否主要对应 high-frequency oscillation content in the fast process。

如果 z 与 x/y 仍高度相关：
    residual axis 定义或 component scaling 需要重新检查，不应过度解释该 3D trajectory。
```

这张图应该替代旧版 P5 `first3 dimred trajectory` 作为主解释 trajectory。旧版 first-three trajectory 只能保留为 diagnostic/QC。

## 9. 后续 P8/P10 要解决的问题和画图设计

P8/P10 xcorr 不使用上面 3D trajectory 的三维作为唯一 source。3D trajectory 是 P5 dynamics 的解释图，不替代完整 dimred efun density space。

后续 P8/P10 应该：

```text
xcorr source = all P5 dimred efun densities
annotation   = P5 strict band-selectivity label for each dimred component
```

因此 P8/P10 top hits 应该 join：

```text
strict_label
theta_effect_z
gamma_effect_z
ripple_effect_z
activity_transform
method
k
component
source_file provenance
```

当前 P8/P10 解释只使用 P5 strict band-selectivity annotation。

### 9.1 P8/P10 核心问题

P8/P10 后续分析要回答五个问题：

1. `dimred BLP efun density` 是否比 `raw BLP efun density` 和 `P2 event density` 更能解释 BOLD efun/deconv_efun。
2. 和 `BOLD efun` 高 xcorr 的 dimred BLP components 是否更偏 `theta_selective`。
3. 和 `BOLD deconv_efun` 高 xcorr 的 dimred BLP components 是否更偏 `ripple_gamma_no_theta / gamma_selective / ripple_selective`。
4. 上述 pattern 在 P8 和 P10 是否一致。
5. 高 xcorr BOLD modes 的 ROI summary / activation map 是否跨 dataset 一致。

这里的重点不是只看 top1/top5 的平均值，而是看不同 topN 下 label composition、enrichment、lag 和 spatial consistency 是否稳定。

关键分层规则：

```text
P8/P10 所有正式统计必须按 BOLD observable 分开生成。

第一层分开：
    global_svd100
    gsvd100_ds
    HP_svd100
    roi_mean

第二层再分：
    efun
    deconv_efun
```

原因：

```text
HP_svd100 更偏 local BOLD observable；
global_svd100 / gsvd100_ds / roi_mean 更偏 global or ROI-summary observable。
如果四种 observable 混在一起统计，theta slow-state 和 RG fast-trigger
pattern 很容易被平均掉，也无法判断某个结论到底来自哪个 BOLD observable。
```

`all observables` 汇总可以保留为 appendix/QC，但不能作为主判断。

默认 topN sensitivity：

```text
top1, top3, top5, top10, top20, top50
```

### 9.2 P8/P10 top-hit annotation table

P8/P10 需要先生成一个统一 top-hit annotation table。每一行是一个 xcorr hit：

```text
pipeline                 # P8 or P10
dataset
BLP condition
BOLD observable
BOLD feature family      # efun or deconv_efun
density_source_family    # raw_efun, dimred_efun, p2_event_density
activity_transform       # abs or adaptive_envelope when applicable
method                   # svd/nmf/mds/umap for dimred_efun
k
component
strict_label
theta_effect_z
gamma_effect_z
ripple_effect_z
BOLD mode index
peak_abs_corr
signed_corr
lag_sec
source_file provenance
```

这个表是所有 P8/P10 图的源头。不要直接从 summary figure cache 判断完整性。

当前实现脚本：

```text
scripts\analyze_p8_p10_strict_band_coupling.py
```

默认输入：

```text
P5 labels:
results\<dataset>_p5_p2_band_event_response_v2_selective_envelope_20260528\component_band_event_response_v2.csv

P8 xcorr:
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\**\xcorr_csplit_standardize_rmsenv_adaptive_peaks__*.csv

P10 xcorr:
E:\DataPons_processed\<dataset>\pipeline10_dimred_xcorr\<run_tag>__<p9_feature>__<p9_method_k>\**\dimred_xcorr_csplit_standardize_rmsenv_adaptive_peaks__*.csv
```

当前补算脚本：

```text
scripts\script_run_standardized_csplit_p8_for_p12.m
scripts\script_run_standardized_csplit_p10_for_p12.m
```

### 9.2.1 P8/P10 主图 0：BOLD observable 分层 top-density trace plate

这组图现在作为 P8/P10 的第一类主图，而不是临时 preview。它的目的不是做最终统计，
而是先直接检查：

```text
每个 BOLD observable 里，
BOLD efun / BOLD deconv_efun 分别和哪些 BLP density/component 最像？

这些相似性在 whole trace 上是不是看得见？
efun 是否更像 slow/session-wise state variable？
deconv_efun 是否更像 fast intrinsic-trigger / perturbation-like activity？
```

图的基本单位：

```text
one figure = one dataset x one pipeline x one density class

rows    = BOLD observable x BOLD feature family
          BOLD feature family = efun, deconv_efun

columns = top unique BLP density/component hits
          default = top3
```

这里必须使用 `unique_density` selection，而不是简单的 top xcorr pair。
原因是同一个 BLP density/component 可能同时匹配多个 BOLD modes；如果直接取
top xcorr pair，top1/top2 可能只是同一条 BLP trace 的重复。正式主图规则是：

```text
对每个 BOLD observable x efun/deconv x density class：
    按 peak_abs_corr 排序
    每个 BLP density/component 只保留一次
    对该 BLP density/component 保留它最匹配的 BOLD mode
    取 top3 unique BLP density/component
```

默认 density class 分开画：

```text
raw_efun_density
dimred_efun_density
event_density      # optional/reference，不作为当前 trace 主图第一优先级
```

正式 trace 口径：

```text
raw/dimred BLP density trace:
    unsmoothed

BOLD trace:
    efun_real or deconv_efun_real
    unsmoothed
    按 xcorr signed_corr 符号翻转，使相似性方向一致

normalization:
    whole-trace robust z
    zoom 图也沿用 whole-trace normalization，不在局部窗口重新缩放

x axis:
    full trace = concatenated session time, minutes
    zoom trace = time within fixed selected window, minutes

QC overlays:
    grey shaded spans = xcorr masked samples
    dashed vertical lines = session borders
```

图分两类：

```text
full_trace_by_observable:
    看 whole-trace 结构，尤其是 session-wise slow state change。

fixed_zoom_by_observable:
    看局部细节，但窗口不能按局部相关挑选。
    默认固定窗口包括：
        raw/dimred density high-activity window
        BOLD feature high-activity window
        middle/session-representative window
```

正式实现脚本：

```text
scripts\plot_bold_density_trace_by_dataset_observable.py
```

推荐命令模板：

```powershell
wsl.exe -d Ubuntu-22.04 --cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished `
  /home/kdshao/anaconda3/bin/python scripts/plot_bold_density_trace_by_dataset_observable.py `
    --dataset <dataset> `
    --pipeline P8 `
    --bold-feature-filter real `
    --density-classes raw_efun_density dimred_efun_density `
    --results-dir results/<dataset>_bold_density_trace_by_observable_<date> `
    --window-sec 120 `
    --top-hits-per-group 3 `
    --selection-mode unique_density `
    --align-xcorr-sign
```

P10 使用同一脚本，只改：

```text
--pipeline P10
```

推荐输出命名：

```text
p8_<dataset>_raw_efun_density_unique_density_top3_full_trace_by_observable.png
p8_<dataset>_raw_efun_density_unique_density_top3_fixed_zoom_by_observable.png
p8_<dataset>_dimred_efun_density_unique_density_top3_full_trace_by_observable.png
p8_<dataset>_dimred_efun_density_unique_density_top3_fixed_zoom_by_observable.png

p8_<dataset>_raw_efun_density_unique_density_top3_selected_rows.csv
p8_<dataset>_dimred_efun_density_unique_density_top3_selected_rows.csv
p8_<dataset>_unique_density_top3_trace_records.csv
```

解读优先级：

```text
1. 先看 dimred_efun_density full trace：
   BOLD efun top hits 是否更常落在 theta_selective 或 slow-looking components。
   BOLD deconv_efun top hits 是否更常落在 RG group 或 fast-looking components。

2. 再看 raw_efun_density full trace：
   BOLD efun 是否对应 slow raw modes / session-wise state changes。
   BOLD deconv_efun 是否对应 fast raw modes / transient activity。

3. zoom 图只用于解释局部形状，不单独作为全局证据。
   全局证据仍来自 full trace、topN statistics、cross-dataset replication。
```

E10gb1 当前验证过的 preview：

```text
results\e10gb1_bold_density_trace_by_observable_unique_density_top3_preview_20260529\
```

### 9.3 图 1：density-source family competition

目的：

```text
dimred BLP efun density 是否比 raw efun density / P2 event density 更强？
```

图形：

```text
x-axis  = topN
y-axis  = fraction of topN hits or best-hit presence
color   = density_source_family
facet   = pipeline x BOLD observable x BOLD feature family
```

分别画：

```text
P8 efun
P8 deconv_efun
P10 efun
P10 deconv_efun
```

解读：

```text
如果 dimred_efun 在多个 topN、多个 dataset 里稳定占优，
说明降维后的 BLP efun density 提供了比 raw/event density 更有用的 BOLD coupling source。
```

### 9.4 图 2：topN strict-label composition

目的：

```text
BOLD efun 和 BOLD deconv_efun 的 top dimred hits 是否偏向不同 P5 strict labels？
```

图形：

```text
x-axis  = topN
y-axis  = label fraction
color   = strict_label
facet   = pipeline x BOLD observable x BOLD feature family
```

其中 strict labels 至少包括：

```text
theta_selective
ripple_gamma_no_theta
gamma_selective
ripple_selective
mixed_or_partial
inactive
```

理想 pattern：

```text
BOLD efun top hits:
    theta_selective fraction higher

BOLD deconv_efun top hits:
    ripple_gamma_no_theta / gamma_selective / ripple_selective fraction higher
```

### 9.5 图 3：label enrichment vs available background

不能只看 top hits 的 label fraction，因为某些 label 本身在 P5 里更多。必须画 enrichment：

```text
enrichment(label, topN) =
    fraction(label among topN dimred hits)
    /
    fraction(label among all available dimred components for that dataset/condition/transform)
```

图形：

```text
x-axis  = topN
y-axis  = enrichment ratio
color   = strict_label
facet   = pipeline x BOLD observable x BOLD feature family
```

解读：

```text
enrichment > 1:
    该 label 在 high-xcorr hits 中过度出现。

enrichment ~ 1:
    只是因为该 label 本来就多。

enrichment < 1:
    high-xcorr hits 避开该 label。
```

### 9.6 图 4：theta-vs-ripple-gamma effect scatter for top hits

目的：

```text
不用只看离散 label，而是看 top hits 在连续 P2 band-effect plane 上的位置。
```

图形：

```text
x-axis = theta_effect_z
y-axis = max(gamma_effect_z, ripple_effect_z)
color  = BOLD feature family: efun / deconv_efun
shape  = pipeline: P8 / P10
size   = peak_abs_corr
facet  = dataset or BOLD observable
```

解读：

```text
efun hits 如果偏右下/右侧：
    更像 theta slow/state variable coupling。

deconv_efun hits 如果偏左上/上侧：
    更像 ripple-gamma fast/intrinsic-trigger coupling。
```

### 9.7 图 5：method-k hit map

目的：

```text
看哪些 dimred method 和 component number 产生稳定、可解释的 high-xcorr hits。
```

图形：

```text
rows    = method
columns = k03:k08
color   = fraction of topN dimred hits with target strict-label group
text    = n_hits / n_available or best component index
facet   = pipeline x BOLD observable x BOLD feature family x activity_transform
```

建议 target strict-label groups：

```text
theta group = theta_selective
RG group    = ripple_gamma_no_theta + gamma_selective + ripple_selective
mixed group = mixed_or_partial
inactive    = inactive
```

### 9.8 图 6：lag distribution by strict label

目的：

```text
检查 high xcorr 的 lag 是否和 slow-state / fast-trigger 叙事一致。
```

图形：

```text
x-axis  = lag_sec
y-axis  = density or jittered hits
color   = strict_label group
facet   = pipeline x BOLD observable x BOLD feature family
```

解读：

```text
theta_selective + efun:
    可能更慢、更宽，lag 分布可较宽。

RG group + deconv_efun:
    如果是 intrinsic-trigger-like coupling，lag 方向和范围应更集中。
```

### 9.9 图 7：BOLD ROI / activation map consistency by P5 label

目的：

```text
高 xcorr BOLD modes 是否对应跨 dataset 一致的 spatial pattern？
```

做法：

1. 对每个 BOLD mode，找它最强 coupling 的 P5 dimred component label。
2. 按 winner label 分组：

```text
theta-linked BOLD modes
RG-linked BOLD modes
mixed-linked BOLD modes
inactive-linked BOLD modes
```

3. 分别比较 ROI summary / activation map 的 cross-dataset consistency。

图形：

```text
ROI profile correlation heatmap by label group
activation-map similarity heatmap by label group
dataset x ROI summary bars for representative winners
```

解读：

```text
如果 theta-linked BOLD efun maps 跨 dataset 一致，
且 RG-linked deconv_efun maps 也跨 dataset 一致，
则 P5 band-selective subprocess 与 BOLD spatial dynamics 的联系更可信。
```

### 9.10 P8/P10 最终判断

一个 dataset 或参数组合更值得保留，如果它同时满足：

```text
1. P5 有 strict theta + RG candidates。
2. P8/P10 top hits 中 dimred_efun density 优于 raw/event density。
3. BOLD efun hits enriched for theta_selective。
4. BOLD deconv_efun hits enriched for RG group。
5. lag distribution 不乱。
6. 对应 BOLD ROI / activation map 有跨 dataset 一致性。
```

## 10. 2026-05-28 result update: 7-dataset P5-first summary

Analysis target:

```text
condition = complex_split_projected_vlambda_standardize
method    = svd, nmf, mds, umap
k         = 03:08
score     = P2 band-event effect_z vs non-event baseline
```

`both` means one method-k contains both:

```text
1. at least one strict theta_selective component
2. at least one strict ripple_gamma_no_theta / gamma_selective / ripple_selective component
```

Current summary:

```text
dataset   abs both   adaptive_envelope both   interpretation
e10gb1    23/24      24/24                    anchor; strong two-subprocess split
e10fV1    18/24      18/24                    partial but substantial replication
e10gh1    24/24      24/24                    strongest replication of E10gb1
e10gw1    18/24      18/24                    partial but substantial replication
f12m01    14/24      16/24                    weaker but present
k13m17     0/24       0/24                    outlier under current strict rules
k13m23    10/24      11/24                    weak/partial replication
```

Strict label counts by dataset:

```text
e10gb1
  abs:               theta=27, RG-no-theta=34, mixed=38, inactive=33
  adaptive_envelope: theta=28, RG-no-theta=45, mixed=33, inactive=26

e10fV1
  abs:               theta=19, RG-no-theta=25, ripple=17, mixed=49, inactive=22
  adaptive_envelope: theta=19, RG-no-theta=33, ripple=16, mixed=46, inactive=18

e10gh1
  abs:               theta=28, RG-no-theta=36, ripple=1, mixed=34, inactive=33
  adaptive_envelope: theta=28, RG-no-theta=39, mixed=44, inactive=21

e10gw1
  abs:               theta=36, RG-no-theta=23, ripple=3, mixed=38, inactive=32
  adaptive_envelope: theta=36, RG-no-theta=30, mixed=40, inactive=26

f12m01
  abs:               theta=20, RG-no-theta=18, mixed=43, inactive=51
  adaptive_envelope: theta=21, RG-no-theta=27, mixed=50, inactive=34

k13m17
  abs:               theta=0, RG-no-theta=14, gamma=4, mixed=58, inactive=56
  adaptive_envelope: theta=0, RG-no-theta=33, gamma=5, mixed=55, inactive=39

k13m23
  abs:               theta=14, RG-no-theta=13, mixed=47, inactive=58
  adaptive_envelope: theta=14, RG-no-theta=20, mixed=52, inactive=46
```

Interpretation:

```text
E10gb1 and E10gH1 show a very strong P5-level two-subprocess split.
E10fV1 and E10gW1 show substantial but not complete replication.
f12m01 and k13m23 show weaker/partial replication.
k13m17 does not show strict theta_selective components under the current
thresholds, even though it has ripple/gamma-like components. This should be
treated as a dataset-specific diagnostic target before using its P8/P10 results
to argue for or against the slow-state / fast-trigger narrative.
```

Generated result folders:

```text
results\e10gb1_p5_p2_band_event_response_v2_selective_envelope_20260528\
results\e10fV1_p5_p2_band_event_response_v2_selective_envelope_20260528\
results\e10gh1_p5_p2_band_event_response_v2_selective_envelope_20260528\
results\e10gw1_p5_p2_band_event_response_v2_selective_envelope_20260528\
results\f12m01_p5_p2_band_event_response_v2_selective_envelope_20260528\
results\k13m17_p5_p2_band_event_response_v2_selective_envelope_20260528\
results\k13m23_p5_p2_band_event_response_v2_selective_envelope_20260528\
```

## 11. 2026-05-29 implementation update: P5/P8/P10 current run

P5 band-selectivity one-dataset runner:

```text
scripts\run_pipeline5_band_selectivity.py
```

Current 7-dataset batch runner:

```text
scripts\run_current7_p5_band_selectivity.ps1
```

Current 7-dataset P5 wrapper status:

```text
dataset   wrapper steps   failed steps
e10gb1    5               0
e10fV1    5               0
e10gh1    5               0
e10gw1    5               0
f12m01    5               0
k13m17    5               0
k13m23    5               0
```

Important QC note:

```text
k13m17 has no strict theta/RG method-k pair under the current strict rule.
The diagnostic-plate and paired-window scripts now treat this as a valid
no_strict_theta_rg_pair QC outcome, not as a runtime failure.
```

P8/P10 strict-band coupling implementation:

```text
scripts\analyze_p8_p10_strict_band_coupling.py
```

P8/P10 observable-separated trace-plate implementation:

```text
scripts\plot_bold_density_trace_by_dataset_observable.py
```

Current accepted design:

```text
selection-mode        = unique_density
top-hits-per-group    = 3
bold-feature-filter   = real
density classes       = raw_efun_density, dimred_efun_density
trace normalization   = whole-trace robust z
BOLD sign handling    = xcorr signed-corr aligned
mask visualization    = grey shaded masked samples
session visualization = dashed session borders
zoom policy           = fixed summary windows, not local-correlation-selected windows
```

Observable-stratified P8/P10 statistics are now part of the same implementation.
The original all-observable CSV/PNG outputs remain available as appendix/QC, but
the formal statistical interpretation should use the by-observable outputs:

```text
density_class_topN_fraction_by_observable.csv
dimred_strict_label_topN_fraction_by_observable.csv
dimred_strict_label_enrichment_by_observable.csv
method_k_strict_label_hit_counts_by_observable.csv

01_density_source_competition_by_observable\
02_dimred_strict_label_composition_by_observable\
03_dimred_strict_label_enrichment_by_observable\
04_dimred_effect_scatter_by_observable\
05_method_k_strict_label_hit_maps_by_observable\
06_lag_distribution_by_observable\
```

E10gb1 accepted preview:

```text
results\e10gb1_bold_density_trace_by_observable_unique_density_top3_preview_20260529\
```

Current output:

```text
tables:
results\pipeline8_10_strict_band_coupling_current\

figures:
E:\DataPons_processed\summary_figures\pipeline8_10_strict_band_coupling_current\
```

Current P8/P10 strict-band run summary:

```text
datasets                  = e10gb1, e10fV1, e10gh1, e10gw1, f12m01, k13m17, k13m23
P8 save tag               = xcorr_csplit_standardize_rmsenv_adaptive
P10 save tag              = dimred_xcorr_csplit_standardize_rmsenv_adaptive
P5 label transform        = adaptive_envelope
total xcorr hit rows      = 885248
dimred xcorr hit rows     = 817024
P5 label background rows  = 924
dimred hits missing label = 0
```

Remaining missing P8/P10 contexts:

```text
P8  k13m23  pv_gsvd100     global_svd100   missing_root
P8  k13m23  pv_gsvd100_ds  gsvd100_ds      missing_root
P10 k13m23  pv_gsvd100     global_svd100   missing_peaks
P10 k13m23  pv_gsvd100_ds  gsvd100_ds      missing_peaks
```

Targeted backfill attempt:

```text
scripts\script_run_k13m23_missing_standardized_csplit_p8_p10_for_p12.m
```

Outcome:

```text
P8: no current-best P7 candidate runs were found for K13m23 global_svd100/gsvd100_ds.
P10: no P9 dimred candidates were found for K13m23 global_svd100/gsvd100_ds.

Therefore this is not a P8/P10 xcorr cache problem. It is an upstream
BOLD/P7/P9 availability issue for those two K13m23 BOLD observable branches.
```

## 11.1 2026-05-30 run update: seven-dataset by-observable P8/P10 analysis

按照当前 SOP，七个 dataset 的 P8/P10 strict-band coupling analysis 已重新生成，
并且所有正式统计都增加了 BOLD-observable 分层。

Run command equivalent:

```text
scripts\analyze_p8_p10_strict_band_coupling.py
    --processed-root /mnt/e/DataPons_processed
    --datasets e10gb1 e10fV1 e10gh1 e10gw1 f12m01 k13m17 k13m23
```

Current tables:

```text
results\pipeline8_10_strict_band_coupling_current\
```

Current figures:

```text
E:\DataPons_processed\summary_figures\pipeline8_10_strict_band_coupling_current\
```

Run summary:

```text
total xcorr hit rows      = 885248
dimred xcorr hit rows     = 817024
P5 label background rows  = 924
dimred hits missing label = 0
```

By-observable output tables:

```text
density_class_topN_fraction_by_observable.csv
dimred_strict_label_topN_fraction_by_observable.csv
dimred_strict_label_enrichment_by_observable.csv
method_k_strict_label_hit_counts_by_observable.csv
```

By-observable output figure folders:

```text
01_density_source_competition_by_observable\
02_dimred_strict_label_composition_by_observable\
03_dimred_strict_label_enrichment_by_observable\
04_dimred_effect_scatter_by_observable\
05_method_k_strict_label_hit_maps_by_observable\
06_lag_distribution_by_observable\
```

Expected figure counts were verified:

```text
01 density competition by observable        = 16
02 strict-label composition by observable   = 16
03 strict-label enrichment by observable    = 16
04 effect scatter by observable             = 16
05 method-k strict-label hit maps           = 80
06 lag distribution by observable           = 16
```

The observable-separated trace plates were also generated for every available
dataset x P8/P10 context:

```text
E:\DataPons_processed\summary_figures\pipeline8_10_strict_band_coupling_current\00_observable_trace_plate\
```

Each dataset x pipeline folder contains:

```text
raw_efun_density unique-density top3 selected rows
raw_efun_density unique-density top3 full trace
raw_efun_density unique-density top3 fixed zoom
dimred_efun_density unique-density top3 selected rows
dimred_efun_density unique-density top3 full trace
dimred_efun_density unique-density top3 fixed zoom
trace_records.csv
```

Trace-plate coverage:

```text
e10gb1   P8/P10  full 4-observable coverage
e10fV1   P8/P10  full 4-observable coverage
e10gh1   P8/P10  full 4-observable coverage
e10gw1   P8/P10  full 4-observable coverage
f12m01   P8/P10  full 4-observable coverage
k13m17   P8/P10  full 4-observable coverage
k13m23   P8/P10  HP_svd100 and roi_mean only
```

Remaining missing contexts are unchanged and are upstream BOLD availability
issues, not P8/P10 analysis-script failures:

```text
P8  k13m23  global_svd100  missing_root
P8  k13m23  gsvd100_ds     missing_root
P10 k13m23  global_svd100  missing_peaks
P10 k13m23  gsvd100_ds     missing_peaks
```

2026-06-01 update: the K13m23 global_svd100 and gsvd100_ds gaps above have
now been backfilled through the BOLD side of the standardized complex-split
workflow. Treat the 2026-05-30 missing-context note as historical provenance,
not as the current audit state.

Current K13m23 BOLD backfill state:

```text
global_svd100  -> P7/P9/P8/P10 available
gsvd100_ds     -> P7/P9/P8/P10 available
HP_svd100      -> already available
roi_mean       -> already available
```

Naming note:

```text
pv_gsvd100     = global_svd100 BOLD branch
pv_gsvd100_ds  = gsvd100_ds BOLD branch
pv_hp100       = HP_svd100 BOLD branch
pv_roi         = roi_mean BOLD branch
```

This `pv_*` run-tag naming is still confusing and should be cleaned up in a
future provenance pass, but it is not a current computational blocker.

K13m23 BOLD model confidence notes:

```text
global_svd100  pat40 current best val = 0.145979, best epoch = 21
gsvd100_ds     pat40 current best val = 0.103137, best epoch = 109
```

A CPU probe for K13m23 global_svd100 with stronger regularization
(`reg = 0.01`) was worse:

```text
best val = 0.330054, best epoch = 104
```

So the early global_svd100 best checkpoint is not explained away by simply
using too little regularization. Keep global_svd100 as usable but lower
confidence than gsvd100_ds for K13m23 until a broader BOLD hyperparameter
comparison is run.

## 11.2 2026-06-01 run update: K13m23-complete seven-dataset P8/P10 refresh

After the K13m23 global_svd100 and gsvd100_ds BOLD branches were backfilled,
the current P8/P10 strict-band coupling and cross-session summaries were
regenerated with all seven standardized complex-split datasets:

```text
e10gb1
e10fV1
e10gh1
e10gw1
f12m01
k13m17
k13m23
```

Strict-band coupling:

```text
script:
scripts\analyze_p8_p10_strict_band_coupling.py

tables:
results\pipeline8_10_strict_band_coupling_current\

figures:
E:\DataPons_processed\summary_figures\pipeline8_10_strict_band_coupling_current\

total xcorr hit rows = 953344
```

P8 cross-session consistency:

```text
script:
scripts\summarize_pipeline8_cross_session_consistency.py

tables:
results\pipeline8_cross_session_consistency_current\

figures:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\cross_session_consistency\

readable hit rows = 176800
figure count      = 74
```

P10 cross-session consistency:

```text
script:
scripts\summarize_pipeline10_cross_session_consistency.py

tables:
results\pipeline10_cross_session_consistency_current\

figures:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p10\cross_session_consistency\

readable hit rows = 1615488
figure count      = 14
```

The conservative P5-stable-method-k-only P8/P10 coupling summary was also
regenerated. This analysis excludes K13m17 by default because K13m17 remains
visually and electrophysiologically suspicious at the P2/P5 level.

Stable method-k set:

```text
mds_k04
mds_k05
mds_k06
mds_k07
mds_k08
nmf_k04
umap_k04
umap_k05
umap_k06
umap_k08
```

Stable-method output:

```text
tables/figures:
results\pipeline8_10_strict_band_coupling_p5_stable_method_k_20260601\

summary copy:
E:\DataPons_processed\summary_figures\pipeline8_10_strict_band_coupling_p5_stable_method_k_20260601\
```

K13m23 observable-stratified trace plates were regenerated after the missing
BOLD branches were filled:

```text
tables/figures:
results\k13m23_bold_density_trace_by_observable_current\

summary copy:
E:\DataPons_processed\summary_figures\k13m23_bold_density_trace_by_observable_current\
```

The K13m23 trace refresh includes both all-method-k and P5-stable-method-k
views:

```text
P8  all method-k       raw + dimred density
P10 all method-k       raw + dimred density
P8  P5-stable method-k dimred density only
P10 P5-stable method-k dimred density only
```

Interpretation rule after this refresh:

```text
1. Use all-method-k P8/P10 figures to inspect the full search space.
2. Use P5-stable-method-k-only figures as the conservative main claim check.
3. Use K13m17-excluded summaries for cross-dataset claims unless explicitly
   discussing K13m17 as a questionable dataset.
4. Keep BOLD observable separated: global_svd100/gsvd100_ds/HP_svd100/roi_mean
   should not be collapsed in the first-pass interpretation.
```

## 11.3 2026-06-02 update: ROI profile consistency review

After the P5-stable P8/P10 coupling summaries were regenerated, ROI profile
consistency was moved into the second-stage review.  This step asks:

```text
For high-coupling BOLD modes/components, is the ROI profile spatially similar
across datasets?
```

Implementation:

```text
export numeric ROI vectors:
scripts\export_p8_p10_roi_profile_consistency_sources.m

exact-match ROI consistency:
scripts\plot_p8_p10_roi_profile_consistency.py

aggregated ROI consistency:
scripts\plot_p8_p10_roi_profile_consistency_aggregated.py
```

The aggregated script was added because exact component/density names are too
strict for P10.  The formal second-stage interpretation should use the
aggregated view first, then use exact-match rows as QC/appendix.

Current combined ROI vector tables:

```text
results\pipeline_roi_profile_consistency_current\p8_roi_profiles_long_combined_20260602.csv
rows = 87841

results\pipeline_roi_profile_consistency_current\p10_roi_profiles_long_combined_20260602.csv
rows = 830209
```

Current figures:

```text
exact P8:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\roi_profile_consistency\

exact P10:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p10\roi_profile_consistency\

aggregated P5-stable P8/P10:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_roi_profile_consistency_aggregated_20260602\p5_stable\
```

Current aggregated P5-stable output:

```text
P8 aggregated groups       = 80
P8 consistency rows        = 60
P10 aggregated groups      = 560
P10 consistency rows       = 420
summary heatmaps           = 5
ROI profile example plates = 32
```

Important ROI-mapping limitations:

```text
k13m17:
    P8 ROI export has no current P7 candidate.
    P10 ROI export is available after K13.m17 ROI-directory alias support.
    k13m17 remains a suspicious dataset and should not drive the main claim.

k13m23:
    HP_svd100 ROI export is available.
    global_svd100 and gsvd100_ds are ROI-mapping blocked because the current
    BOLD_POST files do not provide a usable source-to-voxel mapping for
    global_svd100### / gsvd100_ds### source labels.
    roi_mean is ROI-mapping blocked because the source label LC_CGn_mean is
    not present in the current K13.m23 roiTs label list.
    HP_svd100 is a local/single-region observable, so it is not suitable for
    multi-ROI Pearson profile consistency.
```

Current ROI consistency result:

```text
P8:
    ROI profile consistency is strongest and clearest.
    gsvd100_ds efun is the most spatially repeatable group.
    Top P8 P5-stable examples are gsvd100_ds efun with mds_k04-k08:
        mean pairwise ROI corr ~= 0.68
    global_svd100 efun is also high:
        mean pairwise ROI corr ~= 0.62-0.66
    deconv_efun is still repeatable but slightly weaker:
        mean pairwise ROI corr ~= 0.58-0.66 depending on observable/method-k

P10:
    ROI profile consistency is moderate rather than strong.
    gsvd100_ds/global_svd100 are clearly better than roi_mean.
    gsvd100_ds efun with BOLD P9 mds_k08 and P5-stable LFP method-k gives the
    best current P10 ROI consistency:
        mean pairwise ROI corr ~= 0.57-0.58
    gsvd100_ds deconv_efun top groups are lower:
        mean pairwise ROI corr ~= 0.49-0.51
    roi_mean is low:
        mean pairwise ROI corr ~= 0.29-0.30
```

Interpretation:

```text
ROI summary supports that the BOLD modes involved in top BLP coupling have
some cross-dataset spatial structure, especially for P8 gsvd100_ds/global_svd100
efun.

It does not by itself prove the slow-state / intrinsic-trigger story.
The main mechanistic evidence remains:
    1. P5 theta vs ripple-gamma subprocess gate
    2. P8/P10 xcorr label enrichment and raw-density competition
    3. whole-trace behavior
    4. ROI profile consistency as a second-stage spatial sanity check
```

### 11.4 2026-06-03 update: ROI consistency must be subprocess-stratified

The previous `aggregated ROI consistency` view is now treated as QC only.
It averages over all dimred density top-hit components inside the same
method-k/density source, so theta, ripple-gamma, mixed, and inactive P5
components can be blended into one ROI profile.  That is not the right
scientific question.

Formal ROI consistency should ask:

```text
For BOLD modes/components that couple strongly to a specific P5 subprocess,
is the resulting ROI profile spatially consistent across datasets?
```

Therefore ROI profiles must be exported and plotted with `strict_label_group`
as a first-class grouping variable:

```text
theta
ripple_gamma
mixed
inactive   # QC only; do not use as a main mechanistic claim
```

Principle: keep the three analysis gates separate.

```text
Gate A: P5 subprocess pass
    Unit:
        dataset x BLP condition x dimred method-k
    Question:
        Does this dimred BLP component space contain both a theta-like
        component and a ripple-gamma-like component?
    Meaning:
        The method-k is interpretable at the BLP subprocess level.
    Non-meaning:
        It does not guarantee that every subprocess will enter P8/P10 top
        coupling hits for every BOLD observable.

Gate B: P8/P10 subprocess coupling
    Unit:
        dataset x BOLD observable x BOLD efun family x dimred density hit
    Question:
        Among actual top xcorr hits, which P5 subprocess does the BOLD feature
        couple to most strongly?
    Meaning:
        This is the BOLD-LFP coupling evidence.
    Non-meaning:
        It is not a spatial consistency test yet.

Gate C: ROI profile consistency
    Unit:
        dataset x BOLD observable x BOLD efun family x P5 subprocess x method-k
    Question:
        For BOLD modes/components that couple to the same P5 subprocess, are
        their ROI profiles similar across datasets?
    Meaning:
        This is a second-stage spatial sanity check on the coupling result.
    Non-meaning:
        It cannot rescue a missing P5/P8/P10 result, and it should not average
        across different subprocess labels.
```

This is why P5-pass and ROI availability are not identical:

```text
P5-pass means:
    the method-k contains candidate theta and ripple-gamma components.

ROI availability additionally requires:
    1. the subprocess-specific component appears in P8/P10 top xcorr hits;
    2. the corresponding BOLD mode/component has a usable ROI mapping;
    3. the dataset has a current P7/P9 context discoverable by the exporter;
    4. enough datasets have the same analysis key to support a cross-dataset
       ROI profile correlation.
```

Main figures should not use available-case denominators silently.  The plotting
rule is:

```text
Every main overview heatmap must either:
    a. enforce a fixed minimum n_datasets threshold, or
    b. display n_datasets / denominator explicitly in the figure.

Current default:
    plot_min_datasets = 5

Future after K13m23 is connected:
    use plot_min_datasets = 6 for the six-dataset main ROI figure,
    or show 5/6 and 6/6 explicitly if comparing availability.
```

K13m23 status for this ROI step:

```text
K13m23 has P8 strict xcorr rows and current BOLD_POST files.
It was not included in the current P8 by-subprocess ROI export because the
exporter did not discover K13m23 through the current P7 candidate registry.
This is an implementation/registry issue, not a scientific exclusion.
The current P8 ROI by-subprocess figure is therefore a 5 non-K13 dataset
temporary main figure, not a final six-dataset ROI consistency result.
```

Implementation update:

```text
scripts\export_p8_p10_roi_profile_consistency_sources.m
    group_by_strict_subprocess = true

scripts\plot_p8_p10_roi_profile_consistency_by_subprocess.py
```

The exporter now reads:

```text
results\pipeline8_10_strict_band_coupling_current\p8_p10_strict_band_hits_long.csv
```

and writes subprocess-specific ROI vectors:

```text
results\pipeline_roi_profile_consistency_current\p8_roi_profiles_by_subprocess_long.csv
results\pipeline_roi_profile_consistency_current\p10_roi_profiles_by_subprocess_long.csv
```

New provenance columns:

```text
strict_label_group
strict_label
subprocess_rank
subprocess_peak_abs_corr_sum
subprocess_peak_abs_corr_fraction
subprocess_n_hit_rows
subprocess_density_indices
```

The key interpretation change:

```text
Do not average all xcorr hits before asking about ROI profile consistency.
First split by the P5 subprocess that generated the large xcorr hit, then
compare ROI profiles within the same subprocess group.
```

Current P8 status:

```text
P8 subprocess ROI export completed.
P8 P5-stable by-subprocess ROI summary:
    subprocess aggregated groups = 312
    consistency rows             = 208

workspace figures:
results\pipeline_roi_profile_consistency_current\figures_by_subprocess_p8_only\p5_stable\

summary figures:
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_roi_profile_consistency_by_subprocess_20260603_p8_only\p5_stable\
```

Important plotting rule:

```text
Main heatmaps must use a fixed availability threshold.
Current P8 main figures use plot_min_datasets = 5.
Rows with n_datasets < 5 remain in the CSV for QC but are not averaged into
the main overview heatmaps.
```

This rule prevents variable-denominator means, where a high value from only
two datasets could be visually mixed with true cross-dataset rows.

Current P8 observation, using `n_datasets >= 5` as the main-read filter:

```text
ripple_gamma:
    strongest current spatial consistency.
    top P8 examples are mostly gsvd100_ds/global_svd100 efun with
    mds_k04-k08 and umap_k04-k06/k08.
    top mean pairwise ROI corr ~= 0.69.

theta:
    also spatially consistent, especially gsvd100_ds efun with mds/umap.
    top mean pairwise ROI corr ~= 0.67.

mixed:
    present and sometimes high, but less clean for the mechanistic claim.

inactive:
    can show high ROI corr, but mostly with only two datasets in the current
    export, so it is QC-only and should not be used as a main conclusion.
```

Current P10 status:

```text
P10 subprocess ROI export is not yet complete.
The naive exporter rebuilds each P9 context one by one and is too slow for the
full 7-dataset run.  It needs optimization before P10 by-subprocess ROI
consistency can be treated as a current result.
```

## 11.5 2026-06-03 update: ROI subprocess separation check

After the first P8 by-subprocess ROI consistency figures, the next necessary
question is not only whether each subprocess has cross-dataset spatial
consistency.  We must also ask whether different subprocesses are spatially
separable at the ROI-summary level.

Formal question:

```text
Is corr(theta-linked ROI, ripple-gamma-linked ROI) lower than:
    corr(theta-linked ROI across datasets), and
    corr(ripple-gamma-linked ROI across datasets)?

If the cross-subprocess ROI correlation is similar to the within-subprocess
cross-dataset correlation, then ROI summary does not support spatial separation
of the two subprocesses.
```

Implementation:

```text
scripts\analyze_p8_roi_subprocess_separation.py
```

Input:

```text
results\pipeline_roi_profile_consistency_current\p8_roi_profiles_by_subprocess_long.csv
```

Outputs:

```text
results\pipeline_roi_profile_consistency_current\p8_roi_subprocess_separation_20260603\
    p8_roi_subprocess_pairwise_correlations.csv
    p8_roi_subprocess_pairwise_summary.csv
    p8_roi_subprocess_pairwise_summary_by_feature.csv
    p8_roi_subprocess_pairwise_by_observable.csv
    p8_roi_subprocess_pairwise_by_observable_method_k.csv
    p8_theta_rg_roi_separation_by_method_k.csv

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_roi_subprocess_separation_20260603\
    01_p8_roi_subprocess_confusion_cross_dataset__efun.png
    01_p8_roi_subprocess_confusion_cross_dataset__deconv_efun.png
    02_p8_roi_cross_subprocess_same_dataset__efun.png
    02_p8_roi_cross_subprocess_same_dataset__deconv_efun.png
    03_p8_theta_rg_roi_separation_by_observable.png

    by_observable_feature\
        01_p8_roi_subprocess__<observable>__<feature>_confusion_cross_dataset.png
        02_p8_roi_subprocess__<observable>__<feature>_same_dataset.png

    by_observable_feature_method_k\
        p8_theta_rg_roi_separation_method_k__<observable>__<feature>.png
```

Important rule:

```text
Do not average efun and deconv_efun in the formal ROI separation figure.
They have different meanings:
    efun        = BOLD state/eigenfunction feature
    deconv_efun = deconvolved/residual perturbation-like feature

Any combined efun+deconv plot is QC-only/legacy and should not be used for
the main interpretation.

Also do not collapse BOLD observable or dimred method-k in the formal readout.
The formal read order is:
    1. BOLD observable: global_svd100, gsvd100_ds, roi_mean, and HP_svd100
       when ROI mapping is available.
    2. BOLD feature family: efun vs deconv_efun.
    3. dimred method-k: e.g. nmf_k04, mds_k04-k08, umap_k04/k05/k06/k08.

The method-k heatmap has three panels:
    within theta/RG
    theta vs RG cross-dataset
    separation = within - cross

Only a clearly positive separation panel would support ROI-level spatial
separation.  Near-zero separation means theta-linked and RG-linked BOLD ROI
profiles are spatially similar for that observable and method-k.
```

Definitions:

```text
within-subprocess cross-dataset correlation:
    same observable, same efun/deconv feature family, same density method-k,
    same strict_label_group, different datasets.

cross-subprocess cross-dataset correlation:
    same observable, same efun/deconv feature family, same density method-k,
    different strict_label_group, different datasets.

cross-subprocess same-dataset correlation:
    same dataset and same analysis key, but different strict_label_group.
```

Current P8 result, using the current five non-K13 ROI-available datasets and
keeping efun/deconv_efun separated:

```text
efun:
    within theta/RG mean ROI corr:
        0.608
    theta vs ripple-gamma cross-dataset ROI corr:
        0.603
    separation = within - cross:
        +0.005
    theta vs ripple-gamma same-dataset ROI corr:
        0.914

deconv_efun:
    within theta/RG mean ROI corr:
        0.585
    theta vs ripple-gamma cross-dataset ROI corr:
        0.588
    separation = within - cross:
        -0.003
    theta vs ripple-gamma same-dataset ROI corr:
        0.928
```

Observable-level separation is also near zero:

```text
global_svd100 efun:        within - thetaRG_cross = +0.006
global_svd100 deconv_efun: within - thetaRG_cross = -0.004
gsvd100_ds efun:           within - thetaRG_cross = +0.003
gsvd100_ds deconv_efun:    within - thetaRG_cross = -0.001
roi_mean efun:             within - thetaRG_cross = +0.006
roi_mean deconv_efun:      within - thetaRG_cross = -0.004
```

Interpretation:

```text
P8 ROI summary currently supports cross-dataset spatial repeatability of the
BOLD modes involved in BLP coupling.

It does not support ROI-level spatial separation between theta-linked and
ripple-gamma-linked subprocesses.

Therefore, the two-subprocess interpretation should continue to rely on:
    1. P5 band-event subprocess separation,
    2. P8/P10 label enrichment and raw-density competition,
    3. raw efun timescale and whole-trace behavior,
not on ROI summary alone.

ROI summary can still be used as a spatial repeatability/QC layer, but not as
the primary evidence that theta and ripple-gamma subprocesses occupy different
BOLD ROI spaces.
```

## 11.6 2026-06-04 update: P8 raw_csplit_q070 efun-vs-deconv ROI footprint check

After the subprocess-stratified ROI check, add one narrower P8 ROI readout that
does not average over all density sources:

```text
raw_csplit_q070 x BOLD observable x efun_real/deconv_real
```

Formal question:

```text
When the BOLD feature is selected by top xcorr against raw complex-split BLP
eigenfunction density, do BOLD efun and BOLD deconv_efun select the same BOLD
ROI/mode footprint, or different footprints?
```

This check is intentionally narrower than the P5-label ROI analysis:

```text
It does not test theta-vs-ripple-gamma spatial separation directly.
It tests whether raw-density coupling separates BOLD state-like efun from
deconvolved/residual perturbation-like deconv_efun.
```

Implementation:

```text
scripts\plot_p8_raw_density_efun_vs_deconv_roi_similarity.py
scripts\plot_p8_raw_csplit_focus_panel.py
```

Inputs:

```text
P8 top xcorr CSV:
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\feature\efun\xcorr_top__raw_csplit_q070__efun_real.csv
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\feature\deconv_efun\xcorr_top__raw_csplit_q070__deconv_real.csv

P7 all-mode ROI vectors:
results\pipeline_roi_profile_consistency_current\p7_intrinsic_bold_efun_roi_profiles_long.csv
```

Current outputs:

```text
results\pipeline_roi_profile_consistency_current\p8_raw_density_efun_vs_deconv_roi_similarity_20260603\
    p8_raw_density_efun_vs_deconv_roi_summary.csv
    p8_raw_density_efun_vs_deconv_roi_pairwise.csv
    p8_raw_density_efun_vs_deconv_mode_index_overlap_summary.csv

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_raw_density_efun_vs_deconv_roi_similarity_20260603\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_raw_csplit_q070_focus_20260603\
    p8_raw_csplit_q070_efun_vs_deconv_roi_focus_panel.png
```

Current data coverage after the 2026-06-04 e10gw1 backfill:

```text
datasets included = e10gb1, e10fV1, e10gh1, e10gw1, f12m01
raw top rows read = 300
missing raw top files = 0
```

The e10gw1 gap was specifically missing standard-name P8 raw-density top CSVs
for `xcorr_top__raw_abs_q070__...` and `xcorr_top__raw_csplit_q070__...`.
It was backfilled by running a raw-only P8 xcorr with:

```text
cfg_name = E10gW1
observable_modes = global_svd100, gsvd100_ds, roi_mean
density_source_kinds = raw_abs_density, raw_complex_split_density
xcorr_save_tag = xcorr_rawfill
figures/maps/ROI summaries disabled
```

The generated `xcorr_rawfill_*__raw*_q070*.csv` files were then copied to the
standard `xcorr_*__raw*_q070*.csv` names only where the standard targets were
missing.  Existing P8 dimred/adaptive outputs were not overwritten.

Current P8 `raw_csplit_q070` ROI-footprint readout:

```text
global_svd100:
    same-dataset efun-vs-deconv ROI corr  ~= 0.63
    cross-dataset efun within corr        ~= 0.48
    cross-dataset deconv within corr      ~= 0.40
    selected BOLD mode-index overlap      ~= 0.45

gsvd100_ds:
    same-dataset efun-vs-deconv ROI corr  ~= 0.93
    cross-dataset efun within corr        ~= 0.44
    cross-dataset deconv within corr      ~= 0.41
    selected BOLD mode-index overlap      ~= 0.80

roi_mean:
    same-dataset efun-vs-deconv ROI corr  ~= 0.40
    cross-dataset efun within corr        ~= 0.18
    cross-dataset deconv within corr      ~= 0.51
    selected BOLD mode-index overlap      ~= 0.00
```

Interpretation:

```text
global_svd100 and gsvd100_ds:
    efun and deconv_efun often select overlapping or highly similar BOLD ROI
    subspaces.  These observables support a stable BOLD subspace, but not a
    clean efun-vs-deconv ROI split.

roi_mean:
    efun and deconv_efun selected mode indices do not overlap, and deconv_efun
    has stronger cross-dataset ROI repeatability than efun.  This is currently
    the strongest P8 raw-density ROI evidence that BOLD state-like efun and
    deconv/residual perturbation-like efun can occupy different ROI footprints.
```

Use this as a second-stage spatial sanity check only.  The main two-subprocess
claim still comes from:

```text
P5 band-event subprocess gate
P8/P10 label enrichment and raw-density competition
raw efun timescale and whole-trace behavior
```

## 12. Plotting/QC policy: no trace display clipping

For new-dataset visual QC, raw BLP traces, BOLD observable traces, and P5 component-activity traces should be inspected without hard display clipping. Large excursions should remain visible, with y-limits expanding to the un-clipped trace. This matters for K13m17-like cases where theta-event windows can look saturated in older clipped figures.

This is a figure/provenance policy only. It does not change P2 event detection, P5 dimred reduction, P5 density construction, or P8/P10 xcorr values.

Canonical archive:

```text
docs/plotting_no_trace_clipping_archive_2026-05-28.md
```

For K13m17, the P2 official top-window figures and P5 component-activity top30 diagnostics were regenerated after this change. Pre-2026-05-28 raw/trace browsing figures should be treated as stale for visual-amplitude interpretation.

## 13. New dataset reading order

For any new standardized complex-split dataset, do not begin from a large P8/P10
summary panel.  Use the gate order below.

### Step 1. Mainline completeness

Confirm that the dataset has the current mainline artifacts:

```text
P1/P2:
    bandpass BLP and P2 theta/gamma/ripple event detection.

P3/P4:
    standardized complex-split BLP EDMD/ResKoopNet output.

P5:
    standardized complex-split dimred results and P5 band-event subprocess
    labels for abs and adaptive_envelope activity.

P7:
    current-best BOLD_POST for each BOLD observable, keeping raw/standardized
    provenance separate.

P8/P10:
    xcorr top CSVs for raw density, dimred density, and the current save tags
    required by the SOP.
```

If a step is missing, backfill that step before interpreting downstream figures.
P8/P10 summary-figure cache is never evidence of completeness.

### Step 2. Dataset QC

Look first at raw event/data quality:

```text
P2 top theta/gamma/ripple windows
P5 component top-window diagnostic plates
```

Classify the dataset:

```text
usable
suspicious_but_usable
exclude_from_cross_dataset_conclusion
```

For suspicious data, inspect clipping/saturation without display clipping before
using the dataset in cross-dataset conclusions.

### Step 3. P5 subprocess gate

Use only P2-band strict labels for the main gate:

```text
strict selectivity plane
strict two-subprocess candidate map
method-k diagnostic plates
```

The dataset passes the P5 gate when at least one method-k has both:

```text
theta_selective component
ripple_gamma_no_theta / gamma_selective / ripple_selective component
```

Stable datasets should show multiple passing method-k values, not a single
fragile component pair.

### Step 4. P8/P10 coupling gate

Only after P5 passes, inspect P8/P10 by BOLD observable and efun family:

```text
global_svd100, gsvd100_ds, HP_svd100, roi_mean
efun vs deconv_efun
top1/top3/top5/top10/top20 sensitivity
```

The main coupling questions are:

```text
Does deconv_efun preferentially couple to ripple-gamma / fast BLP density?
Does efun show slow/session-state behavior through raw efun density and
whole-trace evidence?
Does this pattern repeat across observables and datasets?
```

Do not collapse BOLD observable or efun/deconv_efun in the first readout.

### Step 5. Spatial sanity check

Use ROI only after the coupling pattern is clear:

```text
P7 intrinsic all-mode BOLD ROI confusion
P8/P10 selected-mode overlays on the P7 ROI matrix
P8 by-subprocess ROI consistency
P8 raw_csplit_q070 efun-vs-deconv ROI footprint check
```

Interpretation rule:

```text
High ROI repeatability means the selected BOLD modes are spatially consistent.
It does not by itself prove theta/RG spatial separation.

Spatial separation requires within-subprocess ROI correlation to exceed
cross-subprocess ROI correlation after keeping observable, efun/deconv, and
method-k separate.
```

Current conservative working conclusion:

```text
P5 can often find two BLP subprocesses.
deconv_efun coupling is the cleaner route to the fast ripple-gamma / trigger
interpretation.
efun coupling is less clean for theta and should be supported by raw efun
timescale plus whole-trace slow-state behavior.
ROI summary is a QC/spatial-repeatability layer, not the primary evidence.
```

### Step 6. Figure layout policy

All first-pass summary figures should be generated in a wide, ultrawide-monitor
friendly layout. The default target is browsing and side-by-side comparison,
not a portrait report page.

Current plotting policy:

```text
dataset panels:
    put datasets in one horizontal row when possible

dataset-pair matrices:
    put dataset pairs on the x-axis and method-k on the y-axis

trace figures:
    use wide full-trace and zoom-trace panels
    keep session borders and masked regions visible

ROI heatmaps:
    preserve the anatomical P7 ROI order
    do not sort ROIs by effect unless the filename explicitly says it is a
    top-ROI plot
    for all-ROI-by-dataset panels, use per-dataset color scales by default
    because the first-pass goal is within-dataset ROI pattern inspection, not
    cross-dataset absolute amplitude comparison
    use a shared color scale only when the filename or command explicitly says
    shared-scale / absolute-amplitude comparison

P5 diagnostic window sheets:
    prefer more columns and fewer rows for top-window examples
```

This is a figure/provenance policy only. It does not change P2 event
detection, P5 labels, P8/P10 xcorr values, selected BOLD mode indices, or ROI
profile statistics.
