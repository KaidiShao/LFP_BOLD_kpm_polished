# Pipeline12 standardized complex-split consistency check

日期：2026-05-25

状态：新增设计。第一版只聚焦一个 BLP condition：

```text
complex_split_projected_vlambda_standardize
```

对应 activity / density 口径：

```text
RMS-envelope adaptive activity
threshold ratio q070
```

## 1. 为什么要做 P12

现在的问题已经不是单纯找最大 xcorr，而是判断一个新 dataset 是否复现
E10gb1 里看到的结构：

```text
theta-like slow/state subprocess
ripple-gamma-like intrinsic-trigger subprocess
```

为了减少版本和 condition 混乱，P12 第一版只看 standardized complex split。
这样每个新 dataset 的问题变成：

```text
在 standardized complex split 下，
dimred/raw/event density 谁更能解释 BOLD efun/deconv_efun？

高 xcorr 的 dimred BLP components 是 theta-like 还是 ripple-gamma-like？

高 xcorr 的 raw BLP efun density 来自 slow modes 还是 fast modes？
```

## 2. P12 要回答的核心问题

### 2.1 Density source competition

比较三类 density source：

```text
event density
raw BLP efun density
dimred BLP efun density
```

主要问题：

```text
dimred BLP efun density x BOLD efun/deconv_efun
是否比 event density 更高？

dimred BLP efun density
是否比 raw BLP efun density 更高、更稳定、或更可解释？
```

这一步必须分开看：

```text
BOLD efun
BOLD deconv_efun
```

也必须分开看 BOLD observable：

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

### 2.2 Dimred label interpretation

对进入 topN 的 dimred BLP density，接 P5 component label：

```text
theta-like
ripple-gamma no-pure-theta
mixed / partial / other
label missing
```

主要问题：

```text
和 BOLD efun 相似的 dimred BLP component
是否更偏 theta/state-like？

和 BOLD deconv_efun 相似的 dimred BLP component
是否更偏 ripple-gamma/intrinsic-trigger-like？
```

### 2.3 Raw efun timescale distribution

对进入 topN 的 raw BLP efun density，接 raw efun index 和 eigenvalue-derived
timescale metadata。

主要问题：

```text
和 BOLD efun 相似的 raw density
是否来自更 slow/state-like modes？

和 BOLD deconv_efun 相似的 raw density
是否来自更 fast/activity-like modes？
```

这一步不能只看 top10。P12 默认看：

```text
top1 / top3 / top5 / top10 / top20 / top50
```

### 2.4 P8 vs P10 consistency

P8 和 P10 都要做同样检查。

P8 的问题是：

```text
BOLD P7 efun/deconv_efun target
是否和 BLP density source 强耦合？
```

P10 的问题是：

```text
BOLD P9 dimred efun/deconv_efun target
是否和 BLP density source 强耦合？
```

如果 P8 和 P10 都给出相似结论，这个 dataset 和 E10gb1 的一致性更强。

### 2.5 BOLD observable-specific interpretation

P12 不再把所有 BOLD observable 混成一个总图。

当前解释规则是：

```text
global/state-like BOLD efun
    更可能对应 slow theta/state subprocess

HP/local BOLD efun
    不一定是 slow variable；
    可能直接对应 local fast / ripple-gamma-like activity

BOLD deconv_efun
    更常对应 residual / intrinsic perturbation / trigger，
    但仍需要按 observable 单独验证
```

## 3. P12 输入

P12 不从 P4 summary 直接计算结论。它消费已经存在的 P5/P8/P10 输出，并
审计缺失项。

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

P4 summary MAT 只用于 provenance/audit。若 summary MAT 不包含 `efuns`，
且同目录没有 `*_outputs_*.mat` full chunks，则 P12 会标记为：

```text
blocked_missing_full_p4_outputs
```

## 4. P12 输出

默认输出在：

```text
results\pipeline12_standardized_csplit_consistency\<dataset>\
```

主要表：

```text
availability.csv
missing_requirements.csv
p12_top_rows.csv
p12_density_class_best_score.csv
p12_density_class_topN_membership.csv
p12_dimred_top_rows_with_labels.csv
p12_raw_top_rows_with_timescale.csv
p12_method_k_best_score.csv
summary.md
```

主要图：

```text
figures\availability.png
figures\density_class_best_score__<pipeline>__<feature_family>.png
figures\density_class_topN_membership.png
figures\dimred_label_fraction_top50.png
figures\method_k_best_score__<pipeline>__<feature_family>.png
figures\raw_timescale_ecdf_by_topN.png
```

如果没有足够数据，图会明确显示 `No rows available` 或 `metadata missing`，
不会用 legacy/current nonstandard 结果假装 standardized result。

## 5. 新 dataset 与 E10gb1 一致性的判据

一个新 dataset 和 E10gb1 更一致，需要同时满足：

```text
1. standardized csplit P5 能找到 theta component 和 ripple-gamma component
2. dimred density 在 P8/P10 中相对 event density 有竞争力
3. 高 xcorr dimred density 的 label 随 BOLD-side target 分化
4. raw density 的 topN timescale 在 efun vs deconv_efun 间有分化
5. 以上 pattern 按 BOLD observable 分开后仍可解释
```

最强的结果不是数值完全相同，而是 pattern 相同：

```text
global/state-like BOLD efun
    对应 theta/slow-state density

HP/local 或 deconv_efun
    对应 ripple-gamma/intrinsic-trigger density
```

## 6. 当前 E10gH1 注意事项

用户提供的 E10gH1 standardized csplit P4 summary 路径是：

```text
E:\autodl_results_new\e10gh1\mlp\outputs\mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gh1_seed1234_projected_vlambda_complex_split\e10gh1_low50_high250_g2_complex_split_single_Python_resdmd_Layer_100_Ndict_713_summary.mat
```

这个文件需要先确认是否含有 full `efuns`。如果只是 summary，而没有
`*_outputs_*.mat` chunks，则不能直接进入 P5/P8/P10 standardized csplit
检查，需要先导出 full EDMD outputs。

## 7. E10gH1 current run log: 2026-05-25

本轮已把 E10gH1 的 standardized complex-split P12 所需数据补齐，当前审查结果为：

```text
results\pipeline12_standardized_csplit_consistency\e10gh1\summary.md
Missing requirements: 0
```

已补齐/确认的内容：

```text
P4 full EDMD chunks:
E:\DataPons_processed\derived_autodl_results_standardize\e10gh1\mlp\outputs\mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gh1_seed1234_projected_vlambda_complex_split\
882 *_outputs_*.mat chunks

P5 standardized csplit RMS-envelope adaptive:
raw density: 1 MAT
dimred density: SVD/NMF/MDS/UMAP k03:k08, 24/24 complete
peak statistics: 24/24 complete
dimred component labels: 264 labels loaded by P12

P8 standardized csplit xcorr:
104 peak CSV files
feature families covered: efun and deconv_efun

P9 current-best deconv_real:
32 MAT files newly generated for current P7 best checkpoints
features/methods: deconv_real x SVD/NMF/MDS/UMAP x k05/k08 x 4 BOLD observables

P10 standardized csplit xcorr:
1664 peak CSV files
feature families covered: efun and deconv_efun
```

一个重要版本处理：

```text
旧的 E10gH1 P9 deconv_real results 来自 stale/non-current-best P7 run。
这些结果没有删除，已移动到：
E:\DataPons_processed\e10gh1\pipeline9_bold_eigenfunction_reduction_legacy_stale_deconv_pre_current_best_20260525\

当前 P9/P10 deconv_real 结果已经重新绑定到 current-best P7 checkpoint。
```

P12 审查脚本同步更新：

```text
scripts\run_pipeline12_standardized_csplit_consistency.py
```

更新点：

```text
1. 增加 --p4-full-output-dir，用于识别 summary MAT 外部导出的 full EDMD chunks。
2. P12 availability 现在显式检查 P8/P10 是否同时覆盖 efun 和 deconv_efun。
3. 当前 E10gH1 P12 输出目录里的 legacy incomplete placeholder figures 已移入 figures\legacy_incomplete_placeholders_20260525。
```

## 8. E10gb1 vs E10gH1 same-mouth P12 comparison: 2026-05-25

为了判断 E10gH1 是否复现 E10gb1，已把 E10gb1 也重新跑成同一套 P12 口径：

```text
results\pipeline12_standardized_csplit_consistency\e10gb1\summary.md
results\pipeline12_standardized_csplit_consistency\e10gh1\summary.md
```

两者当前都是：

```text
Missing requirements: 0
```

补齐过程中发现 E10gb1 原本缺 current-best P10 deconv coverage：

```text
missing before fix:
P10 deconv_efun only had pv_gsvd100 and pv_hp100
missing pv_gsvd100_ds and pv_roi
```

处理方式：

```text
旧 E10gb1 P9 deconv_real gsvd100_ds/roi_mean results 来自 stale/non-current-best P7。
这些结果没有删除，已移动到：
E:\DataPons_processed\e10gb1\pipeline9_bold_eigenfunction_reduction_legacy_stale_deconv_pre_current_best_20260525\

随后为 current-best P7 重新生成：
P9 deconv_real: gsvd100_ds and roi_mean, SVD/NMF/MDS/UMAP k05/k08
P10 standardized csplit deconv xcorr: gsvd100_ds and roi_mean, numeric only, no figures
```

P12 审查脚本也因此增强：

```text
P12 now checks per-pipeline/per-feature observable coverage:
P8 efun
P8 deconv_efun
P10 efun
P10 deconv_efun

Each must cover:
pv_gsvd100
pv_gsvd100_ds
pv_hp100
pv_roi
```

当前同口径比较结论：

```text
1. 高层结论部分一致：
   event density 基本不进入 top10；
   raw/dimred BLP efun density 都明显比 event density 更有竞争力。

2. 但不是 E10gb1 的简单复现：
   E10gb1 的 P10 deconv top10 强烈偏 dimred density；
   E10gH1 的 P10 deconv top10 更偏 raw density。

3. raw timescale split 也不完全一致：
   E10gb1: P8 efun top5/top10 明显偏 long-timescale mode，deconv 明显偏 fast mode。
   E10gH1: deconv 仍偏 fast，但 efun 只在少数 top hits 中出现 long-timescale mode，
            top10/top50 后整体更偏 fast。

4. dimred label 方向不同：
   E10gb1 top dimred xcorr labels mostly theta/partial_or_inactive under the current primary-label table。
   E10gH1 top dimred xcorr labels contain many ripple_selective components。

5. 因此目前更准确的表述是：
   E10gH1 supports the general density-coupling idea,
   but it does not fully reproduce the E10gb1 slow-variable / intrinsic-trigger split.
```

## 9. TopN handling update: 2026-05-25

Top50 会把很多次级 hit 混进来，不适合作为第一判断依据。P12 现在默认只画：

```text
top1 / top3 / top5 / top10 / top20
```

对应脚本：

```text
scripts\run_pipeline12_standardized_csplit_consistency.py
scripts\compare_pipeline12_standardized_csplit_topn.py
```

P12 单 dataset 输出新增：

```text
figures\dimred_label_fraction_top01.png
figures\dimred_label_fraction_top03.png
figures\dimred_label_fraction_top05.png
figures\dimred_label_fraction_top10.png
figures\dimred_label_fraction_top20.png
figures\dimred_label_fraction_by_topN__<pipeline>__<feature_family>.png
figures\raw_timescale_summary_by_topN.png
```

跨 dataset topN 对比输出：

```text
results\pipeline12_standardized_csplit_consistency\topn_comparison_e10gb1_e10gh1\
```

更推荐作为主查看入口的是 matched-panel 版本：

```text
results\pipeline12_standardized_csplit_consistency\matched_topn_panels_e10gb1_e10gh1\
```

这个目录里的每张图都使用同一套布局和坐标口径：

```text
rows = datasets: e10gb1, e10gh1
columns = P8 efun, P8 deconv_efun, P10 efun, P10 deconv_efun
topN = top1, top3, top5, top10, top20
```

当前 matched-panel 图：

```text
01_density_class_fraction_by_topN_matched.png
02_dimred_label_fraction_by_topN_matched.png
03_raw_timescale_median_by_topN_matched.png
04_raw_timescale_slow_fraction_by_topN_matched.png
05_best_score_by_density_class_matched.png
```

重要口径修正：

```text
raw-timescale topN 现在按每个 observable context 内的 rank_within_group 取 topN，
不再把所有 context 混在一起全局排序后取 topN。
```

当前 E10gb1 vs E10gH1 的 topN 观察：

```text
1. event density 在 top1/top3/top5/top10/top20 中仍然基本不出现。

2. E10gb1:
   P10 deconv_efun top3/top5/top10 强烈偏 dimred density。
   P8 efun top1/top3/top5 的 raw slow-mode fraction 很高，
   但随着 topN 增大到 top20，会被更多 fast/secondary hits 稀释。

3. E10gH1:
   P10 deconv_efun top1/top3/top5/top10 更偏 raw density。
   deconv_efun 的 raw hits 仍然几乎全是 fast mode。
   efun 也有少量 slow-mode hit，但不如 E10gb1 稳定。

4. 因此目前建议第一判断看 top3/top5/top10。
   top1 太受单个 mode 影响；
   top20 可以作为鲁棒性检查；
   top50 仅作为探索性附录，不作为主结论。
```
