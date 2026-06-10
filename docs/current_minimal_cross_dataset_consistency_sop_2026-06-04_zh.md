# 当前 cross-dataset consistency 最小 SOP

> 2026-06-10 更新：本文件已降级为旧版 provenance。
> 当前执行版请看：
> `docs/current_minimal_cross_dataset_consistency_sop_2026-06-10_zh.md`

本文件是执行版。历史叙事、被降级的旧想法、完整证据链和数值记录见：

```text
docs/current_minimal_cross_dataset_consistency_rationale_2026-06-04_zh.md
```

## 1. 当前 working hypothesis

第一阶段只验证一个最小、可重复的问题：

```text
在 standardized complex-split BLP 条件下，
roi_mean deconv_efun 的 top xcorr hits
是否比 roi_mean efun 更偏向 RG-like fast subprocess？
```

这里的 RG-like fast subprocess 包括：

```text
ripple_gamma_no_theta
gamma_selective
ripple_selective
```

当前可写成的叙事：

```text
roi_mean deconv_efun shows pooled RG-like enrichment relative to roi_mean efun.
```

暂时不要写成：

```text
every matched context stably prefers deconv_efun over efun
RG-no-theta has an independent ROI network
theta-slow and RG-fast are already spatially separated
```

## 2. 主线输入

第一阶段固定这些条件，不把所有历史版本混在一起：

```text
BLP observable:
    complex_split_projected_vlambda_standardize

BLP activity transform:
    adaptive_envelope 为主
    abs 作为 QC / sensitivity check

P5 component label:
    只用 P2 theta/gamma/ripple band-event strict label
    不用 dominant label 作为主结论

BOLD observable:
    roi_mean

BOLD feature:
    efun
    deconv_efun

P8/P10 density source:
    raw_csplit_q070
    dimred csplit method-k density
```

## 3. 每个新 dataset 的四个 gate

### Gate 1: Dataset QC

目的：确认这个 dataset 能进入主分析。

最低检查：

```text
P2 theta/gamma/ripple event detection 存在
P5 standardized complex-split reduction/density 存在
P8/P10 roi_mean xcorr 存在
P7 roi_mean ROI profile 存在；如果缺，只跳过 ROI sanity，不影响 coupling 主结论
```

可疑数据，例如 clipping/saturation 明显的数据，先标记为 QC concern。

### Gate 2: P5 subprocess gate

目的：确认 BLP dimred component space 里是否能拆出两个可解释 subprocess。

通过标准：

```text
至少一个 method-k 同时有：
    theta_selective component
    RG-like component:
        ripple_gamma_no_theta / gamma_selective / ripple_selective
```

主图：

```text
strict label counts: abs vs adaptive_envelope
strict theta-vs-ripple/gamma effect scatter
strict label grid by method-k
two-subprocess candidate map
diagnostic plate for selected theta/RG pair
```

解释：

```text
P5 pass 只说明 BLP 侧可以拆出 theta-like 与 RG-like subprocess。
P5 pass 本身还不能说明 BOLD 已经耦合到 RG-like subprocess。
```

P5 density independence check：

在 P5 pass 以后，还应该检查 theta density 与 RG density 本身是否高度相关。
这个检查回答的是：

```text
theta_selective dimred density
和
ripple_gamma_no_theta dimred density

是否只是同一个 density process 的不同标签？
```

当前实现脚本：

```text
scripts/analyze_p5_theta_rg_density_correlation.py
```

输入为 P8/P10 实际使用的 thresholded dimred density：

```text
E:\DataPons_processed\<dataset>\
    pipeline5_dimred_thresholded_density\
        complex_split_projected_vlambda_standardize_rmsenv_adaptive\
            <method>_k<kk>\mat\*.mat
```

主指标：

```text
canonical theta/RG pair:
    每个 method-k 里 theta_effect_z 最强的 theta_selective component
    与 RG effect 最强的 ripple_gamma_no_theta component

density_corr:
    直接对两个 thresholded density traces 算 Pearson correlation

session_demeaned_density_corr:
    先在每个 session 内去均值，再算 correlation
    用来降低 shared session-wise state offset 的影响
```

当前结果目录：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    p5_theta_rg_density_correlation/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        p5_theta_rg_density_correlation\
```

当前整体结果：

```text
canonical pairs:
    n = 354
    median full density corr              = -0.011
    median session-demeaned density corr  = -0.033
    |session-demeaned corr| < 0.2 fraction = 0.81
    |session-demeaned corr| < 0.5 fraction = 1.00

all theta/RG pairs:
    n = 1339
    median full density corr              = -0.007
    median session-demeaned density corr  = -0.033
    |session-demeaned corr| < 0.2 fraction = 0.87
    |session-demeaned corr| < 0.5 fraction = 1.00
```

当前解释：

```text
在 E10 和多数 F12 数据中，theta density 与 RG-no-theta density 大多低相关，
甚至略负相关。

这支持 P5 层面的 subprocess split：
    theta-like density 与 RG-like density 不是同一个 density trace 的重命名。

例外：
    f12m05、k13m18、k13m23 的 canonical pair correlation 偏高一些；
    k13m17/k13m19/k13m20/k13m21 因 theta strict component 不稳定，
    不能纳入这个 theta/RG pair correlation 主结论。
```

### Gate 3: P8/P10 coupling gate

目的：回答主问题：

```text
roi_mean deconv_efun 是否比 roi_mean efun 更偏 RG-like fast subprocess？
```

主图：

```text
roi_mean x efun/deconv_efun x topN strict-label composition
roi_mean deconv-vs-efun RG enrichment by topN
roi_mean x efun/deconv_efun density-source competition
roi_mean x efun/deconv_efun whole-trace examples
raw efun timescale distribution by topN
```

主指标：

```text
pooled RG-like fraction:
    topN hits 中 RG-like label 的比例

delta:
    deconv RG-like fraction - efun RG-like fraction

matched-unit win rate:
    在 matched context 中，
    deconv RG-like fraction > efun RG-like fraction 的比例
```

判读规则：

```text
pooled delta > 0, matched-unit win rate 不强:
    只能写 pooled enrichment / aggregate trend

pooled delta > 0, matched-unit win rate 稳定 > 0.5:
    可以写 cross-context stable deconv preference
```

当前 working result：

```text
结果目录：
    results/p8_p10_roi_mean_deconv_vs_efun_rg_enrichment_20260604/

P8 top10:
    efun RG-like fraction   ~= 0.414
    deconv RG-like fraction ~= 0.610
    delta                   ~= +0.196

P10 top10:
    efun RG-like fraction   ~= 0.526
    deconv RG-like fraction ~= 0.712
    delta                   ~= +0.186

限制：
    matched-unit win rate 在 top5/top10 不够强。
```

### Gate 4: ROI spatial sanity

目的：只做 spatial sanity，不作为主证据。

ROI 现在回答两个问题：

```text
1. deconv-selected BOLD modes 是否有可重复 ROI footprint？
2. efun-selected 和 deconv-selected BOLD ROI profiles 是否明显不同？
```

主图：

```text
roi_mean RG-no-theta all ROI by dataset
roi_mean RG-no-theta dataset-pair ROI consistency
roi_mean RG-no-theta efun-vs-deconv ROI consistency
roi_mean subprocess ROI confusion matrix
```

ROI 图规范：

```text
all ROI by dataset:
    ROI 顺序保持 P7 anatomical order
    每个 dataset 单独 color scale

dataset-pair consistency:
    dataset pair 放横轴
    method-k 放纵轴

subprocess ROI confusion matrix:
    必须分开 efun 和 deconv_efun
    第一阶段主看 roi_mean
    主看 cross-dataset matrix
    same-dataset matrix 只作为 QC
```

RG-no-theta ROI 特殊性的判据：

```text
支持 spatial special:
    RG-within cross-dataset ROI corr 高
    RG-vs-theta/mixed/inactive cross-dataset ROI corr 低

不支持 spatial special:
    RG-within 和 RG-vs-other 差不多
    或所有 subprocess 都落在 shared high-consistency BOLD ROI subspace
```

当前 ROI result：

```text
结果目录：
    results/p8_roi_rg_special_confusion_20260605/

当前数据范围：
    e10gb1, e10fV1, e10gh1, e10gw1, f12m01
    K13 standardized csplit 暂时还没纳入 P8 ROI subprocess profile CSV

roi_mean x deconv_efun:
    RG-within cross-dataset ROI corr mean      ~= 0.56
    theta-within cross-dataset ROI corr mean   ~= 0.61
    theta-vs-RG cross-dataset ROI corr mean    ~= 0.59
    same-dataset theta-vs-RG ROI corr mean     ~= 0.89

当前解释：
    ROI confusion matrix 不能证明 RG-no-theta 有独立 ROI footprint。
    ROI 目前只支持 spatial sanity，不支持 RG-no-theta spatial special。
```

2026-06-09 target footprint 更新：

用户当前更关心的 ROI 问题不是泛泛地问 RG-no-theta 是否特殊，而是问三类
selected BOLD ROI footprint 是否能支持 slow/transient 双过程叙事：

```text
1. efun x theta_selective
       是否像 emergent slow component / slow state 的 spatial footprint？

2. efun x ripple_gamma_no_theta
       是否和 efun x theta_selective 对应不同 ROI？

3. deconv_efun x ripple_gamma_no_theta
       是否像 transient intrinsic perturbation 的 spatial footprint？
```

正式检验脚本：

```text
scripts/plot_roi_mean_sop_target_footprints.py
```

输入：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_rg_no_theta_formal_roi_tests/
        roi_mean_rg_no_theta_selected_roi_profiles_long__top5.csv
        roi_mean_rg_no_theta_selected_roi_profiles_long__top10.csv
        roi_mean_rg_no_theta_selected_roi_profiles_long__top20.csv
```

输出：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_rg_no_theta_formal_roi_tests/
        target_footprint_tests/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        roi_mean_target_footprint_tests\
```

主图：

```text
01_target_roi_footprints_by_dataset__p8__top10.png
01_target_roi_footprints_by_dataset__p10__top10.png
02_target_roi_footprint_contrasts__p8__top10.png
02_target_roi_footprint_contrasts__p10__top10.png
03_target_roi_similarity_summary_by_topN.png
```

当前 top10 结果：

```text
baseline all-mode cross-dataset median ROI corr:
    0.099

baseline same sorted-position cross-dataset median ROI corr:
    0.175

P8:
    efun theta within cross-dataset median      = 0.211
    efun RG within cross-dataset median         = 0.373
    deconv RG within cross-dataset median       = 0.493
    efun theta vs efun RG same-dataset median   = 0.931
    efun theta vs deconv RG same-dataset median = 0.871
    efun RG vs deconv RG same-dataset median    = 0.919

P10:
    efun theta within cross-dataset median      = 0.159
    efun RG within cross-dataset median         = 0.167
    deconv RG within cross-dataset median       = 0.261
    efun theta vs efun RG same-dataset median   = 0.985
    efun theta vs deconv RG same-dataset median = 0.667
    efun RG vs deconv RG same-dataset median    = 0.726
```

当前解释：

```text
deconv RG footprint 的 cross-dataset ROI consistency 高于 all-mode baseline，
尤其 P8 明显。

但是 efun x theta 与 efun x RG 在同一 dataset 内 ROI profile 极其相似：
P8 median ~= 0.93，P10 median ~= 0.99。

因此 ROI 层面目前不能支持：
    efun x theta slow footprint
    和
    efun x RG fast/RG footprint
    是两个空间上清楚分离的 network。

更稳妥的写法是：
    slow/transient 区分主要由 xcorr feature type、subprocess label、
    density class 和 timescale/trace 证据支持；
    ROI 只能说明 selected BOLD modes 有一定可重复 spatial footprint，
    不能作为 theta/RG 空间分离的主证据。
```

2026-06-09 signed ROI 更新：

上面的 target footprint 结果使用的是 `mean_abs` magnitude 口径。
这个口径会丢掉 sign / polarity，因此只能说明：

```text
magnitude footprint 没有清楚分开。
```

不能说明：

```text
signed spatial footprint 也没有分开。
```

因此新增四种 ROI value mode：

```text
mean_abs       # magnitude baseline
real_mean      # signed real part
positive_real  # positive real part only
negative_real  # negative real part magnitude
```

P7 direct exporter 已支持：

```text
scripts/export_p7_roi_mean_profiles_direct_from_bold_post.py
    --roi-value-mode mean_abs / real_mean / positive_real / negative_real
```

当前 focused test 脚本：

```text
scripts/analyze_roi_mean_signed_target_modes.py
```

这个脚本只回答两个问题：

```text
1. selected BOLD mode index overlap:
       efun x theta_selective
       vs
       deconv_efun x ripple_gamma_no_theta

2. signed ROI / top ROI set:
       mean_abs / real_mean / positive_real / negative_real
       下 top ROI set 是否更分开？
```

当前结果目录：

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    roi_mean_signed_target_mode_tests/

E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
    standardized_csplit_k03_k16_all_current_20260607\
        roi_mean_signed_target_mode_tests\
```

当前 top10 结果：

```text
P8 selected P7 BOLD sorted-mode overlap:
    exact Jaccard median            = 0.211
    same/adjacent fraction median   = 0.714

P8 ROI profile correlation, efun-theta vs deconv-RG:
    mean_abs       median = 0.881
    real_mean      median = 0.642
    positive_real  median = 0.612
    negative_real  median = 0.754

P8 top10 ROI-set Jaccard:
    mean_abs       median = 0.538
    real_mean      median = 0.333
    positive_real  median = 0.429
    negative_real  median = 0.429

P10 ROI profile correlation, efun-theta vs deconv-RG:
    mean_abs       median = 0.632
    real_mean      median = 0.482
    positive_real  median = 0.512
    negative_real  median = 0.519

P10 top10 ROI-set Jaccard:
    mean_abs       median = 0.333
    real_mean      median = 0.250
    positive_real  median = 0.250
    negative_real  median = 0.250
```

当前解释：

```text
signed / positive / negative ROI 口径确实比 mean_abs 更能拉开
efun-theta 与 deconv-RG。

所以之前 magnitude ROI 不分开，不应该被解释成
theta/RG spatial footprint 完全不分。

但是 P8 中 selected BOLD sorted mode 经常落在相同或相邻 mode：
same/adjacent fraction median ~= 0.71。

因此更准确的说法是：
    efun-theta 与 deconv-RG 不是完全不同的 BOLD eigenspace；
    但它们在 signed ROI polarity / top ROI set 上比 magnitude profile 更可分。

P10 的 selected unit overlap 不作为主解释，
因为 efun 与 deconv_efun 属于不同 P9 reduction spaces，
component index 不能直接当作同一个 BOLD basis index 比较。
```

## 4. 最小命令入口

### 单 dataset P5 gate

```powershell
wsl.exe -d Ubuntu-22.04 bash -lc "cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished && /home/kdshao/anaconda3/bin/python scripts/analyze_e10gb1_p5_dimred_band_event_response_v2.py --dataset <dataset> --summary-only"
```

### 多 dataset minimal consistency

先 dry run：

```powershell
wsl.exe -d Ubuntu-22.04 bash -lc "cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished && /home/kdshao/anaconda3/bin/python scripts/run_minimal_new_dataset_consistency.py --datasets <dataset1> <dataset2> --tag <tag> --dry-run"
```

再正式运行：

```powershell
wsl.exe -d Ubuntu-22.04 bash -lc "cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished && /home/kdshao/anaconda3/bin/python scripts/run_minimal_new_dataset_consistency.py --datasets <dataset1> <dataset2> --tag <tag>"
```

结果默认放在：

```text
results/<tag>/
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\<tag>\
```

## 5. 缺数据时的规则

```text
P5 缺：
    先补 P5，否则不能解释 P8/P10 subprocess coupling。

P8/P10 缺：
    记录为 coupling unavailable。
    不用 summary_figures cache 代替原始结果。

P7 ROI profile 缺：
    ROI sanity 跳过或补 export。
    不影响 P5/P8/P10 coupling 主结论。

P4 缺：
    minimal pipeline 不自动补。
    必须回到 P4 主训练/主导出流程。
```

## 6. 最小成功标准

一个新 dataset 可以纳入当前主结论，需要满足：

```text
1. Dataset QC 没有严重 exclude 标记。
2. P5 至少一个 method-k 有 theta + RG-like pair。
3. P8/P10 roi_mean deconv_efun top hits 中 RG-like label 明显出现。
4. pooled RG-like fraction:
       deconv_efun > efun
5. 如果要写成 stable deconv preference:
       matched-unit win rate 也应 > 0.5。
```

如果 ROI 不支持 spatial separation：

```text
不推翻 deconv/RG coupling。
只说明 spatial footprint 证据不足。
```
