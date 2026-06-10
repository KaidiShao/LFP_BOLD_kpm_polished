# 当前 cross-dataset consistency rationale 与证据记录

日期：2026-06-04

这份文档保存当前叙事形成过程、历史叙事降级原因、具体数值和图路径。
真正的新 dataset 执行入口见短版 SOP：

```text
docs/current_minimal_cross_dataset_consistency_sop_2026-06-04_zh.md
```

更早的旧文档仍然保留为历史记录：

```text
docs/new_dataset_cross_dataset_consistency_sop_2026-05-28.md
```

## 1. 当前 working hypothesis

当前主叙事先收窄到一条最稳的线：

```text
standardized complex-split BLP Koopman density 中的
ripple/gamma-like, no-pure-theta fast subprocess

preferentially couples to

BOLD deconv_efun / deconvolved intrinsic BOLD perturbation branch.
```

中文表述：

```text
BOLD deconv_efun 更像是在捕捉与 ripple/gamma 相关的 fast intrinsic-trigger
或 perturbation-like BLP subprocess。
```

当前最保守可写的生物学解释是：

```text
ripple/gamma-triggered network reorganization
```

可以作为工作假设联系到：

```text
memory consolidation / replay-associated intrinsic state perturbation
```

但现阶段不能写成因果结论，只能写成 cross-dataset consistency 支持的 working hypothesis。

## 2. 当前最小主线口径

第一阶段只看最稳的条件，不再把所有条件混在一起：

```text
BLP side:
    complex_split_projected_vlambda_standardize

BLP activity:
    adaptive_envelope 为主
    abs 作为对照 / QC

P5 label:
    只用 P2 theta/gamma/ripple band-event strict label
    不用 dominant label 作为主结论

BOLD side:
    第一阶段只看 roi_mean
    global_svd100 / gsvd100_ds / HP_svd100 暂时作为第二阶段展开

P8/P10 xcorr:
    分开看 efun 和 deconv_efun
    当前主检验优先看 deconv_efun

ROI:
    ROI summary 是 spatial sanity check
    不是主证据
```

第一阶段不要再问“所有方法、所有 observable、所有 density 谁最好”。当前只问：

```text
在更多 dataset 里，roi_mean deconv_efun 的 top xcorr 是否不仅稳定偏向
ripple_gamma_no_theta / gamma_selective / ripple_selective 这类 fast subprocess，
而且比 roi_mean efun 更偏向这些 RG-like fast subprocess？
```

## 3. 需要验证的核心问题

### 3.1 P5 层面

先证明每个 dataset 的 BLP dimred space 里能不能拆出可解释 subprocess。

核心问题：

```text
这个 dataset 是否有至少一个 method-k 同时包含：

1. theta_selective component
2. ripple_gamma_no_theta / gamma_selective / ripple_selective component
```

解释：

```text
P5 pass = 这个 dataset 的 BLP dimred space 可以支持 two-subprocess 解释。
P5 fail = 不应该继续把这个 dataset 的 P8/P10 解释成 subprocess coupling。
```

当前主叙事最关心的是第二类 fast component 是否稳定存在：

```text
ripple_gamma_no_theta / gamma_selective / ripple_selective
```

theta_selective 仍然记录，但现在不强行把 BOLD efun 解释成 theta。

### 3.2 P8/P10 coupling 层面

P5 pass 后再看 P8/P10。

核心问题：

```text
roi_mean deconv_efun top xcorr hits 是否比 roi_mean efun 更偏向 fast subprocess label？
```

需要分开：

```text
efun
deconv_efun
```

当前预期：

```text
deconv_efun:
    相比 efun，更偏 ripple_gamma_no_theta / gamma_selective / ripple_selective
    更支持 fast intrinsic-trigger 叙事

efun:
    可能更像 slow/session-state branch
    但不应该直接等同于 theta_selective
```

### 3.3 raw efun density 层面

raw BLP efun density 仍然重要，因为它可以揭示 time scale。

要看：

```text
和 BOLD efun top xcorr 的 raw BLP efun timescale
和 BOLD deconv_efun top xcorr 的 raw BLP efun timescale
是否不同？
```

当前预期：

```text
BOLD efun:
    更可能对应 slow/session-wise state change

BOLD deconv_efun:
    更可能对应 fast perturbation / intrinsic trigger
```

这个结论需要 whole trace 支持，不能只靠 zoom window。

### 3.4 ROI 层面

ROI 现在只回答 spatial sanity：

```text
deconv_efun top xcorr selected BOLD modes 是否有跨 dataset 类似的 ROI footprint？
```

ROI 不再作为证明 theta/RG 空间分离的主证据。

当前 ROI 图应该先看：

```text
roi_mean × RG-no-theta × deconv_efun
all ROI by dataset
dataset-pair ROI profile consistency
subprocess ROI confusion matrix
```

all-ROI-by-dataset 图默认每个 dataset 单独 color scale，因为这张图第一目的
是看每个 dataset 内部 ROI pattern，而不是比较绝对幅度。

ROI confusion matrix 的目的不是预设证明 RG-no-theta 一定特殊，而是检验：

```text
RG-no-theta selected BOLD modes 是否在 ROI footprint 上形成一个可重复、
并且区别于 theta/mixed/inactive 的 spatial group？
```

正式判据：

```text
支持 RG-no-theta spatial special:
    within RG-no-theta cross-dataset ROI correlation 高
    RG-no-theta vs theta/mixed/inactive cross-dataset ROI correlation 低
    same-dataset 的跨 subprocess 相关不能单独作为支持证据

不支持 RG-no-theta spatial special:
    within RG-no-theta 和 RG-vs-other 差不多
    或者所有 subprocess 都落在同一片 high-consistency BOLD ROI subspace
```

因此 ROI confusion matrix 是 spatial sanity / falsification 图：
它可以支持“RG-no-theta 有可重复 ROI footprint”，也可以排除
“RG-no-theta 在 ROI 层面明显特殊”这个说法。

## 4. 之前叙事的当前状态

### 4.1 Event density alone

旧叙事：

```text
P2 theta/gamma/ripple event density 本身可能解释 BOLD coupling。
```

当前状态：

```text
降级为 control / baseline。
```

理由：

```text
P8/P10 top xcorr 更常落在 raw BLP efun density 或 dimred BLP efun density，
而不是 event density。
```

### 4.2 Pure ripple component

旧叙事：

```text
寻找 pure ripple selective component。
```

当前状态：

```text
被 ripple-gamma-no-theta fast subprocess 取代。
```

理由：

```text
ripple 和 gamma 经常共同出现。
sharp-wave-ripple 也可能含有 theta/gamma 成分。
要求 pure ripple 太严格，反而不符合数据结构。
```

当前保留：

```text
ripple_gamma_no_theta
gamma_selective
ripple_selective
```

### 4.3 完整 theta-slow + RG-fast two-branch model

旧叙事：

```text
theta-like BLP subprocess -> BOLD efun
ripple-gamma-like BLP subprocess -> BOLD deconv_efun
```

当前状态：

```text
部分保留。
```

保留：

```text
ripple-gamma-like BLP subprocess -> BOLD deconv_efun
```

降级：

```text
theta-like BLP subprocess -> BOLD efun
```

理由：

```text
deconv_efun 与 RG-like label 的关系更稳定。
efun 的 dimred label composition 更 mixed，不能直接写成 theta-specific。
```

### 4.4 ROI spatial separation

旧叙事：

```text
theta-linked 和 RG-linked BOLD modes 应该有不同 ROI footprint。
```

当前状态：

```text
降级为 spatial sanity check。
```

理由：

```text
BOLD intrinsic ROI space 本身就有强 shared subspace。
很多 selected modes 可能只是落在同一片高一致性 BOLD subspace。
```

当前 ROI 可支持：

```text
selected BOLD modes 有跨 dataset spatial repeatability。
```

当前 ROI 不能单独支持：

```text
theta/RG spatially separated。
```

### 4.5 历史叙事降级的证据链

这一节的目的不是说历史叙事“完全错”，而是说明它们不能再作为第一阶段主结论。
每个被降级的叙事都要对应明确的分析和图，避免后续又回到“什么都可能”的状态。

#### A. event density alone 不够

被检验的问题：

```text
P2 event density 本身是不是 BOLD coupling 的主要来源？
```

用过/应该继续保留的图：

```text
P8/P10 density-source competition 图：
    event_density vs raw_efun_density vs dimred_efun_density

P8/P10 topN density-class composition 图：
    看 top hits 里 event_density 占比是否真的高

whole-trace examples：
    看 event density、raw/dimred efun density 和 BOLD feature 的实际时间序列
```

降级理由：

```text
当前多数据结果里，strong top hits 更常来自 raw BLP efun density 或
dimred BLP efun density。P2 event density 仍然是必要 control，但不能单独解释
BOLD efun/deconv_efun coupling。
```

#### B. pure ripple component 过窄

被检验的问题：

```text
是否应该只寻找 ripple-only component？
```

用过/应该继续保留的图：

```text
P5 P2-band strict label grid by method-k
P5 strict theta-vs-ripple/gamma effect scatter
P5 two-subprocess candidate map
P5 selected component peri-event mean by band
P5 paired top-window diagnostic plate
```

降级理由：

```text
真实 P2 band event 里 gamma 和 ripple 经常共同出现。
strict pure-ripple component 不是稳定主结构。

更稳的 fast subprocess 是：
    ripple_gamma_no_theta / gamma_selective / ripple_selective
```

#### C. theta-slow -> BOLD efun 目前不够稳

被检验的问题：

```text
BOLD efun top hits 是否稳定偏 theta_selective dimred BLP component？
```

用过/应该继续保留的图：

```text
P8/P10 efun vs deconv_efun strict-label composition by topN
P8/P10 deconv-vs-efun RG enrichment by topN
P8/P10 density-source competition, split by efun/deconv_efun
raw efun timescale ECDF / histogram, split by efun/deconv_efun
whole-trace examples, split by efun/deconv_efun
```

降级理由：

```text
deconv_efun 与 RG-like fast label 的关系更稳定。
更严格地说，deconv_efun 的 RG-like fraction 需要高于 efun；
否则只能说 deconv 有 RG hits，不能说 deconv 比 efun 更像 fast process。
efun 的 label composition 更 mixed；它可能仍然反映 slow/session state，
但不能直接写成 theta-specific branch。
```

当前只保留：

```text
efun may reflect slow/session-state dynamics.
```

暂时不写：

```text
efun = theta subprocess.
```

#### D. ROI spatial separation 目前不能证明 two-subprocess 空间分离

被检验的问题：

```text
theta-selected BOLD modes 和 RG-selected BOLD modes 是否有不同 ROI footprint？
```

用过/应该继续保留的图：

```text
P7 intrinsic BOLD efun ROI confusion matrix
P8 selected BOLD mode index overlay on P7 ROI matrix
P8 subprocess ROI confusion / theta-vs-RG ROI separation
P8 roi_mean RG-no-theta all-ROI-by-dataset
P8 roi_mean RG-no-theta dataset-pair ROI consistency
```

降级理由：

```text
BOLD ROI space 自身存在高一致性 shared subspace。
selected modes 可能只是从这个 shared BOLD subspace 里挑出了和 BLP coupling
较强的 mode，而不是形成两个空间上明显分离的 ROI footprint。
```

当前 ROI 问题改成两个更稳的 sanity check：

```text
1. deconv-selected BOLD modes 是否有可重复 ROI footprint？
2. efun-selected 和 deconv-selected BOLD modes 的 ROI profile 是否明显不同？
```

## 5. 每个新 dataset 要画的最小图

### 5.1 Dataset QC

用于判断 dataset 是否可进入主分析：

```text
P2 theta top windows
P2 gamma top windows
P2 ripple top windows
P5 paired subprocess diagnostic plates
P5 paired top-window sheets
```

判断标签：

```text
usable
suspicious_but_usable
exclude_from_cross_dataset_conclusion
```

### 5.2 P5 subprocess gate

必看图：

```text
00_strict_label_counts_abs_vs_adaptive_envelope.png
04_adaptive_envelope_strict_theta_vs_ripple_gamma_effect_scatter.png
02_adaptive_envelope_strict_label_grid_by_method_k.png
03_adaptive_envelope_strict_label_composition_by_method_k.png
05_adaptive_envelope_strict_two_subprocess_candidate_map.png
```

解释图：

```text
method-k diagnostic plate
paired top-window sheet
```

通过标准：

```text
至少一个 method-k 有 theta component + fast RG component。
多个 method-k pass 更可靠。
```

### 5.3 P8/P10 coupling gate

第一阶段最小主图：

```text
roi_mean × efun/deconv_efun × topN strict-label composition
roi_mean × deconv-vs-efun RG enrichment by topN
roi_mean × efun/deconv_efun × density-source competition
roi_mean × efun/deconv_efun whole-trace examples
raw efun timescale distribution by topN
```

当前最关心：

```text
roi_mean deconv_efun 是否比 roi_mean efun 更稳定偏 RG-like fast subprocess。
主指标：
    deconv RG-like fraction - efun RG-like fraction
    matched units 中 deconv RG-like fraction > efun RG-like fraction 的比例

判读规则：
    如果 pooled fraction 支持 deconv > efun，但 matched-unit win rate 不高，
    只能写成 pooled enrichment / aggregate trend。
    只有当 matched-unit win rate 也稳定高于 0.5，才写成 cross-context stable preference。
```

2026-06-04 当前已得到的 working result：

```text
输入/输出：
    results/p8_p10_roi_mean_deconv_vs_efun_rg_enrichment_20260604/
    E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
        p8_p10_roi_mean_deconv_vs_efun_rg_enrichment_20260604\

核心观察：
    roi_mean deconv_efun 的 top xcorr hits 在 pooled 统计上
    比 roi_mean efun 更偏 RG-like / ripple_gamma_no_theta fast subprocess。

关键数值：
    P8 top10:
        efun RG-like fraction   ≈ 0.414
        deconv RG-like fraction ≈ 0.610
        delta                   ≈ +0.196

    P10 top10:
        efun RG-like fraction   ≈ 0.526
        deconv RG-like fraction ≈ 0.712
        delta                   ≈ +0.186

限制：
    matched-unit win rate 在 top5/top10 不是很强。
    因此当前正式表述应为：
        roi_mean deconv_efun shows pooled RG-like enrichment relative to efun。

    暂时不要写成：
        every matched context stably prefers deconv_efun over efun。
```

### 5.4 ROI sanity check

第一阶段最小主图：

```text
P8 roi_mean RG-no-theta all ROI by dataset
P8 roi_mean RG-no-theta dataset-pair consistency
P8 roi_mean RG-no-theta efun-vs-deconv ROI consistency
P8 roi_mean raw_csplit_q070 efun-vs-deconv ROI consistency
P8 roi_mean subprocess ROI confusion matrix
```

画图规范：

```text
all ROI by dataset:
    ROI 顺序保持 P7 anatomical order
    每个 dataset 单独 color scale

dataset-pair consistency:
    dataset pair 放横轴
    method-k 放纵轴

efun-vs-deconv ROI consistency:
    same dataset efun vs deconv ROI correlation
    cross dataset efun-within ROI correlation
    cross dataset deconv-within ROI correlation
    cross dataset efun-vs-deconv ROI correlation
    同时画 selected BOLD mode index overlap

subprocess ROI confusion matrix:
    rows/columns = strict subprocess label
        theta
        ripple_gamma_no_theta
        mixed_or_partial
        inactive
    cell value = selected BOLD ROI profile correlation
    必须分开画 efun 和 deconv_efun
    必须分开画 observable，第一阶段主看 roi_mean
    主看 cross-dataset matrix；same-dataset matrix 只能作为 QC
    主判据是 RG-within 是否高于 RG-vs-other
```

2026-06-05 当前可用 P8 ROI profile 检查：

```text
输入：
    results/pipeline_roi_profile_consistency_current/p8_roi_profiles_by_subprocess_long.csv

输出：
    results/p8_roi_rg_special_confusion_20260605/
    E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
        p8_roi_rg_special_confusion_20260605\

数据范围：
    e10gb1, e10fV1, e10gh1, e10gw1, f12m01
    K13 standardized csplit 数据暂时还没有纳入这个 P8 ROI subprocess profile CSV。

roi_mean × deconv_efun 的当前结果：
    RG-within cross-dataset ROI corr mean ≈ 0.56
    theta-within cross-dataset ROI corr mean ≈ 0.61
    theta-vs-RG cross-dataset ROI corr mean ≈ 0.59
    same-dataset theta-vs-RG ROI corr mean ≈ 0.89

解释：
    当前 ROI confusion matrix 不能证明 RG-no-theta 有独立 ROI footprint。
    更像是 selected BOLD modes 落在 shared/high-consistency BOLD ROI subspace。
    因此 ROI 层面目前只能作为 spatial sanity check，
    不能作为 RG-no-theta 特殊性的主证据。
```

## 6. Cross-dataset 主判断

最终不要堆很多图，只回答五个问题：

```text
1. 哪些 dataset 通过 P5 subprocess gate？
2. 哪些 method-k 在多个 dataset 中稳定 pass？
3. roi_mean deconv_efun top hits 是否比 roi_mean efun 更稳定偏 RG-like fast subprocess？
4. 这些 deconv-selected BOLD ROI profiles 是否有跨 dataset spatial repeatability？
5. efun-selected 和 deconv-selected BOLD ROI profiles 是否真的不同？
```

当前主结论模板：

```text
In standardized complex-split BLP Koopman space, most usable datasets contain
fast ripple/gamma-like no-pure-theta dimred components.  In P8/P10 coupling,
BOLD deconv_efun, especially under roi_mean readout, preferentially selects
these fast subprocess densities more strongly than BOLD efun does.  This supports a working hypothesis that
ripple/gamma-associated BLP subprocesses may act as intrinsic triggers for BOLD
network reorganization.
```

中文模板：

```text
在 standardized complex-split BLP Koopman space 中，多数可用 dataset 能找到
ripple/gamma-like、非 pure-theta 的 fast dimred component。P8/P10 中，BOLD
deconv_efun，尤其是 roi_mean readout 下，相比 efun 更倾向于和这类 fast subprocess
density 发生强 xcorr。因此当前数据支持一个工作假设：ripple/gamma-associated
BLP subprocess 可能对应 intrinsic trigger，并参与 BOLD network reorganization。
```

## 7. Minimal pipeline

这个 minimal pipeline 不重新定义 P1-P10。它只在当前主线结果已经存在时，
生成当前叙事需要的最小分析。

### 7.1 单个 dataset 的 P5 gate

WSL/conda Python：

```bash
cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished
/home/kdshao/anaconda3/bin/python scripts/run_pipeline5_band_selectivity.py \
  --dataset <dataset>
```

快速检查输入是否齐：

```bash
/home/kdshao/anaconda3/bin/python scripts/run_pipeline5_band_selectivity.py \
  --dataset <dataset> \
  --dry-run
```

### 7.2 多 dataset 的最小 consistency run

使用：

```bash
/home/kdshao/anaconda3/bin/python scripts/run_minimal_new_dataset_consistency.py \
  --datasets <dataset1> <dataset2> ... \
  --tag <analysis_tag>
```

这个 runner 会按当前 SOP 执行：

```text
1. P5 band-selectivity gate
2. P8/P10 strict-band coupling annotation
3. P8/P10 trace examples
4. roi_mean RG-no-theta ROI summary, if ROI profile CSV is available
5. P8 roi_mean efun-vs-deconv ROI consistency, if ROI profile CSV is available
```

输出默认写到：

```text
results/<analysis_tag>/
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\<analysis_tag>\
```

### 7.3 缺数据时的规则

如果某个步骤缺输入：

```text
P5 缺：
    先补 P5，不能解释 P8/P10 subprocess coupling。

P8/P10 缺：
    记录为 coupling unavailable，不用 summary_figures cache 代替。

P7 ROI profile 缺：
    ROI sanity check 跳过或补 export，不影响 P5/P8/P10 主 coupling 结论。

P4 缺：
    不能由 minimal pipeline 自动补，必须回到 P4 主计算。
```

## 8. 当前最小 pipeline 的成功标准

一个新 dataset 可以纳入主结论，需要满足：

```text
1. Dataset QC 不是 exclude。
2. P5 至少一个 method-k 有 theta + fast RG pair。
3. P8/P10 roi_mean deconv_efun top hits 中，RG-like fast label 明显出现，
   并且 pooled RG-like fraction 高于 roi_mean efun。
4. 如果要写成“稳定更偏 RG”，matched-unit deconv>efun win rate 也应高于 0.5。
5. raw efun / trace 结果不与 fast-trigger 解释明显矛盾。
```

如果只满足 1 和 3，但 P5 theta pair 不稳定：

```text
可以纳入 RG/deconv fast-branch 分析；
不要用于支持完整 theta-slow + RG-fast two-branch model。
```

如果 ROI 不稳定：

```text
不推翻 RG/deconv coupling；
只说明 spatial footprint 证据不足。
```
