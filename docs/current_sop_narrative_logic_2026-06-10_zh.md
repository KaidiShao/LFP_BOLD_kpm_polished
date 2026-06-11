# 当前 SOP 叙事逻辑

日期：2026-06-10

本文件记录当前 cross-dataset SOP 背后的叙事链。它不是执行脚本，也不是完整图清单；它回答的是：

```text
我们为什么按这个顺序看 P5、P8/P10、timescale 和 ROI？
每一层证据分别支持什么？
哪些结果会加强或削弱当前 working hypothesis？
```

相关文档：

```text
notes/koopman_lfp_bold_framework.md
notes/koopman_event_narrative.md
docs/current_minimal_cross_dataset_consistency_sop_2026-06-10_zh.md
docs/current_sop_figure_candidate_review_2026-06-10_zh.md
docs/roi_real_mean_mainline_rationale_2026-06-10_zh.md
```

---

## 1. 当前 working hypothesis

当前最小叙事是：

```text
standardized complex-split BLP dimred eigenfunction density
可以拆出两个可解释的 LFP/BLP subprocess：

1. theta-like slow/state subprocess
2. ripple/gamma-no-theta fast subprocess

roi_mean BOLD efun 更像 global slow state readout。
roi_mean BOLD deconv_efun 更像 fast intrinsic perturbation readout。

如果 deconv_efun 比 efun 更偏向 ripple/gamma-no-theta subprocess，
并且这个结果跨 dataset、topN、P8/P10 都相对稳定，
那么它支持：

fast hippocampal ripple/gamma-like subprocess
可能作为 intrinsic perturbation / trigger
影响 global BOLD pattern。
```

可以进一步作为一个工作解释，联系到：

```text
ripple-triggered system reorganization
或 memory-consolidation 相关的 fast event-triggered network perturbation
```

但目前不要写成因果结论。当前证据只能支持 coupling / readout / perturbation-compatible interpretation。

---

## 2. 为什么第一步先看 P5，而不是直接看 BOLD coupling

P8/P10 的 xcorr 搜索空间很大。如果直接问哪个 density source 和 BOLD efun/deconv_efun 最相关，很容易把几个不同问题混在一起：

```text
这个 LFP component 本身是否可解释？
它是不是只是普通 event density 的替代？
theta / ripple / gamma 是否真的能被 dimred space 区分？
BOLD efun 和 deconv_efun 的区别是否只是数值搜索造成的？
```

所以第一步必须在 LFP/BLP 侧先过 P5 subprocess gate。

P5 gate 的关键问题是：

```text
只看 standardized complex-split BLP dimred space，
是否能在同一个 dataset x method-k 中同时找到：

1. theta_selective component
2. ripple_gamma_no_theta / RG-like component
```

这里的 label 来自 P2 theta/gamma/ripple band-event windows，并且已经和 non-event baseline 比较过。因此它不是简单地看 event overlap，而是在问：

```text
这个 dimred component 在 theta / gamma / ripple event 附近
是否相对 baseline 有更强 activity？
```

如果 P5 过不了，后续 BOLD coupling 即使有相关，也很难解释成两个 subprocess。

---

## 3. P5 pass 之后，还需要看 theta/RG density correlation

P5 pair-pass 证明的是：

```text
同一个 method-k 中能找到 theta-like component 和 RG-like component。
```

但这还不够。还需要证明这两个 component 的 activity density 不是同一条 signal 的两个名字。

因此需要看：

```text
corr(theta-like density, RG-like density)
```

解释口径：

```text
如果 theta/RG density correlation 低或中等：
    支持它们是两个不同 subprocess。

如果 theta/RG density correlation 很高：
    说明 P5 label 虽然不同，但实际 time course 可能没有充分分离。
    这个 dataset 或 method-k 应该标成 weak / ambiguous。
```

这一层是在 LFP 内部证明 subprocess separability，不依赖 BOLD。

---

## 4. 为什么要比较 event density、raw efun density 和 dimred efun density

当前叙事不是：

```text
传统 event density 本身解释 BOLD。
```

而是：

```text
BLP eigenfunction activity，尤其是 dimred subprocess density，
比 simple event_density 更接近 BOLD dynamical readout。
```

所以 P8/P10 必须做 density-source competition：

```text
event_density
raw BLP efun density
dimred BLP efun density
```

支持当前叙事的结果应该是：

```text
raw/dimred BLP efun density > simple event_density
```

更强的结果是：

```text
dimred BLP efun density > raw BLP efun density > event_density
```

这个结果的含义是：

```text
P5 dimred subprocess representation 不只是把传统 event 重新包装了一遍，
而是提供了更接近 BOLD coupling target 的 latent representation。
```

如果 event_density 和 dimred density 差不多，甚至 event_density 更强，那么当前 P5 subprocess 叙事就会变弱。

---

## 5. BOLD efun 与 deconv_efun 的功能区分

当前最重要的 BOLD 口径先固定为：

```text
BOLD observable = roi_mean
```

解释：

```text
roi_mean 反映更 global 的 BOLD pattern。
BOLD efun 是这个 global pattern 的 Koopman latent state readout。
BOLD deconv_efun 更像从 BOLD efun 中提取出的 modal innovation / residual trigger。
```

因此主问题不是简单地问：

```text
哪个 LFP density 和 BOLD 最相关？
```

而是问：

```text
efun 和 deconv_efun 是否偏向不同的 LFP subprocess？
```

当前希望看到的模式是：

```text
BOLD efun:
    更偏 slow/state-like signal。

BOLD deconv_efun:
    更偏 ripple_gamma_no_theta / RG-like fast subprocess。
```

如果 deconv_efun 比 efun 更稳定地偏 RG-like density，那么可以支持：

```text
RG-like fast subprocess 更像 trigger / intrinsic perturbation，
而不是普通 slow state readout。
```

---

## 6. 如果 BOLD efun 不明显偏 theta，怎么办

theta-like subprocess 与 BOLD efun 的关系可能不会像 deconv-RG 那么干净。

如果：

```text
BOLD efun top hits 不明显偏 theta-like dimred density
```

这不一定推翻 slow-state 叙事。因为 BOLD efun 可能更接近：

```text
session-wise slow state
global drift / latent state
raw BLP slow eigenfunction density
```

而不一定正好落在 strict theta_selective dimred component 上。

因此需要另一条证据链：

```text
1. raw BLP efun density 与 BOLD efun/deconv_efun 的 xcorr
2. top raw BLP efun 的 timescale distribution
3. P8 selected BOLD mode eigenvalue / timescale
```

支持 slow-vs-fast readout 的结果应该是：

```text
efun selected BOLD modes:
    更慢
    mode index 更靠前
    eigenvalue 更接近单位圆

deconv_efun selected BOLD modes:
    更快
    mode index 更靠后
    eigenvalue 离单位圆更远
```

这条证据不要求 efun 必须严格偏 theta。它证明的是：

```text
efun branch 更像 slow/global state readout，
deconv branch 更像 fast intrinsic perturbation readout。
```

---

## 7. 为什么还要看 ROI，而且为什么用 signed ROI

时间相似性只能说明：

```text
某个 LFP density 和某个 BOLD efun/deconv_efun 在时间上相关。
```

但 current hypothesis 还需要空间解释：

```text
slow process 和 fast perturbation 是否对应不同的 BOLD spatial footprint？
```

之前的 mean_abs / magnitude ROI 口径只能回答：

```text
某个 ROI 是否参与这个 BOLD mode。
```

它会抹掉正负方向，因此不适合证明 spatial footprint 是否不同。当前 ROI 主线应该用：

```text
real_mean signed ROI profile
top positive ROI set
top negative ROI set
selected BOLD mode overlap
```

需要比较的四种 target 是：

```text
efun x theta-like density
efun x RG-like density
deconv_efun x theta-like density
deconv_efun x RG-like density
```

最关键的问题是：

```text
efun x theta-like
和
deconv_efun x RG-like

是否比其他组合更不同？
```

支持当前叙事的 ROI 结果应该是：

```text
efun-theta 与 deconv-RG 的 signed ROI profile correlation 更低；
top positive/negative ROI set overlap 更低；
selected BOLD mode overlap 更低或落在不同 mode group。
```

如果 ROI 仍然很像，那么要谨慎表述：

```text
ROI 层面目前支持的是 efun/deconv branch 的 temporal/timescale 区分，
但还不能证明 slow/RG 是两个空间分离的 BOLD network。
```

---

## 8. 为什么 P10 和 topN robustness 不能省

P8 是直接 BLP density 与 BOLD efun/deconv_efun 的 xcorr。
P10 进一步经过 BOLD-side reduction / reconstruction，因此可能受以下因素影响：

```text
BOLD dimred component number
BOLD ResKoopNet dictionary size
BOLD feature construction
topN cutoff
```

因此 P10 的结论不能只看 top10，也不能只看一个 component number。

当前至少需要：

```text
top3
top5
top10
top20
```

如果结论只在 top10 成立，但 top3/top5 不成立，说明它可能是较宽松的 pooled trend。

如果 top3/top5/top10/top20 都一致，说明结果更稳。

后续如果 P10 不稳定，需要把这些变量作为正式 sweep：

```text
BOLD dimred component number
BOLD ResKoopNet dictionary size
P10 topN cutoff
```

---

## 9. 支持、削弱和反驳当前叙事的结果

强支持：

```text
1. 多数 dataset 通过 P5 theta/RG pair gate。
2. theta/RG density correlation 低或中等。
3. dimred density > raw density > event density。
4. deconv_efun 比 efun 更偏 RG-no-theta。
5. efun selected BOLD modes 更慢、更接近单位圆；
   deconv selected modes 更快。
6. real_mean signed ROI 显示 efun-theta 与 deconv-RG
   有不同 spatial footprint 或不同 selected mode group。
7. P8 与 P10、top3/top5/top10/top20 趋势一致。
```

中等支持：

```text
1. P5 pair gate 成立，但 theta/RG density correlation 中等偏高。
2. deconv RG enrichment 只在 pooled/top10 成立，matched context 不稳定。
3. density-source competition 显示 raw/dimred > event，
   但 dimred 不明显强于 raw。
4. timescale 支持 efun/deconv 区分，但 ROI 不明显分开。
```

削弱：

```text
1. 多数 dataset 过不了 P5 pair gate。
2. theta/RG density 高相关，像同一 subprocess。
3. event_density 与 dimred density 一样强或更强。
4. deconv_efun 不比 efun 更偏 RG。
5. efun/deconv selected BOLD mode timescale 没区别。
6. real_mean ROI 和 selected mode overlap 完全不能区分 efun/deconv。
```

---

## 10. 当前最稳的写法

目前最稳妥的表述是：

```text
在通过 standardized complex-split P5 theta/RG pair gate 的数据中，
BLP dimred eigenfunction density 可以提供两个可解释的 subprocess：
theta-like slow/state subprocess 和 ripple/gamma-no-theta fast subprocess。

P8/P10 density-source competition 显示，
BLP efun-derived density，尤其是 dimred density，
比 simple event density 更能解释 roi_mean BOLD efun/deconv_efun coupling。

roi_mean BOLD deconv_efun 相比 efun 更偏向 RG-no-theta subprocess，
而 P8 selected BOLD mode eigenvalue/timescale 显示 efun branch 更慢、
deconv branch 更快。

这支持一个 working model：
BOLD efun 更像 global slow state readout，
BOLD deconv_efun 更像 fast intrinsic perturbation readout，
后者可能与 ripple/gamma-triggered network reorganization 有关。

ROI 层面需要使用 signed real_mean profile、top positive/negative ROI set
和 selected mode overlap 做 spatial sanity。
在 ROI 证据稳定前，不应把 theta/RG 写成已证明的两个空间分离网络。
```

