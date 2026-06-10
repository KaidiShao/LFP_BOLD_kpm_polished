# K13 跨数据叙事整合笔记

日期：2026-06-04

这份文档是中文整理版，目的不是再堆图，而是把目前 K13m18-K13m21 以及之前 E10gb1/E10gH1 相关分析里出现过的叙事线捋清楚。

当前最重要的变化是：

```text
现在最稳的主线不是完整的 “theta slow state + ripple-gamma fast trigger” 双边叙事。

现在最稳的是：

BOLD deconv_efun 和 RG-no-pure-theta / ripple-like BLP dimred components 的关系。
```

也就是说，当前可以先把故事收窄到：

```text
deconvolved BOLD eigenfunction branch preferentially couples to
ripple/gamma-like Koopman BLP density components.
```

中文理解：

```text
BOLD deconv_efun 更像是在捕捉和 ripple/gamma 相关的 fast intrinsic-trigger 或 perturbation-like BLP subprocess。
```

这个叙事已经挺有用，可以进一步联系到：

```text
ripple-trigger network reorganization
```

也可以作为工作假说联系到：

```text
memory consolidation / replay-associated intrinsic state perturbation
```

但目前还不能写成因果结论，只能写成 working hypothesis。

## 0. 历史叙事比较

这一节专门整理我们之前可能有过的叙事线。重点不是谁对谁错，而是把每条线现在应该放在什么位置说清楚。

### 叙事 A：simple event density 可以解释 BOLD coupling

最早的直觉：

```text
theta/gamma/ripple event 本身的 density，可能就能解释 BOLD efun/deconv efun。
```

后来看到的结果：

```text
P8/P10 top xcorr hits 大多数不是 event density，
而是 raw BLP efun density 或 dimred BLP efun density。
```

现在的状态：

```text
降级为 background / baseline。
```

怎么理解：

```text
P2 event density 仍然有用，因为它定义 theta/gamma/ripple event，
也用于 P5 band selectivity label。

但它不是当前最有解释力的 BOLD coupling source。
```

保留价值：

```text
event density 是对照组。
如果 Koopman density 比 event density 更强，
说明 P5/P8/P10 找到的不是简单 event occurrence，
而是 latent BLP activity / subprocess。
```

### 叙事 B：找 pure ripple selective component

早期想法：

```text
最好找到只在 ripple 上 active、不在 theta/gamma 上 active 的 component。
```

后来修正：

```text
ripple 和 gamma 经常共同出现。
sharp-wave-ripple 也可能含有 theta/gamma 成分。
因此要求 pure ripple 太严格。
```

现在的状态：

```text
被 RG-no-pure-theta 叙事替代。
```

当前更合理的 label：

```text
ripple_gamma_no_theta
ripple_selective
gamma_selective
```

核心思想：

```text
对 fast subprocess 来说，gamma coactivity 是可以接受的。
真正需要避免的是 pure-theta 主导。
```

所以当前主叙事不是 pure ripple，而是：

```text
RG-like / no-pure-theta fast subprocess
```

### 叙事 C：standardized complex-split 是更适合找 subprocess 的 BLP observable

早期观察：

```text
E10gb1 里 standardized complex-split 以后 fast components 更明显，
更容易看到 ripple/gamma-like component。
```

后来结果：

```text
standardized complex-split 确实成为当前 P5/P8/P10 新数据主线。
K13m18-K13m21 也都是按 standardized csplit 跑的。
```

现在的状态：

```text
保留为当前主线 condition。
```

注意：

```text
standardized csplit 不保证每个 dataset 都有 strict theta component。
但它目前仍然是最适合检查 RG/deconv 叙事的 BLP condition。
```

### 叙事 D：完整 two-subprocess 模型

之前最漂亮的故事：

```text
theta-like slow/state BLP subprocess
    -> BOLD efun

ripple-gamma-like fast/intrinsic-trigger BLP subprocess
    -> BOLD deconv_efun
```

这个故事的来源：

```text
E10gb1 里 P5 能看到 theta-like 和 ripple-gamma-like dimred components。
trace 图也显示 theta-like 更慢、更宽，RG-like 更尖、更快。
P8/P10 里 deconv_efun 对 RG-like density 的关系比较强。
```

K13 加入后的变化：

```text
RG/deconv 这半边保留下来了。
theta/efun 这半边在 K13m18-K13m21 里不稳定。
```

现在的状态：

```text
拆成两条：

1. RG/deconv:
   当前主线，可以继续讲。

2. theta/efun:
   暂时降级为待验证假说。
```

现在不能再强讲：

```text
所有 dataset 都有稳定 theta + RG pair。
```

但可以讲：

```text
在部分 dataset 中，P5 可以拆出 theta-like 和 RG-like components。
跨 K13 最稳定的是 RG-like component 与 deconv_efun 的关系。
```

### 叙事 E：BOLD efun 对应 slow/session-wise state，deconv_efun 对应 fast activity

之前 whole-trace 和 zoom trace 给出的直觉：

```text
BOLD efun 和 raw BLP density 的强 coupling 常常像 session-wise slow state change。
BOLD deconv_efun 和 raw BLP density 的强 coupling 更像 fast activity / perturbation。
```

现在的状态：

```text
保留，但需要分开处理。
```

deconv_efun 部分：

```text
和 RG-like dimred label 结果一致。
可以支持 fast intrinsic-trigger 叙事。
```

efun 部分：

```text
还不能直接等同于 theta。
需要 raw efun timescale distribution 和 whole-trace 证据继续支持。
```

所以现在应该说：

```text
deconv_efun = fast / perturbation-like branch，这条比较稳。

efun = slow/state-like branch，这条可能成立，
但 theta-specific 解释还不稳。
```

### 叙事 F：dimred component label 可以解释 P8/P10 top coupling

之前目标：

```text
给每个 P5 dimred component 打 label，
然后看 P8/P10 top xcorr hit 到底偏哪个 label。
```

现在结果：

```text
这个思路是对的，但不能把所有 label 混成一个大结论。
```

当前最有用的部分：

```text
P10 deconv_efun top dimred hits 明显偏 RG-no-theta / ripple-like。
```

当前不稳定的部分：

```text
P8 比 P10 弱。
BOLD efun 的 dimred label composition 比较 mixed。
theta label 在 K13m19-K13m21 不稳定。
```

现在的状态：

```text
保留为主分析方法。
但主读法收窄为：

deconv_efun 是否偏 RG-like label？
```

### 叙事 G：ROI 可以证明两个 subprocess 有不同空间 footprint

之前希望：

```text
theta-linked BOLD ROI profile 和 RG-linked BOLD ROI profile 不同。
```

后来看到：

```text
很多情况下，同一 BOLD observable 内的 ROI profiles 本来就有很强的 shared subspace。
特别是 gsvd100_ds，经常 efun/deconv 都落在同一个稳定 BOLD subspace。
```

现在的状态：

```text
降级为 spatial sanity check / QC。
```

ROI 现在能支持什么：

```text
某些 selected BOLD modes 是否有跨 dataset spatial repeatability。
```

ROI 现在不能单独支持什么：

```text
theta 和 RG subprocess 占据不同 BOLD ROI space。
```

目前最有意思的 ROI 结果：

```text
raw_csplit_q070 x roi_mean:
    deconv-selected ROI profile 跨 dataset 更一致；
    efun-selected ROI profile 不太一致；
    efun vs deconv 的 overlap/correlation 比较低。
```

所以 ROI 叙事现在应该写成：

```text
roi_mean 可能显示 deconv-related raw-density coupling 有更一致的 BOLD ROI footprint。
```

而不是：

```text
ROI 已经证明 theta/RG spatially separated。
```

### 叙事 H：BOLD observable 的选择本身很重要

之前可能会把 BOLD observable 混在一起看：

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

现在已经明确：

```text
这些不能混在一起解释。
```

当前理解：

```text
gsvd100_ds:
    很稳定，但可能太 global。
    容易显示 shared BOLD subspace。

global_svd100:
    中间状态。

roi_mean:
    更可能显示 efun/deconv 分离，
    尤其在 raw_csplit_q070 ROI footprint check 里。

HP_svd100:
    local observable。
    不能和 global ROI-vector correlation 混在一起。
```

现在的状态：

```text
保留为正式分层变量。
```

后续任何 P8/P10 结论都应该先分 BOLD observable，再讨论是否能合并。

### 叙事 I：ripple-trigger network reorganization / memory consolidation

这是目前最有生物学意义的一条解释线。

它来自：

```text
deconv_efun -> RG-no-pure-theta / ripple-like BLP subprocess
```

可以理解成：

```text
ripple/gamma-associated BLP subprocess
可能标记或驱动 deconvolved BOLD perturbation。
```

进一步的工作假说：

```text
这种 perturbation 可能对应 ripple-trigger network reorganization。
在 resting/spontaneous data 中，它可能和 replay / memory consolidation
相关的 intrinsic state perturbation 有关。
```

现在的状态：

```text
可以作为讨论方向和工作假说。
不能写成已经证明 memory consolidation。
```

推荐措辞：

```text
These results are consistent with a ripple-trigger network reorganization
hypothesis, potentially related to memory-consolidation-like intrinsic
state perturbations.
```

不推荐措辞：

```text
This proves memory consolidation.
```

### 历史叙事总表

| 历史叙事 | 现在状态 | 怎么处理 |
|---|---|---|
| event density 解释 BOLD | 降级 | 作为 baseline / control |
| pure ripple selective component | 替换 | 改成 RG-no-pure-theta |
| standardized csplit 更适合 fast component | 保留 | 当前主线 condition |
| theta slow + RG fast 双 subprocess | 部分保留 | RG/deconv 保留，theta/efun 降级 |
| efun slow state vs deconv fast activity | 部分保留 | deconv fast 更稳，efun slow 待补 |
| dimred label 解释 top coupling | 保留 | 重点看 deconv 是否偏 RG |
| ROI 证明 spatial separation | 降级 | 作为 spatial QC，不作主证据 |
| BOLD observable 可以混合看 | 修正 | 必须分 observable |
| ripple-trigger network reorganization | 保留为假说 | 可讨论，不写成因果证明 |

一句话总结：

```text
之前很多叙事都不是完全错，而是层级不同。

现在主线应该收窄到：
deconv_efun -> RG-no-pure-theta / ripple-like BLP subprocess。

theta/efun、ROI spatial separation、memory consolidation 都先作为后续扩展层，
不要和当前主证据混在一起。
```

## 1. 当前最稳的结论

最稳的结论是：

```text
BOLD deconv_efun 的 top coupling 富集在 RG-no-pure-theta / ripple-like BLP dimred components。
```

四个新 K13 数据里，P10 deconv_efun 特别明显：

```text
P10 deconv_efun top dimred hits:

top1:
    RG-no-theta 80.4%
    ripple       8.8%

top5:
    RG-no-theta 75.4%
    ripple      15.1%

top10:
    RG-no-theta 69.1%
    ripple      15.6%
```

这说明：

```text
deconv_efun 这条 BOLD 分支，比普通 BOLD efun 更干净地对应 fast / ripple-gamma-like BLP subprocess。
```

目前可以把它作为主线。

## 2. 之前的双 subprocess 叙事要降级

之前我们比较想讲的是：

```text
theta-like slow/state BLP subprocess
    -> BOLD efun

ripple-gamma-like fast/intrinsic-trigger BLP subprocess
    -> BOLD deconv_efun
```

这个叙事在 E10gb1 上很漂亮，但加入 K13m18-K13m21 后，需要更保守。

K13m18-K13m21 的 P5 strict gate 结果是：

```text
k13m18:
    有 strict theta + RG pair
    5 个 method-k 通过

k13m19:
    没有 strict theta + RG pair

k13m20:
    没有 strict theta + RG pair

k13m21:
    没有 strict theta + RG pair
```

重要解释：

```text
这不说明 RG/deconv 结论不成立。

它说明 “theta + RG 两边都跨 K13 稳定复现” 这个说法目前太强。
```

所以现在应该改成：

```text
RG/deconv 是主线。
theta/efun 暂时只是次级假说。
```

## 3. 当前叙事分层

### 第一层：可以作为主结论

```text
Fast RG-like BLP subprocess couples to BOLD deconv_efun.
```

证据：

```text
P10 deconv_efun top dimred hits 强烈偏 RG-no-theta / ripple-like labels。
P8 方向也类似，只是更弱、更 mixed。
```

推荐写法：

```text
The deconvolved BOLD eigenfunction branch preferentially couples to
ripple/gamma-like BLP subprocesses.
```

中文解释：

```text
BOLD deconv_efun 可能对应一种 fast intrinsic perturbation，
而这种 perturbation 和 ripple/gamma-like BLP Koopman density 有关。
```

### 第二层：可以作为技术支持

```text
BLP Koopman/eigenfunction density 比 simple event density 更有解释力。
```

证据：

```text
P8 和 P10 的 top xcorr hits 基本都偏 raw/dimred BLP efun density，
而不是 P2 event density。
```

这点很重要，因为它说明：

```text
有用的不是简单的 event occurrence，
而是 Koopman-derived latent activity / density representation。
```

换句话说，P5/P8/P10 的价值不是在重复 event detection，而是在找 latent BLP subprocess。

### 第三层：保留为假说，不作为当前结论

```text
BOLD efun = theta / slow-state variable
```

这个现在还不能作为强结论。

原因：

```text
K13m18-K13m21 里，strict theta 在多数数据中不稳定。
BOLD efun top dimred hits 也不干净，经常是 mixed / inactive / RG-like。
```

所以现在只能说：

```text
BOLD efun 可能和 slow/session-state-like BLP dynamics 有关，
但 theta-specific interpretation 还没有在 K13 子集里稳定成立。
```

后面如果要继续支持这条，需要看：

```text
raw BLP efun timescale distribution
whole-trace examples
session-wise state-change evidence
theta-active 的 softer label
```

### 第四层：ROI 只能作为空间 sanity check

ROI 结果现在不能作为两个 subprocess 空间分离的主证据。

目前情况是：

```text
ROI 结论强烈依赖 BOLD observable。
```

K13m18-K13m21 的 raw_csplit_q070 ROI 结果：

```text
gsvd100_ds:
    非常稳定，但是 efun 和 deconv 往往落在同一个稳定 BOLD subspace。
    所以它支持 stable BOLD subspace，
    但不支持 efun/deconv spatial separation。

roi_mean:
    更有意思。
    deconv-selected ROI profiles 跨 dataset 更一致。
    efun-selected ROI profiles 不太一致。
    efun vs deconv 的 same-dataset ROI corr 也低。
```

所以：

```text
ROI 可以用来做空间重复性 / QC。
但不能单独证明 theta 和 RG subprocess 占据不同 BOLD ROI space。
```

目前比较合理的说法是：

```text
roi_mean 下，deconv-related raw-density coupling 可能有更一致的 BOLD ROI footprint。
```

## 4. 所有 K13 后续应该怎么放在一起看

不要再把所有 K13 直接塞进一堆 summary 图里。

下一步应该先做一个 K13 dataset-status matrix。

每一行是一个 dataset：

```text
k13m17
k13m18
k13m19
k13m20
k13m21
k13m23
以后更多 K13
```

每一列是一个 gate 或现象：

```text
data QC status
P5 strict theta 是否存在
P5 RG-like component 是否存在
P5 strict theta+RG pair 是否通过
P8 deconv RG enrichment
P10 deconv RG enrichment
raw/dimred density preference
raw_csplit_q070 roi_mean deconv ROI consistency
HP/local caveat
notes
```

这样可以避免一个问题：

```text
某个 dataset theta gate 不过，
不等于它不能支持 RG/deconv 这条主线。
```

现在最应该看的 K13 问题是：

```text
1. 哪些 K13 数据是 usable？
2. 哪些 K13 有 RG-like BLP components？
3. deconv_efun -> RG enrichment 是否跨 K13 重复？
4. 这个 enrichment 是否 P10 比 P8 更强？
5. roi_mean 是否显示 deconv-related ROI footprint 更稳定？
```

这个读法比强行问 “theta + RG 是否每个 dataset 都成对出现” 更适合当前数据。

## 5. 建议保留的主图

### 图 1：K13 dataset-status matrix

目的：

```text
先告诉自己哪些 dataset 能用于什么结论。
```

这是最重要的总览图。

### 图 2：P5 RG-like component presence map

目的：

```text
看 RG-like BLP components 是否跨 K13 存在。
```

行：

```text
dataset
```

列：

```text
method-k
```

颜色：

```text
RG-like strict label present / absent
```

theta 可以作为附加标记，但不要让 theta 变成当前主图的中心。

### 图 3：P8/P10 deconv strict-label composition

目的：

```text
看 deconv_efun 的 top dimred-density hits 是否偏 RG-no-theta / ripple-like。
```

面板：

```text
P8 deconv_efun
P10 deconv_efun
```

分：

```text
top1
top5
top10
```

这是当前最重要的 evidence figure。

### 图 4：Density-class competition

目的：

```text
证明 BLP efun-derived density 比 event density 更有用。
```

分类：

```text
event density
raw BLP efun density
dimred BLP efun density
```

这张图支持方法学价值。

### 图 5：raw_csplit_q070 ROI footprint check

目的：

```text
看 efun 和 deconv 是否选中了不同 BOLD ROI/mode footprint。
```

重点看：

```text
roi_mean
```

原因：

```text
gsvd100_ds 往往太稳定，容易把 efun/deconv 都压到同一个 global subspace。
roi_mean 反而更可能显示 efun/deconv 分离。
```

## 6. 当前最保守、最清楚的整合结论

可以这样写：

```text
Across the currently analyzed standardized complex-split K13 datasets,
the most reproducible cross-modal finding is not a full theta-vs-ripple
two-process split.

Instead, the robust signal is that BOLD deconv_efun, especially after
BOLD-side dimensional reduction in P10, preferentially couples to
ripple/gamma-like BLP Koopman density components that are not pure-theta.

This supports a fast intrinsic-trigger interpretation:
ripple/gamma-associated BLP subprocesses may drive or mark deconvolved
BOLD perturbations, consistent with a ripple-trigger network reorganization
hypothesis.

The theta/slow-state BOLD efun branch remains plausible but is not yet
established in the new K13 subset.
```

中文意思：

```text
目前跨 K13 最稳定的发现，不是完整的 theta-vs-ripple 双 subprocess 分离。

更稳定的是：
BOLD deconv_efun，尤其是 P10 里的 BOLD-side dimred deconv_efun，
更偏向和 RG-no-pure-theta / ripple-like BLP Koopman density components 耦合。

这支持一个 fast intrinsic-trigger 的解释：
ripple/gamma-associated BLP subprocess 可能标记或驱动 deconvolved BOLD perturbation，
从而对应 ripple-trigger network reorganization。

theta/slow-state BOLD efun 仍然可能成立，但现在还不是 K13 子集里的强结论。
```

## 7. 下一步最值得做的事

1. 做 all-K13 dataset-status matrix。

2. 把 k13m17 和 k13m23 加进来，但标清 caveat：

```text
k13m17:
    theta event/raw signal QC 可疑

k13m23:
    需要确认 current P8/P10 coverage 和 BOLD observable availability
```

3. 做 all-K13 P8/P10 deconv label-composition figure。

4. 做 all-K13 density-class competition figure。

5. theta 先不要硬讲，后面单独用 softer theta-active label + raw timescale + whole-trace 去看。

## 8. 一句话版本

如果你现在很累，只记这一句：

```text
当前最稳的故事是 deconv_efun -> RG-no-pure-theta / ripple-like BLP subprocess。

这可以支持 ripple-trigger network reorganization 的工作假说。

theta/efun 那条线先不要作为主结论，等后面用 softer label 和 raw timescale 再看。
```
