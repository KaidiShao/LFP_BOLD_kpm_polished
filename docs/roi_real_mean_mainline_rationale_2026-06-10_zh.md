# ROI 主线口径：为什么当前只用 real_mean

日期：2026-06-10

本文件记录当前 cross-dataset SOP 中 ROI summary 的口径决定：

```text
ROI mainline = real_mean
mean_abs     = legacy / QC only
positive/negative ROI = 从 real_mean 派生的展示，不作为独立 ROI condition
```

对应主 SOP：

```text
docs/current_minimal_cross_dataset_consistency_sop_2026-06-10_zh.md
```

## 1. 当前 ROI 要回答的问题

ROI 层面现在不是用来单独证明 theta/RG 是两个空间分离网络，而是用来回答：

```text
P8 top xcorr 选出来的 BOLD modes
是否有稳定、可解释的 signed spatial footprint？

efun selected modes 和 deconv_efun selected modes
是否对应不同的 BOLD spatial readout？
```

因此 ROI 图必须保留正负方向。否则我们只能看到 involvement strength，
看不到空间 readout 的 polarity / signed pattern。

## 2. 为什么不用 mean_abs 做主口径

`mean_abs` 的含义是：

```text
每个 ROI 内取 mode value 的绝对值再平均。
```

它回答的是：

```text
这个 BOLD mode 在该 ROI 上参与强不强？
```

但它不回答：

```text
这个 ROI 在该 BOLD mode 中是正向还是负向？
不同 ROI 之间是否构成 signed spatial pattern？
```

之前的结果显示，`mean_abs` 会把很多可能有解释意义的正负差异折叠掉。
这会导致：

```text
theta/RG target 的 ROI profile 看起来更像；
efun/deconv 的差异也容易被压扁；
最后只能得到“magnitude footprint 没有明显分开”的弱结论。
```

所以当前结论是：

```text
mean_abs 只能作为 legacy / QC baseline。
它不能作为当前主叙事的 ROI 证据口径。
```

## 3. 为什么 positive_real / negative_real 不作为独立主口径

`positive_real` 和 `negative_real` 是有用的展示方式，但不适合作为独立主分析条件。

原因是：

```text
它们其实是 real_mean signed vector 的两个派生视角；
如果把它们当成和 real_mean 并列的 condition，
会把同一个 signed spatial pattern 人为拆成多套结果，
让 ROI 结论变得更碎、更难比较。
```

当前使用方式应该是：

```text
先生成 real_mean ROI profile。
再从同一个 real_mean vector 中提取：
    top positive ROI
    top negative ROI
    signed profile correlation
    signed top-ROI overlap
```

也就是说：

```text
positive / negative 是 real_mean 的可视化和解释标签，
不是新的 ROI value mode。
```

## 4. 为什么 real_mean 是当前主线

`real_mean` 的含义是：

```text
每个 ROI 内取 BOLD mode real part 的 signed mean。
```

它保留了三个当前最需要的信息：

```text
1. ROI 参与方向：positive vs negative
2. ROI spatial pattern：不同 ROI 的 signed contrast
3. across-dataset comparison：可以直接做 signed profile correlation
```

对于当前叙事，它最适合检验：

```text
efun×theta-like selected modes
和 deconv_efun×RG-like selected modes
是否落在不同或相似的 signed ROI readout 上。
```

也适合检验更稳的版本：

```text
efun selected modes vs deconv_efun selected modes
是否有 feature-family 层面的 signed ROI difference。
```

## 5. 之前结果带来的具体判断

之前做过的 ROI/topN/control 分析给出的方向是：

```text
1. mean_abs ROI 太粗，容易抹平差异。

2. signed ROI 比 mean_abs 更合理，
   但 ROI 差异主要体现为 efun vs deconv_efun feature-family difference，
   而不是同一个 feature family 内 theta vs RG 的强空间分离。

3. efun_theta 和 efun_RG 的 ROI profile 可以很像；
   deconv_theta 和 deconv_RG 的 ROI profile 也可以很像。

4. 因此 ROI 不能作为“theta/RG 是两个完全不同空间网络”的主证据。

5. 更稳的主证据仍然是：
   P5 subprocess gate
   P8/P10 deconv RG enrichment
   P8 selected BOLD mode eigenvalue/timescale

6. ROI 的角色是：
   检查 selected BOLD modes 是否有 signed spatial footprint，
   并判断 efun/deconv 两类 readout 的 spatial interpretation 是否合理。
```

## 6. 当前正式图

ROI 主线图全部使用 `real_mean`。

必画：

```text
1. real_mean all-ROI-by-dataset
   - 保留原始 anatomical ROI 顺序
   - dataset 分面
   - 每个 dataset 可用自己的 symmetric signed color scale
   - 横向 layout，适合宽屏查看

2. real_mean pairwise ROI profile correlation
   - targets:
       efun_theta
       efun_RG
       deconv_theta
       deconv_RG
   - topN:
       top3, top5, top10, top20

3. real_mean top positive / top negative ROI labels
   - 从同一个 real_mean vector 排序得到
   - 不重新定义 ROI value mode

4. selected BOLD mode overlap
   - exact mode overlap
   - adjacent / nearby sorted-mode overlap
   - 用于判断 ROI 相似是否只是选到了同一批 BOLD modes
```

可选：

```text
mean_abs legacy comparison
```

但它只用于说明为什么不用 magnitude 口径，不进入主结论图组。

## 7. 实现口径

P7 ROI profile 导出应使用：

```text
--roi-value-mode real_mean
```

后续 ROI 脚本应默认读取 / 生成 `real_mean`：

```text
scripts/export_p7_roi_mean_profiles_direct_from_bold_post.py
scripts/analyze_roi_mean_signed_target_modes.py
scripts/analyze_roi_mean_signed_target_pairwise.py
```

如果需要 top positive / negative ROI：

```text
不要重新导出 positive_real / negative_real 作为主结果；
直接从 real_mean profile 中排序：
    largest positive values
    most negative values
```

## 8. 推荐写法

可以写：

```text
ROI analyses were performed on signed real_mean ROI profiles, preserving
the polarity of BOLD Koopman modes.  Magnitude-only mean_abs profiles were
treated as legacy/QC summaries because they collapse positive and negative
spatial structure.
```

中文：

```text
ROI 分析使用 signed real_mean profile，以保留 BOLD Koopman mode 的空间极性。
mean_abs 只作为旧版/QC 口径，因为它会把正负空间结构折叠成 magnitude，
不适合作为当前 efun/deconv spatial readout 叙事的主证据。
```

当前结论应写成：

```text
real_mean ROI profile 用于支持或约束 efun/deconv readout 的空间解释；
它目前不是 theta/RG spatial separation 的主证据。
```
