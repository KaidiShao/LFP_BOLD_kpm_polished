# Residual EDMD Internal Perturbation Test

这个文件夹整理了当前这一轮 `EDMD output -> timescale -> internal perturbation` 的测试材料，重点是比较不同 internal perturbation 定义和不同 eigenvalue 来源对热图结果的影响。

## 目录内容

- `script_visualize_edmd_outputs_test_residual_edmd.m`
  - 当前分析脚本的快照版本。
  - 已改成默认把图直接保存到当前文件夹。
- `01_postprocess_main.*`
  - 基础 EDMD 可视化。
- `02_timescale_diagnostics.*`
  - timescale 诊断图。
- `deconv_*.png / deconv_*.fig`
  - 8 个 internal perturbation 条件下的对照图。
- `03_deconv_condition_similarity.*`
  - 8 个条件的热图相似度比较。

## 八个 Condition

当前比较是三维组合：

1. `method`
   - `koopman_residual`
     - 定义为 `u_t = phi_t - lambda * phi_{t-1}`。
     - 这是最直接贴合离散 Koopman 递推关系的定义。
   - `wiener`
     - 把 `phi_t` 视为指数核 `lambda^(k-1)` 卷积后的结果，再做 Wiener deconvolution。
     - 更偏信号处理视角。

2. `lambda_source`
   - `edmd`
     - 直接使用 EDMD 输出的 eigenvalue。
   - `empirical_complex`
     - 不借 EDMD 相位，直接从复数 eigenfunction 序列做一步复数回归拟合：
       - `lambda_hat = sum(conj(phi_{t-1}) * phi_t) / sum(|phi_{t-1}|^2)`

3. `normalization`
   - `global`
     - 用整段时间序列做列归一化，再显示当前窗口。
     - 当前默认在归一化时排除 `t=1`，避免起点边界效应主导颜色。
   - `window`
     - 只在当前可视化窗口内逐列归一化。

因此 8 个条件分别是：

- `residual | edmd | global`
- `residual | empirical complex | global`
- `wiener | edmd | global`
- `wiener | empirical complex | global`
- `residual | edmd | window`
- `residual | empirical complex | window`
- `wiener | edmd | window`
- `wiener | empirical complex | window`

## Similarity 图怎么读

`03_deconv_condition_similarity.png` 目前比较的是两张 `ALL` 热图：

- `ALL |u|`
- `ALL Re(u)`

每种热图都给出两类相似度：

- `corr`
  - 把整张热图拉平成向量后做 Pearson correlation。
  - 越接近 `1`，说明图案结构越相似。
- `mean |diff|`
  - 两张热图逐点取绝对差，再求平均。
  - 越接近 `0`，说明数值上越接近。

## 当前 Similarity 的主要结论

从当前保存的相似度图来看，有两个比较稳定的现象：

1. `global` 和 `window` 各自内部都比较相似。
   - 也就是说，在同一种 normalization 里，`residual/wiener` 和 `edmd/empirical complex` 的结果差异不算大。

2. 真正更明显的差别主要来自 `global` 和 `window` 两种 normalization 本身。
   - 也就是说，归一化方式带来的视觉差异，大于 `method` 或 `lambda_source` 带来的差异。

因此，这一轮 comparison 更像是在说明：

- `method` 和 `lambda_source` 会影响结果，但影响相对较小；
- `global` vs `window` 会更明显地改变热图外观。

## 为什么最后选择 `residual + EDMD lambda` 作为主分析

最终把 `residual + EDMD lambda` 作为主分析，主要基于理论自洽性，而不是因为别的组合完全错误。

理由如下：

1. `residual` 最贴近离散 Koopman 递推定义。
   - 它直接对应：
     - `phi_t = lambda * phi_{t-1} + u_t`
     - `u_t = phi_t - lambda * phi_{t-1}`

2. `EDMD lambda` 与当前 eigenfunction 来自同一个 EDMD 算子估计。
   - 因此它和 `phi` 是同一套模型内部的成对量，定义最自洽。

3. `empirical_complex` 适合作为 robustness check。
   - 它是很有价值的敏感性分析，用来检查“如果完全从数据本身拟合 complex lambda，结果会不会变很多”。
   - 但从主分析角度，它更像 supplementary comparison，而不是第一定义。

## 当前建议

- 主分析：
  - `residual | edmd`
- 归一化：
  - `global` 用作主结果比较
  - `window` 用作局部时窗展示
- 其余组合：
  - 作为 sensitivity / supplementary comparison 保留

## 备注

- 如果后续重新运行本文件夹里的脚本，图会继续直接保存到这个文件夹。
- 如果只想保留 PNG，可以把脚本里的 `viz_cfg.save_fig = false;`。
