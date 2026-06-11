# 用 Koopman Eigenfunction 定义海马瞬态事件：从频段标签到内在时间尺度子过程

## 核心观点

这条叙事的中心句是：

> **Frequency bands describe the oscillatory signature of hippocampal events, whereas Koopman eigenfunction clusters describe the intrinsic dynamical subprocesses that generate those signatures. Classical events such as theta, gamma, ripple, and SPW-R are recurrent trajectory motifs in a low-dimensional space of Koopman subprocesses.**

换句话说，这个解释不是要推翻 theta、gamma、ripple 等传统频段标签，而是把这些标签放进一个更完整的 Koopman spectral representation 中。传统频段主要描述事件的振荡内容，而 Koopman eigenfunctions 进一步提供事件的内在时间尺度、持续性、衰减和状态转移结构。

---

## 1. 为什么单纯用频段定义事件不够

经典做法通常将海马 LFP 瞬态事件定义为特定频段能量的增强，例如 theta、gamma、ripple 或 sharp wave-ripple。这个定义有清晰的生理和文献基础，但它把事件简化成了人为指定频段的 power feature。

实际数据中有几个现象说明这种定义不充分：

1. gamma 和 ripple 经常共同出现，因此它们未必是两个互斥事件类别；
2. isolated ripple 和 SPW-R 共享高频 ripple 成分，但 SPW-R 还包含 slow sharp-wave 成分；
3. 一个 event 可以同时包含多个频率成分；
4. transient event 可能更像一个 state transition 或 trajectory motif，而不是某个时间点的 spectral peak。

因此，频段标签描述了事件的 oscillatory signature，但没有完整定义产生该 signature 的 dynamical subprocess。

![Multi-region Koopman event definition](koopman_events_multiregion.png)

---

## 2. NMF 的贡献与局限

Band-free spectral decomposition，例如 NMF，是对固定频段定义的重要推进。NMF 将 spectrogram 分解为若干 data-driven spectral profiles：

\[
S(f,t) \approx \sum_{k=1}^{K} W_k(f)H_k(t).
\]

其中，\(W_k(f)\) 表示 spectral profile，\(H_k(t)\) 表示该 profile 随时间的贡献。

NMF 的优点是：它不需要预先指定 theta、gamma、ripple 等频段，可以发现 mixed-frequency patterns。例如 SPW-R 可以被识别为 low-frequency sharp wave 与 high-frequency ripple 的复合事件。

但 NMF component 仍主要是 **spectral factor**，不是 **dynamical state variable**。它告诉我们：

> 这个 transient 的频谱形状是什么？

但它不自然回答：

> 这个 transient 属于哪个内在动力学子过程？它如何持续、衰减、转移，或者与其他子过程共同激活？

所以，Koopman 的位置不是否定 NMF，而是进一步推进：

\[
\text{frequency band}
\rightarrow
\text{spectral profile}
\rightarrow
\text{dynamical subprocess}.
\]

---

## 3. Koopman eigenvalue 同时编码频率与时间尺度

在连续时间 Koopman generator 表示中，一个 eigenfunction \(\phi_i\) 沿轨迹满足：

\[
\phi_i(x_t) \sim e^{\nu_i t},
\qquad
\nu_i = \alpha_i + i\omega_i.
\]

其中：

- \(\omega_i = \operatorname{Im}(\nu_i)\)：振荡频率；
- \(\alpha_i = \operatorname{Re}(\nu_i)\)：增长、衰减、持续性和内在时间尺度。

在离散时间 EDMD/Koopman 表示中，特征值通常写作：

\[
\mu_i = \rho_i e^{i\theta_i}.
\]

如果采样间隔为 \(\Delta t\)，则可以转换为连续时间指数：

\[
\nu_i = \frac{\log \mu_i}{\Delta t}
= \frac{\log \rho_i}{\Delta t} + i\frac{\theta_i}{\Delta t}.
\]

因此：

\[
\alpha_i = \frac{\log \rho_i}{\Delta t},
\qquad
\omega_i = \frac{\theta_i}{\Delta t}.
\]

衰减时间尺度可以定义为：

\[
\tau_i^{\mathrm{decay}} = -\frac{1}{\alpha_i}
= -\frac{\Delta t}{\log \rho_i},
\qquad \alpha_i < 0.
\]

所以，传统频段主要对应 \(\omega_i\)，而 subprocess 更依赖 \(\alpha_i\) 或 \(\rho_i\) 所刻画的 persistence / decay / timescale。

这就是这个解释优雅的地方：

> 频段理论不是错的；它只是看了 Koopman eigenvalue 的一部分。

---

## 4. Subprocess 应该是 eigenfunction timescale cluster，而不是单个 eigenfunction

真实数据中不应把单个 eigenfunction 直接解释成一个生物事件。更稳妥的定义是：

> **一个 Koopman subprocess 是一组具有相似 eigenvalue timescale、相似 activation pattern、并参与相似 event trajectory motifs 的 eigenfunctions。**

形式上，可以定义一个 eigenfunction cluster：

\[
\mathcal P_r = \{i : \tau_i \text{ similar},\; \mathrm{activation}_i \text{ similar},\; \mathrm{trajectory}_i \text{ similar}\}.
\]

对应的 subprocess activity 可以写作：

\[
S_r(t)=\left(\sum_{i\in\mathcal P_r} w_i |\phi_i(x_t)|^2\right)^{1/2}.
\]

这样，dimension-reduced eigenfunction 不是任意 UMAP 可视化，而是高维 eigenfunctions 中冗余 timescale 信息的低维 summary。

![Low-dimensional Koopman event state space](koopman_lowdim_event_states.png)

---

## 5. Event 是低维 Koopman subprocess space 中的短轨迹 motif

最推荐的事件定义是：

> **A transient event is a recurrent short trajectory motif in a low-dimensional Koopman subprocess space.**

设低维 subprocess coordinate 为：

\[
z(t) = [S_1(t), S_2(t), \ldots, S_R(t)].
\]

那么一个 event 不是单点 \(z(t)\)，而是短轨迹片段：

\[
z(t:t+L).
\]

这样可以自然解释传统事件之间的关系：

| 传统事件 | Koopman 解释 | 动力学含义 |
|---|---|---|
| Theta | slow subprocess motif | 慢状态或慢 regime axis |
| Gamma | fast subprocess activation | fast transition / burst axis |
| Ripple | high-frequency fast subprocess | 强 fast activation |
| Isolated ripple | fast subprocess only | 缺少 slow sharp-wave component |
| SPW-R | slow + fast compound motif | slow sharp-wave subprocess 与 fast ripple subprocess 协同激活 |

因此，SPW-R 与 isolated ripple 并不是完全无关的两个事件。它们共享 fast ripple subprocess，但 SPW-R 额外包含 slow sharp-wave subprocess。

---

## 6. 与 neural mass model 的对应

Bi-state neural mass model 是这个解释的 proof-of-concept。在这个模型中，ground truth 包含：

- metastable State 0；
- metastable State 1；
- state transition；
- transition 过程中的 fast oscillation。

Koopman eigenfunctions 可以分别捕捉这些结构：metastable states、transition path 和 fast oscillations。这说明 Koopman eigenfunctions 不是任意 feature，而是能够在已知动力系统中分离 state、transition 和 oscillatory subprocess 的坐标。

![Bi-state neural mass Koopman eigenfunctions](koopman_neural_mass_model.png)

对应关系可以写作：

| Neural mass model | Hippocampal LFP interpretation |
|---|---|
| Metastable state | baseline / theta-like slow regime |
| Transition eigenfunction | sharp-wave / event onset subprocess |
| Fast oscillatory eigenfunction | gamma/ripple fast subprocess |
| Trajectory in eigenfunction space | transient event motif |
| State transition + oscillation | SPW-R / gamma-ripple compound event |

---

## 7. 与后续 LFP--BOLD 课题的关系

这个 event definition 是后续 LFP--BOLD 课题的前置问题。首先需要证明：

\[
\text{LFP Koopman subprocess density}
\rightarrow
\text{BOLD modal innovation}.
\]

将高采样率 LFP 映射到 BOLD 时间尺度时，可以定义三种 density：

1. frequency-defined event density；
2. raw LFP eigenfunction activation density；
3. dimension-reduced LFP eigenfunction subprocess density。

如果第三种 density 最能预测 BOLD modal innovation，则说明真正与全局 BOLD network reorganization 相关的不是单一频段 event，也不是高维 eigenfunction 噪声，而是低维 Koopman-defined hippocampal subprocess。

---

## 8. 可验证预测

这条叙事产生几个清楚的预测：

1. Koopman eigenfunction clusters 应显示不同 decay / persistence / autocorrelation timescales；
2. manual event labels 在低维 Koopman subprocess space 中不应完全互斥，而应表现为 overlapping trajectory motifs；
3. gamma 和 ripple 应共享 fast subprocess axis；
4. SPW-R 应表现为 slow subprocess 与 fast subprocess 的复合轨迹；
5. dimension-reduced Koopman subprocess density 应比 frequency-band density 或 NMF component density 更好地预测 BOLD modal innovation；
6. 在 neural mass model 等有 ground truth 的系统中，Koopman eigenfunctions 应能分离 state、transition 和 oscillation。

---

## 9. 推荐写法

> Classical hippocampal events are usually defined by power in predefined frequency bands, such as theta, gamma and ripple. This convention captures the oscillatory frequency of an event, but it does not fully define its underlying dynamical subprocess. In Koopman theory, oscillatory frequency corresponds to the imaginary part of the continuous-time Koopman exponent, whereas persistence, decay and intrinsic timescale are encoded by the real part, or equivalently by the modulus of the discrete-time eigenvalue. Thus, frequency-band labels can be viewed as projections of a richer modal description. We therefore define events not by single frequency bands, but as recurrent short trajectory motifs in a low-dimensional space of Koopman eigenfunction clusters. Under this view, theta reflects a slow subprocess, gamma and ripple reflect fast subprocesses, and sharp-wave ripples correspond to compound motifs involving coordinated activation of slow and fast subprocesses.

---

## 10. Source notes

- Besserve, Safavi, Schölkopf, and Logothetis, *The complex spectral structure of transient LFPs reveals subtle aspects of network coordination across scales and structures*. This work motivates band-free decomposition of transient LFPs and shows that SPW-R is a compound low-frequency and high-frequency event.
- Shao, *Annual Meeting 2025: Peering into Intrinsic Brain Dynamics with the Koopman Operator*. Slides 7--8 motivate Koopman eigenfunctions as event-related state coordinates and low-dimensional summaries of redundant intrinsic timescale information.
- Bi-state neural mass model example from the annual meeting material. The model illustrates how Koopman eigenfunctions can separately capture metastable states, transitions and fast oscillations.
