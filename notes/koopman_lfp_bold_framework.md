# 非线性 local-global latent system、Koopman modal signatures 与两个 LFP-BOLD 课题的统一框架

**版本：2026-06-10**

这份笔记总结当前对两个相关课题的统一数学理解：

1. **Event-related LFP subprocess -> global BOLD pattern / BOLD eigenfunction residual**：局部 LFP subprocess 如何作为 intrinsic perturbation，触发或预测全局 BOLD modal innovation。
2. **Non-causal LFP-BOLD IR / pre-activation**：为什么 local BOLD 相对于 local LFP 会出现 apparent pre-activation，以及如何用 closed-loop local-global latent dynamics 和 Koopman modal representation 解释。

核心思想是：

> LFP、local BOLD、global BOLD 不应被简单看作 feed-forward input-output 关系，而应看作一个 local-global nonlinear latent system 的不同观测读出。Koopman eigenfunctions 将这个复杂非线性系统转化为可解释的 modal variables，使 global gating state、pre-activation modes、global ROI map 和 loop backbone 成为可操作的 phenomenological modal signatures。

---

## 1. 最底层模型：分块的非线性 local-global latent system

令 latent state 分为 local 和 global 两部分：

$$
x_t =
\begin{bmatrix}
x_t^{(l)} \\
x_t^{(g)}
\end{bmatrix}.
$$

其中：

- $x_t^{(l)}$：local latent state，例如 hippocampal local neural state、local vascular state、local circuit excitability；
- $x_t^{(g)}$：global latent state / internal environment / global order-parameter-like state，例如 global brain state、neuromodulatory state、brainstem-thalamic-cortical state、global vascular/arousal state；
- $u_t$：local LFP subprocess / local neural action。

一个 generic nonlinear model 可以写成：

$$
x_{t+1}^{(l)} = F_l\!\left(x_t^{(l)}, x_t^{(g)}, u_t\right),
$$

$$
x_{t+1}^{(g)} = F_g\!\left(x_t^{(g)}, x_t^{(l)}, u_t\right),
$$

或简写为：

$$
x_{t+1}=F(x_t,u_t).
$$

这个模型的关键不是假设 LFP 是外部刺激，而是把 $u_t$ 视为 local subsystem 的 action-like subprocess。它可以暂时作为 input-like variable 写入模型，但在 closed-loop 视角下，它本身也是系统内部生成的。

---

## 2. BOLD 是 nonlinear readout，不是机制本身

local BOLD 和 global BOLD 都是 latent state 的观测读出：

$$
y_t^{(lB)} = h_l\!\left(x_t^{(l)}, x_t^{(g)}\right),
$$

$$
y_t^{(gB)} = h_g\!\left(x_t^{(l)}, x_t^{(g)}\right).
$$

因此，local BOLD 不是单纯由 local LFP 经过 HRF 产生的 output。它可以同时包含：

$$
y_t^{(lB)} = \text{local component} + \text{global component}.
$$

这点非常关键：如果 local BOLD 提前读出了某个 global latent state，那么 local BOLD 相对于未来的 local LFP subprocess 就可能出现 apparent pre-activation。

---

## 3. Global-to-local 路径：global state gates local LFP subprocess

当前对 pre-activation 的最核心解释是：

$$
x_t^{(g)} \rightarrow u_{\mathrm{LFP},t+\Delta}.
$$

更具体地，可以写成：

$$
x_t^{(g)} \rightarrow x_{t+\Delta}^{(l)} \rightarrow u_{\mathrm{LFP},t+\Delta}.
$$

意思是：

> 一个 global latent state / order-parameter-like variable 先出现，然后调节 local hippocampal condition，使某个 local LFP subprocess 更容易发生。

如果要构成完整 closed loop，还可以有 local-to-global 的反向一段：

$$
u_{\mathrm{LFP},t} \rightarrow x_{t+\Delta}^{(g)}.
$$

因此整体闭环可以概括为：

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}} \rightarrow x^{(g)}.
$$

但对于 **pre-activation** 的解释，最关键的是前半段：

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}}.
$$

---

## 4. 为什么 local BOLD-LFP 会产生 apparent non-causal IR？

传统 feed-forward HRF 模型假设：

$$
y_{\mathrm{localBOLD}}(t)
=
\sum_{\tau\ge 0} h(\tau) u_{\mathrm{LFP}}(t-\tau)+\epsilon(t).
$$

这个模型默认 LFP 是 exogenous input，BOLD 是 passive output。因此，如果 estimated IR 在 $\tau<0$ 有 pre-activation，就看起来像 non-causal。

在 local-global latent model 中，情况不同。某个 global state $x_t^{(g)}$ 可以同时：

1. 被 local BOLD 提前读出：

$$
x_t^{(g)} \rightarrow y_t^{(lB)},
$$

2. gate 未来的 local LFP subprocess：

$$
x_t^{(g)} \rightarrow u_{\mathrm{LFP},t+\Delta}.
$$

于是观测上可能出现：

$$
y_{\mathrm{localBOLD},t}
\quad \text{早于} \quad
u_{\mathrm{LFP},t+\Delta}.
$$

如果这时仍用 feed-forward LFP-to-BOLD deconvolution，就会得到负时滞的 pre-activation。它并不意味着：

$$
y_{\mathrm{localBOLD}} \rightarrow u_{\mathrm{LFP}},
$$

而是说明：

$$
x^{(g)} \rightarrow
\begin{cases}
y_{\mathrm{localBOLD}},\\
u_{\mathrm{LFP}}.
\end{cases}
$$

因此，apparent non-causality 是 closed-loop endogenous system 被 open-loop feed-forward model 误读后的现象。

---

## 5. 从非线性系统到增广自治系统

如果 $u_t$ 是系统内部生成的 local neural action，而不是外部输入，那么可以写作：

$$
u_t = \pi\!\left(x_t^{(g)}, x_t^{(l)}, u_{t-1},u_{t-2},\ldots\right)+\epsilon_t.
$$

把 $u_t$ 及其历史增广到状态中：

$$
s_t = \begin{bmatrix}
x_t \\
u_t \\
u_{t-1} \\
\vdots
\end{bmatrix},
$$

整个 closed-loop system 可以写成 autonomous form：

$$
s_{t+1}=F_{\mathrm{cl}}(s_t).
$$

相应的观测读出可写为：

$$
u_t = C_u s_t,
$$

$$
y_t^{(lB)} = C_l s_t,
$$

$$
y_t^{(gB)} = C_g s_t.
$$

这一步非常重要：它说明 local LFP subprocess、local BOLD、global BOLD 都可以被看作同一个 autonomous closed-loop system 的不同 readouts。

---

## 6. Koopman modal representation

对增广自治系统：

$$
s_{t+1}=F_{\mathrm{cl}}(s_t),
$$

Koopman eigenfunction coordinates 记为 $a_t$，则在 modal coordinate 中：

$$
a_{t+1}=\Lambda a_t.
$$

观测量可以近似线性读出：

$$
u_{\mathrm{LFP},t}=C_u a_t,
$$

$$
y_{\mathrm{localBOLD},t}=C_l a_t,
$$

$$
y_{\mathrm{globalBOLD},t}=C_g a_t.
$$

Koopman 的好处在于：原始 latent dynamics 和 observation functions 都可以是非线性的，但在 lifted/modal space 中，dominant modes、phase relationships 和 spatial readouts 可以被线性或近似线性地分析。

---

## 7. Pre-activation mode 的定义

对于某个 LFP-BOLD anchor pair $(p,q)$，其中：

- $p$：某个 local LFP band/subprocess，例如 PGO、theta、ripple；
- $q$：某个 local BOLD ROI 或 local BOLD readout；

如果 empirical cross-kernel $h_{q,p}(\tau)$ 在 $\tau<0$ 出现 pre-activation，则可以在 Koopman modes 中寻找贡献这个 negative-lag component 的 modes。

定义：

$$
\mathcal I_{\mathrm{pre}}(q,p)
=
\{i: \text{mode } i \text{ contributes to the negative-lag component of pair }(q,p)\}.
$$

这些 modes 就是该 anchor pair 的 **pre-activation modes**。

它们的 biological interpretation 是：

> 这些 modes 最自然对应一个 global latent gating state / order-parameter-like variable。这个 global state 先被 BOLD readout 观测到，然后 gate 未来的 local LFP subprocess。

重要边界：pre-activation modes 不是 causal modes，而是 pair-conditioned modal signatures。

---

## 8. Phase lead 的 Koopman 解释

对某个 mode $i$，假设其 continuous-time eigenvalue 为：

$$
\nu_i = \alpha_i + \mathrm{i}\omega_i.
$$

在 LFP 和 BOLD readout 上的 complex weights 可写为：

$$
C_u(i)=|C_u(i)|e^{\mathrm{i}\phi_{u,i}},
$$

$$
C_l(i)=|C_l(i)|e^{\mathrm{i}\phi_{l,i}}.
$$

如果 local BOLD readout 的 phase 领先 LFP readout，则该 mode 会在 LFP-BOLD cross-kernel 中贡献 apparent pre-activation。对应的相位时间差可写作：

$$
\Delta t_{lB-u,i}
=
\frac{\operatorname{wrap}(\phi_{l,i}-\phi_{u,i})}{\omega_i},
$$

具体符号取决于 phase convention。直观上：

$$
\text{BOLD phase leads LFP phase}
\quad \Rightarrow \quad
\text{negative-lag contribution in empirical IR}.
$$

---

## 9. Global ROI map 的含义

对某个 anchor pair $(q,p)$ 的 pre-activation modes，global ROI map 可定义为：

$$
M_{q,p}(r)
=
\sum_{i\in\mathcal I_{\mathrm{pre}}(q,p)} s_i |C_g(r,i)|,
$$

其中：

- $r$：global BOLD ROI；
- $C_g(r,i)$：mode $i$ 在 ROI $r$ 上的 global BOLD readout weight；
- $s_i$：mode $i$ 对 pair-specific pre-activation 的贡献强度，例如 negative-lag area、mode amplitude 或 pairwise mode weight。

解释：

> global ROI map 是 pre-activation modes 在全脑 BOLD ROI 空间中的 spatial footprint。

在 closed-loop / slaving-principle 解释中，它进一步对应：

> gate local LFP subprocess 的 global latent state / order-parameter-like variable 在 BOLD 空间中的可观测足迹。

它不是“哪些 BOLD ROI 反向驱动 LFP”的 causal map，而是“哪些 ROI 共同表达了这个 global gating state”的 map。

---

## 10. Loop backbone 的含义

对同一组 pre-activation modes，每个 ROI 的 global BOLD readout 不仅有 amplitude，还有 phase：

$$
C_g(r,i)=|C_g(r,i)|e^{\mathrm{i}\phi_{r,i}}.
$$

global ROI map 使用：

$$
|C_g(r,i)|,
$$

而 loop backbone 使用：

$$
\phi_{r,i}.
$$

因此：

> global ROI map 描述这个 global state 在哪里表达；loop backbone 描述这个 global state 如何在 ROI 空间中按相位顺序表达。

可以把每个 mode 的 phase order 转成一个 phase-defined adjacency：

$$
A_i(r_1,r_2)=1
\quad \text{if ROI } r_1 \text{ phase-leads ROI } r_2
\quad \text{within mode } i.
$$

然后对 pre-activation modes 聚合：

$$
A_{q,p}^{\mathrm{backbone}}(r_1,r_2)
=
\sum_{i\in\mathcal I_{\mathrm{pre}}(q,p)} s_i A_i(r_1,r_2).
$$

loop backbone 的解释是：

> pre-activation modes 所定义的 global gating state，在全脑 ROI 空间中的 phase-ordered expression / phase portrait。

它是 candidate biological representation / phenomenological modal signature，不是解剖 feedback loop 本身，也不是 interventional causality。

---

## 11. Haken / slaving principle 的对应关系

这个框架可以与 Haken 的 slaving principle 对应起来：

$$
\text{distributed brain regions}
\Rightarrow
\text{global order parameter}
\Rightarrow
\text{local fast subprocess}.
$$

对应到当前问题：

$$
\text{distributed ROI-level activity}
\Rightarrow
x^{(g)}
\Rightarrow
u_{\mathrm{LFP}}.
$$

这里：

- $x^{(g)}$ 是 global order-parameter-like state / slaving variable；
- global ROI map 是 $x^{(g)}$ 的 spatial BOLD footprint；
- loop backbone 是 $x^{(g)}$ 在 ROI 空间中的 phase portrait；
- local LFP subprocess 是被 $x^{(g)}$ gate / enslave 的 fast local process。

这也形成一种 circular causality：

$$
\text{distributed global state}
\rightarrow
\text{local LFP subprocess}
\rightarrow
\text{global state update}.
$$

---

## 12. 两个课题在同一框架中的位置

### 12.1 课题一：Event-related LFP subprocess -> global BOLD pattern / BOLD eigenfunction residual

这个课题对应 local-to-global 路径：

$$
u_{\mathrm{LFP}} \rightarrow x^{(g)} \rightarrow y_{\mathrm{globalBOLD}}.
$$

实际操作中，先对 BOLD 单独建立 Koopman model：

$$
b_{k+1}=\Lambda_B b_k + \eta_k^B,
$$

其中 BOLD modal innovation / BOLD eigenfunction residual 为：

$$
\eta_k^B=b_{k+1}-\Lambda_B b_k.
$$

然后将高采样率 LFP subprocess 压到 BOLD 时间尺度，得到 density：

$$
D_r^{\mathrm{LFP}}[k]
=
\text{density / occupancy of LFP subprocess } r
\text{ in BOLD window } k.
$$

核心模型是：

$$
\eta_{m,k}^B
=
\sum_r\sum_{\ell}
\beta_{m,r,\ell}D_r^{\mathrm{LFP}}[k-\ell]
+
\epsilon_{m,k}.
$$

如果 $D_r^{\mathrm{LFP}}$ 能预测未来的 $\eta^B$，则说明：

> local LFP subprocess acts as an intrinsic perturbation that triggers or predicts global BOLD modal innovations.

注意：BOLD residual 本身不是 LFP subprocess activation；BOLD residual 是 BOLD modal dynamics 的 input-like deviation。只有当它被 LFP subprocess density 预测时，才支持 LFP subprocess 作为 intrinsic perturbation 的解释。

---

### 12.2 课题二：Non-causal LFP-BOLD IR / pre-activation

这个课题对应 global-to-local 路径：

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}}.
$$

具体表现为：

1. global latent state $x^{(g)}$ 先出现；
2. local BOLD 读出了 $x^{(g)}$ 的 local/global mixture；
3. 同一个 $x^{(g)}$ gate 未来的 local LFP subprocess；
4. 因此用 LFP-to-local-BOLD feed-forward IR 会看到 apparent pre-activation。

Koopman 分析中：

- anchor pair $(q,p)$ 选出 pre-activation modes $\mathcal I_{\mathrm{pre}}(q,p)$；
- global ROI map $M_{q,p}(r)$ 表示这些 modes 的 spatial footprint；
- loop backbone $A_{q,p}^{\mathrm{backbone}}$ 表示这些 modes 在 ROI 空间中的 phase-ordered expression。

这个课题说明：

> BOLD pre-activation is not reverse neurovascular causality; it is a phenomenological modal signature of a global gating state that precedes and organizes local LFP subprocesses.

---

## 13. 两个课题合起来的 closed-loop interpretation

两个课题分别对应闭环的两半：

| 路径 | 对应课题 | 解释 |
|---|---|---|
| $u_{\mathrm{LFP}} \rightarrow x^{(g)}$ | Event-related LFP subprocess -> BOLD innovation | local subprocess perturbs / updates global BOLD modal state |
| $x^{(g)} \rightarrow u_{\mathrm{LFP}}$ | Non-causal IR / pre-activation | global gating state precedes and gates local LFP subprocess |

合起来就是：

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}} \rightarrow x^{(g)}.
$$

也就是：

> global order-parameter-like state gates local hippocampal subprocess, and local subprocesses in turn perturb/update global BOLD network state.

这个 closed-loop local-global self-organization 是两个课题共享的核心机制性解释。

---

## 14. Phenomenological 边界与写作措辞

为了避免 overclaim，需要明确：

- pre-activation modes 不是 causal modes；
- global ROI map 不是 causal feedback map；
- loop backbone 不是 anatomical feedback loop；
- 它们都是 phenomenological / operational Koopman modal signatures。

但它们也不是任意的 correlation pattern，因为它们是 conditioned on pair-specific pre-activation modes：

$$
\mathcal I_{\mathrm{pre}}(q,p).
$$

推荐表述：

> The pre-activation ROI map and loop backbone are phenomenological modal signatures derived from the subset of Koopman modes that contribute to the negative-lag component of a given LFP-BOLD cross-kernel. The ROI map describes where the corresponding global order-parameter-like state is expressed in BOLD space, whereas the loop backbone describes how this spatial footprint is phase-ordered across ROIs. We interpret the backbone as a candidate BOLD-observable representation of a latent feedback/slaving process, not as the feedback mechanism itself.

中文表述：

> pre-activation ROI map 和 loop backbone 是由解释某个 LFP-BOLD pair negative-lag cross-kernel 的 Koopman modes 提取出来的现象学模态表征。ROI map 描述这些 modes 对应的 global order-parameter-like state 在全脑 BOLD 空间中的表达位置；loop backbone 描述这个空间足迹在 ROI 之间的相位有序表达。因此，backbone 可以被解释为 latent feedback/slaving process 的候选 BOLD 可观测表征，而不是反馈机制本身。

---

## 15. 一句话总结

> 一个分布式 global latent state / order-parameter-like variable 先于 local LFP subprocess 出现，并 gate 未来 local LFP subprocess。local BOLD 因为同时读出 local 与 global latent components，所以可以在 LFP 之前显示 pre-activation。Koopman pre-activation modes 捕捉这种相位领先的 shared latent dynamics；global ROI map 是这些 modes 对应的 global gating state 的空间足迹；loop backbone 是这个空间足迹在 ROI 间的相位有序表达。另一个 event-related 课题则沿着相反方向说明 local LFP subprocess 可以作为 intrinsic perturbation 触发或预测 global BOLD modal innovations。两个课题合起来描述了 local-global spontaneous brain dynamics 的 closed-loop self-organization。
