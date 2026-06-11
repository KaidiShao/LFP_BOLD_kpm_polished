# A unified framework for nonlinear local-global latent dynamics, Koopman modal signatures, and two LFP-BOLD projects

**Version: 2026-06-10**

This note summarizes the current mathematical framework connecting two related projects:

1. **Event-related LFP subprocess -> global BOLD pattern / BOLD eigenfunction residual**: how local LFP subprocesses can act as intrinsic perturbations that trigger or predict global BOLD modal innovations.
2. **Non-causal LFP-BOLD impulse responses / pre-activation**: why local BOLD can show apparent pre-activation relative to local LFP, and how this can be interpreted using closed-loop local-global latent dynamics and a Koopman modal representation.

The core idea is:

> LFP, local BOLD, and global BOLD should not be treated as a simple feed-forward input-output chain. They are better understood as different observable readouts of a nonlinear local-global latent system. Koopman eigenfunctions transform this complex nonlinear system into interpretable modal variables, making the global gating state, pre-activation modes, global ROI maps, and loop backbones operational phenomenological modal signatures.

---

## 1. Underlying model: a nonlinear local-global latent system

Let the latent state be decomposed into local and global components:

$$
x_t =
\begin{bmatrix}
x_t^{(l)} \\
x_t^{(g)}
\end{bmatrix}.
$$

Here:

- $x_t^{(l)}$ is a **local latent state**, such as hippocampal local neural state, local vascular state, or local circuit excitability.
- $x_t^{(g)}$ is a **global latent state / internal environment / global order-parameter-like state**, such as global brain state, neuromodulatory state, brainstem-thalamic-cortical state, global vascular state, or arousal state.
- $u_t$ is a **local LFP subprocess / local neural action**.

A generic nonlinear model can be written as

$$
x_{t+1}^{(l)} = F_l\!\left(x_t^{(l)}, x_t^{(g)}, u_t\right),
$$

$$
x_{t+1}^{(g)} = F_g\!\left(x_t^{(g)}, x_t^{(l)}, u_t\right),
$$

or more compactly:

$$
x_{t+1}=F(x_t,u_t).
$$

The key point is not that LFP is an external stimulus. Rather, $u_t$ is treated as an action-like subprocess of the local subsystem. It can be written as an input-like variable for modeling convenience, but in the closed-loop view it is generated internally by the system.

---

## 2. BOLD is a nonlinear readout, not the mechanism itself

Local BOLD and global BOLD are observable readouts of the latent state:

$$
y_t^{(lB)} = h_l\!\left(x_t^{(l)}, x_t^{(g)}\right),
$$

$$
y_t^{(gB)} = h_g\!\left(x_t^{(l)}, x_t^{(g)}\right).
$$

Therefore, local BOLD is not simply the output of a local HRF driven by local LFP. It may contain both a local and a global component:

$$
y_t^{(lB)} = \text{local component} + \text{global component}.
$$

This is crucial: if local BOLD reads out a global latent state before that state gates a future local LFP subprocess, local BOLD can show apparent pre-activation relative to LFP.

---

## 3. Global-to-local path: global state gates local LFP subprocesses

The central explanation for pre-activation is:

$$
x_t^{(g)} \rightarrow u_{\mathrm{LFP},t+\Delta}.
$$

More explicitly:

$$
x_t^{(g)} \rightarrow x_{t+\Delta}^{(l)} \rightarrow u_{\mathrm{LFP},t+\Delta}.
$$

This means:

> A global latent state / order-parameter-like variable appears first, then modulates the local hippocampal condition, making a local LFP subprocess more likely to occur.

To form a complete closed loop, there may also be a local-to-global update:

$$
u_{\mathrm{LFP},t} \rightarrow x_{t+\Delta}^{(g)}.
$$

The overall loop can therefore be summarized as

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}} \rightarrow x^{(g)}.
$$

However, for explaining **pre-activation**, the most important part is the first half:

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}}.
$$

---

## 4. Why local BOLD-LFP can produce an apparent non-causal IR

The standard feed-forward HRF model assumes:

$$
y_{\mathrm{localBOLD}}(t)
=
\sum_{\tau\ge 0} h(\tau) u_{\mathrm{LFP}}(t-\tau)+\epsilon(t).
$$

This treats LFP as an exogenous input and BOLD as a passive output. Under this assumption, if the estimated impulse response has negative-lag pre-activation, it appears non-causal.

In the local-global latent model, a global state $x_t^{(g)}$ can simultaneously:

1. be read out early by local BOLD:

$$
x_t^{(g)} \rightarrow y_t^{(lB)},
$$

2. gate a future local LFP subprocess:

$$
x_t^{(g)} \rightarrow u_{\mathrm{LFP},t+\Delta}.
$$

Thus the observations may show

$$
y_{\mathrm{localBOLD},t}
\quad \text{preceding} \quad
u_{\mathrm{LFP},t+\Delta}.
$$

If one still applies a feed-forward LFP-to-BOLD deconvolution model, the resulting cross-kernel can show negative-lag pre-activation. This does **not** imply

$$
y_{\mathrm{localBOLD}} \rightarrow u_{\mathrm{LFP}}.
$$

Instead, it reflects the shared influence of a latent global state:

$$
x^{(g)} \rightarrow
\begin{cases}
y_{\mathrm{localBOLD}},\\
u_{\mathrm{LFP}}.
\end{cases}
$$

The apparent non-causality therefore arises when a closed-loop endogenous system is interpreted using an open-loop feed-forward model.

---

## 5. From nonlinear system to augmented autonomous system

If $u_t$ is generated internally, it can be written as

$$
u_t = \pi\!\left(x_t^{(g)}, x_t^{(l)}, u_{t-1},u_{t-2},\ldots\right)+\epsilon_t.
$$

By augmenting the state with $u_t$ and its history,

$$
s_t = \begin{bmatrix}
x_t \\
u_t \\
u_{t-1} \\
\vdots
\end{bmatrix},
$$

the closed-loop system can be written in autonomous form:

$$
s_{t+1}=F_{\mathrm{cl}}(s_t).
$$

The readouts can be written as

$$
u_t = C_u s_t,
$$

$$
y_t^{(lB)} = C_l s_t,
$$

$$
y_t^{(gB)} = C_g s_t.
$$

This step is important because it shows that the local LFP subprocess, local BOLD, and global BOLD can all be treated as different readouts of the same autonomous closed-loop system.

---

## 6. Koopman modal representation

For the augmented autonomous system

$$
s_{t+1}=F_{\mathrm{cl}}(s_t),
$$

let $a_t$ denote Koopman eigenfunction coordinates. In modal coordinates:

$$
a_{t+1}=\Lambda a_t.
$$

The observables can be approximated as linear readouts:

$$
u_{\mathrm{LFP},t}=C_u a_t,
$$

$$
y_{\mathrm{localBOLD},t}=C_l a_t,
$$

$$
y_{\mathrm{globalBOLD},t}=C_g a_t.
$$

The value of the Koopman representation is that the original latent dynamics and observation functions may be nonlinear, but in lifted/modal space the dominant modes, phase relationships, and spatial readouts can be analyzed linearly or approximately linearly.

---

## 7. Definition of pre-activation modes

For an LFP-BOLD anchor pair $(p,q)$:

- $p$ denotes a local LFP band/subprocess, such as PGO, theta, or ripple.
- $q$ denotes a local BOLD ROI or local BOLD readout.

If the empirical cross-kernel $h_{q,p}(\tau)$ shows pre-activation at $\tau<0$, we can identify Koopman modes that contribute to this negative-lag component.

Define:

$$
\mathcal I_{\mathrm{pre}}(q,p)
=
\{i: \text{mode } i \text{ contributes to the negative-lag component of pair }(q,p)\}.
$$

These are the **pre-activation modes** for the anchor pair.

Biological interpretation:

> These modes most naturally correspond to a global latent gating state / order-parameter-like variable. This global state is read out by BOLD before it gates future local LFP subprocesses.

Important boundary:

> Pre-activation modes are not causal modes. They are pair-conditioned modal signatures.

---

## 8. Koopman phase-lead explanation

For mode $i$, let the continuous-time eigenvalue be

$$
\nu_i = \alpha_i + \mathrm{i}\omega_i.
$$

The complex readout weights for LFP and local BOLD can be written as

$$
C_u(i)=|C_u(i)|e^{\mathrm{i}\phi_{u,i}},
$$

$$
C_l(i)=|C_l(i)|e^{\mathrm{i}\phi_{l,i}}.
$$

If the BOLD readout phase leads the LFP readout phase, this mode contributes to apparent pre-activation in the LFP-BOLD cross-kernel. A phase-derived time difference can be written as

$$
\Delta t_{lB-u,i}
=
\frac{\operatorname{wrap}(\phi_{l,i}-\phi_{u,i})}{\omega_i},
$$

with the exact sign depending on the phase convention. Conceptually:

$$
\text{BOLD phase leads LFP phase}
\quad \Rightarrow \quad
\text{negative-lag contribution in empirical IR}.
$$

---

## 9. Meaning of the global ROI map

For the pre-activation modes of an anchor pair $(q,p)$, define the global ROI map as

$$
M_{q,p}(r)
=
\sum_{i\in\mathcal I_{\mathrm{pre}}(q,p)} s_i |C_g(r,i)|,
$$

where:

- $r$ is a global BOLD ROI.
- $C_g(r,i)$ is the global BOLD readout weight of mode $i$ at ROI $r$.
- $s_i$ is the strength with which mode $i$ contributes to the pair-specific pre-activation, for example negative-lag area, mode amplitude, or pairwise mode weight.

Interpretation:

> The global ROI map is the spatial footprint of the pre-activation modes in whole-brain BOLD ROI space.

In the closed-loop / slaving-principle interpretation, it further corresponds to:

> the BOLD-observable footprint of the global latent state / order-parameter-like variable that gates the local LFP subprocess.

It is **not** a causal map of BOLD ROIs driving LFP. It is a map of which ROIs jointly express the global gating state.

---

## 10. Meaning of the loop backbone

For the same pre-activation modes, each ROI readout has both amplitude and phase:

$$
C_g(r,i)=|C_g(r,i)|e^{\mathrm{i}\phi_{r,i}}.
$$

The global ROI map uses the amplitude:

$$
|C_g(r,i)|,
$$

whereas the loop backbone uses the phase:

$$
\phi_{r,i}.
$$

Thus:

> The global ROI map describes where the global state is expressed. The loop backbone describes how that global state is phase-expressed across ROIs.

For each mode, one may define a phase-based adjacency:

$$
A_i(r_1,r_2)=1
\quad \text{if ROI } r_1 \text{ phase-leads ROI } r_2
\quad \text{within mode } i.
$$

Aggregating over pre-activation modes gives

$$
A_{q,p}^{\mathrm{backbone}}(r_1,r_2)
=
\sum_{i\in\mathcal I_{\mathrm{pre}}(q,p)} s_i A_i(r_1,r_2).
$$

The loop backbone is interpreted as:

> the phase-ordered expression, or phase portrait, of the global gating state in ROI space.

It is a candidate biological representation / phenomenological modal signature. It is not the anatomical feedback loop itself, and it does not by itself establish interventional causality.

---

## 11. Relation to Haken's slaving principle

This framework can be related to Haken's slaving principle:

$$
\text{distributed brain regions}
\Rightarrow
\text{global order parameter}
\Rightarrow
\text{local fast subprocess}.
$$

In the current setting:

$$
\text{distributed ROI-level activity}
\Rightarrow
x^{(g)}
\Rightarrow
u_{\mathrm{LFP}}.
$$

Here:

- $x^{(g)}$ is a global order-parameter-like state / slaving variable.
- The global ROI map is the spatial BOLD footprint of $x^{(g)}$.
- The loop backbone is the phase portrait of $x^{(g)}$ across ROIs.
- The local LFP subprocess is the fast local process gated or enslaved by $x^{(g)}$.

This also forms a kind of circular causality:

$$
\text{distributed global state}
\rightarrow
\text{local LFP subprocess}
\rightarrow
\text{global state update}.
$$

---

## 12. Relationship to the two projects

### 12.1 Project 1: Event-related LFP subprocess -> global BOLD pattern / BOLD eigenfunction residual

This project corresponds to the local-to-global path:

$$
u_{\mathrm{LFP}} \rightarrow x^{(g)} \rightarrow y_{\mathrm{globalBOLD}}.
$$

Operationally, one first builds a BOLD-only Koopman model:

$$
b_{k+1}=\Lambda_B b_k + \eta_k^B,
$$

where the BOLD modal innovation / BOLD eigenfunction residual is

$$
\eta_k^B=b_{k+1}-\Lambda_B b_k.
$$

The high-rate LFP subprocesses are then converted to the BOLD time scale as density variables:

$$
D_r^{\mathrm{LFP}}[k]
=
\text{density / occupancy of LFP subprocess } r
\text{ in BOLD window } k.
$$

The core predictive model is

$$
\eta_{m,k}^B
=
\sum_r\sum_{\ell}
\beta_{m,r,\ell}D_r^{\mathrm{LFP}}[k-\ell]
+
\epsilon_{m,k}.
$$

If $D_r^{\mathrm{LFP}}$ predicts future $\eta^B$, then:

> the local LFP subprocess acts as an intrinsic perturbation that triggers or predicts global BOLD modal innovations.

Important boundary:

> The BOLD residual itself is not the LFP subprocess activation. It is an input-like deviation of BOLD modal dynamics. It supports the LFP-subprocess-as-intrinsic-perturbation interpretation only if it is predicted by LFP subprocess density.

---

### 12.2 Project 2: Non-causal LFP-BOLD IR / pre-activation

This project corresponds to the global-to-local path:

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}}.
$$

The proposed sequence is:

1. A global latent state $x^{(g)}$ appears first.
2. Local BOLD reads out a local/global mixture that includes $x^{(g)}$.
3. The same $x^{(g)}$ gates a future local LFP subprocess.
4. Therefore a feed-forward LFP-to-local-BOLD IR can show apparent pre-activation.

In the Koopman analysis:

- The anchor pair $(q,p)$ selects the pre-activation modes $\mathcal I_{\mathrm{pre}}(q,p)$.
- The global ROI map $M_{q,p}(r)$ is the spatial footprint of these modes.
- The loop backbone $A_{q,p}^{\mathrm{backbone}}$ is their phase-ordered expression across ROIs.

This project argues:

> BOLD pre-activation is not reverse neurovascular causality. It is a phenomenological modal signature of a global gating state that precedes and organizes local LFP subprocesses.

---

## 13. Integrated closed-loop interpretation

The two projects correspond to two halves of the same closed loop:

| Path | Project | Interpretation |
|---|---|---|
| $u_{\mathrm{LFP}} \rightarrow x^{(g)}$ | Event-related LFP subprocess -> BOLD innovation | local subprocess perturbs / updates global BOLD modal state |
| $x^{(g)} \rightarrow u_{\mathrm{LFP}}$ | Non-causal IR / pre-activation | global gating state precedes and gates local LFP subprocess |

Together:

$$
x^{(g)} \rightarrow u_{\mathrm{LFP}} \rightarrow x^{(g)}.
$$

In words:

> A global order-parameter-like state gates local hippocampal subprocesses, and local subprocesses in turn perturb or update the global BOLD network state.

This closed-loop local-global self-organization is the shared mechanistic interpretation of the two projects.

---

## 14. Interpretation boundaries

The framework should be stated carefully:

- A **pre-activation mode** is not a causal mode; it is a pair-conditioned modal signature.
- A **global ROI map** is not a causal feedback map; it is the spatial footprint of pre-activation modes.
- A **loop backbone** is not an anatomical feedback loop; it is the phase-ordered expression of a global gating state in ROI space.
- The whole interpretation is phenomenological but mechanistically informative.

A concise manuscript-friendly formulation is:

> Pre-activation modes identify a global order-parameter-like state associated with an LFP-BOLD anchor pair. The global ROI map describes where this state is expressed in BOLD space, while the loop backbone describes how that state is phase-expressed across ROIs. These objects are phenomenological modal signatures of a latent feedback/slaving process, rather than direct anatomical or interventional feedback pathways.
