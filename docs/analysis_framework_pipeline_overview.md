# Analysis Framework Pipeline Overview

This note is a meeting-ready draft of the full analysis framework. The goal is
to show how the BLP and BOLD analyses fit into one ResKoopNet / Koopman
dynamics workflow, and to separate the scientific questions from the many
implementation conditions.

## One-Slide Pipeline Diagram

```mermaid
flowchart LR
    D["Datasets<br/>E10gb1, E10fV1, E10gH1, F12m01"]

    D --> BLP_RAW["Raw BLP sessions"]
    D --> BOLD_RAW["Raw BOLD ROI time series"]

    subgraph PRE["Preprocessing"]
        direction TB

        subgraph BLP_PRE["BLP preprocessing"]
            direction TB
            BLP_SPEC["Region-mean spectrograms<br/>abs and complex"]
            BLP_OBS["ResKoopNet dictionaries<br/>abs, complex_split"]
            BLP_EVENTS["Bandpass events"]
            BLP_DENS["Event density<br/>2 s bins"]
            BLP_STATES["Consensus states"]
            BLP_DIV["State-diversity windows"]
            BLP_SPEC --> BLP_OBS
            BLP_EVENTS --> BLP_DENS
            BLP_EVENTS --> BLP_STATES
            BLP_STATES --> BLP_DIV
        end

        subgraph BOLD_PRE["BOLD preprocessing"]
            direction TB
            BOLD_OBS["BOLD observables<br/>roi_mean, eleHP, HP, svd,<br/>HP_svd100, global_svd100,<br/>slow_band_power"]
            BOLD_QC["Pre-ResKoopNet QC"]
            BOLD_OBS --> BOLD_QC
        end
    end

    BLP_RAW --> BLP_SPEC
    BLP_RAW --> BLP_EVENTS
    BOLD_RAW --> BOLD_OBS

    subgraph RK["ResKoopNet / Koopman learning"]
        direction TB

        BLP_MLP["BLP MLP ResKoopNet<br/>observable: abs or complex_split<br/>residual: projected_kv or projected_vlambda"]
        BLP_CNN["BLP CNN ResKoopNet<br/>optional comparison branch"]
        BOLD_MLP["BOLD MLP ResKoopNet<br/>7 BOLD observable modes<br/>x projected_kv / projected_vlambda"]
    end

    BLP_OBS --> BLP_MLP
    BLP_SPEC --> BLP_CNN
    BOLD_OBS --> BOLD_MLP

    subgraph POST["Postprocessing and interpretation"]
        direction TB

        subgraph BLP_POST["BLP postprocessing"]
            direction TB
            BLP_TOP["Top state-diversity windows"]
            BLP_DR["Eigenfunction dimension reduction<br/>no UMAP branch"]
            BLP_SPIKE["Spike-residual correlation"]
        end

        subgraph BOLD_POST["BOLD postprocessing"]
            direction TB
            BOLD_EFUN["BOLD eigenfunction postprocessing"]
            BOLD_DECONV["Deconvolved eigenfunctions"]
            BOLD_TS["Timescale diagnostics"]
            BOLD_ACT["Mode activation maps"]
        end

        subgraph XMODAL["BLP-BOLD coupling"]
            direction TB
            XCORR["Lagged cross-correlation<br/>BOLD eigenfunctions/deconv<br/>vs BLP density signals"]
            LAG["Lag interpretation<br/>positive lag = BLP density leads BOLD"]
        end
    end

    BLP_MLP --> BLP_TOP
    BLP_MLP --> BLP_DR
    BLP_MLP --> BLP_SPIKE
    BLP_CNN -. optional .-> BLP_TOP

    BOLD_MLP --> BOLD_EFUN
    BOLD_EFUN --> BOLD_DECONV
    BOLD_EFUN --> BOLD_TS
    BOLD_EFUN --> BOLD_ACT

    BLP_DENS --> XCORR
    BLP_DIV --> XCORR
    BLP_DR --> XCORR
    BOLD_DECONV --> XCORR
    BOLD_EFUN --> XCORR
    XCORR --> LAG

    subgraph QUESTIONS["Scientific readouts"]
        direction TB
        Q1["Which neural states/events organize BLP dynamics?"]
        Q2["Which Koopman modes align with spiking?"]
        Q3["Which BOLD spatial modes are interpretable?"]
        Q4["Do BLP dynamics lead or follow BOLD dynamics?"]
        Q5["Which model/observable condition is most robust?"]
    end

    BLP_TOP --> Q1
    BLP_SPIKE --> Q2
    BOLD_ACT --> Q3
    LAG --> Q4
    BLP_DR --> Q5
    BOLD_TS --> Q5
```

## Compact Version For A Slide

Use this if the full diagram is too dense for one slide.

```mermaid
flowchart LR
    A["Raw data<br/>BLP + BOLD"] --> B["Preprocessing<br/>BLP spectrograms, events, states<br/>BOLD observables and QC"]
    B --> C["ResKoopNet<br/>BLP MLP/CNN conditions<br/>BOLD MLP conditions"]
    C --> D["Within-modality interpretation<br/>BLP top windows, DR, spike corr<br/>BOLD postprocess, deconv, activation maps"]
    D --> E["Cross-modal analysis<br/>BOLD eigenfunctions vs BLP densities<br/>session-aware lagged correlation"]
    E --> F["Scientific decisions<br/>primary condition, robustness checks,<br/>biological interpretation"]
```

## Condition Design To Explain

The conditions should be framed as controlled perturbations of the modeling
choice, not as a brute-force list.

| Axis | Conditions | Interpretation |
|---|---|---|
| Model family | MLP, CNN | MLP is the primary vector-observable model. CNN is an optional branch to test whether local spectrogram structure helps. |
| BLP observable | abs, complex_split | abs focuses on power magnitude. complex_split keeps real and imaginary spectrogram information. |
| Residual form | projected_kv, projected_vlambda | Two residual definitions for Koopman-consistent reconstruction and dynamics. |
| BOLD observable | roi_mean, eleHP, HP, svd, HP_svd100, global_svd100, slow_band_power | Tests spatial scale, targeted regions, and dimensionality reduction choices. |
| Postprocessing target | state diversity, eigenfunction DR, spike correlation, BOLD activation, BLP-BOLD correlation | Tests whether learned dynamics map onto state structure, spiking, spatial BOLD maps, and cross-modal lagged coupling. |

## Suggested Figure Set For The Discussion

Keep the meeting centered on the framework. Show representative figures rather
than every condition.

| Figure | Purpose |
|---|---|
| Pipeline overview | Establish the full analysis framework. |
| Condition matrix | Explain why each condition exists and which are primary vs robustness branches. |
| BLP top state-diversity window | Show how Koopman dynamics align with interpretable BLP state structure. |
| Eigenfunction DR/state-space plot | Show whether states separate in learned coordinates. |
| Spike-residual correlation heatmap or top-bar plot | Show whether residual dynamics relate to spiking. |
| BOLD activation map | Show whether BOLD Koopman modes have spatial interpretation. |
| BLP-BOLD lagged cross-correlation plot | Show directionality and timing of cross-modal coupling. |

## Proposed Meeting Questions

- Which BLP condition should be the primary analysis: abs or complex_split?
- Which residual form should be primary: projected_kv or projected_vlambda?
- Should CNN remain a supplemental comparison branch?
- Which BOLD observable modes are biologically most meaningful?
- Should BLP-BOLD coupling use event density, state diversity, eigenfunction density, or all as robustness checks?
- How should we define and report lag direction in the BLP-BOLD correlation analysis?
