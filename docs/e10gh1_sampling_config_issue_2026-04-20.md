# E10.gH1 Sampling And Config Issue

Date: 2026-04-20

## Summary

`E10.gH1` currently needs a manual config review before its MLP/postprocessing results are treated as final.

The issue is not a tiny floating-point difference. Raw `E10.gH1` sessions mix acquisition-derived time bases across data products:

- `blp.dx` is effectively consistent across sessions, about `0.001515 s`.
- `roiTs{1}.dx` has two values: `2 s` for sessions `1-14`, and `1 s` for sessions `15-25`.
- `Spkt.dx` has two values: `0.025 s` for sessions `1-16, 20`, and `0.005 s` for sessions `17-19, 21-25`.

The current config file should be checked manually:

`configs/cfg_E10gH1.m`

At the time this note was written, its included sessions were:

```matlab
cfg.sessions(1).session_id = 1:5;          % polar
cfg.sessions(2).session_id = [6:14, 21:25]; % spont
```

The `21:25` block is the main problem for spike correlation because it has `Spkt.dx = 0.005 s`, while the earlier included sessions have `Spkt.dx = 0.025 s`.

## Reference Config Evidence

The original reference file is:

`D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10gH1.m`

Relevant lines:

```matlab
% SPONTANEOUS(EXPS) - Poor neuro-data [17:19 21:25]
GRP.spont.exps = [6:14];
GRP.polar.exps = [1:5];
```

So the original session grouping appears to exclude `21:25` from the good spontaneous neuro-data block.

## Raw Sampling Audit

An audit CSV was written here:

`E:\DataPons_processed\e10gh1\postprocessing\e10gh1_raw_sampling_audit.csv`

Important rows:

```text
sessions 1-5:   roiTs dx=2s, roi_n=128, spkt dx=0.025s
sessions 6-14:  roiTs dx=2s, roi_n=300, spkt dx=0.025s
sessions 15-16: roiTs dx=1s, roi_n=600, spkt dx=0.025s
sessions 21-25: roiTs dx=1s, roi_n=600, spkt dx=0.005s
```

The audit included `15:16` because an earlier inspected config state included them. The current config should still be reviewed because it includes `21:25`.

## Processed Dictionary Evidence

The processed dictionary used by the current MLP outputs is:

`E:\DataPons_processed\e10gh1\reskoopnet_dictionary\e10gh1_low50_high250_g2_abs_single.mat`

Its `session_dx` values are all BLP-like, about `0.001515 s`. This means the EDMD/MLP time axis itself is based on the BLP sampling grid. The problematic mixed `0.025 s` versus `0.005 s` time base appears when aligning the full-time deconvolved residual to spike bins.

## Affected Outputs

Treat current `e10gh1` full-time spike-residual correlation outputs as not final:

`E:\DataPons_processed\e10gh1\postprocessing\mlp_fulltime_spike_residual_correlation`

Also treat current `e10gh1` MLP/postprocessing outputs as tied to the old session selection until the config is corrected and the upstream inputs are regenerated:

`E:\DataPons_processed\e10gh1\postprocessing\mlp_top_state_diversity_postprocessing`

`E:\autodl_results\e10gh1\mlp\outputs`

If the config changes the included sessions, rerunning only the final plotting/postprocessing stage is not enough. The dictionary/observables and MLP outputs need to correspond to the corrected session set.

## Manual Config Decision Needed

Do not edit this automatically. The config should be changed manually after deciding the intended `E10.gH1` session policy.

Conservative/reference-matching option:

```matlab
cfg.sessions(1).session_id = 1:5;   % polar
cfg.sessions(2).session_id = 6:14;  % spont
```

Alternative, only if `21:25` must be analyzed:

- Keep `21:25` as a separate group.
- Do not pool its spike correlation together with `1:14`.
- Report it separately because it has different `Spkt.dx` and was flagged as poor neuro-data in the reference config.

## Rerun Plan After Manual Config Fix

After `cfg_E10gH1.m` is manually corrected, rerun the `e10gh1` chain from the first stage that depends on session selection.

At minimum, regenerate:

1. `e10gh1` observables / ResKoopNet dictionary inputs.
2. `e10gh1` MLP EDMD outputs under `E:\autodl_results\e10gh1\mlp\outputs`.
3. `e10gh1` consensus state diversity windows if their inputs depend on the corrected session set.
4. `e10gh1` top-window postprocessing plots.
5. `e10gh1` full-time spike-residual correlation.

The other datasets and `F12m01` are not part of this issue.
