# E10.gH1 Sampling And Session Selection Note

Date: 2026-04-20

## Summary

`E10.gH1` should not use sessions `[6:16, 21:25]` as one pooled spontaneous block for spike-residual analyses.

The original session definition file:

`D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10gH1.m`

defines:

| group | source sessions | note |
| --- | ---: | --- |
| `polar` | `1:5` | visual stimulation |
| `spont` | `6:14` | spontaneous with fMRI |
| `vspont` | `31:35` | spontaneous block collected during visual-stim protocol |

The source file also comments:

```text
SPONTANEOUS(EXPS) - Poor neuro-data [17:19 21:25]
GRP.spont.exps = [6:14]
```

The previous local config included `[6:16, 21:25]` for spontaneous sessions. That mixed sessions that are outside the source `GRP.spont.exps`, including sessions `21:25`, which are explicitly flagged as poor neuro-data.

## Sampling Audit

Audit CSV:

`E:\DataPons_processed\e10gh1\postprocessing\e10gh1_raw_sampling_audit.csv`

For the sessions previously included by `cfg_E10gH1.m`:

| sessions | `blp.dx` | `roiTs{1}.dx` | `Spkt.dx` | comment |
| --- | ---: | ---: | ---: | --- |
| `1:5` | about `0.001515 s` | `2 s` | `0.025 s` | source `polar` |
| `6:14` | about `0.001515 s` | `2 s` | `0.025 s` | source `spont` |
| `15:16` | about `0.001515 s` | `1 s` | `0.025 s` | not in source `GRP.spont.exps` |
| `21:25` | about `0.001515 s` | `1 s` | `0.005 s` | source marks as poor neuro-data |

Additional check over sessions `1:25`:

| sessions | `roiTs{1}.dx` | `Spkt.dx` |
| --- | ---: | ---: |
| `1:14` | `2 s` | `0.025 s` |
| `15:16` | `1 s` | `0.025 s` |
| `17:19` | `1 s` | `0.005 s` |
| `20` | `1 s` | `0.025 s` |
| `21:25` | `1 s` | `0.005 s` |

`blp.dx` is effectively uniform across these sessions; the large mismatch is in `roiTs` TR and `spkt` bin size.

## Consequence

The prior full-time `e10gh1` spike-residual correlation result should not be treated as final, because pooled spike correlation mixed `Spkt.dx = 0.025 s` and `Spkt.dx = 0.005 s` sessions. The `0.005 s` sessions contribute five times as many spike bins per second, which changes pooled weighting and time smoothing.

The MLP/EDMD dictionary metadata itself uses the BLP time base, so the immediate issue is not BLP chunk alignment. The issue is session inclusion and pooled spike comparison.

## Config Decision

`cfg_E10gH1.m` should use the source groups:

```matlab
polar = 1:5
spont = 6:14
```

Sessions outside those source groups should be excluded from the main config. Sessions `21:25` should not be pooled with `1:14` for spike-residual correlation. If they are analyzed later, they should be handled as a separate explicit group with their own sampling interpretation.

## Rerun Implication

Changing `cfg_E10gH1.m` changes the intended session span for `e10gh1`. Existing processed `e10gh1` artifacts generated with the old config may no longer match the new config, including:

- raw BLP cache
- spectrograms
- event detection and density
- consensus states
- consensus-state diversity windows
- reskoopnet dictionaries
- MLP/EDMD outputs
- postprocessing and spike-residual correlation outputs

For final `e10gh1` results, regenerate the upstream observable and MLP/EDMD outputs with the corrected config before re-running the top-window postprocessing pipeline.
