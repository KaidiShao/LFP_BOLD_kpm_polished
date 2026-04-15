# Configuration Establishment Summary (2026-04-14)

## Scope

This note records the current establishment status of the newly added dataset configuration files:

- `configs/cfg_E10aw1.m`
- `configs/cfg_E10bv1.m`
- `configs/cfg_E10eA1.m`
- `configs/cfg_E10fV1.m`
- `configs/cfg_E10gH1.m`
- `configs/cfg_E10gb1.m`

The main comparison sources used in this round were:

- `D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys`
- `D:\Onedrive\experimental_event_data\ripple_monkeys`

## Establishment Results

### `cfg_E10aw1.m`

- Current reference status: only found under `forNikos`
- Active channel map: `{'pl', 'hp', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', 'pl'}`
- Active region labels: `{'hp', 'pl'}`
- Active sessions:
  - `spont = 1:25`
- Current interpretation:
  - Single-source config
  - No `polar` session is currently defined

### `cfg_E10bv1.m`

- Current reference status: only found under `forNikos`
- Active channel map: `{'pl', 'hp', 'pl', 'pl', 'hp', 'hp', 'hp', 'pl', 'pl'}`
- Active region labels: `{'hp', 'pl'}`
- Active sessions:
  - `spont = 1:25`
- Current interpretation:
  - Single-source config
  - No `polar` session is currently defined

### `cfg_E10eA1.m`

- Current reference status:
  - found under `forNikos`
  - found under `ripple_monkeys`
- Reference mismatch:
  - `forNikos` channel map uses `pl/hp`
  - `ripple_monkeys` channel map uses `pl/nos/sr`
- Active channel map in cfg:
  - `{'pl', 'nos', 'nos', 'sr', 'sr', 'sr', 'nos', 'nos', 'pl'}`
- Active region labels in cfg:
  - `{'sr', 'pl'}`
- Active sessions in cfg:
  - `spont = 1:10`
  - `vstim = 11:30`
- Current interpretation:
  - The cfg now uses the more conservative Michel / `ripple_monkeys` channel labeling
  - The session split remains the conservative `1:10` and `11:30`
  - No `polar` session is currently defined

### `cfg_E10fV1.m`

- Current reference status: only found under `forNikos`
- Active channel map: `{'lgn', 'lgn', 'lgn', 'st', 'lgn', 'st', 'lgn', 'hp', 'hp', 'hp', 'hp', 'pl', 'hp', 'hp', 'pl'}`
- Active region labels: `{'hp', 'pl'}`
- Active sessions:
  - `polar = 1:5`
  - `spont = 6:25`
- Current interpretation:
  - Single-source config
  - Explicit `polar` block is available and should be treated as first-priority follow-up data

### `cfg_E10gH1.m`

- Current reference status: only found under `forNikos`
- Active channel map: `{'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'st', 'st', 'st', 'st', 'pl', 'pl', 'pl', 'hp', 'hp', 'hp', 'hp', 'hp'}`
- Active region labels: `{'hp', 'pl'}`
- Active sessions:
  - `polar = 1:5`
  - `spont = [6:16, 21:25]`
- Current interpretation:
  - Single-source config
  - Explicit `polar` block is available and should be treated as first-priority follow-up data

### `cfg_E10gb1.m`

- Current reference status:
  - found under `forNikos`
  - found under `ripple_monkeys`
- Reference mismatch:
  - `forNikos` channel map: `{'lgn', 'lgn', 'st', 'lgn', 'lgn', 'lgn', 'st', 'lgn', 'pl', 'hp', 'pl', 'hp', 'hp', 'hp', 'pl', 'hp'}`
  - `ripple_monkeys` channel map: `{'lgn', 'lgn', 'st', 'lgn', 'lgn', 'lgn', 'st', 'st', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', 'pl'}`
- Active channel map in cfg:
  - follows the `ripple_monkeys` version
- Active region labels in cfg:
  - `{'hp', 'pl'}`
- Active sessions in cfg:
  - `polar = 1:5`
  - `spont = [6:7, 9:13, 20:25]`
  - `vspont = [31:35]`, but currently excluded
- Current interpretation:
  - This cfg is intentionally kept as the validation baseline
  - It preserves the broader manually annotated spont subset used in prior runs
  - More conservative `ripple_monkeys` spont subset for later comparison:
    - `[6:7, 10:11, 20:24]`

## Working Conclusions

### 1. Current cfg roles are now separated more clearly

- `E10eA1` is now treated conservatively on channel labeling
- `E10gb1` is kept as a validation-baseline cfg rather than a strict conservative cfg
- `E10aw1`, `E10bv1`, `E10fV1`, and `E10gH1` remain single-source configs from `forNikos`

### 2. The next priority should be datasets with explicit `polar` sessions

Reason:

- `polar` is already explicitly encoded in the current cfg files
- these datasets are the cleanest next step for comparing behavior across the new pipeline
- they avoid mixing the question of session-type selection with the question of conservative relabeling

## Next-Step Priority List: datasets with `polar` sessions

### Priority 1

- `configs/cfg_E10fV1.m`
  - `polar = 1:5`
  - `spont = 6:25`
  - region labels: `{'hp', 'pl'}`

- `configs/cfg_E10gH1.m`
  - `polar = 1:5`
  - `spont = [6:16, 21:25]`
  - region labels: `{'hp', 'pl'}`

- `configs/cfg_E10gb1.m`
  - `polar = 1:5`
  - validation-baseline `spont = [6:7, 9:13, 20:25]`
  - conservative reference `spont = [6:7, 10:11, 20:24]`
  - region labels: `{'hp', 'pl'}`

### Not first priority for the `polar` pass

- `configs/cfg_E10aw1.m`
  - only `spont`

- `configs/cfg_E10bv1.m`
  - only `spont`

- `configs/cfg_E10eA1.m`
  - `spont + vstim`
  - no explicit `polar` block in the current cfg

## Local Data Availability Check On This Machine

At the time of this note, local raw-data folders under `E:\DataPons` were observed as:

- present:
  - `E10.aW1`
  - `E10.gb1`
- not currently found:
  - `E10.bv1`
  - `E10.eA1`
  - `E10.fV1`
  - `E10.gH1`

This means the conceptual next priority is still the `polar` datasets above, but on the current machine the immediately runnable one appears to be `E10gb1`.
