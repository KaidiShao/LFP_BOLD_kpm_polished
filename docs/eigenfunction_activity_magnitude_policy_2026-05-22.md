# Eigenfunction Activity Magnitude Policy

Date: 2026-05-22

## Core Conclusion

For the current mainline analysis, eigenfunction activity should be defined by
magnitude/envelope, not by unanchored signed polarity.

In practical terms:

- Use `abs(phi)`, envelope, rectified activity, RMS, or power-like summaries to
  decide whether an eigenfunction/component is active.
- Do not interpret the positive or negative direction of `real(phi)` as a
  biological direction unless a dedicated sign/phase anchoring step has been
  applied and saved in metadata.
- Treat signed real/imag outputs as diagnostic views by default.

This makes the current scientific question:

```text
Does this eigenfunction/component activity magnitude increase near theta/ripple
events, and is that activity selective?
```

not:

```text
Does the positive direction of this eigenfunction mean theta/ripple activation?
```

## Why This Is Necessary

EDMD eigenfunctions have arbitrary scale. For real eigenfunctions, this includes
an arbitrary sign: `phi` and `-phi` are equivalent. For complex eigenfunctions,
the ambiguity is larger: `phi * exp(i*theta)` is an equivalent phase rotation.

Therefore, without explicit anchoring:

- The sign of `real(phi)` is not stable.
- The phase of complex `phi` is not biologically interpretable.
- A positive threshold on an unanchored signed signal can flip if the
  eigenfunction sign/phase flips.

## Important Distinction: Observable Branch vs Eigenfunction Transform

`abs_projected_vlambda` is an observable/source-branch label. It means the P4
input observable was constructed from an absolute-value BLP feature.

It does **not** guarantee that downstream EDMD eigenfunctions are real, positive,
or polarity-stable.

These are different concepts:

```text
abs_projected_vlambda        = observable/input branch label
abs(phi)                     = downstream eigenfunction magnitude
abs(component) / envelope    = downstream dimred component activity magnitude
```

The first does not solve polarity. The second and third avoid the polarity
question by defining activity as magnitude.

## Complex Eigenfunctions

For complex eigenfunctions, the mainline should use magnitude/envelope:

```text
activity(t) = abs(phi(t))
```

or a smoothed envelope/power-like feature derived from it.

This is not a compromise; it is the natural phase-invariant activity definition
for complex modes. Signed `real(phi)` and `imag(phi)` can be useful diagnostics,
but they should not drive current event selectivity conclusions.

## Real Eigenfunctions

For real eigenfunctions, the mainline can still use magnitude activity:

```text
activity(t) = abs(real(phi(t)))
```

or rectified/envelope/RMS activity.

This answers whether the eigenfunction is active around an event. It does not
answer whether positive or negative excursions have different biological
meaning. If that question becomes important, it should be handled as a separate
polarity-anchored analysis.

## Dimension-Reduced Components

Dimred component signs are also generally arbitrary for methods such as SVD,
MDS, and UMAP. NMF is more naturally nonnegative, but preprocessing and
component construction can still make the exact interpretation method-specific.

Mainline dimred activity should therefore use:

```text
abs(component(t))
```

or an envelope/energy feature before thresholding or event detection.

This is especially important for P5 dimred thresholded density. A threshold on
raw signed component values can detect opposite events after an arbitrary sign
flip. A threshold on magnitude/envelope detects activity independent of sign.

## Residual Correlations

Residual signals do not need polarity anchoring for the current coupling use
case.

For residual-vs-density or residual-vs-BLP comparisons:

- Rank coupling strength by `abs(correlation)`.
- Keep signed correlation for diagnostics.
- Do not interpret the sign as biological polarity unless both signals have
  explicit, compatible anchoring.

Thus `r = -1` and `r = +1` both indicate strong coupling in the current
strength-consistency analysis.

## Pipeline Implications

### P5

Mainline P5 event activity/selectivity should use magnitude/envelope features:

- raw eigenfunction density: `abs(phi)` or envelope-derived activity;
- dimred component density: `abs(component)` or envelope-derived activity;
- peak statistics: event-vs-baseline tests on magnitude/envelope activity;
- selectivity: compare event-family activity magnitudes, not signed polarity.

Signed real/imag outputs remain diagnostic unless explicitly anchored.

For raw eigenfunction density, the mainline density columns should also preserve
mode identity and timescale provenance:

```text
raw_efun_index
source_mode_index
lambda_discrete
lambda_continuous_real
lambda_continuous_imag
timescale_sec
frequency_hz
activity_transform
activity_window_policy
```

This is required because P8/P10 raw-density xcorr should be able to answer
whether strong coupling is concentrated in fast eigenfunction modes, slow modes,
or a broad mix.

### P5 Peak Statistics

Peak statistics should be computed on activity-magnitude signals, not on raw
signed eigenfunctions/components.

Recommended flow:

```text
raw eigenfunction/component
-> activity transform: magnitude or rms_envelope
-> event-window and baseline-window summary
-> event-vs-baseline statistics
-> selectivity and cross-session consistency
```

The statistic should answer whether activity magnitude increases near an event:

```text
peak(activity) near event
mean(activity) near event
event activity - baseline activity
```

It should not use an unanchored signed positive threshold such as:

```text
peak(real(phi))
mean(real(phi))
positive crossing of signed component
```

unless that signal has an explicit polarity/phase anchor.

The current `maxabs` peak-state branch is direction-safe in spirit, but future
mainline naming should make the activity transform explicit. Prefer names such
as:

```text
pipeline5_eigenfunction_activity_by_state
pipeline5_eigenfunction_peaks_by_state_activity_abs
pipeline5_eigenfunction_peaks_by_state_activity_rms_timescale
```

Every peak-statistics output should save the activity metadata described below,
including the transform, window policy, and source eigenfunction/component MAT.

### P5 Dimred Component Process Labels

P5 should expose a derived label table that describes what each dimred LFP
component appears to represent before that component is used as a P8/P10 density
source.

The label table should reuse the existing event-family selectivity logic as the
starting point, but run it on activity-magnitude peak statistics:

```text
component(t)
-> activity transform: abs(component) or rms_envelope(component)
-> event-vs-baseline peak statistics
-> theta/ripple/gamma family selectivity
-> per-component process label
```

Mainline label outputs should therefore be based on `abs(component)` or
envelope/RMS activity, not signed component values and not unanchored polarity.
Current `maxabs` selectivity summaries may be used as prototypes while this
branch is being implemented, but they should be marked transitional.

Recommended output columns:

```text
dataset
condition
method
k
method_tag
component_idx
density_name_for_p8_p10
density_index
theta_active
theta_selective
ripple_active
ripple_selective
gamma_active
gamma_selective
active_family_set
selective_family_set
primary_process_label
label_confidence
source_peak_stats_file
lfp_activity_transform
lfp_activity_window_policy
```

Recommended initial label vocabulary:

```text
theta_selective_similar
theta_selective_unequal
ripple_selective_similar
ripple_selective_unequal
theta_ripple_joint
mixed_theta_ripple
pan_event
partial_or_inactive
nonselective
unlabeled
```

These labels are interpretation metadata for P8/P10. They should be joined onto
xcorr top-hit and consistency tables, but they should not replace the numeric
coupling strength, lag, or ROI/profile consistency scores.

### P5 Consensus Trajectory

Consensus-state trajectory QC should also have an activity-space version.

The current signed dimred trajectory answers how samples move in arbitrary
component coordinates. That is useful as a diagnostic, but the axes can flip for
methods with arbitrary sign or phase. The mainline scientific trajectory should
instead show activity magnitude:

```text
dimred temporal components
-> activity transform: abs(component) or rms_envelope(component)
-> optional normalization
-> trajectory in activity-component space
-> color by consensus state or event family
```

Recommended figure families:

```text
activity_trajectory_first3
activity_trajectory_selected
```

`activity_trajectory_first3` is a simple diagnostic based on the first three
activity components.

`activity_trajectory_selected` is the preferred scientific view. It should use
peak statistics/selectivity to choose theta/ripple-active components and then
plot the corresponding activity dimensions.

Required trajectory metadata:

```text
trajectory_space = activity_magnitude
activity_transform
activity_window_policy
selected_component_idx
selection_event_family
selection_score
source_peak_statistics
source_dimred_result
```

Signed first-three trajectory plots should remain available only as
`first3_signed_diagnostic` or `unanchored_signed_diagnostic`.

## Envelope Implementation Policy

The activity transform should be implemented as a reusable preprocessing step
before thresholding, density estimation, peak statistics, or cross-modal
coupling.

Recommended function shape:

```matlab
[A, meta] = compute_eigenfunction_activity(X, dt, params, mode_meta)
```

where `X` is `[time x mode_or_component]` and `A` is the corresponding
activity-magnitude signal.

Supported transforms:

```text
magnitude       A(t) = abs(X(t))
abs_smooth      A(t) = movmean(abs(X(t)), window)
rms_envelope    A(t) = sqrt(movmean(abs(X(t)).^2, window))
hilbert         A(t) = abs(hilbert(real(X(t))))  [optional/diagnostic]
```

Default should remain simple and robust:

```text
activity_transform = magnitude
```

For envelope-style analyses, the preferred mainline option is:

```text
activity_transform = rms_envelope
```

`rms_envelope` is preferable to a Hilbert envelope as a default because
eigenfunctions and dimred components are not guaranteed to be narrow-band
oscillatory signals.

### Session-Aware Smoothing

Any moving-window envelope must be computed within session boundaries. Do not
let `movmean`, RMS, or Hilbert edge effects cross session boundaries.

Implementation rule:

```text
for each session:
    compute activity/envelope on that session only
concatenate session-wise activity outputs
```

This prevents artificial activity at session joins.

### Eigenvalue-Based Envelope Windows

A fixed envelope window such as `0.1 sec` or `0.2 sec` can be useful as a
fallback, but a more principled option is to choose the smoothing window from
the eigenvalue-derived timescale.

The key point is that EDMD eigenvalues are often discrete-time eigenvalues. For
decay-rate interpretation, convert to continuous time first:

```text
lambda_c = log(lambda_d) / dt
decay_rate = -real(lambda_c)
tau = 1 / decay_rate
```

Equivalently, when only the magnitude is needed:

```text
tau = -dt / log(abs(lambda_d))
```

for stable modes with `0 < abs(lambda_d) < 1`.

Then choose the RMS/envelope window as:

```text
window_sec_j = clamp(alpha * tau_j, min_window_sec, max_window_sec)
```

Recommended initial settings:

```text
alpha = 0.25 to 0.5
min_window_sec = 0.03
max_window_sec = 1.0 to 2.0
```

The clamp is important because eigenvalues close to the unit circle can imply
very long timescales, which would oversmooth event-scale activity.

For edge cases:

- if `abs(lambda_d) <= 0`, use `min_window_sec` or mark the timescale invalid;
- if `abs(lambda_d) >= 1`, cap with `max_window_sec` and record that the mode is
  non-decaying/unstable by the discrete eigenvalue criterion;
- if `dt` is missing, fall back to a fixed sample window and mark the provenance.

### Dimred Component Timescales

A dimred component usually mixes multiple eigenfunctions, so it does not have a
single native eigenvalue. Its envelope window should be derived from the
eigenfunction timescales and the component loading weights.

Recommended policy:

```text
tau_component = weighted_median(tau_modes, abs(component_loading))
```

The same metadata should also include mean/median summaries so P8/P10 can
compare raw efun top-hit timescales against dimred efun/component timescales:

```text
component_timescale_median_sec
component_timescale_mean_sec
component_timescale_weighted_median_sec
component_timescale_weighted_mean_sec
component_timescale_iqr_sec
component_timescale_source
```

Alternative summaries that may be useful for sensitivity checks:

```text
weighted_geometric_mean_tau
weighted_mean_log_tau
top_loading_mode_tau
```

Use `abs(component_loading)` or squared loading magnitude as weights. Save which
weight transform was used.

For methods where the loading matrix is not directly available or not
interpretable, fall back to:

```text
component_window_sec = fixed_default_window_sec
```

and record that the component timescale was not resolved.

### Metadata To Save

Every activity/envelope-derived output should save enough metadata to make the
activity definition auditable:

```text
activity_transform
activity_window_policy
activity_smooth_sec or window_sec_by_mode
session_aware
dt
lambda_discrete
lambda_continuous_method
timescale_sec_discrete_log
frequency_hz_discrete_angle
timescale_sec_bilinear
frequency_hz_bilinear
timescale_sec_preferred
window_alpha
min_window_sec
max_window_sec
dimred_timescale_summary_method
dimred_loading_weight_transform
input_was_complex
```

This metadata is required before P11 treats an envelope-derived product as a
current mainline source.

### P6

P6 heatmaps should be interpreted as:

- `abs` views: activity magnitude, mainline-interpretable;
- `real` views: unanchored signed diagnostics.

For deconvolved/residual eigenfunctions:

- `abs(u)` or envelope is the mainline activity magnitude;
- `real(u)` is diagnostic;
- residual cross-correlation strength should be ranked by `abs(r)`.

### P8/P10

Cross-modal coupling should prioritize phase/sign-invariant density sources:

- `efun_abs`;
- `deconv_abs`;
- dimred density from `abs(component)` or envelope;
- event-family density based on magnitude/envelope activity.

Signed `efun_real`, `deconv_real`, `raw_real`, and `raw_imag` outputs may be
kept as optional diagnostics, but they should not be the default conclusion
source for cross-session consistency.

For raw LFP efun density sources, P8/P10 top-hit and peak tables should retain
the raw eigenfunction index and eigenvalue-derived timescale metadata. The
default consistency summaries should include a numeric-first timescale
concentration check:

```text
top raw efun hits
-> join raw_efun_index to lambda/timescale metadata
-> plot top xcorr raw_efun_index distribution
-> plot normalized raw_efun_index quantile distribution
-> join dimred density hits to component mean/median timescale metadata
-> compare raw-hit timescales with dimred component timescales
-> summarize weighted median/mean timescale using peak_abs_corr weights
-> report fast/intermediate/slow or within-run quantile concentration
```

This check should be interpreted separately from density-source preference and
from dimred component subprocess labels.

### P11

P11 should report the current mainline as magnitude-based:

```text
mainline = eigenfunction/component activity magnitude
diagnostic = unanchored signed real/imag polarity
```

Completeness and summary views should label signed products as diagnostic unless
their provenance includes explicit polarity/phase anchoring metadata.

## Implementation Status, 2026-05-22

Implemented for P5 thresholded density:

- `build_blp_eigenfunction_reduction_params.m` defaults
  `density_value_transform='abs'`,
  `lfp_activity_transform='abs_magnitude'`, and
  `lfp_activity_window_policy='samplewise_abs_no_envelope'`.
- `get_thresholded_density.m` writes `mode_metadata` with raw efun index,
  eigenvalue, preferred discrete-log decay timescale/frequency, legacy bilinear
  continuous-time decay/frequency, activity transform, and threshold policy.
- `get_dimred_thresholded_density.m` writes
  `component_timescale_metadata` with weighted component timescale/frequency
  summaries and top source raw efun indices. These weighted component
  timescales now prefer the discrete-eigenvalue definition and keep bilinear
  values only as fallback/provenance.
- `scripts/script_run_current_p5_activity_density_unblocked.m` reruns the
  current unblocked P5 density grid one dataset at a time.
- P11 marks older P5 density MATs as stale when they predate the current
  activity-metadata code.

Still pending:

- activity-magnitude peak-statistics regeneration;
- dimred component process labels from the regenerated peak-statistics stage;
- selected-component 3D consensus-state trajectory export;
- optional RMS/envelope implementation using eigenvalue-aware windows.

## Optional Future Branch: Polarity/Phase Anchoring

If a future question requires signed biological interpretation, add a separate
anchoring branch before using positive/negative thresholds.

Possible anchoring policies:

- event-anchor: choose sign/phase so the mean response near a target event
  family is positive;
- reference-correlation anchor: choose sign/phase so correlation with a chosen
  reference signal is positive;
- loading-anchor: choose sign/phase so the largest loading or top-loading sum is
  positive;
- complex phase anchor: rotate `phi` by a phase factor so a chosen reference
  projection is real-positive.

Any anchored branch must save metadata:

```text
anchor_type
anchor_reference
sign_flip or phase_factor
anchor_score
source_signal
```

Until such metadata exists, signed polarity should not be interpreted as a
mainline biological direction.
