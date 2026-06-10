# K13 cross-dataset narrative integration note

Date: 2026-06-04

This note organizes the current interpretation after the K13m18-K13m21
standardized complex-split P5/P8/P10 summary, and connects it with the older
E10gb1/E10gH1-centered narrative.

## 1. Current strongest narrative

The strongest current result is:

```text
BOLD deconv_efun top coupling is enriched for RG-no-pure-theta / ripple-like
BLP dimred components.
```

This is especially clear in the four new K13 datasets in P10:

```text
P10 deconv_efun top dimred hits:
    top1:  RG-no-theta 80.4%, ripple 8.8%
    top5:  RG-no-theta 75.4%, ripple 15.1%
    top10: RG-no-theta 69.1%, ripple 15.6%
```

Interpretation:

```text
BOLD deconv_efun is the cleaner BOLD-side feature for fast ripple/gamma-like
BLP activity.
```

This can support a biological hypothesis phrased cautiously:

```text
The deconvolved/residual BOLD eigenfunction component may index a fast
ripple-trigger-related network reorganization process.
```

A possible neuroscience framing is:

```text
ripple-trigger network reorganization during spontaneous/resting activity,
potentially related to memory consolidation or replay-associated intrinsic
state perturbations.
```

This should be presented as a working hypothesis, not as a finished causal
claim.

## 2. What changed after adding K13m18-K13m21

The older working narrative had two linked subprocesses:

```text
theta-like slow/state BLP subprocess
    -> expected to relate to BOLD efun

ripple-gamma-like fast/intrinsic-trigger BLP subprocess
    -> expected to relate to BOLD deconv_efun
```

The new K13 subset changes the balance:

```text
The RG/deconv side is supported.
The theta/efun side is not yet supported by strict P5 labels in most K13 data.
```

P5 strict gate in K13m18-K13m21:

```text
k13m18: strict theta + RG pair present in 5 method-k settings.
k13m19: no strict theta + RG pair.
k13m20: no strict theta + RG pair.
k13m21: no strict theta + RG pair.
```

Important interpretation:

```text
The lack of strict theta in k13m19-k13m21 does not invalidate the RG/deconv
result.  It means the two-sided theta-vs-RG narrative is currently too strong
for the new K13 subset.
```

## 3. Current narrative tiers

### Tier 1: Keep as main narrative

```text
Fast RG-like BLP components couple to BOLD deconv_efun.
```

Evidence:

```text
P10 deconv_efun top dimred hits are strongly enriched for RG-no-pure-theta and
ripple-selective labels.

P8 also supports the same direction, though weaker and with more mixed labels.
```

Suggested wording:

```text
The BOLD deconvolved eigenfunction branch preferentially couples to
ripple/gamma-like BLP subprocesses rather than to simple event density.
```

### Tier 2: Keep as technical support

```text
BLP Koopman/eigenfunction densities outperform simple event density.
```

Evidence:

```text
In both P8 and P10, top xcorr hits overwhelmingly come from raw or dimred BLP
eigenfunction densities, not from P2 event density.
```

Meaning:

```text
The useful BLP signal is not just event occurrence.  It is a Koopman-derived
latent activity/density representation.
```

### Tier 3: Keep as hypothesis, not as conclusion

```text
BOLD efun corresponds to theta/slow state.
```

Current status:

```text
Not cleanly supported by K13m18-K13m21 strict dimred labels.
Efun hits often include mixed, inactive, or RG-like labels.
```

What is still needed:

```text
raw BLP efun timescale distribution
whole-trace examples
session-wise state-change evidence
possibly softer theta-active labels
```

Suggested wording:

```text
BOLD efun may reflect slower state-like BLP dynamics, but the theta-specific
interpretation is not yet stable across K13 under the current strict label.
```

### Tier 4: Use as QC/spatial support only

```text
ROI summary / activation maps
```

Current status:

```text
ROI evidence depends strongly on BOLD observable.
```

K13m18-K13m21 raw_csplit_q070 ROI result:

```text
gsvd100_ds:
    very stable ROI subspace, but efun/deconv often collapse onto the same
    BOLD spatial subspace.

roi_mean:
    more interesting separation.
    deconv-selected ROI profiles are cross-dataset consistent.
    efun-selected ROI profiles are weak or inconsistent.
```

Interpretation:

```text
ROI is useful for spatial sanity checking and for identifying stable BOLD
subspaces.  It should not be the primary proof of two subprocesses.
```

## 4. How to put all K13 datasets together

The next K13-level summary should not merge all figures into one large plot.
It should make a compact dataset-status matrix first.

Rows:

```text
k13m17
k13m18
k13m19
k13m20
k13m21
k13m23
future K13 datasets
```

Columns:

```text
data QC status
P5 strict theta present
P5 RG-like component present
P5 strict theta+RG pair pass
P8 deconv RG enrichment
P10 deconv RG enrichment
raw/dimred density preference
raw_csplit_q070 roi_mean deconv ROI consistency
HP/local caveat
notes
```

The first K13-level readout should ask:

```text
1. Which K13 datasets have usable standardized csplit P5/P8/P10 outputs?
2. Is RG-like BLP activity present even when strict theta is absent?
3. Is deconv_efun -> RG enrichment repeated across K13?
4. Is this stronger in P10 than P8?
5. Does roi_mean show a deconv-related ROI footprint more consistently than
   efun?
```

This makes the current story less confusing because the primary axis becomes:

```text
deconv/RG consistency across datasets
```

instead of forcing every dataset into:

```text
theta-vs-RG two-subprocess separation
```

## 5. Recommended figure set for all-K13 integration

Only a few figures should be promoted to the main K13 story.

### Figure 1. K13 dataset status matrix

Purpose:

```text
Show which K13 datasets are usable and which gates each dataset passes.
```

This should be the first figure, because it prevents incomplete or suspicious
datasets from being interpreted as biological failures.

### Figure 2. P5 RG-like component presence map

Purpose:

```text
Show that RG-like BLP components are present across K13, even when strict theta
is weak.
```

Rows:

```text
dataset
```

Columns:

```text
method-k
```

Color:

```text
RG-like strict label present / absent
```

Theta can be shown as a separate annotation, but should not dominate the
current story.

### Figure 3. P8/P10 deconv strict-label composition

Purpose:

```text
Show that deconv_efun top dimred-density hits prefer RG-no-theta / ripple-like
labels.
```

Panels:

```text
P8 deconv_efun
P10 deconv_efun
```

Split by:

```text
top1, top5, top10
```

This is the current main evidence figure.

### Figure 4. Density-class competition

Purpose:

```text
Show that BLP efun-derived densities are more informative than event density.
```

Panels:

```text
P8 efun
P8 deconv_efun
P10 efun
P10 deconv_efun
```

Classes:

```text
event density
raw BLP efun density
dimred BLP efun density
```

### Figure 5. raw_csplit_q070 ROI footprint check

Purpose:

```text
Ask whether efun and deconv select the same BOLD ROI/mode footprint.
```

Main observable to highlight:

```text
roi_mean
```

Reason:

```text
gsvd100_ds is often dominated by a shared stable BOLD subspace.
roi_mean shows clearer efun/deconv separation in the current K13 subset.
```

## 6. Conservative integrated conclusion

Current best integrated wording:

```text
Across the currently analyzed standardized complex-split K13 datasets, the
most reproducible cross-modal finding is not a full theta-vs-ripple two-process
split.  Instead, the robust signal is that BOLD deconv_efun, especially after
BOLD-side dimensional reduction in P10, preferentially couples to
ripple/gamma-like BLP Koopman density components that are not pure-theta.

This supports a fast intrinsic-trigger interpretation: ripple/gamma-associated
BLP subprocesses may drive or mark deconvolved BOLD perturbations, consistent
with a ripple-trigger network reorganization hypothesis.  The theta/slow-state
BOLD efun branch remains plausible but is not yet established in the new K13
subset and should be treated as a separate hypothesis requiring raw-timescale
and whole-trace evidence.
```

## 7. Immediate next actions

1. Build the all-K13 dataset-status matrix.

2. Add `k13m17` and `k13m23` to the same matrix, but mark their caveats:

```text
k13m17: suspicious theta-event/raw-signal QC.
k13m23: check current P8/P10 coverage and BOLD observable availability.
```

3. Generate an all-K13 P8/P10 deconv label-composition figure.

4. Generate an all-K13 density-class competition figure.

5. Keep theta as a secondary hypothesis until a softer theta-active analysis
and raw-timescale/whole-trace evidence are reviewed.
