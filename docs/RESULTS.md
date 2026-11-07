# Results

What was measured, what it means, and what it does not show. Hypotheses are marked as
such. Measures are defined in [METHOD.md §6](METHOD.md#6-models-and-measures); the
result files are in `results/`.

## 1. Real data: does acceleration help on UTD?

Published feature matrices of Zhang et al. (2017): 9 211 five-second windows of 20
participants, 349 (3.79 %) labelled as artifacts by at least two of three experts.
EDA-only detector (24 features) against EDA + acceleration magnitude (28 features),
leave one participant out. Command `python -m edasynth realdata --stored --raw`,
output `results/realdata.json`.

| Model | AUPRC, EDA only | AUPRC, + acceleration | Δ | 95 % CI of Δ |
|---|---|---|---|---|
| Logistic regression | 0.4486 | 0.4754 | +0.0269 | −0.0022 to +0.0580 |
| Random forest, 500 trees | 0.4641 | 0.4900 | +0.0259 | −0.0206 to +0.0555 |
| Logistic regression, features recomputed from the raw signals (G7) | 0.4564 | 0.4818 | +0.0253 | −0.0013 to +0.0572 |

- Per participant (17 with artifacts) the median Δ is +0.004; Wilcoxon *p* = 0.76.
- Acceleration shuffled within participants: mean Δ −0.001 over 20 shuffles.
- Guessing scores 0.038.

**Interpretation.** The data neither show a gain nor show that there is none: gains
from slightly below zero to about +0.06 are compatible with them. The result does not
say whether acceleration helps for real contact disturbances, because one labelling
criterion refers to movement, so the labels are not independent of it.

## 2. Is the generator trustworthy?

All checks on UTD with the final version 1.4.0 (Linux, 28 October 2026):

| Check | Result |
|---|---|
| G1 unit tests | PASS |
| G2 parameter recovery, cvxEDA's response shape | PASS: τ₀ error 0.0 s, rate error 0.144, KS *p* = 0.37 |
| G2, varied response shapes | PASS: τ₀ error 0.15 s, rate error 0.125, KS *p* = 0.33 |
| G3 realism, clean windows | PASS: median AUC 0.59 (10 repetitions 0.50–0.66) |
| G4 realism, artifact windows | PASS: median AUC 0.555 (0.45–0.73; 2 of 10 above 0.70) |
| G5 no link / shifted acceleration | no gain anywhere; small losses: Δ<sub>A</sub> −0.032, Δ<sub>A4</sub> −0.002, Δ<sub>B</sub> −0.002 (accepted exception) |
| G6 positive controls | PASS: C4 Δ<sub>B</sub> +0.008 (+0.006 to +0.010); mixture *c* = 1 Δ<sub>B</sub> +0.050 (+0.040 to +0.060) |
| G7 real-data consistency | PASS (§1, last row) |

**How the generator got there.** The first run on the raw files (version 1.1.2) failed
G3 (AUC 0.98) and G4 (0.95). Each failure led to a documented correction, never to a
changed threshold:

| Version | Correction | Failure that led to it |
|---|---|---|
| 1.2.0 | noise from blocks of the decomposition residual | G3: the noise estimate was 0, synthetic signals too smooth |
| 1.2.0 | artifact sizes relative to the signal level; S18 excluded as physiology donor; response sizes from clean windows only | sizes measured at one level added to another; one participant with almost no EDA; impossible response heights |
| 1.3.0 | artifacts use the measured course of real episodes | G4: synthetic artifacts too smooth (AUC 0.86) |
| 1.3.0 | control gain for C4; tonic follows the donor's own stages | G6: no C4 artifacts; negative clean signals |
| 1.4.0 | contact state for level shifts; realism judged against a generator that never saw the compared participant; varied response shapes in G2; delays up to 10 s | an error found in the 1.3.0 grid (below), and three circular checks |

**A withdrawn run.** The main grid was first computed with 1.3.0 and its summary was
seen. Checking an odd result revealed an error: overlapping artifacts multiplied, and
about 20 % of the samples in some cells exceeded twice the clean signal's 99th
percentile. That grid is not used; everything below comes from 1.4.0. The changes made
after seeing it are listed in [PREREGISTRATION.md](PREREGISTRATION.md).

**Caveat.** Realism of artifact windows is borderline: single repetitions reach 0.73,
and a split-half version fails in one direction (0.76), mostly because of S10, whose
recording has a coarser quantisation step (0.006 µS against about 0.001). The
leave-one-out scheme and the median over 10 repetitions were chosen after seeing this.

## 3. The synthetic map

Version 1.4.0, tag `prereg-synth-v3`. 300 cells (226 computed), 50 replicates of 20
synthetic participants each (11 300 rows, macOS, 28–30 October 2026). Outputs
`grid-summary-utd.csv`, `compatibility*.csv`, `sensitivity-utd.csv`.

**3.1 Where the accelerometer is useful** (mean Δ ≥ 0.01 and lower bound > 0), cells out
of 60 per composition:

| Composition | Δ<sub>A</sub> (96 features) | Δ<sub>A4</sub> (4 features) | Δ<sub>B</sub> (4 coupling features) |
|---|---|---|---|
| C1 peaks | 2 | 22 | 40 |
| C2 drops | 0 | 12 | 28 |
| C3 level shifts | 0 | 15 | 25 |
| C5 signal loss | 0 | 0 | 0 |
| mixture | 0 | 16 | 28 |
| **all 300** | **2** | **65** | **121** |

**3.2 The decisive parameter.** Mean Δ in the mixture by *c* (12 cells each):

| *c* | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| Δ<sub>A</sub> | −0.027 | −0.034 | −0.033 | −0.029 | −0.020 |
| Δ<sub>A4</sub> | −0.001 | −0.002 | +0.002 | +0.009 | +0.019 |
| Δ<sub>B</sub> | −0.001 | +0.002 | +0.009 | +0.018 | +0.032 |
| useful cells A4 / B | 0 / 0 | 0 / 0 | 0 / 6 | 4 / 10 | 12 / 12 |

Each parameter varied alone around the reference cell changes the mean Δ<sub>B</sub> by:
*c* 0.049, π<sub>r</sub> 0.012, δ<sub>max</sub> 0.006, κ 0.002. Mixture means of
Δ<sub>B</sub>: δ<sub>max</sub> 1 s +0.015, 5 s +0.012, 10 s +0.008; π<sub>r</sub> 0
+0.017, 0.5 +0.007.

**3.3 Signal loss (C5).** EDA alone reaches AUPRC 0.98–1.00, and every acceleration set
lowers it (Δ<sub>B</sub> about −0.025; Δ<sub>A4</sub> down to −0.058 at *c* = 1).

**3.4 Cost.** At 95 % specificity, the share of movement-evoked genuine responses that
fall in flagged windows (π<sub>r</sub> = 0.5) rises from 0.256 (EDA only) to 0.285
(+FB). This measure also counts windows that contain a real artifact, so it overstates
the cost.

**3.5 Compatibility with the real gain** (Δ = 0.027 within the central 95 % of a cell's
replicates): 9 of 60 mixture cells, one at *c* = 0.75 and eight at *c* = 1. With noisy
training labels (κ ≈ 0.62): 8 of 60, all at *c* = 1. No cell with *c* ≤ 0.5.

**3.6 Variants** (half-widths of the intervals ≤ 0.004):

| Variant (mixture unless noted) | *c* | Δ<sub>A4</sub> | Δ<sub>B</sub> |
|---|---|---|---|
| 15 % prevalence | 0 / 0.5 / 1 | −0.001 / +0.002 / +0.033 | +0.001 / +0.007 / +0.027 |
| window positive only if mostly corrupted | 0 / 0.5 / 1 | −0.003 / −0.005 / −0.005 | −0.001 / −0.001 / +0.001 |
| V2 multiplicative, C2 | 0 / 0.5 / 1 | −0.001 / +0.004 / +0.020 | −0.001 / +0.019 / +0.049 |
| V2 multiplicative, C3 | 0 / 0.5 / 1 | −0.003 / +0.002 / +0.027 | −0.001 / +0.020 / +0.059 |
| C4 only, control gain | 1 | +0.002 | +0.009 |

**3.7 Power of a laboratory study** (share of 50 simulated studies whose bootstrap
interval of Δ excludes zero; three compatible cells at the quartiles of Δ<sub>B</sub>,
all *c* = 1; standard error up to 0.07):

| Cell (mixture, *c* = 1) | Features | 5 | 10 | 15 | 20 participants |
|---|---|---|---|---|---|
| κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0.5 | FB / FA4 | 0.08 / 0.06 | 0.38 / 0.14 | 0.64 / 0.18 | 0.72 / 0.34 |
| κ = 0, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0.5 | FB / FA4 | 0.12 / 0.02 | 0.34 / 0.08 | 0.44 / 0.16 | 0.62 / 0.34 |
| κ = 1, δ<sub>max</sub> = 5 s, π<sub>r</sub> = 0 | FB / FA4 | 0.14 / 0.00 | 0.62 / 0.04 | 0.94 / 0.18 | 1.00 / 0.26 |

## 4. Checks on data the generator never saw

Planned before the runs (PREREGISTRATION.md, last section); command
`python scripts/external_checks.py edabe|wesad <file>`; outputs `results/external-*.json`.

**EDABE** (43 Shimmer recordings, other experts, no acceleration; 53 590 windows, 19.1 %
with an artifact). Features z-scored within each recording, as planned, because levels
differ between devices.

| Check | Result | Criterion |
|---|---|---|
| E1 artifact windows, real vs synthetic | AUC 0.678 (0.651–0.718); unscaled 0.744 | ≤ 0.70: PASS |
| E2 clean windows | AUC 0.546 (0.529–0.583); unscaled 0.914 | ≤ 0.70: PASS |
| E3 AUPRC on EDABE: trained on synthetic data | 0.556 (0.537–0.583) | ≥ (b) − 0.05: PASS |
| (b) trained on real UTD / (c) on EDABE itself / guessing | 0.537 / 0.611 / 0.191 | |

**WESAD** (W1: movement of 15 E4 wrist recordings in place of UTD movement; mixture,
κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0; 50 replicates of 15 participants). WESAD
has fewer movement events (median 0.76 against 1.22 per minute).

| *c* | Δ<sub>B</sub> WESAD / UTD | Δ<sub>A4</sub> WESAD / UTD |
|---|---|---|
| 0 | +0.000 / −0.001 | −0.002 / −0.002 |
| 0.5 | **+0.014** / **+0.018** | **+0.011** / +0.004 |
| 1 | **+0.043** / **+0.048** | **+0.055** / **+0.027** |

Bold = useful. The usefulness of FB is the same in all three cells: PASS.

## 5. What it means

- **One quantity decides the answer:** the share *c* of artifacts that movement
  actually causes. It matters about four times more than any other link parameter, and
  below 0.5 the accelerometer is never useful.
- **How acceleration is used matters.** Four features that test whether EDA changes
  right after a movement beat four summaries of movement strength (121 vs 65 useful
  cells) and are more stable across devices. The 96 features of Zhang et al. lose
  almost everywhere; with about 370 artifact windows per replicate a 120-input linear
  model probably overfits (hypothesis).
- **Acceleration should not be used to flag signal loss,** which EDA finds by itself.
- **The real UTD gain is reproduced only with a strong link (*c* ≥ 0.75),** if the
  generator is right. Two cautions: the labels partly depend on movement, and with
  less movement (WESAD) the same gain would correspond to a smaller *c*.
- **The model transfers.** Trained only on synthetic data, a detector works on another
  device's data as well as one trained on the real recordings.
- **Next step:** measure *c* in a laboratory study, about 15–20 participants for the
  coupling features.

## 6. What these results do not show

- How often movement disturbs the electrode in real life: the map is conditional on *c*.
- That synthetic data are realistic in every respect: the realism checks use a linear
  classifier on 24 features and pass close to their limit for artifacts. The synthetic
  tonic contains a slow oscillation that real recordings do not show (figure F1).
- Synthetic artifacts are easier to find than real ones (EDA-only AUPRC about 0.65 vs
  0.45), so the same Δ may mean something different at a different baseline.
- Power assumes UTD-like recordings and labels; another device or protocol may need
  more participants.
