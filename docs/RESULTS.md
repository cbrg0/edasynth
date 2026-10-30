# Results

This page reports what has been measured so far and separates four things: the
measured numbers, their interpretation, possible explanations that are still
hypotheses, and conclusions that the data do not allow. Measures such as AUPRC are
explained in [METHOD.md §6](METHOD.md#6-models-and-performance-measures).

| Part | Data | State |
|---|---|---|
| [1. Real-data comparison](#1-real-data-comparison) | UTD, published features | **final**: run and verified |
| [2. Development run](#2-development-run-on-aww-not-a-result) | AWW, one person | run; checks the code, **not a result** |
| [3. Checks on the raw UTD files](#3-first-run-on-the-raw-utd-files-version-112) | UTD raw recordings, versions 1.1.2 to 1.3.0 | run; with 1.3.0 all checks pass except G5, which fails by small losses (accepted exception) |
| [4. Synthetic map](#4-synthetic-map) | UTD-based semi-synthetic data | **run with 1.4.0** (the 1.3.0 run was withdrawn, §3.8) |
| [5. Checks on independent data](#5-checks-on-independent-data) | EDABE (43 recordings), WESAD (15) | run; all four checks pass |

---

## 1. Real-data comparison

**Question.** On the real UTD recordings, does a detector that also sees the
acceleration magnitude find the expert-labelled artifacts better than a detector
that sees EDA only?

**Data and procedure.** The feature matrices published by Zhang et al. (2017):
9 211 five-second windows from 20 participants, of which 349 (3.79 %) were labelled
as artifacts by at least two of three experts. Detector A uses the 24 EDA features;
detector B adds the mean, SD, maximum and minimum of the acceleration magnitude
(28 features). Both are evaluated with leave-one-participant-out. Command:
`python -m edasynth realdata --stored`, settings in `configs/default.yaml`
(section `realdata`), output `results/realdata.json`.

### 1.1 Measured

| Model | AUPRC, EDA only | AUPRC, EDA + acceleration | Δ AUPRC | 95 % CI of Δ (participant bootstrap) | AUROC, EDA only → with acceleration |
|---|---|---|---|---|---|
| Logistic regression | 0.4486 | 0.4754 | +0.0269 | −0.0022 to +0.0580 | 0.923 → 0.921 |
| Random forest, 500 trees | 0.4641 | 0.4900 | +0.0259 | −0.0206 to +0.0555 | 0.928 → 0.933 |
| Logistic regression, features recomputed from the raw signals (check G7) | 0.4564 | 0.4818 | +0.0253 | −0.0013 to +0.0572 | — |

The random-forest row was obtained on Linux. On macOS (Apple Silicon, Python 3.12) the
same command gave 0.4643 → 0.4903 (Δ +0.0261, CI −0.0210 to +0.0561, Wilcoxon
*p* = 0.569). The difference in the fourth decimal comes from the numerical
libraries; the logistic-regression values were identical on both systems.

The last row repeats the comparison with the 24 EDA features and 4 acceleration
features computed by this code from the raw PhysioNet recordings, z-scored within
each participant like the published matrices (`python -m edasynth realdata --raw`,
run on macOS on 23 October 2026). Both AUPRC values lie within 0.01 of the
published-feature values and Δ has the same sign, so check G7 passes. The Wilcoxon
test over participants gives *p* = 0.91.

Additional measurements on the same data:

- **Per participant.** Of the 17 participants with at least one artifact window,
  the median Δ is +0.004 (logistic regression). A Wilcoxon signed-rank test finds no
  consistent direction (*p* = 0.76; random forest *p* = 0.47).
- **Shuffled control.** When the acceleration features are shuffled between
  windows of the same participant, so that movement no longer lines up with the
  EDA, the mean Δ over 20 shuffles is −0.001 (SD 0.001).
- **Single features.** Used alone as a score, the standard deviation of the EDA
  first difference reaches AUPRC 0.502, higher than the 24-feature logistic
  regression (0.449). This comparison was made after seeing the data (post hoc).
- **Reference point.** A detector that guesses scores about 0.038, the artifact
  share.

### 1.2 Interpretation

- Both models score far above chance, so the EDA features carry much information
  about the expert labels.
- Adding acceleration raises AUPRC by about 0.026 to 0.027 in both models, but the
  95 % confidence interval includes zero in both. The data therefore **do not show a
  gain, and they do not show that there is no gain**. Gains between slightly below
  zero and about +0.06 are all compatible with them.
- The shuffled control lies below the observed Δ, which suggests that whatever
  information the acceleration adds depends on its timing relative to the EDA.
  This is a weak indication, not a significance test of the observed gain.
- The per-participant results vary in sign, so a pooled gain, if real, is not shared
  by every participant.

### 1.3 Possible explanations (hypotheses, not tested here)

- Many movements may not disturb the electrode, so acceleration often signals
  nothing (Zhang et al., 2017).
- Summary statistics of acceleration magnitude may be the wrong kind of information;
  features describing whether EDA and movement change together (FB, see
  [METHOD.md §5](METHOD.md#5-features)) may do better.
- Because one labelling criterion refers to movement, part of any gain may reflect
  agreement with how the labels were made rather than better detection of contact
  disturbances.

### 1.4 What these data cannot show

- Whether the accelerometer helps for real contact disturbances: the labels are not
  an independent ground truth.
- Whether the result holds for other devices, body sites, protocols or daily life.
- Whether a better model or better acceleration features would change the
  conclusion. The single-feature observation above suggests that the model
  comparison itself is not yet optimal.

### 1.5 Differences from an earlier analysis

An earlier, unpublished analysis of the same data by the same author reported
AUPRC 0.449 → 0.474 (Δ = +0.025) for logistic regression and Δ = 0.054 for the
random forest. With the settings documented here, logistic regression gives
0.4486 → 0.4754 (Δ = +0.027). Regularised logistic regression has a single optimum,
so the difference must come from the solver: with its default tolerance, `lbfgs`
reproduces the earlier values in this code, while `newton-cholesky` converges in a
few iterations to a slightly different, better-converged solution. The random-forest value 0.054 could **not be reproduced**; this code
gives +0.026, and the settings of the earlier run were not recorded.

---

## 2. Development run on AWW (not a result)

**Purpose.** AWW contains recordings of a single person. It was used to run every
step of the pipeline end to end (`make smoke`, settings in `configs/smoke.yaml`)
while the raw UTD files were not available. Because all "participants" are days of
the same person and the grid is tiny (8 cells, 3 replicates), **none of these
numbers is a finding about the accelerometer**. They are reported because they show
the state of the generator.

### 2.1 What was run

All eight steps (estimate, validate, detectors, pilot, grid, analyze, power,
figures) completed without errors, in about 50 minutes on 2 CPU cores with 1.1.x and 63 minutes with 1.2.0. Outputs are in
`results/smoke/` and `figures/smoke/` after running `make smoke`.

### 2.2 Validation checks on AWW

The first column is the run with version 1.1.x. The second is version 1.2.0
(residual noise, relative artifact sizes, response sizes from clean windows),
25 October 2026, Linux, 2 cores.

| Check | 1.1.x | 1.2.0 | What it means |
|---|---|---|---|
| G2 parameter recovery | fail (response sizes, KS *p* = 0.023) | fail (median τ<sub>0</sub> error 1.0 s; rate error 5 %; response sizes KS *p* < 0.001) | See below. |
| G3 realism, clean windows | fail (AUC 0.87) | pass (AUC 0.59, limit 0.70) | A linear classifier separates real from synthetic clean AWW windows only slightly better than chance. |
| G4 realism, artifact windows | fail (AUC 0.99) | fail (AUC 0.93) | Synthetic artifacts still look clearly different from the labelled AWW artifacts. |
| G5 negative controls | fail | fail | No link: Δ<sub>A</sub> −0.039 (−0.058 to −0.019); Δ<sub>A4</sub> −0.005 and Δ<sub>B</sub> −0.002, both intervals contain 0. Shifted acceleration: all three intervals contain 0 (Δ<sub>A</sub> −0.037, Δ<sub>A4</sub> +0.002, Δ<sub>B</sub> +0.002), but Δ<sub>A</sub> is larger than 0.01 in size, so it fails. |
| G6 positive controls | fail | fail | The walking control again produced no artifacts (*G* = 0.0003). With every artifact caused by movement, Δ<sub>B</sub> = +0.017 (−0.020 to +0.054) and Δ<sub>A</sub> = +0.003. |

**Why G2 now fails on AWW.** With response sizes taken from clean windows only,
the AWW responses are small (median height about 0.09 µS, close to the 0.05 µS
threshold), and the residual noise is no longer zero. In a diagnostic run of four
synthetic participants generated with τ<sub>0</sub> = 2.0 s, the estimate was
2.0, 3.0, 4.0 and 4.0 s; it was right for the participant with the smallest
residual (SD 0.002 µS) and wrong for those with larger residuals (SD 0.005 to
0.016 µS). A wrong τ<sub>0</sub> changes the measured heights, which explains the
KS result. So τ<sub>0</sub> cannot be recovered when responses are small compared
with the noise (finding of a diagnostic run, not a preregistered check). The smoke
configuration searches τ<sub>0</sub> only at 2, 3 and 4 s, which makes every error
at least 1 s. Whether this happens on UTD is not known; with 1.1.2, G2 passed there.

An intermediate run of 1.2.0, before the response-size fix, gave G2 pass, G3 AUC
0.56 and G4 AUC 0.92; figure F2 of that run showed a synthetic response of about
11 µS, which led to the fix ([METHOD.md §3.4](METHOD.md#34-responses-and-response-rates)).

### 2.3 What the development run shows

- The pipeline runs, writes every output and can be resumed.
- With 1.1.x the generator was not realistic. Version 1.2.0 makes the clean AWW
  windows pass the realism check; the artifact windows still fail it. The class
  rules put most AWW episodes into the level-shift class (514 of 799 with 1.2.0),
  which is one candidate cause ([STATUS.md](STATUS.md), open decision 3). Figure F1
  of the 1.2.0 run shows it: the real day has short downward spikes, while the
  synthetic participant climbs in steps from about 0.95 to 1.35 µS through added
  level shifts.
- Inspecting the figures of this run found two errors in the code, which were fixed
  before this version (see [CHANGELOG.md](../CHANGELOG.md)).

Whether version 1.2.0 passes the checks on UTD is not known until they are run there. A pass on AWW is weak evidence: AWW has one person and the smoke settings use few replicates.

---

## 3. First run on the raw UTD files (version 1.1.2)

**What was run.** `make before-tag` on macOS (Apple Silicon, 4 jobs), 23 October
2026, version 1.1.2: download and verification, all tests, inspection, alignment,
real-data comparison, parameter estimation, checks G1–G7 and the pilot. Every step
completed. The checks were run before the preregistration tag, as planned, so that
the generator can still be corrected.

### 3.1 Measured: data and estimation

- **Inspection.** All 20 records have 8 stage annotations and no stage error. The
  header gives the EDA unit as "NU"; the medians range from 0.03 to 7.8 per
  participant, which is plausible for µS. Two participants have very low levels:
  S18 (median 0.033, maximum 0.076) and S19 (median 0.070).
- **Alignment.** 20 of 20 participants accepted, correlation 1.0000, no offset.
- **Noise.** σ = 0 for all 20 participants. The quantisation step is about
  0.001 µS (0.006 µS for S10).
- **τ<sub>0</sub>.** 12 of 20 participants chose a value at the edge of the range
  2–4 s: 9 at 4.0 s, 3 at 2.0 s.
- **Responses.** Rates between 0.05 and 4.0 per minute; four participants (S06,
  S11, S18, S19) have 6 or fewer detected responses.
- **Artifact episodes.** 208 episodes: C3 147 (71 %), C2 43 (21 %), C1 14 (7 %),
  C5 2, C4 2. The walking gain *G* is −0.0009, practically zero.

### 3.2 Measured: validation checks

| Check | Outcome | Numbers |
|---|---|---|
| G1 unit tests | PASS | all tests passed |
| G2 parameter recovery | PASS | median τ<sub>0</sub> error 0.0 s; median relative rate error 0.094; KS test of response sizes *p* = 0.41 |
| G3 realism, clean windows | FAIL | AUC 0.981 (limit 0.70) |
| G4 realism, artifact windows | FAIL | AUC 0.946 |
| G5 no link (*c* = 0) | FAIL | Δ<sub>A</sub> −0.050 (−0.061 to −0.039); Δ<sub>A4</sub> −0.004 (−0.009 to −0.000); Δ<sub>B</sub> −0.002 (−0.005 to +0.001), passes |
| G5 acceleration shifted | FAIL | Δ<sub>A</sub> −0.052; Δ<sub>A4</sub> −0.006 (−0.008 to −0.004); Δ<sub>B</sub> −0.003 (−0.006 to +0.000), passes |
| G6 C4 only | FAIL | no positive windows: with *G* ≈ 0 the class produces no visible artifact |
| G6 mixture, *c* = 1 | passes | Δ<sub>A</sub> −0.007 (−0.022 to +0.009); Δ<sub>B</sub> +0.045 (+0.031 to +0.060) |
| G7 real-data consistency | PASS | §1.1, last row |

**Pilot.** The largest SD of Δ between replicates over the pilot cells is 0.019,
which would need *R* = 55 replicates for a half-width of 0.005, so the cap *R* = 50
applies.

### 3.3 Interpretation

- G2 passes: the estimation recovers the physiology it is given.
- G3 and G4 fail by a wide margin: the synthetic windows are easy to tell from the
  real ones, so a map made with this generator would not describe UTD-like data.
- G5 fails only through losses of FA and FA4. Without a link no method gains,
  which is what the check is for; the 96 acceleration features add noise. FB
  passes both negative controls, and in the full-link control FB gains while FA
  does not.
- G6 fails because the walking class cannot be generated when *G* ≈ 0.

### 3.4 Causes found (diagnosis, not yet tested by a rerun)

- **Noise.** With σ = 0 the generator added no noise, so the synthetic signal
  was smoother than the real one. The features that differ most between real and
  synthetic clean windows are those of fast changes: the Kolmogorov–Smirnov distance
  is 0.83 for the maximum and the minimum of the second difference, 0.77 for its SD
  and 0.72–0.73 for the maximum of the first difference and of the finest Haar level. *Changed in 1.2.0*: the noise is made of blocks
  of the donor's decomposition residual ([METHOD.md §3.2](METHOD.md#32-noise-and-quantisation)).
- **Artifact sizes in µS.** Episode sizes measured at one participant's level were
  added to donors with other levels (medians from 0.03 to 7.8 µS). *Changed in
  1.2.0*: sizes are stored relative to the level before the episode and scaled by
  the local synthetic level ([METHOD.md §3.6](METHOD.md#36-artifact-episodes-and-their-classes)).
- **Donor S18** has almost no measurable EDA. *Changed in 1.2.0*: recordings with
  more than half of the samples below 0.05 µS are not used as physiology donors
  ([METHOD.md §4.1](METHOD.md#41-donors)), and artifact episodes that start below
  0.05 µS are not used as templates.
- **G5 rule.** Kept unchanged; how FA and FA4 gains are read is stated in
  [METHOD.md §7](METHOD.md#7-validation-checks-g1g7).

The rerun with version 1.2.0 is reported in §3.5.

### 3.5 Rerun with version 1.2.0

**What was run.** `estimate`, `validate` and `pilot` on macOS, 25 October 2026,
version 1.2.0, after `make test` (66 passed, 1 skipped). The first `pilot` run
reused the replicates of 1.1.2 because of a fault in the resume logic (fixed in
1.2.1, [STATUS.md](STATUS.md) known problem 11); the numbers below come from a
second run after the old pilot file was deleted.

**Estimation.** S18 was excluded as a physiology donor. A response height larger
than the range of the recording was reported only for S18 (5 heights); with 1.1.2
there had been 25 such heights in 14 participants. No artifact episode started
below 0.05 µS, and the class counts are unchanged (C3 147, C2 43, C1 14, C5 2,
C4 2). τ<sub>0</sub> is at the edge of its range for the same 12 participants.

| Check | 1.1.2 | 1.2.0 |
|---|---|---|
| G1 unit tests | PASS | PASS |
| G2 parameter recovery | PASS (τ<sub>0</sub> 0.0 s, rate 0.094, KS *p* 0.41) | PASS (τ<sub>0</sub> 0.0 s, rate 0.158, KS *p* 0.16) |
| G3 realism, clean windows | FAIL (AUC 0.981) | **PASS (AUC 0.602, limit 0.70)** |
| G4 realism, artifact windows | FAIL (AUC 0.946) | FAIL (AUC 0.855) |
| G5, no link | FAIL | FAIL: Δ<sub>A</sub> −0.035 (−0.039 to −0.031); Δ<sub>A4</sub> −0.003 (−0.005 to −0.000); Δ<sub>B</sub> −0.001 (−0.003 to +0.001), passes |
| G5, acceleration shifted | FAIL | FAIL: Δ<sub>A</sub> −0.040 (−0.047 to −0.033); Δ<sub>A4</sub> −0.003 (−0.004 to −0.002); Δ<sub>B</sub> −0.002 (−0.004 to +0.001), passes |
| G6, C4 only | FAIL (no artifacts) | FAIL (no artifacts, *G* = −0.0009) |
| G6, mixture with *c* = 1 | passes through FB | passes: Δ<sub>A</sub> +0.027 (+0.006 to +0.047); Δ<sub>B</sub> +0.068 (+0.046 to +0.090) |
| G7 | PASS | PASS (unchanged) |

**Pilot.** Largest SD of Δ between replicates 0.0172, so *R* = 46
(1.96² × 0.0172² / 0.005² = 45.5, rounded up).

**Where real and synthetic windows still differ** (largest Kolmogorov–Smirnov
distances over the 24 EDA features):

- clean windows: the level features (mean, minimum, maximum of the EDA: 0.26); the
  fast-change features, which differed most with 1.1.2 (0.83), are now at most 0.15;
- artifact windows: the fast-change features (minimum of the first difference 0.73,
  maximum of the 2-Hz and 4-Hz Haar details 0.70 and 0.69).

**Interpretation.**

- The residual noise removed the main difference between real and synthetic clean
  windows; the clean part of the generator now passes its preregistered realism
  check on UTD.
- The synthetic artifacts still differ from the labelled ones in how fast the signal
  changes. The KS distance does not tell in which direction; that needs a look at
  the distributions before any change.
- G5 fails only through FA and FA4, as decided (no spurious gain, small losses).
  FB passes both negative controls. In the full-link control FB gains more than FA.
- G6 fails only because the C4 control cannot be generated with *G* ≈ 0.

**What these numbers do not show.** Nothing here is a result about the
accelerometer: the checks test the generator. Passing G3 means that a linear
classifier on the 24 EDA features cannot separate clean windows well; other
features could.

### 3.6 Version 1.3.0

**What was run.** `estimate` and `validate` (G2–G6) on Linux (2 cores) on a copy of
the PhysioNet files provided by the user and verified with `check-data`, 25 October
2026. G7 does not depend on the generator and is unchanged. Before the run, a check
of the clean signal built from the bank of 1.2.0 found it below 0 µS in 11 of 20
synthetic participants; the tonic was changed (METHOD.md §4.2).

The bank was written before the 1.3.0 commit existed, so its `gitCommit` field names
the 1.2.1 commit. The files it was built with (estimation, bank, configuration) were
not changed afterwards; the generator was changed once more (the tonic floor) before
`validate` started, and not after it. G2–G4 and G5–G6 were run in two calls
because the first was stopped by a time limit after G4; the version number was raised
to 1.3.0 between the calls, and the merged report carries 1.3.0. The result files of this run are distributed
with the 1.3.0 patch, so that the grid can be run on exactly this bank.

| Check | 1.2.0 | 1.3.0 |
|---|---|---|
| G2 parameter recovery | PASS | PASS (τ<sub>0</sub> error 0.0 s, rate error 0.144, KS *p* = 0.37) |
| G3 realism, clean windows | PASS (AUC 0.602) | PASS (AUC 0.661) |
| G4 realism, artifact windows | FAIL (AUC 0.855) | **PASS (AUC 0.352)** |
| G5, no link | FAIL | FAIL: Δ<sub>A</sub> −0.032 (−0.039 to −0.025); Δ<sub>A4</sub> −0.002 (−0.003 to −0.000); Δ<sub>B</sub> −0.002 (−0.004 to −0.000) |
| G5, acceleration shifted | FAIL | FAIL: Δ<sub>A</sub> −0.039 (−0.044 to −0.033); Δ<sub>A4</sub> −0.002 (−0.004 to +0.001), passes; Δ<sub>B</sub> −0.001 (−0.003 to +0.001), passes |
| G6, C4 only (control gain) | no artifacts | passes: Δ<sub>B</sub> +0.006 (+0.003 to +0.009); Δ<sub>A</sub> −0.019 (−0.023 to −0.016) |
| G6, mixture with *c* = 1 | passes | passes: Δ<sub>B</sub> +0.044 (+0.038 to +0.050); Δ<sub>A</sub> −0.017 (−0.023 to −0.011) |

**Pilot.** Largest SD of Δ between replicates 0.0160, so *R* = 40
(1.96² × 0.0160² / 0.005² = 39.3, rounded up). The pilot ran 40 replicates in
60 minutes on 2 cores.

**Interpretation.**

- With measured courses the synthetic artifact windows can no longer be told from
  the labelled ones by the linear classifier (an AUC below 0.5 in grouped
  cross-validation means no usable separation). The largest KS distance in artifact
  windows fell from 0.73 to 0.22.
- G3 still passes, but closer to its limit (0.66); the level features differ most.
- G5: no feature set gains without a link. All three lose slightly; FB's loss without
  a link (−0.002) is the smallest and its interval reaches −0.0002. This is the
  accepted exception of [PREREGISTRATION.md](PREREGISTRATION.md).
- G6: FB detects both links that exist. FA does not: it loses even when every
  artifact is made from the acceleration (C4) or caused by movement. In the mixture
  control (−0.017) FA is about 0.015 above its no-link loss in the same mixture
  (−0.032), so it carries some of the link, but 96 features with 20 participants cost
  more than they add. (The C4 control uses other artifacts than G5, so its −0.019 is
  not directly comparable with −0.032.)

**What these numbers do not show.** The checks test the generator. They say nothing
yet about when the accelerometer helps; that is the main grid.

### 3.7 Version 1.4.0

**What was run.** On Linux, 28 October 2026, with the bank of 1.3.0 (the estimation
did not change): G2 twice (with cvxEDA's response shape and with varied response
shapes), G3 and G4 with the leave-one-participant-out scheme (10 repetitions) and
the split-half scheme, G5, G6 and the pilot.

| Check | Result |
|---|---|
| G2, cvxEDA's shape | PASS: τ<sub>0</sub> error 0.0 s, rate error 0.144, KS *p* = 0.37 |
| G2, varied shapes (τ<sub>0</sub> ± 25 %, τ<sub>1</sub> 0.35–1.05 s) | PASS: τ<sub>0</sub> error 0.15 s, rate error 0.125, KS *p* = 0.33 |
| G3, leave one out | PASS: median AUC 0.59 (10 repetitions: 0.50–0.66) |
| G4, leave one out | PASS: median AUC 0.555 (10 repetitions: 0.45–0.73; 2 of 10 above 0.70) |
| G3, split half | 0.57 and 0.39 |
| G4, split half | 0.30 and **0.76** |
| G5 | FAIL by small losses, no gain: no link Δ<sub>A</sub> −0.032, Δ<sub>A4</sub> −0.002 (passes), Δ<sub>B</sub> −0.002 (−0.003 to −0.001); shifted Δ<sub>A</sub> −0.033, Δ<sub>A4</sub> −0.004, Δ<sub>B</sub> −0.001 (passes) |
| G6 | PASS: C4 control Δ<sub>B</sub> +0.008 (+0.006 to +0.010); mixture Δ<sub>B</sub> +0.050 (+0.040 to +0.060); Δ<sub>A</sub> negative in both |

**Pilot (1.4.0).** Largest SD of Δ between replicates 0.0206, which would need *R* = 65
(1.96² × 0.0206² / 0.005²), so the cap *R* = 50 applies. 40 replicates took 18 minutes on
2 cores (about 55 s of CPU time per replicate, against about 3 min with 1.3.0).

**Interpretation and honest caveats.**

- Without circularity, the realism of artifact windows is borderline: it passes on
  the median, but single repetitions reach 0.73, and the split-half version fails in
  one direction. In that direction most of the gap comes from S10 (60 of 162 real
  artifact windows of its half, quantisation step 0.006 µS against about 0.001 for the
  others; without S10: 0.61). Synthetic participants can only have the recording
  properties of their donors, so a participant recorded differently is not reproduced.
- The leave-one-out scheme and the use of the median over 10 repetitions were chosen
  after the split-half result and the variation between seeds had been seen. Both
  facts are reported here and in PREREGISTRATION.md.
- G2 with varied response shapes shows that the estimation does not depend on the
  data having exactly cvxEDA's response shape (within the variation tested).

### 3.8 The main grid of 1.3.0 was withdrawn

The main grid (154 cells × 40 replicates) and the variants were computed with 1.3.0
on 25–28 October 2026, and their summary was looked at. While checking an unexpected
result of the 15 % variant (every feature set lost about 0.10 at *c* = 1), an error in
the generator was found: overlapping artifacts multiplied, and in the C3 and mixture
cells about 20 % of the samples exceeded twice the 99th percentile of the clean
signal (METHOD.md §4.3). The results of that grid are therefore not reported as
findings. The grid is run again with 1.4.0. Because its summary had been seen, the
changes of 1.4.0 are listed with their reasons in PREREGISTRATION.md; none of them
was chosen to move the map in a particular direction, and the decision rule and
thresholds are unchanged.

---

## 4. Synthetic map

**What was run.** Version 1.4.0, bank `5d7e3e8ad09e`, tag `prereg-synth-v3`.
The main grid (300 cells, of which 226 are computed and the rest are copies of
equivalent cells; 50 replicates each, 11 300 rows) was computed on macOS from
28 to 30 October 2026. The variants (13 cells × 50) and the label-noise grid
(the 50 mixture cells × 50) were computed on Linux. Logistic regression with
fixed *C* = 1, leave-one-participant-out, 20 synthetic participants per replicate.
`analyze` wrote `grid-summary-utd.csv`, `compatibility-utd.csv`,
`compatibility-labelnoise-utd.csv` and `sensitivity-utd.csv`; `figures` wrote
F1 to F6.

### 4.1 Measured: where the accelerometer is useful

A cell counts as useful when the mean Δ is at least 0.01 and its lower 95 % bound is
above zero ([METHOD.md §8](METHOD.md)). Number of useful cells out of 60 per
artifact composition:

| Composition | Δ<sub>A</sub> (96 features) | Δ<sub>A4</sub> (4 features) | Δ<sub>B</sub> (4 coupling features) |
|---|---|---|---|
| C1 | 2 | 22 | 40 |
| C2 | 0 | 12 | 28 |
| C3 | 0 | 15 | 25 |
| C5 | 0 | 0 | 0 |
| mixture C1+C2+C3+C5 | 0 | 16 | 28 |
| **all 300** | **2** | **65** | **121** |

Mean Δ in the mixture cells by *c*, the share of artifacts caused by movement
(each value is the mean over 12 cells):

| *c* | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| Δ<sub>A</sub> | −0.027 | −0.034 | −0.033 | −0.029 | −0.020 |
| Δ<sub>A4</sub> | −0.001 | −0.002 | +0.002 | +0.009 | +0.019 |
| Δ<sub>B</sub> | −0.001 | +0.002 | +0.009 | +0.018 | +0.032 |
| useful cells Δ<sub>A4</sub> / Δ<sub>B</sub> | 0 / 0 | 0 / 0 | 0 / 6 | 4 / 10 | 12 / 12 |

Other parameters, mixture cells (mean Δ<sub>B</sub>): δ<sub>max</sub> 1 s +0.015,
5 s +0.012, 10 s +0.008; π<sub>r</sub> 0 +0.017, 0.5 +0.007; κ 0 and 1 both +0.012.
The largest mean Δ<sub>B</sub> of the whole grid is +0.074 (C1, *c* = 1, κ = 0,
δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0).

**Sensitivity** (each parameter varied alone around the reference cell; range of the
mean Δ): *c* 0.049 for Δ<sub>B</sub> and 0.029 for Δ<sub>A4</sub>; π<sub>r</sub>
0.012 and 0.002; δ<sub>max</sub> 0.006 and 0.003; κ 0.002 and 0.001.

**C5 (signal loss).** EDA alone finds these artifacts almost perfectly (AUPRC
0.98 to 1.00), and every accelerometer feature set lowers AUPRC: Δ<sub>B</sub>
about −0.025 at every *c*, Δ<sub>A4</sub> from −0.010 at *c* = 0 to −0.058 at *c* = 1.

**Cost** (share of genuine responses ≥ 0.05 µS whose peak falls in a flagged window,
at 95 % specificity). In cells with π<sub>r</sub> = 0.5, movement-evoked responses:
F0 0.256, F0+FA4 0.279, F0+FB 0.285. Spontaneous responses: F0 0.200, F0+FB 0.196.

### 4.2 Measured: compatibility with the real data

A mixture cell is compatible when the real gain on UTD (Δ = 0.027) lies within
the central 95 % of its replicate values of Δ<sub>A4</sub>.

- Exact training labels: **9 of 60** mixture cells are compatible, one with
  *c* = 0.75 (κ = 0, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0) and eight with *c* = 1.
- Training labels made noisy to κ ≈ 0.62: **8 of 60**, all with *c* = 1.
- No cell with *c* ≤ 0.5 is compatible in either version.

### 4.3 Measured: variants

| Variant | *c* | AUPRC F0 | Δ<sub>A</sub> | Δ<sub>A4</sub> | Δ<sub>B</sub> |
|---|---|---|---|---|---|
| prevalence 15 %, mixture | 0 | 0.79 | −0.014 | −0.001 | +0.001 |
| | 0.5 | 0.78 | −0.012 | +0.002 | +0.007 |
| | 1 | 0.73 | +0.019 | +0.033 | +0.027 |
| window positive only if most of it is corrupted, mixture | 0 | 0.51 | −0.036 | −0.003 | −0.001 |
| | 0.5 | 0.53 | −0.049 | −0.005 | −0.001 |
| | 1 | 0.52 | −0.047 | −0.005 | +0.001 |
| C4 (pressure on the electrode), control gain | 1 | 0.97 | −0.018 | +0.002 | +0.009 |
| V2 multiplicative, C2 | 0 | 0.77 | −0.022 | −0.001 | −0.001 |
| | 0.5 | 0.78 | −0.021 | +0.004 | +0.019 |
| | 1 | 0.81 | −0.008 | +0.020 | +0.049 |
| V2 multiplicative, C3 | 0 | 0.62 | −0.031 | −0.003 | −0.001 |
| | 0.5 | 0.62 | −0.034 | +0.002 | +0.020 |
| | 1 | 0.60 | −0.013 | +0.027 | +0.059 |

Half-widths of the 95 % intervals of these means are 0.004 or less.

### 4.4 Measured: power for a laboratory study

Power is the share of 50 replicates in which the participant-bootstrap 95 % interval
of Δ lies above zero, for three compatible mixture cells chosen at the quartiles of
their mean Δ<sub>B</sub> ([METHOD.md §8](METHOD.md)). All three have *c* = 1, because
almost every compatible cell has *c* = 1.

| Cell (all *c* = 1, mixture) | mean Δ<sub>B</sub> | Feature set | 5 | 10 | 15 | 20 participants |
|---|---|---|---|---|---|---|
| κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0.5 | 0.032 | FB | 0.08 | 0.38 | 0.64 | 0.72 |
| | | FA4 | 0.06 | 0.14 | 0.18 | 0.34 |
| κ = 0, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0.5 | 0.035 | FB | 0.12 | 0.34 | 0.44 | 0.62 |
| | | FA4 | 0.02 | 0.08 | 0.16 | 0.34 |
| κ = 1, δ<sub>max</sub> = 5 s, π<sub>r</sub> = 0 | 0.041 | FB | 0.14 | 0.62 | 0.94 | 1.00 |
| | | FA4 | 0.00 | 0.04 | 0.18 | 0.26 |

With 50 replicates each value has a standard error of up to 0.07.

### 4.5 Interpretation

- **The share of movement-caused artifacts decides the answer.** In every
  composition except C5, the gain grows steadily with *c*, and *c* matters about
  four times more than any other link parameter. Below *c* = 0.5 no cell is useful.
- **Coupling features do better than summary features at equal count.** With four
  features each, FB is useful in 121 cells and FA4 in 65. The 96 summary features (FA)
  lose in almost every cell: with about 370 artifact windows per replicate, a
  linear model with 120 inputs overfits (this explanation is a hypothesis; it was not
  tested separately).
- **A longer or more variable delay between movement and artifact weakens the
  gain**, as expected for features computed over lags of 0 to 5 s. With
  δ<sub>max</sub> = 10 s FB is still useful at *c* = 1.
- **Real responses to movement make the task harder**, since a movement then
  predicts a genuine response as well as an artifact (in the mixture cells
  π<sub>r</sub> = 0.5 lowers the mean Δ<sub>B</sub> from 0.017 to 0.007).
- **Where EDA alone is already near perfect (C5), the accelerometer only adds
  errors.**
- **A laboratory study needs about 15–20 participants for FB** in the compatible
  cells (power 0.44 to 0.94 with 15 and 0.62 to 1.00 with 20); FA4 reaches at most 0.34
  with 20, so with four summary features the gain would usually be missed.
- **If the generator is right, the real UTD gain points to a strong link.** The
  real Δ of 0.027 is reproduced only when nearly all artifacts are caused by
  movement. This is a conditional statement about the model, and the caveats below
  limit it.

### 4.6 Caveats

- **The labels.** The real Δ refers to expert labels, and one of the labelling
  criteria refers to movement (§1.3). A link between movement and labels can
  therefore come from the labelling as well as from the electrode. Compatibility at
  *c* = 1 is consistent with both readings.
- **Difficulty differs.** On the synthetic data EDA alone reaches AUPRC of about 0.65
  in the mixture cells, against 0.45 on the real data. The synthetic artifacts are
  easier to find than the real expert-labelled ones, and the same Δ may mean
  something different at a different baseline.
- **Window labels.** With the rule that a window counts as an artifact only when
  most of it is corrupted, the gain disappears even at *c* = 1. The map therefore
  depends on how windows are labelled; the main rule (any corrupted sample) was
  fixed in the preregistration.
- **Cost measure.** The share of lost responses also counts responses whose window
  contains a real artifact, so part of the higher cost with the accelerometer is
  correct flagging of corrupted windows. A cost measure restricted to clean windows
  was not computed; this is planned for the laboratory study.
- **Realism.** The checks of realism pass on the median but are borderline (§3.7).
  The tonic level of the synthetic signals is a line plus a slow sinusoid fitted to
  the donor (METHOD.md §4.2), and in F1 this slow wave is visible where the real
  recordings do not show one.
- **One data source for the generator.** The generator is built from 20 UTD
  participants recorded with one device at one body site. Checks on two independent
  data sets pass (§5), and they also show that the gain of FA4 depends on the movement
  record.
- **Power refers to the synthetic setting.** It assumes UTD-like recordings, the
  window labels of this study and a link as strong as *c* = 1. A laboratory study with
  another device, protocol or labelling can need more participants. The numbers are a
  first orientation for planning, with about 20 participants as the lower end when
  responses to movement are common (π<sub>r</sub> = 0.5).
- **What the map does not say.** It does not measure *c* in the real world. It says
  what gain to expect for each value of *c*, and the laboratory study is designed to
  measure *c* directly.

---

## 5. Checks on independent data

Planned in [EXTERNAL-CHECKS.md](EXTERNAL-CHECKS.md) on 30 October 2026 before they were
run (commits `7e197d8` and `766a2c8`); version 1.4.0 unchanged; Linux; command
`python scripts/external_checks.py edabe|wesad <file>`; outputs
`results/external-edabe.json`, `results/external-wesad.json`.

### 5.1 Measured: EDABE

43 recordings, 53 590 windows of 5 s after removing leading and trailing zeros, 19.1 %
of them with at least one sample labelled as artifact.

| Check | Main version (features z-scored within recording) | Unscaled features | Criterion |
|---|---|---|---|
| E1, artifact windows, real against synthetic | AUC 0.678 (10 repetitions: 0.651–0.718) | 0.744 (0.681–0.781) | ≤ 0.70: **PASS** |
| E2, clean windows | AUC 0.546 (0.529–0.583) | 0.914 (0.885–0.929) | ≤ 0.70: **PASS** |

| E3, AUPRC on all 43 EDABE recordings | |
|---|---|
| (a) detector trained on synthetic data (median of 10) | 0.556 (0.537–0.583) |
| (b) detector trained on the 20 real UTD recordings, expert labels | 0.537 |
| (c) EDABE itself, leave one recording out (reference) | 0.611 |
| guessing (share of artifact windows) | 0.191 |

Criterion (a) ≥ (b) − 0.05: **PASS**.

### 5.2 Measured: WESAD

Movement donors: the 15 WESAD wrist recordings (83–104 min each after cropping to the
protocol). Movement events per minute, median: 0.76 on WESAD against 1.22 on UTD.
Physiology donors: those of UTD. Mixture, κ = 1, δ<sub>max</sub> = 1 s,
π<sub>r</sub> = 0, 50 replicates of 15 synthetic participants per cell; UTD values are
the same cells of the main grid (50 replicates of 20).

| *c* | Δ<sub>B</sub>, WESAD movement | Δ<sub>B</sub>, UTD movement | Δ<sub>A4</sub>, WESAD movement | Δ<sub>A4</sub>, UTD movement |
|---|---|---|---|---|
| 0 | +0.000 (−0.000 to +0.001) | −0.001 (−0.002 to −0.000) | −0.002 (−0.003 to −0.002) | −0.002 (−0.003 to −0.002) |
| 0.5 | **+0.014** (+0.013 to +0.016) | **+0.018** (+0.015 to +0.020) | **+0.011** (+0.010 to +0.013) | +0.004 (+0.003 to +0.005) |
| 1 | **+0.043** (+0.040 to +0.046) | **+0.048** (+0.044 to +0.052) | **+0.055** (+0.052 to +0.058) | **+0.027** (+0.025 to +0.029) |

Bold: useful (mean ≥ 0.01, lower bound > 0). Criterion: the decision for FB is the same
as with UTD movement in all three cells: **PASS**.

### 5.3 Interpretation

- **Synthetic artifacts resemble artifacts marked by other experts on another device**
  within the preregistered limit, though less closely than they resemble the UTD
  artifacts they were built from (E1 0.68 against G4 0.555). The check passes close to
  its limit, and one of 10 repetitions exceeds it (0.72).
- **A detector trained only on synthetic data transfers to EDABE as well as one trained
  on the real UTD recordings** (0.556 against 0.537), and comes within 0.06 of a detector
  trained on EDABE itself. This is the most direct answer available to the concern that
  the synthetic data only reproduce themselves: the generator never saw EDABE.
- **Levels differ between devices.** Without scaling within recording, clean windows of
  EDABE and synthetic data are almost fully separable (AUC 0.91). The comparison of
  shapes therefore needs the scaling, which was fixed in the plan before the run; the
  unscaled numbers are reported so that this choice can be judged.
- **The shape of the map for FB transfers to another wrist device and protocol.** With
  WESAD movement, FB is useless at *c* = 0 and useful at *c* = 0.5 and 1, with gains
  within about 0.005 of those with UTD movement.
- **The gain of the four summary features depends on the movement record.** With WESAD
  movement Δ<sub>A4</sub> is about twice as large (+0.055 against +0.027 at *c* = 1) and
  exceeds Δ<sub>B</sub>. A possible explanation, not tested: WESAD participants sat for
  most of the protocol and moved less often, so any movement is a more specific sign of
  an artifact. One consequence: the compatibility of the real UTD gain with *c* ≥ 0.75
  (§4.2) is a statement about UTD-like movement; with less movement the same real gain
  would correspond to a smaller *c*.
- **What EDABE cannot show.** It has no acceleration, so it says nothing about the gain
  from the accelerometer. Its labels are expert corrections of the EDA alone, which are
  not a ground truth of contact disturbances either.
- **What WESAD cannot show.** It keeps the assumed link between movement and artifacts,
  so it says nothing about *c* on the E4. It has no artifact labels.
