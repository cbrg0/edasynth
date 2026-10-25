# Results

This page reports what has been measured so far and separates four things: the
measured numbers, their interpretation, possible explanations that are still
hypotheses, and conclusions that the data do not allow. Measures such as AUPRC are
explained in [METHOD.md §6](METHOD.md#6-models-and-performance-measures).

| Part | Data | State |
|---|---|---|
| [1. Real-data comparison](#1-real-data-comparison) | UTD, published features | **final**: run and verified |
| [2. Development run](#2-development-run-on-aww-not-a-result) | AWW, one person | run; checks the code, **not a result** |
| [3. First run on the raw UTD files](#3-first-run-on-the-raw-utd-files-version-112) | UTD raw recordings, version 1.1.2 | run; the checks failed, the generator was changed in 1.2.0 |
| [4. Synthetic map](#4-synthetic-map) | UTD-based semi-synthetic data | **not run yet** |

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

Whether these changes make G3 and G4 pass is not known until the checks are run
again with version 1.2.0. The remaining open points (C3 share, the C4 control,
the τ<sub>0</sub> range, the cost of the grid) are listed in [STATUS.md](STATUS.md).

---

## 4. Synthetic map

Not run. The main experiment requires the raw UTD recordings, all validation checks
passing on them, and the preregistration tag
([PREREGISTRATION.md](PREREGISTRATION.md)). Until then, this repository makes no
claim about when the accelerometer helps.
