# Results

This page reports what has been measured so far and separates four things: the
measured numbers, their interpretation, possible explanations that are still
hypotheses, and conclusions that the data do not allow. Measures such as AUPRC are
explained in [METHOD.md §6](METHOD.md#6-models-and-performance-measures).

| Part | Data | State |
|---|---|---|
| [1. Real-data comparison](#1-real-data-comparison) | UTD, published features | **final**: run and verified |
| [2. Development run](#2-development-run-on-aww-not-a-result) | AWW, one person | run; checks the code, **not a result** |
| [3. Synthetic map](#3-synthetic-map) | UTD-based semi-synthetic data | **not run yet** |

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
figures) completed without errors in about 50 minutes on 2 CPU cores. Outputs are in
`results/smoke/` and `figures/smoke/` after running `make smoke`.

### 2.2 Validation checks on AWW

| Check | Outcome | What it means |
|---|---|---|
| G2 parameter recovery | fail | τ<sub>0</sub> and response rates were recovered within limits (median errors 0.0 s and 16 %), but the distribution of response sizes differed (Kolmogorov–Smirnov *p* = 0.023). |
| G3 realism, clean windows | fail | A linear classifier separated real from synthetic clean windows with AUC 0.87 (limit 0.70). |
| G4 realism, artifact windows | fail | AUC 0.99: synthetic artifacts look clearly different from the labelled AWW artifacts. |
| G5 negative controls | fail | With no movement–artifact link, adding 96 acceleration features *lowered* AUPRC (Δ<sub>A</sub> = −0.024, 95 % interval −0.028 to −0.021) and 4 features slightly too (Δ<sub>A4</sub> = −0.005). The coupling features FB passed both negative controls. |
| G6 positive controls | fail | The walking-artifact control produced no artifacts, because the measured walking gain on AWW is almost zero. With every artifact caused by movement, Δ<sub>A</sub> was +0.019 (interval 0.008 to 0.030), but Δ<sub>B</sub> was −0.038. |

### 2.3 What the development run shows

- The pipeline runs, writes every output and can be resumed.
- The generator is not yet realistic enough for the main experiment. Several
  concrete causes were found and are listed with the open design decisions in
  [STATUS.md](STATUS.md): artifact sizes in absolute µS do not fit donors with a low
  signal level, the class rules put most AWW episodes into the level-shift class,
  and every recording chose τ<sub>0</sub> at the lower edge of the search range.
- Inspecting the figures of this run found two errors in the code, which were fixed
  before this version (see [CHANGELOG.md](../CHANGELOG.md)).

Whether the same problems appear on UTD is not known until the checks are run there.

---

## 3. Synthetic map

Not run. The main experiment requires the raw UTD recordings, all validation checks
passing on them, and the preregistration tag
([PREREGISTRATION.md](PREREGISTRATION.md)). Until then, this repository makes no
claim about when the accelerometer helps.
