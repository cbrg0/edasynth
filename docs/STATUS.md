# Implementation status

Last updated: 25 October 2026, version 1.2.0, before the preregistration tag.

This page states for every part of the study whether it exists, whether it has been
tested and whether it has been run on the data it is meant for. It also lists known
problems and the design decisions that are still open.

## How to read the status labels

| Label | Meaning |
|---|---|
| **PASS** | The code was run and the stated check succeeded. |
| **FAIL** | The code was run and the stated check did not succeed. |
| **PARTIAL** | Run and checked on part of the inputs only; the row says which part. |
| **BLOCKED** | Could not be run where the repository was built; the row gives the reason. |
| **NOT RUN** | Could be run, but has not been run yet. |
| **NOT IMPLEMENTED** | No code exists. |
| **NOT VERIFIED** | The code runs, but its correctness has not been checked independently. |

No row is marked PASS on the basis of reading the code. "Development run" means the
end-to-end run on AWW described in [RESULTS.md §2](RESULTS.md#2-development-run-on-aww-not-a-result);
its numbers are not results.

## Environment in which the statuses were obtained

**Linux** (where the code is written): 2 CPU cores, Python 3.13.16, numpy 2.5.3,
scipy 1.18.1, scikit-learn 1.9.1, pandas 3.0.5, joblib 1.6.0, matplotlib 3.11.2.
`cvxopt`, `wfdb` and `PyWavelets` are not installed there and PhysioNet cannot be
reached, so the raw UTD files have never been read on Linux.

**macOS** on Apple Silicon, conda, Python 3.12, with `wfdb` (23 October 2026,
version 1.1.2): `make data` and `check-data` verified every file, `make test-all`
passed (59 tests), and `make before-tag` ran every step up to the pilot on the raw
UTD files ([RESULTS.md §3](RESULTS.md#3-first-run-on-the-raw-utd-files-version-112)).
Version 1.2.0 changes the generator; its checks on UTD have **not been run yet**.

Windows has not been tested.

## Data and readers

| Element | Status | Evidence |
|---|---|---|
| Download script with checksum verification (`make data`) | PASS | run on macOS: both sources downloaded, `check-data` reports all files verified. The first run exposed a bug (the `shasum` check could not read PhysioNet's checksum format and did not stop the script); fixed in 1.1.1 |
| Zhang et al. files at commit `7c4bd46` (43 files) | PASS | `python -m edasynth check-data`: 43 of 43 match. The checksum list was computed from the same copy, so it protects later copies; it does not prove the original download. |
| WFDB format-16 reader with checksum and initial-value checks | PASS | unit tests; all 20 real records read on macOS with 1.1.2, and the readings match the `wfdb` package (`make test-all`) |
| Annotation reader and stage assignment | PASS | unit tests; `make inspect` on macOS: 8 annotations per record, no stage error |
| EDA unit of UTD (header says "NU") | NOT VERIFIED | µS is assumed; the per-participant medians (0.03–7.8) are plausible for µS, which supports but does not prove it |
| UTD label loader | PASS | 9 211 windows, 349 majority artifacts, 20 groups (tests) |
| AWW loader | PASS | tests; used in the development run |

## Method steps

| Step ([METHOD.md](METHOD.md)) | Element | Status | Evidence |
|---|---|---|---|
| — | Real-data comparison on published features, logistic regression | PASS | AUPRC 0.4486 → 0.4754, Δ +0.027, 95 % CI −0.002 to +0.058 ([RESULTS.md §1](RESULTS.md#1-real-data-comparison)); a slow test reproduces it |
| — | Reproduce the earlier random-forest Δ = 0.054 | FAIL | this code gives Δ +0.026 (CI −0.021 to +0.056); the earlier settings were not recorded |
| — | Real-data comparison on features recomputed from the raw signals (G7) | PASS | 0.4564 → 0.4818, Δ +0.025 (macOS, 1.1.2) |
| §2.4 | Label–signal alignment | PASS | 20 of 20 accepted, correlation 1.0000, no offset (macOS, 1.1.2) |
| §3 | Parameter estimation on UTD | PASS (runs to completion, 1.1.2) | σ = 0 for all participants; τ<sub>0</sub> at the edge for 12 of 20; rates 0.05–4.0 per minute. Version 1.2.0 adds residual blocks, episode levels and the donor check: NOT RUN on UTD |
| §3 | Parameter estimation, development run | PASS (runs to completion) | 11 recordings; τ<sub>0</sub> = 2.0 s (lower edge of the range) on all of them; response rates 0.3–3.3 per minute |
| §3.3 | Decomposition solver (numpy) | PARTIAL | recovers known inputs (tests); KKT optimality reported per run; agreement with the authors' cvxEDA code NOT VERIFIED (code not available in the build environment) |
| §3.3 | τ<sub>0</sub> selection by the objective value | PASS | recovers τ<sub>0</sub> = 2.5 s and 4.0 s within 0.5 s on synthetic records (test) |
| §3.4 | Response detection | PASS | recovers onset and height of isolated responses; the rebuilt phasic part matches (tests) |
| §3.2 | Residual blocks for the noise (1.2.0) | PASS | unit tests (clean rest windows, fallback, centring, block length) |
| §3.6 | Artifact class rules | NOT VERIFIED | not yet checked by eye on UTD episodes; 71 % of UTD episodes are C3 |
| §3.6 | Relative artifact sizes (1.2.0) | PASS | unit tests (ratio to the level, floor at 0.05 µS, synthetic size proportional to the local level) |
| §4 | Generator: shapes, labels, donors, prevalence | PASS | unit tests (template heights and durations, label support, C3 transition only, derangement, Poisson counts, noise made of whole residual blocks) |
| §4.1 | Physiology-donor check and pairing (1.2.0) | PASS | unit tests (50 % rule; an excluded donor is never used; nobody paired with themself; balanced use) |
| §4.4 | Movement–artifact link | PASS | tests: *c* = 0 gives no movement-caused artifact; κ scales size with intensity. The delay δ<sub>max</sub> has no separate test. |
| §5 | Feature sets F0, FA, FA4, FB | PARTIAL | unit tests (counts, causality of FB, planted lag, published AWW features reproduced); comparison of the Haar details with PyWavelets skipped (not installed) |
| §6 | LOSO models, AUPRC, bootstrap | PASS | tests; real-data reproduction above |
| §9 | Untrained detectors | PARTIAL | run in the development run; tested on simple planted cases; not compared with the original authors' code |

## Validation checks

| Check | UTD, version 1.1.2 | UTD, version 1.2.0 | Development run (AWW), 1.2.0 |
|---|---|---|---|
| G1 unit tests | PASS (macOS) | PASS on Linux without raw data; `test-all` NOT RUN | same suite |
| G2 parameter recovery | PASS | NOT RUN | FAIL (τ<sub>0</sub> not recovered when responses are small relative to the noise) |
| G3 realism, clean windows | FAIL (AUC 0.98) | NOT RUN | PASS (AUC 0.59; was 0.87) |
| G4 realism, artifact windows | FAIL (AUC 0.95) | NOT RUN | FAIL (AUC 0.93; was 0.99) |
| G5 negative controls | FAIL (losses of FA and FA4; FB passes) | NOT RUN | FAIL (loss of FA; FA4 and FB pass) |
| G6 positive controls | FAIL (C4 control produced no artifacts; the mixture control passes through FB) | NOT RUN | FAIL (C4 control produced no artifacts; no gain in the mixture control) |
| G7 real-data consistency | PASS | unchanged by 1.2.0 | not applicable |

Details: UTD in [RESULTS.md §3](RESULTS.md#3-first-run-on-the-raw-utd-files-version-112),
AWW in [RESULTS.md §2](RESULTS.md#2-development-run-on-aww-not-a-result).

## Experiment and outputs

| Element | UTD | Development run |
|---|---|---|
| Pilot and number of replicates | run with 1.1.2: *R* = 50 (cap); to be repeated with 1.2.0 | PASS (runs to completion) |
| Main grid, 200 cells | NOT RUN (comes after the tag) | PASS (runs to completion: 8 cells × 3 replicates, resumable) |
| Variants and label-noise runs | NOT RUN | NOT RUN |
| Analysis, compatibility, sensitivity | NOT RUN | PASS (runs to completion) |
| Power for the laboratory study | NOT RUN | PASS (runs to completion) |
| Figures F1–F7 | NOT RUN | PASS (produced and inspected; the inspection found two errors that were fixed) |
| Two-hand laboratory study | NOT IMPLEMENTED | — |

## Known problems

1. **Version 1.2.0 not yet checked on UTD.** The changes of 1.2.0 answer the
   failures of the first UTD run, but `estimate`, `validate` and `pilot` have to be
   run again before anything is concluded from them.
2. **Realism failed on UTD with 1.1.2** (G3 AUC 0.98, G4 AUC 0.95). The diagnosed
   causes were changed in 1.2.0 ([RESULTS.md §3.4](RESULTS.md#34-causes-found-diagnosis-not-yet-tested-by-a-rerun)).
3. **Feature scaling differs between real and synthetic data.** The published UTD
   features are z-scored within each participant; the synthetic pipeline only
   standardises over the training data. The real Δ used in the compatibility analysis
   (goal 5) was obtained with the former.
4. **Choice of conditions for the power analysis.** The code uses the first three
   compatible cells in grid order, which favours small values of *c*. A selection
   rule with a stated reason should be fixed before the tag.
5. **Label-noise variant.** It trains on noisy labels and scores against the true
   labels, while on UTD both training and scoring use expert labels. Whether this is
   the right comparison for goal 5 is open.
6. **Optional dependencies.** `pip install -e ".[reference]"` builds `cvxopt` from
   source on some systems and can fail; conda (`environment.yml`) installs a prebuilt
   version. The authors' cvxEDA code must be downloaded by hand
   ([REPRODUCING.md](REPRODUCING.md#6-optional-cross-checks)).
7. **Computing cost of the main grid.** The UTD pilot with 1.1.2 asks for the cap of
   *R* = 50 replicates (largest SD of Δ 0.019). The pilot took
   about 66 s of CPU time per replicate of one cell (40 replicates in 11 minutes on
   4 cores, macOS, 1.1.2). At that rate the main grid (200 cells × 50 replicates)
   needs about 180 CPU hours, about 46 hours on 4 cores, and the label-noise grid
   the same again
   ([REPRODUCING.md §5.3](REPRODUCING.md#53-computing-time)). Whether the run is
   feasible, or the design has to be reduced, has to be decided before the tag.
8. **One physiology donor is used twice.** With S18 excluded, 19 valid donors serve
   20 movement donors, so in every replicate one physiology donor appears in two
   synthetic participants. Their EDA is still drawn independently.
9. **Very large synthetic responses.** Response sizes had been taken from every
   detected response, including those inside artifact windows; in the UTD bank of
   1.1.2, 25 of 1 092 heights exceeded the range of the recording, and the AWW run
   produced a synthetic response of about 11 µS. Fixed in 1.2.0 (sizes from clean
   windows only, [METHOD.md §3.4](METHOD.md#34-responses-and-response-rates)); on
   AWW the largest height is now 0.49 µS. Whether the fix removes all implausible
   sizes on UTD is shown by the warnings of the next `estimate` run.
10. **τ<sub>0</sub> recovery with small responses.** After the fix, G2 fails on AWW:
   with responses close to 0.05 µS and residual noise, τ<sub>0</sub> is not
   recovered ([RESULTS.md §2.2](RESULTS.md#22-validation-checks-on-aww)). If G2
   fails on UTD for the same reason, the estimation of τ<sub>0</sub> (or the way it
   is used) has to be reconsidered before the tag.

## Open design decisions (to settle before the tag `prereg-synth-v1`)

Settled on 25 October 2026, after the first UTD run, before the tag (version 1.2.0):

- **G5 and a significant loss.** The rule is kept: a significant loss also fails
  G5. FA and FA4 gains in the map are read against their own no-link loss
  ([METHOD.md §7](METHOD.md#7-validation-checks-g1g7)).
- **Absolute or relative artifact sizes.** Relative to the level before the
  episode ([METHOD.md §3.6](METHOD.md#36-artifact-episodes-and-their-classes)).
- **Noise model.** Blocks of the decomposition residual instead of Gaussian noise
  with σ, which is 0 on UTD ([METHOD.md §3.2](METHOD.md#32-noise-and-quantisation)).
- **Physiology donors with almost no EDA.** Excluded by the 0.05 µS rule
  ([METHOD.md §4.1](METHOD.md#41-donors)).

Still open:

1. **Positive control C4.** The walking gain measured on UTD is *G* = −0.0009, so
   the C4 control produced no artifacts. A fixed control gain, with a stated reason,
   is needed. It has to be said, when it is chosen, that G6 on UTD had already been
   seen.
2. **τ<sub>0</sub> at the edge of its range.** On UTD 12 of 20 participants chose an
   edge value (9 at 4.0 s, 3 at 2.0 s); on AWW every day chose 2.0 s. Widening the
   range departs from Greco et al. (2016) and needs a stated reason.
3. **Share of level shifts.** The class rules put 71 % of UTD episodes (68 % of AWW
   episodes) into C3. A possible cause is the noise floor of rule 3: with σ = 0 it is
   *q*/√12 ≈ 0.0003 µS, so any shift of at least 0.05 µS counts as a level shift
   (hypothesis). The rules should be checked by eye on UTD episodes.
4. **Coupling features in the positive control.** On UTD (1.1.2) FB gained in the
   full-link control (Δ<sub>B</sub> +0.045), while on AWW it had lost. This points to
   a property of AWW or of the earlier generator rather than a defect of FB; it
   has to be confirmed with 1.2.0.
5. **Feature scaling** and **conditions for the power analysis** (known problems 3
   and 4).
