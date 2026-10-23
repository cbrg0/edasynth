# Implementation status

Last updated: 23 October 2026, version 1.1.0, before the preregistration tag.

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

Linux, 2 CPU cores, Python 3.13.16, numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1,
pandas 3.0.5, joblib 1.6.0, matplotlib 3.11.2. Three optional packages could not be
installed there: `cvxopt` (needed by the reference cvxEDA code), `wfdb` and
`PyWavelets`. PhysioNet could not be reached, so **no raw UTD file was read**. The
package has not been tested on macOS or Windows.

## Data and readers

| Element | Status | Evidence |
|---|---|---|
| Download script with checksum verification (`make data`) | BLOCKED | PhysioNet not reachable from the build environment |
| Zhang et al. files at commit `7c4bd46` (43 files) | PASS | `python -m edasynth check-data`: 43 of 43 match. The checksum list was computed from the same copy, so it protects later copies; it does not prove the original download. |
| WFDB format-16 reader with checksum and initial-value checks | PARTIAL | unit tests on files written by the tests; real `noneeg` files BLOCKED; comparison with the `wfdb` package skipped (not installed) |
| Annotation reader and stage assignment | PARTIAL | unit tests on synthetic annotation files; real annotation files BLOCKED |
| EDA unit of UTD (header says "NU") | NOT VERIFIED | µS is assumed until `make inspect` is run |
| UTD label loader | PASS | 9 211 windows, 349 majority artifacts, 20 groups (tests) |
| AWW loader | PASS | tests; used in the development run |

## Method steps

| Step ([METHOD.md](METHOD.md)) | Element | Status | Evidence |
|---|---|---|---|
| — | Real-data comparison on published features, logistic regression | PASS | AUPRC 0.4486 → 0.4754, Δ +0.027, 95 % CI −0.002 to +0.058 ([RESULTS.md §1](RESULTS.md#1-real-data-comparison)); a slow test reproduces it |
| — | Reproduce the earlier random-forest Δ = 0.054 | FAIL | this code gives Δ +0.026 (CI −0.021 to +0.056); the earlier settings were not recorded |
| — | Real-data comparison on features recomputed from the raw signals (G7) | BLOCKED | needs raw UTD; the command path has never been executed |
| §2.4 | Label–signal alignment | BLOCKED | needs raw UTD. The shift search and the cropping by the found offset are unit-tested. |
| §3 | Parameter estimation on UTD | BLOCKED | needs raw UTD |
| §3 | Parameter estimation, development run | PASS (runs to completion) | 11 recordings; τ<sub>0</sub> = 2.0 s (lower edge of the range) on all of them; response rates 0.3–3.3 per minute |
| §3.3 | Decomposition solver (numpy) | PARTIAL | recovers known inputs (tests); KKT optimality reported per run; agreement with the authors' cvxEDA code NOT VERIFIED (code not available in the build environment) |
| §3.3 | τ<sub>0</sub> selection by the objective value | PASS | recovers τ<sub>0</sub> = 2.5 s and 4.0 s within 0.5 s on synthetic records (test) |
| §3.4 | Response detection | PASS | recovers onset and height of isolated responses; the rebuilt phasic part matches (tests) |
| §3.6 | Artifact class rules | NOT VERIFIED | not yet checked by eye on UTD episodes |
| §4 | Generator: shapes, labels, donors, prevalence | PASS | unit tests (template heights and durations, label support, C3 transition only, derangement, Poisson counts) |
| §4.4 | Movement–artifact link | PASS | tests: *c* = 0 gives no movement-caused artifact; κ scales size with intensity. The delay δ<sub>max</sub> has no separate test. |
| §5 | Feature sets F0, FA, FA4, FB | PARTIAL | unit tests (counts, causality of FB, planted lag, published AWW features reproduced); comparison of the Haar details with PyWavelets skipped (not installed) |
| §6 | LOSO models, AUPRC, bootstrap | PASS | tests; real-data reproduction above |
| §9 | Untrained detectors | PARTIAL | run in the development run; tested on simple planted cases; not compared with the original authors' code |

## Validation checks

| Check | UTD | Development run (AWW) |
|---|---|---|
| G1 unit tests | PASS: all tests pass; 5 are skipped where optional packages or raw data are missing | same suite |
| G2 parameter recovery | BLOCKED | FAIL (sizes of responses: KS *p* = 0.023) |
| G3 realism, clean windows | BLOCKED | FAIL (AUC 0.87, limit 0.70) |
| G4 realism, artifact windows | BLOCKED | FAIL (AUC 0.99) |
| G5 negative controls | BLOCKED | FAIL (adding acceleration features lowered AUPRC without any link) |
| G6 positive controls | BLOCKED | FAIL (C4 control produced no artifacts) |
| G7 real-data consistency | BLOCKED | not applicable |

Details in [RESULTS.md §2](RESULTS.md#2-development-run-on-aww-not-a-result).

## Experiment and outputs

| Element | UTD | Development run |
|---|---|---|
| Pilot and number of replicates | NOT RUN (comes after the checks) | PASS (runs to completion) |
| Main grid, 200 cells | NOT RUN (comes after the tag) | PASS (runs to completion: 8 cells × 3 replicates, resumable) |
| Variants and label-noise runs | NOT RUN | NOT RUN |
| Analysis, compatibility, sensitivity | NOT RUN | PASS (runs to completion) |
| Power for the laboratory study | NOT RUN | PASS (runs to completion) |
| Figures F1–F7 | NOT RUN | PASS (produced and inspected; the inspection found two errors that were fixed) |
| Two-hand laboratory study | NOT IMPLEMENTED | — |

## Known problems

1. **Not run on the intended data.** Everything that needs the raw PhysioNet files
   (inspection, alignment, estimation, checks G2–G7, the grid) has only been run on
   the development data, or not at all.
2. **The development data fail every check.** The causes found so far are listed
   under open decisions 3–6.
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
7. **Computing cost of the main grid.** Extrapolated from the development run, the
   main grid needs roughly 70 CPU hours with 10 replicates and roughly 350 with 50
   ([REPRODUCING.md §5.3](REPRODUCING.md#53-computing-time)). On the development data
   the spread of Δ between replicates (SD about 0.04) would call for more than the cap of
   50 replicates. Whether the run is
   feasible on the available hardware, or the design has to be reduced, has to be
   decided before the tag.

## Open design decisions (to settle before the tag `prereg-synth-v1`)

1. **G5 and a significant loss.** G5 fails when adding acceleration features makes
   detection significantly *worse*. The check was designed to catch a spurious
   *gain*; whether a loss should also fail it must be decided before the checks are
   run on UTD.
2. **Positive control C4.** If the walking gain *G* measured on UTD is also close to
   zero, the C4 control cannot produce artifacts. A fixed control gain would then be
   needed, chosen before seeing G6 on UTD.
3. **τ<sub>0</sub> at the edge of its range.** On AWW every recording chose
   τ<sub>0</sub> = 2.0 s. If this repeats on UTD, the range should be widened.
4. **Absolute or relative artifact sizes.** Artifact sizes are stored in µS. A level
   shift measured on a day with a level of 5–6 µS, added to a donor at about 1 µS,
   pushed a synthetic signal down to 0.2 µS (figure F1 of the development run).
   Either sizes are scaled to the donor's level, or the multiplicative variant V2
   becomes the main variant.
5. **Share of level shifts.** The class rules put 68 % of AWW episodes into C3,
   while the AWW traces mostly show short drops. The rules should be checked by eye
   on UTD episodes before their shares are used.
6. **Coupling features in the positive control.** On AWW the FB features lowered
   AUPRC when every artifact was caused by movement. Whether this is a property of
   AWW or a defect of FB has to be checked on UTD.
