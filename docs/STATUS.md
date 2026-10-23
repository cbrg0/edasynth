# Implementation status

Updated: 2026-10-23, version 1.0.0, before the preregistration tag.

Statuses mean exactly this and nothing more:

| Status | Meaning |
|---|---|
| **PASS** | The code was run and the stated check succeeded. |
| **FAIL** | The code was run and the stated check did not succeed. |
| **PARTIAL** | Run and checked on part of the inputs only (named in the row). |
| **BLOCKED** | Could not be run in the build environment; the reason is named. |
| **NOT RUN** | Could be run, was not run yet. |
| **NOT IMPLEMENTED** | No code exists. |
| **NOT VERIFIED** | Code exists and runs, but its correctness has no independent check yet. |

No row is marked PASS from reading code. "Smoke" means the development run on
AWW (one person, 11 recording days used as donors) with `configs/smoke.yaml`;
its numbers check that the code runs and are **not results**.

## Build environment

Python 3.13.16, numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1, pandas 3.0.5,
joblib 1.6.0, matplotlib 3.11.2; Linux, 2 cores. Not installable there (no
package index access): `cvxopt` with the reference `cvxEDA.py`, `wfdb`, `PyWavelets`. PhysioNet was not
reachable, so no raw UTD file was read.

## Data and readers

| Element | Status | Evidence |
|---|---|---|
| Download script with checksum verification (`make data`) | BLOCKED | PhysioNet not reachable from the build environment |
| Zhang et al. files at commit 7c4bd46, SHA-256 of 43 files | PASS | `edasynth check-data`: 43/43 OK (list computed from the same clone, so this checks integrity of later copies, not the clone itself) |
| WFDB format-16 reader with checksum and initial-value checks | PARTIAL | unit tests on files written by the tests; real `noneeg` files BLOCKED; comparison with the `wfdb` package skipped (not installed) |
| MIT annotation parser, stage assignment | PARTIAL | unit tests on synthetic annotation files; real `.atr` files BLOCKED |
| UTD EDA unit (header says "NU") | NOT VERIFIED | `make inspect` must be run; µS is an assumption |
| UTD expert labels loader | PASS | 9 211 windows, 349 majority artifacts, 20 groups (tests) |
| AWW loader (columns read by name) | PASS | tests; used in the smoke run |

## Method steps

| Step | Element | Status | Evidence |
|---|---|---|---|
| 1 | Real-data baseline, published features (logistic regression, LOSO) | PASS | AUPRC 0.4486 → 0.4754, Δ +0.027, 95 % CI −0.002 to +0.058 (`results/realdata.json`, slow test) |
| 1 | Same with random forest: reproduce the earlier Δ = 0.054 | FAIL | this code gives Δ +0.026, CI −0.021 to +0.056; the earlier settings were not recorded |
| 1 | Real-data analysis on recomputed features (`realdata --raw`) | BLOCKED | needs raw UTD |
| 2 | Label-signal alignment (`make align`) | BLOCKED | needs raw UTD |
| 3 | Parameter estimation on UTD (`make estimate`) | BLOCKED | needs raw UTD |
| 3 | Parameter estimation, smoke (runs to completion) | PASS | 11 recordings; plausibility warnings logged: τ0 = 2.0 s (lower grid edge) on all 11 days; response rates 0.3-3.3 per minute, none above the 10/min warning |
| 3 | cvxEDA numpy solver | PARTIAL | recovers known inputs (tests); KKT optimality reported per run; agreement with the authors' `cvxEDA.py` skipped (not available) |
| 3 | τ0 selection by objective value | PASS | recovers τ0 = 2.5 and 4.0 within 0.5 s on synthetic records (test) |
| 3 | Response detection (peaks of phasic, height ≥ 0.05 µS) | PASS | recovers onset and height of isolated responses (test); smoke rates 0.3-3.3 per minute |
| 4 | Generator, artifact classes C1-C5, labels | PASS | unit tests (shapes, label support, C3 transition only, derangement) |
| 4 | Generator realism on UTD (G3/G4) | BLOCKED | needs the UTD bank |
| 5 | Feature sets F0, FA, FA4, FB | PARTIAL | unit tests; Haar details checked against PyWavelets skipped (not installed); recomputed F0 vs published F0 on UTD BLOCKED |
| 6 | Movement-contact link (c, κ, δ_max, π_r) | PASS | unit tests (c = 0 gives no motion-caused artifacts; κ scales size with movement intensity); the δ_max delay is not tested separately |
| 7 | LOSO models, AUPRC, participant bootstrap | PASS | tests; real-data reproduction above |
| 7 | Untrained detectors (slope, Kleckner, Kocielnik, P_n, SWT) | PARTIAL | run on smoke data; tested on simple planted cases; not compared with the authors' own code |

## Validation gates

| Gate | UTD | Smoke (AWW, not results) |
|---|---|---|
| G1 unit tests | PASS: 52 passed, 5 skipped (reference cvxEDA, wfdb, PyWavelets or raw data missing) | same suite |
| G2 parameter recovery | BLOCKED | FAIL: τ0 error 0.0 s and rate error 0.16 within limits; KS test of response sizes p = 0.023 < 0.05 |
| G3 realism, clean windows (C2ST AUC ≤ 0.70) | BLOCKED | FAIL: AUC 0.87 |
| G4 realism, artifact windows | BLOCKED | FAIL: AUC 0.99 |
| G5 negative controls | BLOCKED | FAIL: with c = 0, adding acceleration features *lowers* AUPRC (Δ_A −0.024, CI −0.028 to −0.021; Δ_A4 −0.005, CI excludes 0); F_B passes both controls |
| G6 positive controls | BLOCKED | FAIL: C4-only control produced no positive windows (C4 gain on AWW ≈ 0.0003). Mixture with c = 1: Δ_A +0.019 (CI 0.008 to 0.030) passes, while Δ_B is −0.038 (CI −0.051 to −0.025) |
| G7 real-data consistency | BLOCKED | not applicable |

Smoke gate failures may differ on UTD and are not fixed by changing thresholds.
Several raise design questions to settle **before** the tag (see below).

## Experiment and outputs

| Element | UTD | Smoke |
|---|---|---|
| Pilot and R | NOT RUN (after gates) | PASS (runs to completion; R = 3 by the smoke cap) |
| Main grid, 200 cells | NOT RUN (after tag) | PASS (runs to completion: 8 cells × 3 replicates, resumable) |
| Variants, label noise | NOT RUN | NOT RUN |
| Analysis, compatibility (goal 5) | NOT RUN | PASS (runs to completion) |
| Sensitivity and power (goal 6) | NOT RUN | PASS (runs to completion) |
| Figures F1-F7 | NOT RUN | PASS (produced and inspected by eye; this inspection found two bugs, see CHANGELOG) |
| Two-hand laboratory study | NOT IMPLEMENTED | — |

## Open questions to settle before the tag `prereg-synth-v1`

1. **G5 with a significant loss.** G5 fails when Δ is significantly *negative*
   (extra uninformative features cost accuracy). The gate exists to catch a
   spurious *gain*. Whether a loss should fail it has to be decided before UTD
   gates are run, not after.
2. **C4 positive control.** If the UTD gain G of the walking oscillation is also
   close to 0, the C4 control cannot produce artifacts. A fixed control gain
   would then be needed; its value has to be chosen before seeing G6 on UTD.
3. **τ0 grid edge.** On AWW every recording chose τ0 = 2.0 s, the lower grid
   edge. If this repeats on UTD, the grid (2-4 s) should be widened before the
   grid run.
4. **Artifact size: absolute or relative.** Templates carry sizes in µS. On AWW,
   level shifts (C3) measured on days with a level of 5-6 µS were added to a
   donor at about 1 µS and pushed the synthetic signal down to 0.2 µS (F1). Either
   sizes are scaled to the donor's level or V2 (relative contact loss) becomes the
   main variant. To decide on UTD before the tag.
5. **C3 share.** The classification rules put 68 % of AWW episodes into C3
   (level shift), while the real AWW traces show mostly short drops. The rules
   should be checked on UTD episodes by eye before their shares are used.
6. **F_B under the positive control.** On AWW, adding the coupling features
   lowered AUPRC when every artifact was caused by movement. Whether this is a
   property of AWW or a defect of F_B has to be checked on UTD.
