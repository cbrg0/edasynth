# Changelog

## 1.0.0 — 2026-10-22
First version of the semi-synthetic test bench, written from scratch.

* Readers for the PhysioNet `noneeg` records (format 16 with checksum verification,
  MIT annotations) and for the expert labels of Zhang et al. (2017).
* Label alignment check, parameter estimation from real recordings, semi-synthetic
  generator with five artifact classes and a parameterised movement-contact link,
  feature sets F0/FA/FA4/FB, LOSO models, untrained detectors, gates G2-G7, grid,
  compatibility and power analyses, figures F1-F7.
* Reproduction of the published-feature real-data baseline
  (AUPRC 0.4486 → 0.4754, logistic regression).
* Carried over from the earlier project repository: MIT licence, CI layout.

### Corrections made during the development (smoke) run, before any UTD result

These were found by tests or by the end-to-end run on AWW and fixed before the
preregistration tag; none of them was made after seeing a synthetic result.

* `estimation.scr_events`: responses are now peaks of the phasic component with
  trough-to-peak height >= 0.05 uS. Counting runs of non-zero driver gave 14-95
  "responses" per minute on AWW because the cvxEDA driver is dense at the optimum
  on noisy, quantised data (checked: the solver's KKT violation was < 1e-4, so the
  density is a property of the solution, not of non-convergence).
* `estimation.decompose_participant`: tau0 is chosen by the smallest cvxEDA
  objective value. The residual norm on clean samples, used first, was almost
  flat in tau0 and pointed to the wrong end of the grid on synthetic records with
  known tau0 (new test `test_tau0_selection_recovers_known_value`); G2 on AWW
  showed the same (median tau0 error 2.0 s).
* `decompose.bateman`: no overflow warning for negative times.
* `estimation.estimate_participant`: the tonic given to the generator is fitted
  to t + r − r_events (signal minus the detected responses). With cvxEDA's tonic
  alone, synthetic AWW records sat several uS below the real ones (seen in F1),
  because the phasic part also holds sub-threshold activity.
* `cli._run_jobs`: resumable results are rewritten with aligned columns; appending
  rows with different columns had shifted values in `grid_metrics` (found in the
  smoke run, new test).
* `figures.f4_map`: one file per artifact composition, constrained layout.
* `figures.f1_real_vs_synthetic`: the real column is the synthetic participant's
  physiology donor (comparable level) and both columns start 1 min before their
  first artifact window.
* Reference cvxEDA: no longer listed as a PyPI dependency (package name not
  verified); `cvxEDA.py` from the authors' repository is used when on the path.
* `Makefile`: `before-tag` and `after-tag` targets so that the preregistration
  tag is pushed between validation and the grid; `clean-results` keeps README files.
* `features.fb`: B4 coherence is 0 where a segment is constant (was NaN).
* Validation writes its report after every gate.
* `estimation.motion_events`: the 0.5-s minimum length is applied to raw runs
  before merging; before, isolated noise samples within 1 s merged into
  spurious events.
* `estimation.noise_sigma`: recordings without clean rest windows fall back to
  all clean windows (source recorded in the bank); episode classification uses
  the noise floor max(sigma, q/sqrt(12)) because sigma can be 0 when noise is
  below the quantisation step.
* `data.stage_type`: AWW resting days ("rest") are now typed as rest.
* `generator`: a donor without clean time of one stage type uses its other rate
  (counted in `Generator.rate_fallbacks`) instead of silently generating no
  responses.
* `estimation.plausibility_warnings`: logs response rates above 10/min, tau0 on
  the grid edge and solver KKT violations above 1e-3 (warnings, not gates).
* `bank.build_bank`: optional parallel estimation (`--jobs`).
* `models`: logistic regression now uses the `newton-cholesky` solver (same L2
  objective as `lbfgs`, ~20x faster here). `lbfgs` stopped at its default
  tolerance before the optimum: the real-data baseline moved from 0.4485 → 0.4739
  (Δ +0.025) to 0.4486 → 0.4754 (Δ +0.027, 95 % CI −0.002 to +0.058). The
  reference values of G7 and goal 5 were updated accordingly before the tag.
* `estimation.quantisation_step`: ignores differences that occur in < 1 % of
  non-zero steps.
