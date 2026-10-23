# Changelog

All notable changes to this project are recorded here. Changes to the method after
the preregistration tag must also be released under a new tag
([docs/PREREGISTRATION.md](docs/PREREGISTRATION.md)).

## 1.1.1 (23 October 2026)

- `scripts/download-data.sh` verifies both data sources with
  `python -m edasynth check-data` and stops when a file is missing or differs.
  PhysioNet's `SHA256SUMS.txt` separates hash and name by one space, which
  `shasum -c` rejects ("no properly formatted SHA checksum lines found"); the
  previous script printed that message and continued without verifying. Found on the
  first run on macOS, where `check-data` then confirmed all files.
- Documentation: the unit tests and `make data` have now run on macOS (Apple Silicon,
  Python 3.12); `shasum` is no longer required.

## 1.1.0 (23 October 2026)

### Naming without underscores

This release renames things that users of 1.0.x refer to; result files and
configuration files written for 1.0.x are not read by 1.1.0.

- All identifiers written by the project use lower camel case (for example
  `deltaMax`, `piR`, `loadConfig`), including configuration keys, keys of the
  parameter bank and result columns (`auprcF0`, `lostMotionF0+FB`).
- Output files and figures use hyphens: `bank-utd.json`, `grid-metrics-utd.csv`,
  `pilot-utd-R.json`, `F4-map-C1.pdf`; `scripts/download-data.sh`;
  `docs/README-RU.md`; module `edasynth.wfdbio` (was `io_wfdb`); test files
  `tests/testEstimation.py` and so on, collected through `python_files` and
  `python_functions` in `pyproject.toml`.
- Figure labels write symbols with proper subscripts (κ, δ<sub>max</sub>,
  π<sub>r</sub>, Δ<sub>A4</sub>).
- The seed of a grid replicate is derived from the condition key, whose text changed
  (for example `pi_r=` became `piR=`). Replicates of 1.1.0 therefore differ from those
  of 1.0.x; the validation checks and the real-data result use other seeds and are
  unchanged (verified by rerunning both).
- `docs/DEVELOPMENT.md` documents the naming rule and its exceptions (Python, library
  interfaces, external data file names).

## 1.0.1 (23 October 2026)

### Documentation

- The documentation was rewritten for readers without prior knowledge of the project.
  New files: `docs/BACKGROUND.md` (the problem and the research gap),
  `docs/RESULTS.md` (what has been measured, kept apart from interpretation),
  `docs/REPRODUCING.md` (installation, data, every step, computing time,
  troubleshooting) and `docs/DEVELOPMENT.md` (architecture and conventions).
  `README.md`, `docs/METHOD.md`, `docs/STATUS.md`, `docs/PREREGISTRATION.md`,
  `docs/GLOSSARY.md`, `docs/REFERENCES.md`, `docs/README-RU.md` and the README files in
  `data/`, `results/` and `figures/` were rewritten.
- Source attributions were checked against the publications. Corrections:
  - The tonic form (linear trend plus a sinusoid with a period of 45–90 s) is that of
    the synthetic test signal of Greco et al. (2016); the source attribution removed
    in 1.0.0 was restored.
  - Greco et al. (2016) chose τ<sub>0</sub> by the residual norm. The selection by the
    objective value used here is now described as a departure from the source.
  - Greco et al. (2016) do not state that the input is z-scored; z-scoring is now
    described as a project decision.
  - Posada-Quintero and Chon (2020) do not state that movement can evoke genuine
    responses. The motivation for movement-evoked responses now cites the example given
    by Taylor et al. (2015) and is labelled an assumption.
  - The 0.05 µS minimum response amplitude is supported by Kleckner et al. (2018) and,
    through Boucsein et al. (2012), by Posada-Quintero and Chon (2020).
  - `docs/REFERENCES.md` now states for every source whether the attributed claim was
    checked against the text.
- Corrected claims: only the parameter bank records the git commit, not every results
  file; the computing time of the main grid was re-estimated from a measurement and
  can reach several hundred CPU hours.

### Code

- The offset found by the label–signal alignment is now applied. Before, a participant
  accepted with a shift was used without it in estimation, validation, the real-data
  comparison on recomputed features and the figures. `Recording.cropped` and
  `results/alignment.csv` (`offsetSamples`, `nWindowsAfterOffset`) were added, and a
  participant is accepted only if enough windows remain after the offset.
- Features recomputed from the raw UTD signals for check G7 are now z-scored within each
  participant, like the published feature matrices they are compared with.
- Module docstrings, configuration comments and Makefile help texts now refer to the
  sections of `docs/METHOD.md`.
- New tests for the alignment offset.

## 1.0.0 (23 October 2026)

First version of the code base, written from scratch.

- Readers for the PhysioNet `noneeg` records (WFDB format 16 with checksum
  verification, MIT annotations) and for the expert labels of Zhang et al. (2017).
- Label alignment check, parameter estimation, semi-synthetic generator with five
  artifact classes and a parameterised movement–artifact link, feature sets F0, FA,
  FA4 and FB, leave-one-participant-out models, untrained reference detectors, checks
  G2–G7, grid, compatibility and power analyses, figures F1–F7.
- Reproduction of the real-data comparison on the published features
  (AUPRC 0.4486 → 0.4754, logistic regression).
- MIT licence and CI layout taken over from the author's earlier repository.

### Corrections made during development, before any UTD result

These were found by tests or by the development run on AWW and fixed before the
preregistration tag; none was made after seeing a synthetic result.

- Responses are detected as peaks of the phasic component with a height of at least
  0.05 µS. Counting runs of non-zero driver produced 14–95 responses per minute on AWW,
  because the cvxEDA driver is dense at the optimum on noisy, quantised data.
- τ<sub>0</sub> is chosen by the smallest cvxEDA objective value; the residual norm did not
  recover τ<sub>0</sub> on synthetic records with known τ<sub>0</sub>.
- The tonic given to the generator is fitted to the signal minus the detected
  responses; with the cvxEDA tonic alone, synthetic AWW records were several µS too
  low.
- Movement events: the 0.5-s minimum length is applied before merging; before, single
  noise samples less than 1 s apart merged into spurious events.
- Noise estimate: recordings without clean rest windows fall back to all clean windows;
  artifact classification uses the noise floor max(σ, *q*/√12).
- AWW resting days are treated as rest.
- A donor without clean time of one stage type uses its other response rate instead of
  generating no responses.
- Plausibility warnings for response rates above 10 per minute, τ<sub>0</sub> on the edge
  of its range and solver convergence.
- Logistic regression uses the `newton-cholesky` solver. With `lbfgs` at its default
  tolerance the real-data result was 0.4485 → 0.4739 (Δ +0.025); the converged solution
  is 0.4486 → 0.4754 (Δ +0.027). The reference values of G7 and goal 5 were updated.
- Results of resumable runs are rewritten with aligned columns; appending rows with
  different columns had shifted values.
- Figure F4 is drawn per artifact composition; figure F1 compares a synthetic
  participant with its physiology donor.
- The reference cvxEDA code is no longer listed as a PyPI dependency; it is used when
  `cvxEDA.py` from the authors' repository is on the Python path.
- `Makefile`: `before-tag` and `after-tag` targets.
- Coherence feature B4 is 0 where a segment is constant; validation writes its report
  after every check.
