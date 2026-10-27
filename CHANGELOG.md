# Changelog

All notable changes to this project are recorded here. Changes to the method after
the preregistration tag must also be released under a new tag
([docs/PREREGISTRATION.md](docs/PREREGISTRATION.md)).

## 1.3.1 (27 October 2026)

First change after the tag `prereg-synth-v1`; released under the tag
`prereg-synth-v2`.

- The label-noise grid (`grid --label-noise`, `make labelnoise`) computes only the
  34 mixture cells. Only these cells enter its analysis (the compatibility with the
  real gain under noisy training labels, goal 5); the 120 single-class cells were
  computed but never analysed. Reason: computing time on the available laptop (about
  two days for the full label-noise grid). No analysis, threshold or reported
  outcome changes. The main grid and the variants were already running with 1.3.0
  and are not affected: the generator, the features and the models are unchanged.

## 1.3.0 (25 October 2026)

Changes before the preregistration tag, made after the UTD checks of 1.2.0 (G4 failed
with AUC 0.86, the C4 control produced no artifacts). From this version on, the UTD
checks were run in the development environment on a copy of the PhysioNet files
(verified with `check-data`).

### Generator

- **Measured courses for C1, C2, C3.** Every episode stores its course relative to
  the level before it (`traceDev`, `traceStep`), and the generator inserts the course
  of a drawn episode, scaled by the local level and the movement scaling. The simple
  shapes are kept for templates without a course. Synthetic artifact windows had a
  median SD of the first difference of 0.005 µS per sample against 0.033 in the real
  ones.
- **Local level of the signal so far.** Relative sizes refer to the signal including
  earlier artifacts; a downward change is limited to 95 % of the signal after
  movement scaling (it reached −10 µS before).
- **Tonic follows the donor's levels.** Each stage runs from the level at its start to
  the end level of the physiology donor's fit of the same stage; before, slopes of
  randomly drawn stages added up and the clean signal fell below 0 µS in 11 of 20
  synthetic participants (down to −6.3 µS). The tonic is floored at 0.05 µS. The tonic
  fits record `durationS`.
- **C4 control gain.** `generator.c4.gain: control` sets the C4 gain so that a
  movement event of median vigour gives a median artifact size of 0.1 µS; the measured
  gain on UTD (−0.0009 µS/g) produced no artifacts.

### Experiment

- Equivalent grid cells (κ and δ<sub>max</sub> at *c* = 0; κ in the C5-only cells)
  are computed once and copied (`sharedWith`): 154 instead of 200 cells.
- Power is computed for the compatible mixture cells at the quartiles of mean
  Δ<sub>B</sub> instead of the first three in grid order.
- `analyze` writes the compatibility with noisy training labels when the label-noise
  grid exists.
- `loadBank` takes the configuration (needed for the C4 control gain).

### Documentation

- `docs/PREREGISTRATION.md` lists the decisions taken after the first UTD checks,
  the accepted G5 exception and what was kept unchanged on purpose.

## 1.2.1 (25 October 2026)

- `pilot` and `grid` write the checksum of the parameter bank (`bankId`) and the
  version into every row and refuse to resume a file made with another bank. Before,
  a rerun of `pilot` after `estimate` with 1.2.0 skipped all 40 replicates because
  the file of 1.1.2 was still there, and printed the old *R* = 50. Found on the
  first UTD run of 1.2.0; a new test covers it.
- Documentation: `docs/RESULTS.md` §3.5 reports the UTD checks with 1.2.0 (G1, G2,
  G3 and G7 pass; G4, G5 and G6 fail) and the pilot (*R* = 46); `docs/STATUS.md`
  updated.

## 1.2.0 (25 October 2026)

Changes to the generator after the first run of the checks on the raw UTD files
(version 1.1.2, [docs/RESULTS.md §3](docs/RESULTS.md#3-first-run-on-the-raw-utd-files-version-112)).
They were made before the preregistration tag; the decisions and the results that
led to them are documented. Banks written by earlier versions are not read; run
`estimate` again.

### Generator and estimation

- **Noise from the decomposition residual.** On UTD the noise estimate σ was 0 for
  every participant (the noise is below the quantisation step), so the generator
  added no noise and the synthetic signal was too smooth (G3 AUC 0.98). The noise is
  now made of 5-s blocks of the physiology donor's residual (signal minus cvxEDA
  tonic and phasic parts) from clean rest windows, centred and drawn with
  replacement. The blocks are stored in the bank (`resid<pid>` in the `.npz`).
  Gaussian noise remains available as `generator.noise: gaussian`.
- **Relative artifact sizes.** Episode sizes are stored relative to the level before
  the episode, divided by max(level, 0.05 µS), and the generator multiplies them by
  the local level of the synthetic signal. Template keys are now `relA`, `relDepth`
  and `relDelta`; episodes record `level` and `levelAfter`. Episodes that start
  below 0.05 µS are left out of the templates and class shares (their number is
  logged and stored as `library.episodesBelowFloor`).
- **Physiology donors.** A recording with more than half of its samples below
  0.05 µS (Kleckner et al., 2018, rule 1; the 50 % limit is a decision) is not used
  as a physiology donor; it remains a movement donor. On UTD this concerns S18.
  `generator.assignDonors` replaces the derangement when donors are excluded; with
  all donors valid it is the same derangement as before.
- **Response sizes only from clean windows.** Sizes were taken from every detected
  response, including those inside artifact windows; on UTD 25 of 1 092 heights
  exceeded the range of the recording (up to 31.7 µS), and on the AWW run a synthetic
  response of about 11 µS appeared in figure F2. Now only responses with onset and
  peak in clean windows supply sizes, and heights above the range of the recording
  are reported as warnings.
- New configuration keys: `estimation.episodes.levelFloor`, `estimation.donor`,
  `estimation.noise.blockS`, `generator.noise`, `generator.levelWindowS`.
- The seeds are unchanged, but the synthetic data differ from 1.1.2 because the
  noise and the sizes are drawn differently.

### Decisions recorded

- G5 keeps its rule: a significant loss without a movement–artifact link also fails
  it. FA and FA4 gains are to be read against their no-link loss
  ([docs/METHOD.md §7](docs/METHOD.md#7-validation-checks-g1g7)).

### Documentation

- `docs/RESULTS.md` §3 reports the first UTD run: inspection, alignment,
  estimation, G1–G7 and the pilot; §1 adds the G7 comparison on recomputed
  features and the random-forest values measured on macOS.
- `docs/RESULTS.md` §2 reports the AWW development run with 1.2.0: G3 passes
  (AUC 0.59), G4 still fails (0.93), and G2 now fails because τ<sub>0</sub> is not
  recovered when responses are small compared with the residual noise.
- `docs/STATUS.md`, `docs/METHOD.md` (§3.2, §3.6, §4.1–4.5, §7, §10, §11) and
  `docs/REPRODUCING.md` (measured computing times) were updated.

## 1.1.2 (23 October 2026)

- WFDB reader: checksums are compared modulo 2<sup>16</sup>. The PhysioNet noneeg
  headers store them as unsigned numbers (57979), while the reader computed the signed
  form of the same value (−7557) and rejected every record. Found by `make test-all`,
  the first run of the reader on the real files; a new test covers headers with
  unsigned checksums.

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
