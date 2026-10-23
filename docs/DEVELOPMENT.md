# Development guide

This page is for people who want to read, change or extend the code. How to run the
study is described in [REPRODUCING.md](REPRODUCING.md), and what each step computes in
[METHOD.md](METHOD.md).

## Architecture

The package `edasynth` is organised by method step. Data flow from left to right:

```
io_wfdb ─┐
         ├─► data ─► alignment ─► estimation ─► bank ─► generator ─► features ─► models
config ──┘                          │                     ▲              │          │
                                decompose                 │              ▼          ▼
                                                    experiment ◄──── detectors   validation
                                                         │
                                                         ▼
                                                      figures        (cli wires everything)
```

| Module | Responsibility | Method section |
|---|---|---|
| `config.py` | loads YAML configs, merges overrides, resolves paths | — |
| `io_wfdb.py` | reads PhysioNet WFDB headers, format-16 signals (with checksum checks) and MIT annotation files | §2.1 |
| `data.py` | `Recording` class; UTD and AWW loaders; protocol stages; expert labels | §2 |
| `alignment.py` | checks that expert windows line up with the raw signal and finds offsets | §2.4 |
| `decompose.py` | Bateman response; cvxEDA decomposition (reference code or built-in solver) | §3.3 |
| `estimation.py` | everything measured on real recordings: noise, responses, tonic, artifact episodes, movement events | §3 |
| `bank.py` | collects the estimates into the parameter bank (`bank_*.json` / `.npz`) | §3 |
| `generator.py` | builds one synthetic participant: clean signal, artifacts, true labels | §4 |
| `features.py` | feature sets F0, FA, FA4, FB | §5 |
| `models.py` | logistic regression and random forest, leave-one-participant-out, metrics, bootstrap | §6 |
| `validation.py` | checks G2–G7 | §7 |
| `experiment.py` | data sets per condition, evaluation, grid, summaries, compatibility, power | §8 |
| `detectors.py` | untrained EDA-only detectors | §9 |
| `realdata.py` | analyses of the real UTD data | [RESULTS.md §1](RESULTS.md#1-real-data-comparison) |
| `figures.py` | figures F1–F7, drawn only from files in `results/` | — |
| `cli.py` | one command per step (`python -m edasynth <command>`) | — |

### Key data structures

- **`Recording`** (`data.py`): one continuous real recording: `eda`, `acc` (n × 3, g),
  optional `temp`, per-sample `stage`, expert labels `labels3` (windows × 3) and
  metadata. `majority()` gives the two-of-three label per window.
- **Bank** (`results/bank_<source>.json` and `.npz`): `phys` (one entry per
  participant: τ<sub>0</sub>, rates, response sizes, tonic fits, σ, *q*), `motion`
  (movement events and quantisation of each movement donor; the acceleration arrays
  are in the `.npz`), `library` (artifact templates by class, class shares, walking
  gain *G*), the full per-participant `estimates`, the configuration path and the git
  commit that produced it.
- **`Condition`** (`generator.py`): one cell of the grid (*c*, κ, δ<sub>max</sub>,
  π<sub>r</sub>, classes, variant, prevalence, label rule). `key()` gives a stable text
  key used to name results and derive seeds.
- **`SynthRecord`** (`generator.py`): one synthetic participant with the observed
  signal `y`, the clean signal `y0`, the label mask, window labels, the list of
  artifacts and responses, and the donors in `meta`.

## Configuration

Every number that can change a result is in `configs/default.yaml`; the code reads it
from the dictionary returned by `config.load_config`. `configs/smoke.yaml` holds the
small settings of the development run and is selected with `--config`. Comments in
the YAML state the origin of each value; [METHOD.md §11](METHOD.md#11-parameter-table)
lists them in one table.

Changing a value in `configs/default.yaml` after the preregistration tag changes the
study. Such a change needs a new tag and a `CHANGELOG.md` entry
([PREREGISTRATION.md](PREREGISTRATION.md)).

## Reproducibility conventions

- Downloaded files in `data/raw/` and `data/external/` are never modified; every step
  writes to `results/` or `figures/`.
- Every random draw comes from a NumPy generator seeded from `configs/*.yaml`.
  The seed of a grid replicate is derived from the CRC32 checksum of the condition
  key and the replicate number, so it does not depend on the order in which jobs run
  or on Python's hash randomisation. The negative control's time shift uses a
  separate random stream, so the shifted and unshifted data sets have identical EDA.
- The same seeds give the same synthetic data. Model results can differ in the last
  digits between machines or library versions.
- Long steps (`pilot`, `grid`) write their results after every batch and skip jobs
  that are already in the output file, so an interrupted run can simply be started
  again.
- The bank records the git commit that produced it. Other result files do not yet;
  record the commit yourself when archiving results.

## Tests

```bash
make test        # unit tests that need no downloaded data (pytest -m "not slow")
make test-all    # all tests, including those that use the downloaded data
```

Tests live in `tests/`. Shared fixtures are in `tests/conftest.py`: `cfg` (the default
configuration) and `toy_bank` (a small, fully synthetic bank of four made-up
participants for exercising the generator and the pipeline). Tests that need the
Zhang et al. repository or the PhysioNet files are marked with `needs_zhang` or
`needs_noneeg` and skip when the data are absent. Cross-checks against optional
packages (`wfdb`, `PyWavelets`, the authors' cvxEDA) skip when those are not
installed. The slow test `test_reproduces_published_feature_baseline` repeats the
real-data comparison.

| File | Covers |
|---|---|
| `test_io_and_data.py` | WFDB reader, annotation reader, label loader, alignment and cropping |
| `test_g1_signal.py` | Bateman shape, templates, labels, donors, Poisson counts, movement–artifact link |
| `test_g1_features.py` | feature counts and order, Haar details, derivatives, FB causality and lag |
| `test_estimation.py` | quantisation, clean windows, noise, episode classes, movement events, τ<sub>0</sub> selection |
| `test_scr.py` | response detection and reconstruction |
| `test_models_and_decompose.py` | leave-one-out isolation, metrics, bootstrap, decomposition, real-data baseline |
| `test_pipeline.py` | data sets, evaluation rows, negative-control shift, detectors, resumable results |

Code style is checked with `ruff check src tests` and formatted with
`ruff format src tests` (settings in `pyproject.toml`). The GitHub Actions workflow
(`.github/workflows/ci.yml`) runs ruff and the data-free tests on Python 3.11, 3.12
and 3.13, and a second job downloads the data and runs all tests.

## Extending the study

**A new feature set.** Add the function to `features.py`, return it from
`experiment.features_for`, and register it in `FEATURE_SETS` (and, for a new gain,
`DELTAS`) in `experiment.py`. Evaluation, summaries and the grid output pick it up
from there. Add a test of its shape and, if it should be causal, of its causality.

**A new artifact class.** Add its measurement and template parameters in
`estimation.classify_episode` and `estimation.template_parameters`, its shape in
`generator.Generator._artifact_signal` (returning the artifact and the mask of
corrupted samples), and the class name to the grid and detector lists in the configs.
Add a test that the template has the requested size and duration and that its label
mask is correct.

**A new untrained detector.** Add a function to `detectors.py` that returns either a
score per window or a 0/1 flag per window and register it in `detectors.run_all`.

**A new command.** Add `cmd_<name>(cfg, args)` to `cli.py` and its parser entry in
`main`; add a Makefile target if it is part of the study.

Any change to the method after the preregistration tag must be recorded in
`CHANGELOG.md` and released under a new tag; results obtained with it are labelled
*post hoc* wherever they are reported.
