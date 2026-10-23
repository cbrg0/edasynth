# Reproducing the study

This guide explains how to install the package, download the data, check the
installation and run every step of the study. What each step computes is explained
in [METHOD.md](METHOD.md); which steps have already been run is listed in
[STATUS.md](STATUS.md).

## 1. Requirements

- Python 3.11 or newer
- `git`, `make`, `curl`, `unzip` and `shasum` (used by the download script)
- about 200 MB of disk space for data and results
- for the full study: several CPU cores and, depending on the number of replicates,
  from tens to several hundred CPU hours (see §5.3)

The code was developed and tested on Linux with Python 3.13. It has not been tested on
macOS or Windows; on Windows, use WSL because the Makefile and the download script
are written for a Unix shell.

## 2. Installation

```bash
git clone https://github.com/cbrg0/edasynth.git
cd edasynth
```

**Option A, conda** (recommended on Apple Silicon, because it installs a prebuilt
`cvxopt`):

```bash
conda env create -f environment.yml
conda activate edasynth
```

**Option B, a plain virtual environment:**

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

`pip install -e ".[dev,reference]"` additionally installs `cvxopt`, `wfdb` and
`PyWavelets` for the optional cross-checks (§6). On some systems `cvxopt` has to be
compiled and the installation fails; the study runs without these packages.

**Check the installation:**

```bash
python -m edasynth --help     # lists the commands
make test                     # unit tests that need no data; about 1.5 min
```

Some tests are reported as skipped. This is expected: they need downloaded data or
optional packages.

## 3. Data

```bash
make data
```

The script `scripts/download-data.sh` downloads:

| Data | Source | Saved in |
|---|---|---|
| PhysioNet *Non-EEG Dataset* 1.0.0 (UTD recordings) | `https://physionet.org/content/noneeg/get-zip/1.0.0/` | `data/raw/noneeg-1.0.0/` |
| Repository of Zhang et al. (2017) at commit `7c4bd467bec9099669ddf3581b00a9a843a2f136` (labels, feature matrices, AWW) | `https://github.com/IdeasLabUT/EDA-Artifact-Detection` | `data/external/eda-artifact-detection/` |

Both are verified: PhysioNet files against the `SHA256SUMS.txt` that ships with them,
the Zhang et al. files against `data/checksums/zhang-eda-artifact-detection-7c4bd46.sha256`.
The script stops if a file does not match, and writes `data/PROVENANCE.txt` with the
download date. `python -m edasynth check-data` repeats the verification later.

The downloaded files are never modified; every step writes to `results/` and
`figures/`.

## 4. Quick checks

| What | Command | Needs | Time |
|---|---|---|---|
| Unit tests | `make test` | nothing | about 1.5 min |
| All tests, including data-dependent ones | `make test-all` | `make data` | a few minutes more |
| Reproduce the real-data comparison | `python -m edasynth realdata --stored` | Zhang et al. files | about 10 min on 2 cores |
| Run every step on the development data | `make smoke` | Zhang et al. files | about 50 min on 2 cores |

The real-data comparison should print, for the logistic regression,
`AUPRC 0.4486 -> 0.4754, delta +0.0269, 95% CI [-0.0022, +0.0580]`
([RESULTS.md §1](RESULTS.md#1-real-data-comparison)).

`make smoke` writes to `results/smoke/` and `figures/smoke/`. These outputs show that
the code runs; they are **not results** (AWW contains one person).

## 5. Running the study

### 5.1 Overview

The study is split into two halves around the preregistration tag
([PREREGISTRATION.md](PREREGISTRATION.md)):

```bash
make before-tag                          # data, tests, data checks, estimation, validation, pilot
# settle the open decisions in docs/STATUS.md, update docs/PREREGISTRATION.md, then:
git tag prereg-synth-v1 && git push origin prereg-synth-v1
make after-tag                           # grid, variants, label noise, analysis, power, figures
```

The number of CPU cores used by the long steps is set with `JOBS`, for example
`make after-tag JOBS=8` (default 4).

### 5.2 Steps, outputs and what to check

| Step | Command | Output in `results/` | What to check |
|---|---|---|---|
| Inspect raw files | `make inspect` | `inspect-utd.json` | EDA range plausible for µS; the `stageError` column empty for every participant |
| Align labels and signals | `make align` | `alignment.csv`; plots in `figures/alignment/` | `accepted` for each participant; the log line `accepted N/20` |
| Real-data comparison | `make realdata` | `realdata.json` | logistic regression as in §4; `G7` entry for the recomputed features |
| Estimate parameters | `make estimate` | `bank-utd.json`, `bank-utd.npz` | `WARNING` lines in the log (response rates, τ<sub>0</sub> at the edge, solver convergence) |
| Validation checks | `make validate` | `validation-utd.json` | `PASS` or `FAIL` for G1–G7 at the end of the log |
| Untrained detectors | `make detectors` | `detectors-by-class.csv` | table printed at the end |
| Pilot | `make pilot` | `pilot-utd.csv`, `pilot-utd-R.json` | the number of replicates *R* |
| Main grid | `make grid` | `grid-metrics-utd.csv` | progress lines `k/N`; can be restarted |
| Variants | `make variants` | `variants-metrics-utd.csv` | — |
| Label noise | `make labelnoise` | `grid-labelnoise-metrics-utd.csv` | — |
| Analysis | `make analyze` | `grid-summary-utd.csv`, `compatibility-utd.csv`, `sensitivity-utd.csv` | number of compatible cells in the log |
| Power | `make power` | `protocol-planning-utd.csv` | — |
| Figures | `make figures` | `figures/F1…F7` (PDF and PNG) | look at F1 first: the real and synthetic recordings should look alike |

`align` must run before `realdata --raw`, `estimate`, `validate` and `figures`,
because these use only the participants accepted by the alignment check, with the
offset it found. The Makefile runs the steps in this order.

### 5.3 Computing time

| Step | Measured | Estimate for UTD |
|---|---|---|
| `estimate` with the built-in solver | about 35 s per τ<sub>0</sub> value for one hour of signal on one core (AWW) | 20 participants × 21 values of τ<sub>0</sub> × about 38 min of signal: roughly 3 h on one core |
| one grid cell, one replicate, default model settings | 66 s on one core with 11 participants (AWW) | roughly 2–3 min with 20 participants |
| main grid | — | 200 cells × *R* replicates × 2–3 min: about 70 h of CPU time for *R* = 10 and about 350 h for *R* = 50 |

The UTD estimates are extrapolations from the development data and may be off by a
factor of two or more. The main grid can be split over several sessions because it
resumes where it stopped. If the cost is too high, reducing the search for the
regularisation strength (`models.logreg.CSelection: fixed`) shortens each replicate
considerably; such a change must be made before the tag.

## 6. Optional cross-checks

| Cross-check | How to enable | Test that runs |
|---|---|---|
| Reader against the `wfdb` package | `pip install wfdb` | `testRecordReaderMatchesWfdb`, `testAnnotationReaderMatchesWfdbOnRealData` |
| Haar details against PyWavelets | `pip install PyWavelets` | `testHaarMatchesPywaveletsIfInstalled` |
| Built-in cvxEDA solver against the authors' code | install `cvxopt`, download `cvxEDA.py` from https://github.com/lciti/cvxEDA into a folder on your Python path | `testNumpyAndReferenceBackendsReachTheSameObjective` |

When the authors' cvxEDA code is importable, `decompose` uses it by default
(`estimation.cvxeda.backend: auto`). Its licence is stated in its own repository; it
is not copied into this one.

## 7. Troubleshooting

| Problem | Cause and remedy |
|---|---|
| `make data` fails while downloading from PhysioNet | The site may be unreachable from your network. Download the ZIP of version 1.0.0 in a browser from the PhysioNet page of the data set, unpack it into `data/raw/noneeg-1.0.0/` and run `make data` again; it will only verify. |
| `shasum: command not found` | Install Perl's `shasum`, or replace `shasum -a 256` by `sha256sum` in the script. |
| `results/alignment.csv not found` | Run `make align` before `estimate`, `validate`, `realdata --raw` or `figures`. |
| `stages could not be parsed` | The annotation file has an unexpected number of entries. `make inspect` shows the codes; the stage list is `utdStages` in `src/edasynth/data.py`. |
| `pip install` fails on `cvxopt` | Install without the `reference` extra, or use conda. |
| A grid run was interrupted | Run the same command again; finished replicates are skipped. |
| Results differ slightly from another machine | Library versions and numerical libraries can change the last digits of model results; the synthetic data themselves are identical for the same seeds. |
