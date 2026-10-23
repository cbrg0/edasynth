# edasynth: does a wrist accelerometer help find motion artifacts in skin-conductance recordings?

A reproducible, semi-synthetic test bench. Real wrist-movement recordings
(20 people, PhysioNet) are combined with synthetic electrodermal activity (EDA) whose
clean part, artifacts and true labels are known exactly. The bench maps **under which
assumed links between movement and electrode contact** adding the accelerometer
improves artifact detection, how large the gain is, what it costs, and what a
laboratory study has to measure to locate the real world on that map.

> **What this repository does not claim.** Synthetic data cannot show that the
> accelerometer helps (or does not help) for real people, because the decisive
> quantity, how often wrist movement actually disturbs the electrode, has not been
> measured in any available data set. Every synthetic result here is conditional:
> "*if* the link is like this, *then* the gain is that". The real-data analysis below
> is the only unconditional result, and it is inconclusive.

---

## Three ways into this repository

| You are | Read first | Then |
|---|---|---|
| **Biologist / physiologist** | [In plain words](#in-plain-words) and [docs/GLOSSARY.md](docs/GLOSSARY.md) | [What the results can and cannot say](#what-the-results-can-and-cannot-say) |
| **Supervisor / reviewer** | [Research design](#research-design) and [Status](#status) | [docs/METHOD.md](docs/METHOD.md), [docs/PREREGISTRATION.md](docs/PREREGISTRATION.md), [docs/STATUS.md](docs/STATUS.md) |
| **Programmer** | [Quick start](#quick-start) and [Repository map](#repository-map) | [Reproducibility rules](#reproducibility-rules), `make help`-style comments in the [Makefile](Makefile) |

---

## In plain words

**EDA (skin conductance)** is how easily a tiny current passes through the skin. Sweat
glands, steered by the sympathetic nervous system, change it, so EDA is used to follow
stress and arousal, increasingly with wristbands.

**The problem.** When the wrist moves, the electrode can slide, be pressed or lose
contact. The recording then shows peaks, drops or steps that are not physiology:
*motion artifacts*. They must be found before the signal is interpreted.

**The obvious idea** is to use the accelerometer that most wristbands already contain:
if the wrist moved, the change is suspicious. But movement does not always disturb the
electrode, and movement can also cause *real* sweat responses. In the only published
test of this idea (Zhang et al., 2017) the gain was tiny, and in our own re-analysis of
the same data it was not distinguishable from zero.

**What we do.** We build artificial recordings where we decide how often movement
causes artifacts. Because we made the artifacts, we know the right answer for every
5-second piece. We then ask a computer model to find the artifacts with and without
the accelerometer and measure the difference, for every assumed link between movement
and artifacts. The result is a map, not a yes/no answer.

**Why it is still useful.** The map shows (1) whether the analysis pipeline is correct
(it must find no gain where we planted none, and a gain where we planted one),
(2) which artifact shapes the usual methods miss, (3) how much the accelerometer could
help and what it costs (real responses during movement wrongly discarded), and
(4) which quantity a laboratory study must measure first and with roughly how many
people.

---

## Research design

**Question.** Does a synchronous wrist accelerometer improve the detection of motion
artifacts in wrist EDA compared with a model that uses EDA only?

**Why a conditional map.** Zhang et al. (2017) offered two explanations for the small
gain: (1) the acceleration features may be unsuited to the task; (2) many movements do
not change electrode contact. These become the two axes of the map:

* **Axis 1, kind of features:** none (F0, 24 EDA features), 96 acceleration statistics
  (F0+FA, Zhang's setting), 4 acceleration statistics (F0+FA4, our real-data analysis),
  4 *coupling* features describing whether EDA and movement change together (F0+FB).
* **Axis 2, movement-contact link:** share of artifacts caused by movement c, dependence
  of artifact size on movement intensity κ, maximum delay δ_max, share of movements
  followed by a real response π_r, and artifact class.

**Goals derived from what synthetic data can and cannot do**

| # | Goal | Output |
|---|---|---|
| 1 | Verify the pipeline: correctness, no leakage, finds a planted effect | gates G1-G7 (`results/validation_utd.json`) |
| 2 | Blind spots of EDA-only methods by artifact class | F7, `results/detectors_by_class.csv` |
| 3 | Conditional map of the accelerometer's benefit | F4, F5, `results/grid_summary_utd.csv` |
| 4 | Cost: real responses during movement flagged as artifacts | F6 |
| 5 | Which conditions are compatible with the small real gain on UTD | `results/compatibility_utd.csv` |
| 6 | What to measure first in the two-hand protocol and how many participants | `results/sensitivity_utd.csv`, `results/protocol_planning_utd.csv` |

Goal 1 gates all others. Full method: [docs/METHOD.md](docs/METHOD.md). Decision
rules fixed in advance: [docs/PREREGISTRATION.md](docs/PREREGISTRATION.md).

---

## Status

Honest status of every element, updated on each run: [docs/STATUS.md](docs/STATUS.md).
Summary at the time of the last commit:

* **Run and verified:** unit tests (52 passed, 5 skipped for missing optional
  packages or data); reproduction of the real-data baseline from the published UTD
  features (see below); end-to-end smoke run of every step on AWW (development
  data, numbers not reportable). In the smoke run all gates G2-G6 **fail**; the
  reasons and the design questions they raise are listed in docs/STATUS.md.
* **Not yet run:** everything that needs the raw UTD signals (alignment, parameter
  estimation, gates G2-G7, the grid). The environment in which this repository was
  built could not reach PhysioNet; `make all` on a normal internet connection runs
  these steps.

### Real-data result reproduced here (published UTD features)

`python -m edasynth realdata --stored` (logistic regression, standardised, balanced
classes, leave-one-participant-out, 9 211 windows, 349 artifacts):

| | EDA only (24) | EDA + \|a\| (28) | Δ AUPRC | 95 % CI (participant bootstrap) | per-participant Wilcoxon |
|---|---|---|---|---|---|
| Logistic regression | 0.4486 | 0.4754 | +0.027 | −0.002 to +0.058 | p = 0.76 |
| Random forest (500 trees) | 0.4641 | 0.4900 | +0.026 | −0.021 to +0.056 | p = 0.47 |

Shuffling the acceleration rows within each participant gives a mean Δ of −0.001
(20 shuffles, logistic regression). The confidence interval includes zero: the data neither show a gain nor
show its absence. An earlier analysis reported Δ = 0.054 for the random forest; with
the settings documented here we obtain 0.026 (the earlier settings were not recorded).

---

## What the results can and cannot say

* A positive cell on the map means: *if* movement disturbs the electrode in that way,
  the accelerometer helps by that much. It does not say that it does so in people.
* A cell "compatible with UTD" means the small real gain could have arisen under that
  condition. Several conditions will usually be compatible; then UTD cannot tell them apart.
* UTD labels were given by experts who looked at the movement trace (one of their
  criteria is "sharp change associated with movement"), so they are not independent of
  the accelerometer. The planned two-hand laboratory protocol is designed to reduce this
  dependence; even there, the still hand is a reference, not a perfect ground truth.
* All numbers inherit the properties of one device (Affectiva Q Sensor) and one group
  (20 students, laboratory). Transfer to other sensors or daily life is not established.

---

## Quick start

Requirements: Python ≥ 3.11, git, curl, unzip, ~200 MB disk.

```bash
git clone https://github.com/cbrg0/edasynth.git
cd edasynth

# Option A (recommended on Apple Silicon): conda-forge has a prebuilt cvxopt
conda env create -f environment.yml && conda activate edasynth
# Option B: plain venv (cvxopt optional; a dependency-free solver is built in)
python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev,reference]"

make data      # downloads PhysioNet noneeg 1.0.0 and Zhang et al. @7c4bd46, verifies SHA-256
make test      # G1 unit tests
make smoke     # ~50 min end-to-end run on AWW (2 cores): proves everything runs (NOT results)
make before-tag  # data, tests, alignment, estimation, gates G1-G7, pilot
# settle the open questions in docs/STATUS.md, then: git tag prereg-synth-v1 && git push --tags
make after-tag   # grid, variants, analysis, power, figures (long: see "Compute budget")
```

Optional cross-check against the authors' cvxEDA (Greco et al., 2016): the reference
code is a single file that is not on PyPI under a name we could verify, so it is not
installed automatically. To use it, download `cvxEDA.py` from
https://github.com/lciti/cvxEDA into a folder on your Python path (with `cvxopt`
installed); `decompose` then uses it and the agreement test runs instead of skipping.

Every step is also a single command, e.g. `python -m edasynth align --plots`
(`python -m edasynth --help` lists them).

### Compute budget

| Step | Approximate cost |
|---|---|
| `estimate` with the reference cvxEDA on the path | not measured |
| `estimate` with the built-in numpy solver | measured on AWW: ~35 s per τ0 value for one hour of signal on one core; UTD (20 × 21 τ0 values, ~38 min each) ≈ 2-3 h, `--jobs` splits it over cores |
| one replicate of one grid cell (4 feature sets, logistic regression, inner C search) | measured on AWW (11 donors, inner 3-fold): ~45 s on one core; UTD with inner 5-fold and 4 C values is expected to take a few minutes |
| full grid: 200 cells × R replicates | R × 200 × ~1 min / JOBS cores (resumable: rerun continues where it stopped) |

---

## Repository map

```
configs/default.yaml        every number that can change a result (frozen by the prereg tag)
configs/smoke.yaml          tiny settings for the development run on AWW
data/                       downloaded data (not committed) + checksum list + README
scripts/download_data.sh    fixed-version download with checksum verification
src/edasynth/
  io_wfdb.py                PhysioNet reader (format 16, checksums) and MIT annotations
  data.py                   UTD and AWW loaders, protocol stages, expert labels
  alignment.py              labels ↔ raw signal check (step 2)
  decompose.py              Bateman response; cvxEDA (reference package or numpy solver)
  estimation.py             noise, quantisation, responses, tonic, artifact episodes, motion events (step 3)
  bank.py                   writes the estimated parameters ("bank") with provenance
  generator.py              semi-synthetic participants, artifact classes, true labels (steps 4, 6)
  features.py               F0, FA, FA4, FB (step 5)
  models.py                 LOSO logistic regression / random forest, AUPRC, bootstrap (step 7)
  detectors.py              untrained EDA-only detectors (goal 2)
  validation.py             gates G2-G7 (step 8)
  experiment.py             data sets, grid, pilot, compatibility, power (step 9)
  realdata.py               real-data baseline (stored and recomputed features)
  figures.py                F1-F7 from results files only (step 10)
  cli.py                    the commands used by the Makefile
tests/                      G1 unit tests; data-dependent tests skip when data are absent
docs/                       METHOD, PREREGISTRATION, STATUS, GLOSSARY
results/, figures/          outputs (not committed; regenerated by make)
```

## Reproducibility rules

* Downloaded files are never modified; every step writes to `results/`.
* Every random operation draws from a seed in the config; reruns are bit-identical.
* Every results file records the git commit that produced it (`bank_*.json`).
* Figures are produced only by code from results files.
* Thresholds of the gates and the decision rule are fixed before the grid by the tag
  `prereg-synth-v1`; changes after it get a new tag and an entry in `CHANGELOG.md`.

## Licences

Code: MIT ([LICENSE](LICENSE)). Data are not redistributed: PhysioNet `noneeg` is
under ODC-By 1.0 (cite Birjandtalab et al., 2016, the data set and PhysioNet);
the Zhang et al. repository is BSD. The optional reference cvxEDA code is not copied into this repository;
check its licence in its own repository before using it.

## How to cite

See [CITATION.cff](CITATION.cff). Please also cite the data and methods this work builds on:

Birjandtalab, J., Cogan, D., Baran Pouyan, M., & Nourani, M. (2016). A non-EEG biosignals dataset for assessment and visualization of neurological status. In *2016 IEEE International Workshop on Signal Processing Systems (SiPS)* (pp. 110–114). IEEE. https://doi.org/10.1109/SiPS.2016.27

Birjandtalab, J., Cogan, D., Baran Pouyan, M., & Nourani, M. (2017). *Non-EEG dataset for assessment of neurological status* (Version 1.0.0) [Data set]. PhysioNet. https://doi.org/10.13026/C26Q2Z

Greco, A., Valenza, G., Lanata, A., Scilingo, E. P., & Citi, L. (2016). cvxEDA: A convex optimization approach to electrodermal activity processing. *IEEE Transactions on Biomedical Engineering, 63*(4), 797–804. https://doi.org/10.1109/TBME.2015.2474131

Zhang, Y., Haghdan, M., & Xu, K. S. (2017). Unsupervised motion artifact detection in wrist-measured electrodermal activity data. In *Proceedings of the 2017 ACM International Symposium on Wearable Computers* (pp. 54–57). Association for Computing Machinery. https://doi.org/10.1145/3123021.3123054

Full reference list with the role of each source: [docs/METHOD.md](docs/METHOD.md) and [docs/REFERENCES.md](docs/REFERENCES.md).
