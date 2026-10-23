# edasynth

**Does a wrist accelerometer help to detect motion artifacts in skin-conductance
recordings, and under which conditions?**

`edasynth` is a reproducible research code base that answers this question with
*semi-synthetic* data: real wrist-movement recordings from 20 people are combined
with computer-generated skin-conductance signals whose artifacts, and therefore
whose correct labels, are known exactly. The code also re-analyses a public data
set in which experts labelled such artifacts.

> **Current state (October 2026).** The code for every step is written and tested.
> One result on real data is reproduced and verified (see [Results](docs/RESULTS.md)).
> The main synthetic experiment has **not been run yet**: it needs the raw PhysioNet
> recordings, and a development run on a second data set (recordings of one person)
> failed all realism and control checks, which raises design questions that must be settled first
> ([Status](docs/STATUS.md)). No synthetic finding is reported in this repository.

---

## The problem in brief

Electrodermal activity (EDA), also called skin conductance, measures how easily a
very small electric current passes through the skin. It changes with sweat-gland
activity, which is controlled by the sympathetic nervous system, so EDA is used to
follow stress and arousal, increasingly with wristbands
(Posada-Quintero & Chon, 2020).

When the wrist moves, the electrodes can slide, press into the skin or lose contact.
The recording then shows peaks, drops or steps that do not come from the body. These
are **motion artifacts**. They have to be found before the signal is interpreted,
and some of them look like genuine skin-conductance responses
(Taylor et al., 2015).

Many wristbands also contain an **accelerometer**, a sensor that records how the
wrist accelerates. Using it to find motion artifacts seems obvious. However, when
Zhang et al. (2017) added accelerometer features to artifact detectors trained on
expert-labelled wrist recordings, the supervised detectors improved by only about
0.4 % on average, and our re-analysis of the same data finds a gain that cannot be
distinguished from zero (see [Results](docs/RESULTS.md)). Two explanations are
possible: the acceleration features may be poorly suited to the task, or many
movements may simply not disturb the electrode (Zhang et al., 2017). The available
data cannot tell these apart: no data set records how often movement actually
disturbs the electrode, and one of the experts' labelling criteria was "a sudden
change correlated with motion", so the labels themselves depend on movement.

A longer introduction for readers new to the field is in
[docs/BACKGROUND.md](docs/BACKGROUND.md).

## What the project does

Because the decisive quantity (how strongly movement is linked to artifacts) is
unknown, the project does not try to give a single yes-or-no answer. It builds a
**conditional map**: the movement–artifact link is set to many different values,
and for each one the gain from adding the accelerometer is measured with a
confidence interval. Every statement on the map has the form *"if movement causes
artifacts in this way, the accelerometer improves detection by this much"*.

The work is organised around six goals, each answering one question:

| Goal | Question | Main output |
|---|---|---|
| 1 | Is the analysis pipeline correct? It must find no gain where none was built in and a gain where one was. | validation report, checks G1–G7 |
| 2 | Which artifact shapes do common EDA-only methods miss? | detection rate by artifact class (figure F7) |
| 3 | How much does the accelerometer help under each assumed link? | the map (figures F4, F5) |
| 4 | What does it cost? Real responses that occur during movement may be discarded as artifacts. | figure F6 |
| 5 | Which assumed links are consistent with the small gain seen on real data? | compatibility table |
| 6 | Which quantity should a laboratory study measure first, and with how many participants? | sensitivity and power tables |

Goal 1 comes first: no other result is interpreted unless the checks pass.

## How it works

```
 Real data (PhysioNet "noneeg", 20 people, 8 Hz)
 wrist acceleration + EDA + expert artifact labels
            │
            ▼
 Measure what the literature does not provide        (docs/METHOD.md §3)
 noise, response rate and size, slow baseline drift,
 shapes of labelled artifacts, movement events
            │
            ▼
 Build synthetic participants                         (§4)
 acceleration of person A (real, unchanged)
 + clean EDA generated with the physiology of person B
 + artifacts placed according to the assumed movement–artifact link
            │
            ▼
 Compute features for every 5-second window           (§5)
 EDA only │ EDA + 96 acc. │ EDA + 4 acc. │ EDA + 4 "coupling" features
            │
            ▼
 Train and test models, leaving one participant out   (§6)
 gain = AUPRC(with acceleration) − AUPRC(EDA only)
            │
            ▼
 Checks G1–G7 (§7)  →  grid of assumed links (§8)  →  map, costs, compatibility, power
```

Taking movement from one person and physiology from another means that any link
between the two in the synthetic data exists only because the generator put it
there, so its strength is known exactly. AUPRC, the main score, is explained in
[docs/METHOD.md §6](docs/METHOD.md#6-models-and-performance-measures).

## Results so far

Only one result is final, and it is on real data:

| Model (leave-one-participant-out) | EDA only (24 features) | EDA + acceleration magnitude (28) | Gain | 95 % CI of the gain |
|---|---|---|---|---|
| Logistic regression | 0.449 | 0.475 | +0.027 | −0.002 to +0.058 |
| Random forest | 0.464 | 0.490 | +0.026 | −0.021 to +0.056 |

AUPRC on the published UTD features (9 211 five-second windows, 349 labelled as
artifact by at least two of three experts). The confidence intervals include zero:
**these data neither show a gain nor show that there is none.** Details, controls
and what may and may not be concluded are in [docs/RESULTS.md](docs/RESULTS.md).

The synthetic experiment has not produced results yet. A development run on a
second data set, recorded from one person, shows that every step executes, and also that the
synthetic data are not yet realistic enough; see the same file.

## Getting started

Requirements: Python 3.11 or newer, `git`, `make`, `curl`, `unzip` and about
200 MB of disk space. The unit tests pass on Linux (Python 3.13) and on macOS with
Apple Silicon (Python 3.12, conda); Windows has not been tested.

```bash
git clone https://github.com/cbrg0/edasynth.git
cd edasynth

# either conda (recommended on Apple Silicon)
conda env create -f environment.yml && conda activate edasynth
# or a plain virtual environment
python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"

make test        # unit tests, about 1.5 min, no data needed
make data        # download both data sets and verify their checksums
```

What to run next depends on the aim:

| Aim | Command | Time (measured or estimated) |
|---|---|---|
| Check that everything runs | `make smoke` | about 50 min on 2 cores (measured) |
| Reproduce the real-data result | `python -m edasynth realdata --stored` | about 10 min on 2 cores (measured) |
| Run the first half of the study (data checks, parameter estimation, validation, pilot) | `make before-tag` | several hours (estimated) |
| Run the main experiment | `make after-tag` | depends on the pilot; see the guide |

The full study is split in two on purpose: the analysis rules must be frozen with a
public git tag *before* the main experiment runs. Step-by-step instructions,
expected outputs and troubleshooting are in
[docs/REPRODUCING.md](docs/REPRODUCING.md).

## Repository layout

```
README.md                  this page
docs/                      all further documentation (see the map below)
configs/default.yaml       every number that can change a result, with its source
configs/smoke.yaml         small settings for the development run
src/edasynth/              the Python package (one module per method step)
tests/                     unit tests; tests that need downloaded data skip without it
scripts/download-data.sh   downloads both data sets at fixed versions, verifies checksums
data/                      downloaded data (not committed) and the checksum list
results/, figures/         generated outputs (not committed)
Makefile                   one target per step of the study (`make help` lists them)
```

## Main limitations

- Synthetic results are conditional. They show what would happen *if* movement
  disturbed the electrode in a given way; they cannot show how often it does so in
  people.
- All parameters come from one device (Affectiva Q Sensor, according to Zhang et
  al., 2017) and one group of 20 healthy participants in a laboratory. Whether the
  conclusions transfer to other sensors or to daily life is not tested.
- The expert labels partly depend on the movement trace, so they cannot serve as an
  independent ground truth for judging the accelerometer.
- The steps that need the raw PhysioNet recordings have not been run; the synthetic
  data have so far failed all realism checks on the development data.

A complete list, including open design decisions, is in
[docs/STATUS.md](docs/STATUS.md).

## Documentation map

| File | Read it to … |
|---|---|
| [docs/BACKGROUND.md](docs/BACKGROUND.md) | understand EDA, accelerometers, motion artifacts and the research gap |
| [docs/METHOD.md](docs/METHOD.md) | see every step, formula, parameter and its source |
| [docs/RESULTS.md](docs/RESULTS.md) | see what has been measured and what it means |
| [docs/STATUS.md](docs/STATUS.md) | see what is implemented, tested, run, or still open |
| [docs/REPRODUCING.md](docs/REPRODUCING.md) | install, download data, run every step, check outputs |
| [docs/PREREGISTRATION.md](docs/PREREGISTRATION.md) | see the analysis rules fixed before the main experiment |
| [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) | work on the code: architecture, tests, conventions |
| [docs/GLOSSARY.md](docs/GLOSSARY.md) | look up a term |
| [docs/REFERENCES.md](docs/REFERENCES.md) | find every source and what it is used for |
| [docs/README-RU.md](docs/README-RU.md) | read a short overview in Russian |
| [CHANGELOG.md](CHANGELOG.md) | see what changed between versions |

## Licence and citation

The code is released under the MIT licence ([LICENSE](LICENSE)). The data are not
redistributed here: PhysioNet `noneeg` is available under the Open Data Commons
Attribution License v1.0, and the repository of Zhang et al. under BSD.

To cite this software, use [CITATION.cff](CITATION.cff). Please also cite the data
and the methods it builds on, in particular Birjandtalab et al. (2016), the
PhysioNet data set (Birjandtalab et al., 2017), Greco et al. (2016) and Zhang et al.
(2017). Full references are in [docs/REFERENCES.md](docs/REFERENCES.md).
