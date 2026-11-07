# edasynth

[![tests](https://github.com/cbrg0/edasynth/actions/workflows/ci.yml/badge.svg)](https://github.com/cbrg0/edasynth/actions/workflows/ci.yml) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cbrg0/edasynth/blob/main/run-in-browser/quickstart.ipynb)

**Can the accelerometer of a wristband help to recognise skin-conductance artifacts
caused by movement?**

A Russian version of this page is in [study-description/README-RU.md](study-description/README-RU.md).

## About the project

Electrodermal activity (EDA) is the change in the electrical conductance of the skin
that sweat glands produce, and it is used to follow emotional arousal and stress.
Wristbands can record it in daily life, but there the signal changes for three
reasons: an emotional reaction, sweating that cools the body, and a poor contact
between the skin and the electrodes. A change of the last kind is called an artifact.
The long-term aim of this project is an algorithm that tells the three sources apart
with the help of other sensors in the same wristband, such as an accelerometer, a skin
thermometer or a heart-rate sensor.

This repository covers the first step. It asks whether the accelerometer helps to find
artifacts and under which conditions it does.

## Result

On public recordings of 20 people whose artifacts were marked by experts (Zhang et
al., 2017), adding acceleration improves the detection of artifacts so little that the
improvement may be zero. Real recordings cannot explain this, because nobody recorded
which artifacts were caused by movement. A computer model of wristband recordings was
therefore built from real material, and in it the cause of every artifact is known.

![Gain from adding acceleration against the share of artifacts caused by movement](result-figures/F0-summary.png)

The figure shows how much the detection improves when acceleration is added (vertical
axis) depending on the share of artifacts that movement causes (horizontal axis). The
accelerometer is useless when movement causes fewer than half of the artifacts. When
movement causes all of them, the most useful are four features that check whether the
skin signal changes right after a movement. The 96 acceleration features used in
earlier work make detection worse in this figure and help in only 2 of the 300
scenarios. The benefit of the accelerometer therefore depends
on a share that no public data set contains, and the next step is a laboratory
experiment that measures it.

## How it was done

1. Real recordings were used to measure how the skin signal behaves, how people move
   and what real artifacts look like.
2. From these measurements the model assembles new recordings in which it inserts the
   artifacts itself and therefore knows their cause.
3. Before the main calculations the model was checked: a classifier told its
   recordings from real ones only weakly, within a limit fixed in advance. The model
   was also checked on recordings from another device (EDABE) and on movements from
   another wristband (WESAD), which it had never seen.
4. The model was run in 300 scenarios, which differ in how movement and artifacts are
   linked. Each scenario was repeated 50 times, so that chance does not decide the
   outcome.

The rules of the analysis were fixed and published before the main results were
obtained. All numbers and their limits are in [study-description/RESULTS.md](study-description/RESULTS.md).

## What is in this repository

| Folder or file | What it contains |
|---|---|
| `study-description/METHOD.md` | the method step by step, with the source of every number |
| `study-description/RESULTS.md` | all results and what they do not show |
| `study-description/PREREGISTRATION.md` | the rules fixed before the results, and every later change with its reason |
| `study-description/REFERENCES.md` | the literature |
| `study-description/README-RU.md` | this page in Russian |
| `result-figures/` | the figures |
| `result-tables/` | the result tables from which the figures are drawn |
| `run-in-browser/` | a notebook that repeats the main results in the browser |
| `program-code/edasynth/` | the program |
| `program-tests/` | automatic tests of the program |
| `data-download-and-external-checks/` | programs that download the data and run the checks on EDABE and WESAD |
| `settings/` | all settings that can change a result |
| `downloaded-data/` | the downloaded data appear here; the folder describes where they come from |
| other files in the main folder | technical files for installation, citation and the licence |

## Running it

The **Open in Colab** button at the top of the page repeats the main results in the
browser in about 15 minutes without any installation. To run everything on your own
computer, install the program with `pip install -e ".[dev]"` and see `make help`.

## Authorship and licence

The research question, the design, the validation criteria and the interpretation are
by Nuray Doskaliyeva. The code was written with the help of an AI assistant (Claude)
and was run and checked by the author. The code is under the MIT licence; it can be
cited with [CITATION.cff](CITATION.cff). The data are not stored in the repository:
PhysioNet `noneeg` (Birjandtalab et al., 2016, 2017), the labels of Zhang et al. (2017),
EDABE (Llanes-Jurado et al., 2023) and WESAD (Schmidt et al., 2018).
