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
The long-term aim of this project is an algorithm that tells these sources apart with
the help of other sensors in the same wristband, such as an accelerometer, a skin
thermometer or a heart-rate sensor.

This repository covers the first step. It asks whether the accelerometer helps to find
artifacts and under which conditions it does.

## Result

Zhang et al. (2017) published recordings of 20 people in which experts marked the
artifacts. On these recordings, a detector that also sees the acceleration finds
artifacts so little better than a detector that sees only the skin signal that the
improvement may be zero. The recordings cannot explain why, because nobody recorded
which artifacts were caused by movement. A computer model of wristband recordings was
therefore built from real material, and in this model the cause of every artifact is
known.

![Gain from adding acceleration against the share of artifacts caused by movement](result-figures/F0-summary.png)

The vertical axis of the figure shows how much the detection of artifacts improves
when acceleration is added. The horizontal axis shows the share of artifacts that
movement causes. The accelerometer is useless when movement causes fewer than half of
the artifacts. When movement causes all of them, the most useful are four features
that check whether the skin signal changes right after a movement. The 96 acceleration
features used in earlier work make the detection worse in this figure and help in only
2 of the 300 scenarios. The benefit of the accelerometer therefore depends on a share
that no public data set contains, and the next step is a laboratory experiment that
measures it.

## How it was done

1. Real recordings were used to measure how the skin signal behaves, how people move
   and what real artifacts look like.
2. From these measurements the model assembles new recordings. It inserts the
   artifacts itself and therefore knows their cause.
3. Before the main calculations the model was checked. A classifier could tell its
   recordings from real ones only weakly, and stayed within a limit fixed in advance.
   The model was also checked on recordings from another device, the EDABE data set,
   and on movements from another wristband, the WESAD data set, neither of which it had
   seen before.
4. The model was run in 300 scenarios, which differ in how movement and artifacts are
   linked. Each scenario was repeated 50 times, so that chance does not decide the
   outcome.

The rules of the analysis were fixed and published before the main results were
obtained. All numbers and their limits are in [study-description/RESULTS.md](study-description/RESULTS.md).

## What is in this repository

| Folder or file | What it contains |
|---|---|
| `study-description/METHOD.md` | The method step by step, with the source of every number. |
| `study-description/RESULTS.md` | All results, and what they do not show. |
| `study-description/PREREGISTRATION.md` | The rules fixed before the results, and every later change with its reason. |
| `study-description/REFERENCES.md` | The literature, with what each source is used for. |
| `study-description/README-RU.md` | This page in Russian. |
| `result-figures/` | The figures of the study. |
| `result-tables/` | The tables with all results, from which the program draws the figures. |
| `run-in-browser/` | A notebook that repeats the main results in a web browser. |
| `program-code/edasynth/` | The program that carries out the study. |
| `program-tests/` | Automatic tests that check the program. |
| `data-download-and-external-checks/` | Programs that download the data and run the checks on EDABE and WESAD. |
| `settings/` | The settings file, which contains every number that can change a result. |
| `downloaded-data/` | The downloaded data are placed here; the description in this folder says where they come from. |
| other files in the main folder | Technical files needed for installing the program, citing it and stating its licence. |

## Running it

The easiest way is the **Open in Colab** button at the top of this page. It opens the
notebook on Google's free computers, so nothing has to be installed. To run the
notebook, choose *Runtime*, then *Run all*. The first part recomputes the main figure
and all checks from the result tables in less than a minute. The second part downloads
the public recordings and repeats the comparison on them in about ten minutes. The
third part runs the automatic tests in about two minutes. These times are approximate,
because they depend on the speed of the computer and of the internet connection. The
full computation of all 300 scenarios is not part of the notebook, because it takes one
to two days on a laptop; its results are stored in `result-tables/`, and the first part
of the notebook works from them.

To run the program on your own computer, you need Python 3.11 or newer. The command
`pip install -e ".[dev]"` installs it, and the command `make help` lists every step of
the study with a short description.

## Authorship and licence

The research question, the design, the validation criteria and the interpretation are
by Nuray Doskaliyeva. The code was written with the help of an AI assistant (Claude)
and was run and checked by the author. The code may be used under the MIT licence and
cited as described in [CITATION.cff](CITATION.cff). The data are not stored in the
repository and come from PhysioNet (Birjandtalab et al., 2016, 2017), from the labels
of Zhang et al. (2017), from EDABE (Llanes-Jurado et al., 2023) and from WESAD (Schmidt
et al., 2018).
