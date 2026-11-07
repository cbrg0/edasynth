# edasynth

**Does a wrist accelerometer help to detect motion artifacts in skin-conductance
(EDA) recordings, and under which conditions?**

Wristbands that measure electrodermal activity (EDA) usually also carry an
accelerometer. Using it to flag motion artifacts seems obvious, yet on real
expert-labelled data the gain is small and uncertain. This project explains why with a
validated semi-synthetic model and says what a laboratory study must measure next.

![Gain from adding acceleration against the share of movement-caused artifacts](figures/F0-summary.png)

## Answer in four points

1. **Real data do not decide it.** On the public UTD recordings (20 people, labels of
   Zhang et al., 2017) adding acceleration raises AUPRC from 0.449 to 0.475
   (Δ +0.027, 95 % CI −0.002 to +0.058).
2. **The gain depends on one unmeasured quantity:** *c*, the share of artifacts that
   movement actually causes. Over 300 simulated scenarios the accelerometer never
   helps when *c* < 0.5; at *c* = 1 four "coupling" features (does EDA change right
   after a movement?) add +0.032 AUPRC. They are useful in 121 of 300 scenarios, four
   acceleration-magnitude features in 65, the 96 features of Zhang et al. in 2.
3. **Where EDA alone suffices (signal loss), acceleration only adds errors.**
4. **The model holds on data it never saw.** Its artifacts pass a realism test against
   EDABE (another device, other experts); a detector trained only on synthetic data
   finds EDABE artifacts as well as one trained on real UTD data (AUPRC 0.556 vs
   0.537); the map keeps its shape with movement from another wristband (WESAD).

So the next step is a laboratory study that measures *c* directly. With the effect
sizes found here it needs about 15–20 participants.

Full numbers and caveats: [docs/RESULTS.md](docs/RESULTS.md).

## How the results were protected from wishful thinking

The analysis rules were fixed in public git tags **before** the main results existed
(`prereg-synth-v1` to `v3`, and a dated plan for the independent checks). Every check
the generator had to pass, every change made after a failure and one withdrawn run are
listed in [docs/PREREGISTRATION.md](docs/PREREGISTRATION.md).

| Check | Result |
|---|---|
| Unit tests (75) | pass |
| Parameter recovery (known physiology re-estimated) | pass |
| Realism, clean / artifact windows (classifier AUC ≤ 0.70) | 0.59 / 0.555 |
| No spurious gain without a movement–artifact link | no gain; small losses (accepted exception) |
| A built-in link is detected | pass |
| Our features reproduce the published real-data baseline | pass |
| EDABE realism, artifact / clean windows | 0.68 / 0.55 |
| Synthetic-trained detector on EDABE vs real-trained | 0.556 vs 0.537 |

## Repository

```
docs/METHOD.md            every step, formula and parameter, with its source
docs/RESULTS.md           what was measured and what it does and does not show
docs/PREREGISTRATION.md   rules fixed before the results, and every later change
docs/REFERENCES.md        sources
src/edasynth/             the code, one module per step
scripts/                  data download; preparation and checks of EDABE and WESAD
tests/                    unit tests
configs/default.yaml      every number that can change a result
results/, figures/        the outputs reported in RESULTS.md
```

## Reproducing

Python 3.11+, about 200 MB of disk.

```bash
git clone https://github.com/cbrg0/edasynth.git && cd edasynth
pip install -e ".[dev]"      # or: conda env create -f environment.yml
make test                    # unit tests, ~2 min
make data                    # downloads PhysioNet "noneeg" and the labels, checks SHA-256
python -m edasynth realdata --stored   # the real-data result, ~10 min
make before-tag              # estimation, checks, pilot: a few hours
make after-tag               # the full grid: about 1–2 days on a laptop
```

`make help` lists every step. Checks on EDABE and WESAD: see the header of
`scripts/external_checks.py`. One grid replicate recomputed on Linux matched the value
computed on macOS to the last digit.

## Authorship

Research question, design decisions, validation criteria and interpretation:
Nuray Doskaliyeva. The code was written with the help of an AI assistant (Claude) and
run and checked by the author.

## Licence and citation

MIT licence ([LICENSE](LICENSE)); cite with [CITATION.cff](CITATION.cff). The data are
not redistributed: PhysioNet `noneeg` (ODC-By 1.0; Birjandtalab et al., 2016, 2017),
labels of Zhang et al. (2017, BSD), EDABE (Llanes-Jurado et al., 2023), WESAD
(Schmidt et al., 2018). A short overview in Russian: [docs/README-RU.md](docs/README-RU.md).
