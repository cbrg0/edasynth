# Results

This file reports what was measured, what it means and what it does not show. Where a
statement is a hypothesis and was not tested, it is called a hypothesis. The measures
used here are explained in [section 6 of METHOD.md](METHOD.md#6-detection-models-and-evaluation-metrics),
and the files with all numbers are in the folder `result-tables/`.

## 1. Baseline on real data

The first comparison uses the feature tables published by Zhang et al. (2017). They
cover 9 211 five-second windows from 20 participants, of which 349, or 3.79 %, were
marked as artifacts by at least two of the three experts. One detector sees only the 24
EDA features. The other detector sees these features and, in addition, the four
features of the acceleration magnitude, 28 features in total. Each detector is trained
on 19 participants and tested on the remaining one, in turn for every participant.

| Detector | AUPRC with EDA only | AUPRC with EDA and acceleration | Gain Δ | 95 % interval of the gain |
|---|---|---|---|---|
| Logistic regression | 0.4486 | 0.4754 | +0.0269 | from −0.0022 to +0.0580 |
| Random forest with 500 trees | 0.4641 | 0.4900 | +0.0259 | from −0.0206 to +0.0555 |
| Logistic regression with features recomputed from the raw signals by this program, which is check G7 | 0.4564 | 0.4818 | +0.0253 | from −0.0013 to +0.0572 |

**Interpretation.** The real data show neither a gain nor its absence:
any gain from slightly below zero to about +0.06 agrees with them. They also cannot
say whether the accelerometer helps against real disturbances of the contact, because
one of the criteria the experts used refers to movement, so the labels are not
independent of it.

## 2. Generator validation

Before the main calculations, the final version of the generator, 1.4.0, had to pass
the checks described in section 7 of METHOD.md. All of them were run on the UTD
recordings.

| Check | Result |
|---|---|
| G1, the automatic tests of the program | passed |
| G2, measuring known physiology again with the response shape that cvxEDA assumes | passed: the error of τ₀ is 0.0 s, the relative error of the response rate is 0.144, and the Kolmogorov–Smirnov test gives *p* = 0.37 |
| G2, the same with response shapes that differ from the assumed one | passed: the error of τ₀ is 0.15 s, the relative error of the rate is 0.125, and *p* = 0.33 |
| G3, telling real clean windows from synthetic ones | passed: the median AUC is 0.59, and the 10 repetitions range from 0.50 to 0.66 |
| G4, telling real artifact windows from synthetic ones | passed: the median AUC is 0.555; the repetitions range from 0.45 to 0.73, and 2 of the 10 lie above the limit of 0.70 |
| G5, no link between movement and artifacts, or acceleration shifted in time | no feature set gains anywhere, but all lose slightly: Δ<sub>A</sub> = −0.032, Δ<sub>A4</sub> = −0.002 and Δ<sub>B</sub> = −0.002; this was accepted as an exception |
| G6, finding a link that was built in | passed: with C4 artifacts only, Δ<sub>B</sub> = +0.008 with an interval from +0.006 to +0.010; with the mixture at *c* = 1, Δ<sub>B</sub> = +0.050 with an interval from +0.040 to +0.060 |
| G7, agreement with the published result on real data | passed, as the last row of the table in section 1 shows |

**Caveat.** The realism of artifact windows is close to its limit. Single
repetitions reach an AUC of 0.73. A version of the check that builds the generator
from one half of the participants and compares it with the other half fails in one
direction, with 0.76. The main reason is participant S10, whose recording has a
coarser quantisation step of 0.006 µS instead of about 0.001 µS. The scheme in which
the compared participant is left out, and the median over 10 repetitions, were chosen
after this was seen.

## 3. Simulation results

The map was computed with version 1.4.0, which is fixed by the tag `prereg-synth-v3`.
It has 300 cells, of which 226 were computed and the rest copied, as explained in
section 8 of METHOD.md. Each computed cell was run 50 times with 20 synthetic
participants, which gives 11 300 runs.

**3.1 Useful cells.** A cell counts as useful when the mean gain
is at least 0.01 and the lower bound of its 95 % interval is above zero. The table gives
the number of useful cells out of the 60 cells of each composition.

| Composition | 96 acceleration features (Δ<sub>A</sub>) | 4 features of movement strength (Δ<sub>A4</sub>) | 4 coupling features (Δ<sub>B</sub>) |
|---|---|---|---|
| C1, peaks | 2 | 22 | 40 |
| C2, drops | 0 | 12 | 28 |
| C3, shifts of level | 0 | 15 | 25 |
| C5, loss of signal | 0 | 0 | 0 |
| mixture of the four | 0 | 16 | 28 |
| **all 300 cells** | **2** | **65** | **121** |

**3.2 Effect of each parameter.** The next table shows the mean gain in
the mixture for each value of *c*, the share of artifacts caused by movement. Each
value is the mean over the 12 mixture cells with that *c*.

| *c* | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| Δ<sub>A</sub> | −0.027 | −0.034 | −0.033 | −0.029 | −0.020 |
| Δ<sub>A4</sub> | −0.001 | −0.002 | +0.002 | +0.009 | +0.019 |
| Δ<sub>B</sub> | −0.001 | +0.002 | +0.009 | +0.018 | +0.032 |
| useful cells for A4 and for B | 0 and 0 | 0 and 0 | 0 and 6 | 4 and 10 | 12 and 12 |

When each parameter is changed alone around the reference cell, the mean Δ<sub>B</sub>
changes by 0.049 for *c*, by 0.012 for π<sub>r</sub>, by 0.006 for δ<sub>max</sub> and
by 0.002 for κ. The share *c* therefore matters about four times more than any other
parameter.

**3.3 Signal loss (C5).** When the artifacts are losses of signal, the detector that
sees only EDA already reaches an AUPRC between 0.98 and 1.00. Every set of acceleration
features makes it worse: Δ<sub>B</sub> is about −0.025, and Δ<sub>A4</sub> falls to
−0.058 at *c* = 1.

**3.4 Cost for genuine responses.** The detectors were set to recognise 95 % of
clean windows as clean. When half of the movements are followed by a genuine response,
the share of such responses that fall into flagged windows rises from 0.256 with EDA
alone to 0.285 with the coupling features added. A flagged window may also contain a
real artifact, so this measure overstates the cost.

**3.5 Compatibility with the real data.** A mixture cell agrees with the real data when the
real gain of 0.027 lies within the middle 95 % of its replicates. This holds for 9 of
the 60 mixture cells, one with *c* = 0.75 and eight with *c* = 1. When the training
labels are made as noisy as the experts' labels, with an agreement of κ ≈ 0.62, it
holds for 8 of the 60 cells, all with *c* = 1. No cell with *c* of 0.5 or less agrees
with the real data.

**3.6 Variants.** The variants of section 4.5 of METHOD.md were computed for up to
three values of *c*. The half-widths of their 95 % intervals are at most 0.004.

| Variant | *c* | Δ<sub>A4</sub> | Δ<sub>B</sub> |
|---|---|---|---|
| 15 % of windows contain artifacts, mixture | 0, 0.5 and 1 | −0.001, +0.002 and +0.033 | +0.001, +0.007 and +0.027 |
| a window counts as an artifact only when most of it is corrupted, mixture | 0, 0.5 and 1 | −0.003, −0.005 and −0.005 | −0.001, −0.001 and +0.001 |
| loss of contact multiplies the signal (V2), drops | 0, 0.5 and 1 | −0.001, +0.004 and +0.020 | −0.001, +0.019 and +0.049 |
| loss of contact multiplies the signal (V2), shifts of level | 0, 0.5 and 1 | −0.003, +0.002 and +0.027 | −0.001, +0.020 and +0.059 |
| only C4 artifacts, with the control gain | 1 | +0.002 | +0.009 |

**3.7 Statistical power.** For three mixture cells that
agree with the real data, all with *c* = 1, 50 laboratory studies of each size were
simulated. The table gives the share of these studies in which the bootstrap interval
of the gain lies above zero, which is called the power. With 50 simulated studies, each
share is uncertain by up to 0.07.

| Cell of the mixture, *c* = 1 | Features | 5 participants | 10 | 15 | 20 |
|---|---|---|---|---|---|
| κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0.5 | FB and FA4 | 0.08 and 0.06 | 0.38 and 0.14 | 0.64 and 0.18 | 0.72 and 0.34 |
| κ = 0, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0.5 | FB and FA4 | 0.12 and 0.02 | 0.34 and 0.08 | 0.44 and 0.16 | 0.62 and 0.34 |
| κ = 1, δ<sub>max</sub> = 5 s, π<sub>r</sub> = 0 | FB and FA4 | 0.14 and 0.00 | 0.62 and 0.04 | 0.94 and 0.18 | 1.00 and 0.26 |

## 4. External validation

These checks were planned before they were run, as the last section of
PREREGISTRATION.md describes.

**EDABE.** This data set has 43 recordings from a Shimmer device, marked by other
experts and without acceleration. It contains 53 590 windows, and 19.1 % of them
contain an artifact. Because different devices give different signal levels, each
feature was scaled within each recording, as planned. The results without this
scaling are given as well.

| Check | Result | Criterion |
|---|---|---|
| E1, telling real artifact windows from synthetic ones | median AUC 0.678, range 0.651 to 0.718; without scaling 0.744 | at most 0.70, passed |
| E2, the same for clean windows | median AUC 0.546, range 0.529 to 0.583; without scaling 0.914 | at most 0.70, passed |
| E3, AUPRC on EDABE of a detector trained only on synthetic data | 0.556, range 0.537 to 0.583 | at least the value of a detector trained on real UTD data minus 0.05, passed |

For comparison, a detector trained on the real UTD recordings reaches 0.537 on EDABE,
a detector trained on EDABE itself reaches 0.611, and guessing reaches 0.191.

**WESAD.** In check W1, the movements of the UTD participants were replaced by the
wrist movements of 15 WESAD participants, recorded with an Empatica E4. The cell used
is the mixture with κ = 1, δ<sub>max</sub> = 1 s and π<sub>r</sub> = 0, run 50 times with
15 participants each. WESAD contains fewer movement events than UTD, with a median of
0.76 instead of 1.22 per minute.

| *c* | Δ<sub>B</sub> with WESAD and with UTD movement | Δ<sub>A4</sub> with WESAD and with UTD movement |
|---|---|---|
| 0 | +0.000 and −0.001 | −0.002 and −0.002 |
| 0.5 | **+0.014** and **+0.018** | **+0.011** and +0.004 |
| 1 | **+0.043** and **+0.048** | **+0.055** and **+0.027** |

Values in bold are useful by the decision rule. The coupling features are useful in
exactly the same cells with both kinds of movement, so the check passed.

## 5. Discussion

- **One quantity decides the answer.** It is the share *c* of artifacts that movement
  actually causes. It matters about four times more than any other parameter of the
  link, and when it is below 0.5 the accelerometer is never useful.
- **How acceleration is used matters.** Four features that test whether EDA changes
  right after a movement are useful in 121 cells, against 65 cells for four summaries
  of movement strength, and their result is more stable across devices. The 96
  features of Zhang et al. lose almost everywhere. A possible explanation, which was
  not tested and is therefore a hypothesis, is that a linear model with 120 inputs
  overfits when each replicate contains only about 370 artifact windows.
- **Acceleration should not be used to find losses of signal,** because EDA finds them
  by itself.
- **The real gain on UTD appears in the model only with a strong link,** with *c* of
  0.75 or more, provided that the generator is right. Two cautions apply. The expert
  labels partly depend on movement, and with less movement, as in WESAD, the same gain
  would correspond to a smaller *c*.
- **The model carries over to another device.** A detector trained only on synthetic
  data works on the data of another device as well as a detector trained on real
  recordings.
- **The next step** is a laboratory study that measures *c*. For the coupling features
  it needs about 15 to 20 participants.

## 6. Limitations

- They do not show how often movement disturbs the electrode in real life. The map
  gives the answer for each assumed value of *c*.
- They do not show that synthetic data are realistic in every respect. The realism
  checks use a linear classifier on 24 features, and for artifacts they pass close to
  their limit. The synthetic tonic level also contains a slow oscillation that real
  recordings do not show, as figure F1 makes visible.
- Synthetic artifacts are easier to find than real ones: the detector that sees only
  EDA reaches an AUPRC of about 0.65 on synthetic data and 0.45 on real data. The same
  gain may therefore mean something different at a different starting level.
- The power calculation assumes recordings and labels like those of UTD. Another
  device or another protocol may need more participants.
