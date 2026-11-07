# Method of the study

This file describes every step of the study, every formula and every number used.
For each number it says where the number comes from. Some numbers are taken from a
cited publication, some were measured on the recordings by the program, and some are
decisions made in this study. The decisions were written into the settings file
`settings/default.yaml` before the main results were obtained, as described in
[PREREGISTRATION.md](PREREGISTRATION.md).

## 1. Question and design of the study

**The question of the study.** Wristbands record the electrodermal activity (EDA) of the skin and,
at the same time, the acceleration of the wrist. The study asks whether a detector of
motion artifacts in the EDA signal finds them better when it also sees the
acceleration than when it sees the EDA alone.

**Why the answer is given under conditions.** The benefit of the accelerometer depends
on how often movement actually disturbs the contact between the skin and the
electrodes, and on how it disturbs it. No existing data set measures this. The study
therefore describes the link between movement and artifacts by a few parameters,
generates recordings for many combinations of their values and measures the benefit
in each combination. Zhang et al. (2017), who found only a small benefit on real
recordings, suggested two possible reasons: the acceleration features may have been
unsuitable, and many movements may not affect the contact at all. The design follows
these two reasons. It compares four sets of features (section 5) and varies five
parameters of the link (section 4.4).

## 2. Recordings, expert labels and independent data sets

**The UTD recordings.** The main data come from the PhysioNet *Non-EEG Dataset for
Assessment of Neurological Status*, version 1.0.0 (Birjandtalab et al., 2016, 2017).
It contains recordings of 20 healthy participants. A wristband recorded the
acceleration of the wrist along three axes in units of g, the skin temperature and the
EDA, all 8 times per second. Each recording follows the same seven stages: relaxation,
physical stress (standing and walking on a treadmill), relaxation, cognitive stress,
relaxation, emotional stress and relaxation. The file headers give the unit of EDA as
"NU", which means that no unit is stated. The values are treated as microsiemens (µS),
because the median values of the participants lie between 0.03 and 7.8, which is
plausible for this unit. Every channel is checked against the checksum stored in its
file header. The program reads these files with its modules `data` and `wfdbio`.

**The expert labels.** Zhang et al. (2017) published artifact labels for these
recordings in their repository IdeasLabUT, which is used here at commit `7c4bd46`.
Three experts looked at every five-second piece of each recording, called a window,
and marked whether it contains an artifact. In this study a window counts as an
artifact when at least two of the three experts marked it; this rule is a decision of
the study. By this rule 349 of the 9 211 windows, or 3.79 %, are artifacts. The
agreement between the three experts, measured by Fleiss' κ, is 0.62. One of the
criteria the experts used was "a sudden change correlated with motion". The labels
therefore depend partly on movement and cannot serve as a fully independent truth
about the cause of an artifact.

**Matching the labels to the signal.** The labels are given per window, so the program
must know which samples of the signal belong to which window. It assumes that window
number *i* covers the samples 40*i* to 40*i* + 39, because five seconds at 8 samples per
second are 40 samples. To check this assumption, the program recomputes for every
window the mean, the standard deviation, the maximum and the minimum of the EDA from
the raw signal and correlates them with the same values in the feature table published
by Zhang et al. A participant is accepted when all four correlations are at least
0.999; this threshold is a decision of the study. All 20 participants were accepted,
and no shift between labels and signal was found. This step is done by the module
`alignment` and the command `align`.

**EDABE and WESAD as independent data.** Two further data sets are used only for the checks described in
section 4 of RESULTS.md. EDABE (Llanes-Jurado et al., 2023) contains 43 recordings made
with a Shimmer3 GSR+ device at 128 samples per second, in which every sample is marked
as artifact or clean; it has no acceleration. WESAD (Schmidt et al., 2018) contains
recordings of 15 participants who wore an Empatica E4 wristband, with EDA at 4 samples
per second, acceleration at 32 samples per second and the stages of the protocol, but
without artifact labels. The scripts `edabe_pack.py` and `wesad_pack.py` in the folder
`data-download-and-external-checks/` prepare these data sets.

## 3. Parameters measured on the recordings

The command `estimate` (module `estimation`) measures how real recordings behave and
stores the results in a file called the bank, `result-tables/bank-utd.json` together
with `bank-utd.npz`. When the generator later builds a recording, it draws values from
this bank. It does not assume a mathematical distribution for them.

**Clean windows.** A window counts as clean when all three experts marked it as clean
and both neighbouring windows are clean as well. The margin of one window on each side
follows rule 4 of Kleckner et al. (2018); its application here is a decision of the
study.

**Noise taken from the remainder of the decomposition.** The usual estimate of measurement noise uses the finest Haar wavelet
coefficients *d*₁ and computes σ = median|*d*₁| / 0.6745 (Donoho & Johnstone, 1994). On
the UTD recordings this estimate is zero for every participant, because the noise is
smaller than the smallest step the device can record. This step, called the
quantisation step *q*, is about 0.001 µS, and 0.006 µS for participant S10. The
generator therefore takes its noise from the part of the real signal that the
decomposition described next does not explain, *e* = *x* − *t* − *r*, where *x* is the
signal, *t* its slow tonic part and *r* its fast responses. This remainder is cut into
five-second blocks inside clean rest windows, each block is shifted to a mean of zero,
and the generator draws blocks at random with replacement. This approach is a decision
of the study.

**Splitting the signal into a slow and a fast part with cvxEDA.** The method cvxEDA
(Greco et al., 2016) describes a clean EDA signal as

$$
y(t) = s(t) + \sum_k a_k\,h(t-t_k) + \varepsilon(t),\qquad
h(t) = \left(e^{-t/\tau_0} - e^{-t/\tau_1}\right)[t\ge 0].
$$

Here *s* is the slow tonic level, each response *k* has the size *a<sub>k</sub>* and
starts at the time *t<sub>k</sub>*, every response has the shape *h*, called the Bateman
function, and ε is the noise. The method finds the parts by solving the optimisation
problem

$$
\min_{p\ge 0,\ \ell,\ d}\ \tfrac12\lVert Kp + B\ell + Cd - y\rVert^2 + \alpha\sum_n p_n + \tfrac{\gamma}{2}\lVert\ell\rVert^2 .
$$

In it, *p* is the driving signal of the sweat glands, *K* turns it into responses of the
Bateman shape, *Bℓ* is a smooth curve made of cubic pieces joined every 10 seconds, and
*Cd* is a constant offset with a linear trend. The constants τ₁ = 0.7 s, τ₀ between 2
and 4 s, α = 8·10⁻⁴ and γ = 10⁻² are taken from Greco et al. (2016). The signal is
scaled to a mean of zero and a standard deviation of one before solving and scaled back
afterwards; this is a decision of the study. The program tries τ₀ in steps of 0.1 s
and keeps the value at which the optimisation problem reaches its smallest value. Greco
et al. chose τ₀ by the smallest difference between the model and the signal instead.
This study departs from them because, on synthetic recordings with a known τ₀, only
the smallest value of the problem recovered the true τ₀; the program test
`testTau0SelectionRecoversKnownValue` shows this. When the authors' code, which needs
the package `cvxopt`, is not installed, a built-in solver (`decompose.numpySolve`) is
used.

**Rate and size of responses.** A response is a peak of the fast part that rises at least 0.05 µS above
its surroundings. This minimum size is the accepted one in the literature (Kleckner et
al., 2018; Boucsein et al., 2012, as summarised by Posada-Quintero & Chon, 2020). The
sizes of responses are measured only on responses whose start and peak both lie in
clean windows, so that artifacts cannot distort them; this is a decision of the study.
The rate of responses is the number of responses *K* divided by the duration *T*,
measured separately for the rest stages and for the stages with physical, cognitive or
emotional load.

**Tonic level.** Every stage longer than 180 seconds is described by the curve
*s*(τ) = *s*₀ + *b*τ + *A<sub>s</sub>* sin(2πτ / *T<sub>t</sub>* + φ), with a period
*T<sub>t</sub>* between 45 and 90 seconds. Greco et al. (2016) used this curve to
generate synthetic tonic levels; using it to describe real stages is a decision of this
study. The curve is fitted to the signal after the detected responses have been
removed.

**Artifact episodes.** An episode is a run of consecutive windows that the experts
marked as artifacts. For each episode the program draws a straight baseline between
the median of the two seconds before the episode and the median of the two seconds
after it. Relative to this baseline it stores the amplitude, the times of the start
and the end, the change of level between before and after, the share of the signal's
power above 0.4 Hz, called *P<sub>n</sub>* (Hossain et al., 2022), and the whole
course of the episode sample by sample. Sizes are stored relative to the level before
the episode, with a lower limit of 0.05 µS for that level, so that the same episode
can later be placed on a signal of a different level; this is a decision of the study.
Each episode is assigned to one class by the following rules, applied in this order,
which are decisions of the study:

1. **C5, loss or saturation of the signal**, if a value is at or below 0.01 µS or at
   or above 60 µS, or if the signal stays flat near zero for at least one second. The
   limit of 60 µS is taken from Kleckner et al. (2018).
2. **Unclassified**, if the amplitude is smaller than 0.05 µS in either direction.
3. **C3, shift of level**, if the level after the episode differs from the level before
   it by more than three times the effective noise σ<sub>eff</sub> and by at least
   0.05 µS.
4. **C4, oscillation while walking**, if at least half of the power lies above 0.4 Hz
   and the person is walking.
5. **C2, drop**, if the amplitude is negative, and otherwise **C1, peak**.

On UTD there are 208 episodes: 147 shifts of level, 43 drops, 14 peaks, 2 losses of
signal and 2 oscillations while walking.

**Movement events.** The program measures how abruptly the wrist moves by the jerk,
*j*[*n*] = *f<sub>s</sub>* ‖**a**[*n*] − **a**[*n*−1]‖ in g per second, where **a** is
the acceleration vector and *f<sub>s</sub>* the number of samples per second. A
movement event is a stretch of at least half a second in which the jerk is higher than
it is in 99 % of the same participant's rest time. Events less than one second apart
are joined into one, and the intensity of an event is its largest jerk. These rules are
decisions of the study.

**Strength of the walking artifact C4.** On UTD, the fast part of EDA follows the fast
part of acceleration with a strength of −0.0009 µS per g, which is far too small to
produce visible artifacts. The class C4 therefore uses a strength chosen by the study,
called the control gain, which makes the median C4 artifact 0.1 µS large. C4 serves as
a positive control, a case in which the right answer is known in advance. It is not a
claim about how large walking artifacts are in reality.

## 4. Generator of synthetic recordings

The generator, in the module `generator`, builds synthetic recordings from the bank.

**4.1 Donors.** Each synthetic participant combines the real acceleration and stages of
one participant, the movement donor, with the measured physiology of another
participant, the physiology donor. Because the two are different people, the synthetic
recording contains no link between movement and EDA except the one that the generator
adds on purpose with a known strength. A recording in which more than half of the EDA
values lie below 0.05 µS is not used as a physiology donor; on UTD this excludes
participant S18. The limit of 0.05 µS is rule 1 of Kleckner et al. (2018), and its use
here is a decision of the study. One run of the generator, called a replicate,
contains 20 synthetic participants.

**4.2 The clean signal.** The generator builds the clean signal at 64 samples per
second and then averages it down to 8 samples per second, the rate of the real device;
this is a decision of the study. The clean signal has four parts:

- The tonic level is made of the physiology donor's own stage curves, joined without
  jumps and kept at or above 0.05 µS.
- Spontaneous responses appear at random times, with the donor's rates, and their
  sizes are drawn from the donor's measured responses.
- Responses evoked by movement follow a movement event with the probability
  π<sub>r</sub>. They start one to five seconds after the movement, which is the usual
  delay of a response after a stimulus (Boucsein et al., 2012, as summarised by
  Posada-Quintero & Chon, 2020). Applying this delay to movement is an assumption.
- The noise consists of the donor's remainder blocks, and the result is rounded to the
  donor's quantisation step *q*.

**4.3 Artifacts.** Artifacts of the classes C1, C2 and C3 copy the measured course of a
real episode of the same class. The course is multiplied by the local level of the
clean synthetic signal, which is the median of the five seconds before the artifact
with a lower limit of 0.05 µS; this is a decision of the study. A shift of level (C3)
changes the state of the contact: the next shift replaces the previous offset instead
of adding to it, so that overlapping artifacts do not pile up. A downward change can
remove at most 95 % of the signal, and the result always stays between 5 % of the clean
signal and 60 µS. A C4 artifact is the control gain multiplied by the acceleration
magnitude minus its moving average over two seconds, and it appears only inside a
movement event. A C5 artifact replaces the signal by 0 µS or by the upper limit of the
device.

**4.4 The link between movement and artifacts.** These five parameters form the axes
of the grid of scenarios:

| Parameter | What it means | Values in the grid |
|---|---|---|
| *c* | the probability that an artifact is caused by a movement | 0, 0.25, 0.5, 0.75 and 1 |
| κ | how strongly the size of an artifact follows the intensity *I* of its movement; the size is multiplied by (*I* / median *I*)<sup>κ</sup>, so κ = 0 means no dependence | 0 and 1 |
| δ<sub>max</sub> | the longest delay between a movement and the artifact it causes; the delay is drawn evenly between zero and δ<sub>max</sub> | 1, 5 and 10 seconds |
| π<sub>r</sub> | the share of movements that are followed by a genuine EDA response | 0 and 0.5 |
| composition | which artifacts appear: only C1, only C2, only C3, only C5, or a mixture of the four in the shares measured on UTD | 5 options |

Artifacts that are not caused by a movement start at random times.

**4.5 How many artifacts, and how windows are labelled.** The generator adds artifacts
until 3.8 % of the windows contain one, which is the share in the UTD labels (Zhang et
al., 2017). A sample counts as corrupted when the artifact moves it by at least
0.05 µS. For a shift of level only the moment of change counts, and the lasting offset
after it does not; this is a decision of the study. A window is labelled as an
artifact when any of its samples is corrupted. A variant of the study labels it only
when more than half of its samples are corrupted.

**Variants of the generator.** Besides the main grid, four variants are computed. In the first, called
V2, a loss of contact multiplies the signal instead of adding to it,
*y* = *y*₀(1 + course), for drops and shifts of level. In the second, 15 % of the windows
contain artifacts instead of 3.8 %. In the third, a window counts as an artifact only
when most of it is corrupted. In the fourth, only C4 artifacts appear.

## 5. Features computed for each window

A detector does not look at the raw signal. It looks at a set of numbers computed for
each five-second window, called features. The module `features` computes them, with
the same code for real and synthetic recordings.

| Set | What it contains | Source |
|---|---|---|
| F0 | 24 features of the EDA signal: the mean, the standard deviation, the maximum and the minimum of the signal, of its first and second differences, and of its Haar wavelet details at three levels | Zhang et al. (2017) |
| FA | the same 24 statistics for the acceleration magnitude and for each of the three axes, 96 features in total | Zhang et al. (2017) |
| FA4 | the mean, the standard deviation, the maximum and the minimum of the acceleration magnitude, which are the first four features of FA | Zhang et al. (2017) |
| FB | four coupling features, described below | a decision of this study |

The coupling features ask whether the EDA changes right after a movement. Let *v* =
*f<sub>s</sub>*|Δ*y*| be the speed at which the EDA changes and *j* the jerk.
**B1** is the largest correlation between the jerk and the later speed of EDA change
over delays from 0 to 5 seconds. **B2** is the delay at which B1 is reached. **B3** is
the product of how unusually fast the EDA changes and how unusually strong the jerk
is, each compared with its median over the preceding 60 seconds. **B4** is the mean
coherence between the acceleration magnitude and the EDA between 0.5 and 3 Hz.
All four use only the current window and the five seconds before it, so a wristband
could compute them in real time.

## 6. Detectors and measures of their quality

The modules `models` and `experiment` train the detectors and measure their quality.

- **The detector** is a logistic regression with L2 regularisation. The features are
  scaled before training, the two classes are weighted so that the rare artifacts
  count as much as the clean windows, and the regularisation constant is *C* = 1, the
  same as in the comparison on real data; these settings are decisions of the study.
  On the real data a random forest with 500 trees is trained as well, to check that
  the result does not depend on the type of model.
- **Leaving one participant out.** The detector is trained on all participants except
  one and tested on that one. This is repeated for every participant, and the
  predictions of all repetitions are pooled. The detector is thus always tested on a
  person it has never seen.
- **AUPRC**, the area under the precision–recall curve, is the main measure of quality.
  It shows how well a detector ranks true artifacts above clean windows. A detector
  that guesses scores the share of artifacts, 0.038 here, and a perfect detector scores
  1. AUPRC is preferred to the more common AUROC because artifacts are rare, and AUROC
  can look good for a detector that still raises many false alarms.
- **The gain** from a set of acceleration features X is
  Δ<sub>X</sub> = AUPRC(F0 + F<sub>X</sub>) − AUPRC(F0), where X is A, A4 or B. It says
  how much better the detector becomes when it sees these features in addition to the
  EDA features.
- **The cost.** A detector set to recognise 95 % of the clean windows as clean still
  flags some windows. The cost is the share of genuine EDA responses whose peak lies in
  a flagged window and would therefore be thrown away. A flagged window may also contain a real
  artifact, and such windows are counted too, so the measure overstates the cost.
- **Uncertainty.** On real data the 95 % interval is obtained by drawing the
  participants at random with replacement 1 000 times. On synthetic data it is the 95 %
  *t*-interval over the replicates.

## 7. Checks of the generator

The command `validate` (module `validation`) runs seven checks before the main
calculations. If a check fails, the generator is corrected; its threshold is never
changed to make it pass.

| Check | What it tests | When it passes |
|---|---|---|
| G1 | the automatic tests of the program | all tests pass |
| G2 | 20 synthetic participants without artifacts are built from known physiology, and the program measures that physiology again; this is repeated with response shapes that differ from the assumed one (τ₀ changed by up to 25 %, τ₁ between 0.35 and 1.05 s) | the median error of τ₀ is at most 0.3 s, the median relative error of the response rate is at most 0.2, and a Kolmogorov–Smirnov test does not distinguish the measured response sizes from the true ones (*p* > 0.05) |
| G3 | whether a classifier can tell real clean windows from synthetic ones; the synthetic participant compared with a real one is built from the other 19 participants, so the generator has never seen the person it is compared with | a logistic regression on the F0 features, tested by five-fold cross-validation grouped by participant, reaches a median AUC of at most 0.70 over 10 repetitions; 0.5 would mean that the two cannot be told apart at all |
| G4 | the same question for windows with artifacts | as for G3 |
| G5 | that no gain appears when there is no link, once with *c* = 0 and once with a full link but with the acceleration shifted by at least 60 seconds | the 95 % interval of the gain contains zero and the mean gain is smaller than 0.01 |
| G6 | that a link built into the data is found, with C4 artifacts only and with the mixture at *c* = 1 | the lower bound of Δ<sub>A</sub> or Δ<sub>B</sub> is above zero |
| G7 | that the features of this study, computed from the raw UTD signals, reproduce the result obtained with the features published by Zhang et al. | the gain has the same sign, and AUPRC differs by no more than 0.10 |

G5 also fails when a set of features loses significantly. In fact every set loses
slightly when there is no link. This was accepted as an exception, as recorded in
PREREGISTRATION.md, and the gains in the map are therefore also compared with that
loss.

## 8. Grid of scenarios, number of replicates and analysis

The commands `pilot`, `grid`, `analyze` and `power` run the experiment.

- **The grid.** The five parameters of section 4.4 give 5 × 2 × 3 × 2 × 5 = 300
  scenarios, also called cells. In some cells a parameter has no effect: κ and
  δ<sub>max</sub> do nothing when *c* = 0, and κ does nothing for C5. Such cells are
  computed once and copied, so 226 cells are actually computed.
- **The number of replicates.** A small preliminary run, the pilot, gives the largest
  standard deviation *s* of the gain. The number of replicates per cell is then
  *R* = ⌈(1.96 *s* / 0.005)²⌉ with an upper limit of 50. Without the limit, this number would make the 95 % interval
  of each cell's mean about ±0.005 wide. These rules are decisions of the study. The pilot
  gave *R* = 50.
- **The decision rule.** A cell counts as useful when the mean gain is at least 0.01
  and the lower bound of its 95 % interval is above zero; this rule is a decision of
  the study. All cells are reported, whatever their result.
- **Agreement with the real data.** A mixture cell agrees with the real data when the
  real gain of 0.027 lies within the middle 95 % of the gains Δ<sub>A4</sub> of its
  replicates. This is repeated with training labels that were flipped at random until
  their agreement with the true labels, measured by Cohen's κ, is about 0.62, the
  agreement between the UTD experts. Using that value as a target is an assumption.
- **Sensitivity.** Starting from a reference cell (*c* = 0.5, κ = 1,
  δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0, mixture), each parameter is changed alone to
  see how much it moves the result.
- **Power.** For three mixture cells that agree with the real data, chosen at the
  quartiles of their mean Δ<sub>B</sub>, 50 laboratory studies with 5, 10, 15 or 20
  participants are simulated. The power is the share of these studies in which the
  bootstrap interval of the gain lies above zero.

## 9. Limitations of the method

- The synthetic results hold under the assumed link and the assumed shapes of
  artifacts. They do not say how often movement disturbs the electrode in real people.
- The generator is built from one device, one protocol and 20 people. The checks on
  two other data sets, described in section 4 of RESULTS.md, support its transfer to
  other devices but do not prove it.
- The expert labels partly depend on movement, so they are not an independent truth
  about the cause of an artifact.
- Artifact sizes are stored relative to the signal level, which assumes that artifacts
  grow with the level of the signal.
- The noise made of remainder blocks keeps its internal correlations only within five
  seconds.
- The features published for UTD were scaled within each participant, while the
  synthetic pipeline scales them over the training set. This difference matters for
  the comparison with the real data.
