# Method

This document describes every step of the study, the formulas behind it and the
origin of every number. Each step names the module in `src/edasynth/` and the
command that runs it. Readers new to EDA may want to start with
[BACKGROUND.md](BACKGROUND.md); terms are also collected in
[GLOSSARY.md](GLOSSARY.md).

Every number is labelled by its origin:

- **[source]**: taken from the cited publication; the reference is given.
- **[measured]**: estimated from real recordings by this code.
- **[decision]**: chosen by this project. All decisions are fixed in
  `configs/default.yaml` before the main experiment (see
  [PREREGISTRATION.md](PREREGISTRATION.md)).

Where a step departs from the source it builds on, the departure is stated.

**Contents.**
[1 Question and design](#1-question-and-design) ·
[2 Data](#2-data) ·
[3 Measuring parameters](#3-measuring-parameters-on-real-recordings) ·
[4 Generator](#4-generating-semi-synthetic-participants) ·
[5 Features](#5-features) ·
[6 Models and measures](#6-models-and-performance-measures) ·
[7 Validation](#7-validation-checks-g1g7) ·
[8 Experiment](#8-experiment-and-analysis) ·
[9 Reference detectors](#9-untrained-reference-detectors) ·
[10 Assumptions](#10-assumptions-and-limitations-of-the-method) ·
[11 Parameter table](#11-parameter-table)

---

## 1. Question and design

**Question.** Does a synchronous wrist accelerometer improve the detection of motion
artifacts in wrist EDA compared with a detector that uses EDA alone?

**Why a conditional answer.** The benefit depends on how movement is linked to
contact disturbances, and no available data set measures that link
([BACKGROUND.md §5](BACKGROUND.md#5-the-gap-this-project-addresses)). The study
therefore treats the link as a set of parameters, generates data for many values of
them and measures the benefit for each. Zhang et al. (2017) suggested two reasons
for the small benefit they observed: unsuitable acceleration features and movements
that do not affect the contact. These become the two axes of the design:

- **Axis 1, kind of acceleration information.** Four feature sets are compared
  (§5): EDA only (F0), EDA plus the 96 acceleration features of Zhang et al. (F0+FA),
  EDA plus 4 acceleration-magnitude features (F0+FA4) and EDA plus 4 features that
  describe whether EDA and movement change together (F0+FB).
- **Axis 2, the movement–artifact link.** Five parameters (§4.4): the share of
  artifacts caused by movement *c*, how strongly artifact size depends on movement
  intensity *κ*, the longest delay between movement and artifact *δ*<sub>max</sub>,
  the share of movements followed by a genuine response *π*<sub>r</sub>, and the
  artifact class.

The six goals listed in the [README](../README.md#what-the-project-does) follow from
this design. Goal 1 (validation) must succeed before the others are interpreted.

## 2. Data

Module `edasynth.data`, readers in `edasynth.wfdbio`; download with `make data`.

### 2.1 UTD recordings

The PhysioNet *Non-EEG Dataset for Assessment of Neurological Status*, version 1.0.0
(Birjandtalab et al., 2016, 2017), contains wrist recordings of 20 healthy
participants. Each `SubjectN_AccTempEDA` record holds three acceleration axes (g),
skin temperature (°C) and EDA, sampled at 8 Hz [source: record headers]. The
protocol has seven stages: relaxation (5 min), physical stress (standing, then
treadmill walking at 1 and 3 mph, about 5 min), relaxation, cognitive stress
(preceded by about 40 s of instructions), relaxation, emotional stress and a final
relaxation [source: PhysioNet data set page]. The annotation file marks the stage
transitions; its codes are not documented, so the stage parser relies on the order
of the annotations and stops with an error if their number is unexpected.

The header gives the EDA unit as "NU" (no unit). Treating the values as µS is an
assumption until `make inspect` has been run on the raw files
([STATUS.md](STATUS.md)).

The reader in `wfdbio.py` decodes the WFDB format 16 and checks every channel
against the initial value and checksum stored in the header, so a corrupted or
misread file stops the pipeline.

### 2.2 Expert labels

Zhang et al. (2017) published labels for UTD in their repository
(IdeasLabUT, n.d.), used here at commit `7c4bd46`. Three experts labelled every
non-overlapping 5-second window (40 samples) as clean or artifact. A window counts
as an artifact here when at least two of the three experts said so; Zhang et al. do
not state a combination rule, so the majority rule is a **[decision]**. The files
contain 9 211 windows, 349 of them artifacts (3.79 %) [measured]. Agreement between
the three experts is Fleiss' κ = 0.62 [measured from the label files].

### 2.3 AWW recordings

*Alan Walks Wales* (AWW) consists of recordings of one person on 11 days, labelled
in the same way by Zhang et al. (2017). Because all days come from one person, AWW
is used only to develop and test the code (`make smoke`); its numbers are never
reported as results.

### 2.4 Linking labels to signals

Module `edasynth.alignment`, command `align`.

Labels and signals were published separately, and no document states which label
belongs to which piece of signal. The working hypothesis is that label group *g*
belongs to `Subject(g+1)` and that window *i* covers samples 40*i* … 40*i*+39, as in
Zhang et al.'s feature script. It is checked per participant: the mean, SD, maximum
and minimum of each window are recomputed from the raw signal, z-scored within the
participant, and correlated with the first four columns of the published feature
matrix. The participant is accepted when all four correlations are at least 0.999
**[decision]** and enough raw windows exist. If the check fails, shifts of −10 to
+10 windows and 0 to 39 samples are searched; the best shift is stored and all later
steps drop that many samples from the start of the record. A participant that still
fails is excluded, with the reason recorded in `results/alignment.csv`.

---

## 3. Measuring parameters on real recordings

Module `edasynth.estimation`, command `estimate`, output
`results/bank-utd.json` and `.npz` (the "bank").

The literature gives the shape of a skin-conductance response, but not the noise
level, response rate, artifact shapes or movement statistics of this particular
device and protocol. These are measured here on the real recordings, and the
generator later resamples the measured values instead of assuming a distribution.

### 3.1 Clean windows

A window is treated as clean when all three experts called it clean and so did both
neighbouring windows. The one-window margin corresponds to the 5-s zone around
invalid data in rule 4 of Kleckner et al. (2018), because labelling errors tend to
occur at the edges of artifacts **[decision]**. Noise, response rates and
decomposition errors are measured only in clean windows.

### 3.2 Noise and quantisation

The white-noise level σ is estimated from the level-1 Haar wavelet details of the
clean windows recorded at rest:

$$
d_k = \frac{x_{2k} - x_{2k+1}}{\sqrt{2}}, \qquad
\hat\sigma = \frac{\operatorname{median}_k |d_k|}{0.6745}.
$$

For white Gaussian noise of standard deviation σ, the details *d<sub>k</sub>* also
have standard deviation σ, and the median absolute value of a normal variable is
0.6745 σ. Slow physiological changes hardly affect *d<sub>k</sub>* because
neighbouring samples are subtracted (Donoho & Johnstone, 1994). A recording without
clean rest windows (some AWW days) falls back to all clean windows; the bank records
which source was used.

The quantisation step *q* is the smallest non-zero absolute difference between
neighbouring samples that occurs in at least 1 % of all non-zero differences
**[decision]**; the 1 % rule ignores rare tiny differences caused by rounding. When
σ is smaller than the quantisation step the estimate can be 0; in that case the
noise floor σ<sub>eff</sub> = max(σ, *q*/√12) is used for classifying artifacts
(§3.6), because rounding alone produces noise with standard deviation *q*/√12.

### 3.3 Separating tonic and phasic parts

**Model of the clean signal.** A clean recording is modelled as

$$
y(t) = s(t) + r(t) + \varepsilon(t), \qquad
r(t) = \sum_k a_k\, h(t - t_k), \qquad
h(t) = \left(e^{-t/\tau_0} - e^{-t/\tau_1}\right)\,[t \ge 0],
$$

where *s* is the tonic level (µS), *r* the phasic part (µS), ε white noise, and each
response *k* starts at time *t<sub>k</sub>* with size *a<sub>k</sub>*. The response
shape *h* (the Bateman function) rises with time constant τ<sub>1</sub> and decays
with the slower time constant τ<sub>0</sub>. It reaches its maximum

$$
h_{\max} = h(t_{\text{peak}}), \qquad
t_{\text{peak}} = \frac{\tau_0\tau_1}{\tau_0-\tau_1}\ln\frac{\tau_0}{\tau_1},
$$

so a single isolated response has visible height *A<sub>k</sub>* = *a<sub>k</sub>*
*h*<sub>max</sub>. τ<sub>1</sub> = 0.7 s and τ<sub>0</sub> ∈ [2, 4] s
[source: Greco et al., 2016]. This is a model: it describes a typical response
shape, it does not claim that every response has this shape.

**cvxEDA.** Greco et al. (2016) estimate the parts by solving a convex optimisation
problem. With the signal z-scored (mean 0, SD 1) and sampled every δ = 1/8 s:

$$
\min_{p \ge 0,\ \ell,\ d}\ \ \tfrac12\,\lVert K p + B\ell + C d - y \rVert^2
+ \alpha \sum_n p_n + \tfrac{\gamma}{2}\lVert \ell \rVert^2 .
$$

Here *p* ≥ 0 is the sudomotor driver (one value per sample, mostly zero), *K* turns
the driver into the phasic part using a discrete form of the Bateman function,
*B*ℓ is a cubic spline with knots every 10 s that represents the tonic level, *Cd*
adds an offset and a linear trend, α penalises the total driver (favouring few
responses) and γ penalises large spline coefficients (favouring a smooth tonic level).
After solving, the parts are multiplied by the SD of the recording and the mean is
added back to the tonic part, so all parts are again in µS. α = 8·10⁻⁴, γ = 10⁻², 10-s knots
[source: Greco et al., 2016]. Greco et al. do not state that the input is z-scored;
z-scoring is a **[decision]** of this project, made so that the same α and γ work
for recordings with different levels.

The authors' reference implementation (lciti, n.d.) needs the `cvxopt` package. The
package also contains an independent solver (`decompose.numpySolve`, accelerated
projected gradient with exact elimination of the tonic variables) that solves the
same problem without extra dependencies and reports how well the optimality
(KKT) conditions are met. Its agreement with the reference implementation has not
been tested yet ([STATUS.md](STATUS.md)).

**Choosing τ<sub>0</sub>.** The decomposition is repeated for τ<sub>0</sub> = 2.0,
2.1, …, 4.0 s (step 0.1 s **[decision]**). Greco et al. (2016) chose the
τ<sub>0</sub> that minimised the ℓ<sub>2</sub> norm of the residual. **This project
departs from that rule** and chooses the τ<sub>0</sub> with the smallest value of the
whole objective above **[decision]**. On synthetic records with known τ<sub>0</sub>
the residual norm was almost the same for all τ<sub>0</sub> and pointed to the wrong
end of the range, while the objective recovered the true value within 0.5 s (test
`testTau0SelectionRecoversKnownValue`). The residual norm is still stored.

### 3.4 Responses and response rates

Responses are taken from the phasic part *r* found by the decomposition: every peak
of *r* whose height above the surrounding trough (its prominence) is at least
0.05 µS counts as a response. 0.05 µS is an accepted minimum response amplitude
[source: Kleckner et al., 2018; Posada-Quintero & Chon, 2020, citing Boucsein et al., 2012].
The onset is the latest lowest point of *r* in the 5 s before the peak, and the size
used by the generator is *a* = *A*/*h*<sub>max</sub>, so that an isolated synthetic
response has the measured height *A* **[decision]**.

An earlier version counted each run of non-zero driver as one response. That rule
was rejected because the optimal driver is dense on noisy, quantised signals and the
runs split into many spurious responses (14–95 per minute on AWW, with the solver
verified to be at its optimum).

The response rate is estimated separately for rest stages (the four relaxation
stages) and load stages (all others) as the maximum-likelihood rate of a Poisson
process, λ = *K*/*T*: *K* responses with onset in clean windows of that stage type
divided by the clean time *T* in seconds [measured]. A rate above 10 per minute
triggers a warning in the log. This threshold is a **[decision]** used only as a
plausibility flag; a published reference range for spontaneous response rates has
not been verified for this project.

### 3.5 Tonic level

The tonic level of each protocol stage is described by a straight line plus one slow
oscillation,

$$
s(\tau) = s_0 + b\,\tau + A_s \sin\!\left(\frac{2\pi\tau}{T_t} + \varphi\right),
$$

with τ the time since the start of the stage. The parameters are fitted by least
squares for every period *T<sub>t</sub>* = 45, 46, …, 90 s, and the best period is
kept. This form and the period range are those of the synthetic tonic signal used by
Greco et al. (2016) to test cvxEDA [source]; applying it to describe real stages is a
**[decision]**. Stages shorter than 180 s (twice the longest period) are skipped.

The fit is made to *s* + *r* − *r*<sub>events</sub>, the signal without the detected
responses, where *r*<sub>events</sub> = Σ *a<sub>k</sub>h*(*t* − *t<sub>k</sub>*)
is rebuilt from the responses of §3.4 **[decision]**. The cvxEDA tonic alone sits too
low, because the phasic part also contains small activity below the 0.05 µS
threshold; with it, synthetic AWW records were several µS below the real ones.

### 3.6 Artifact episodes and their classes

An **episode** is a run of consecutive artifact windows (majority label). For each
episode a local baseline is drawn as a straight line between the median of the 2 s
before and the median of the 2 s after it, and the deviation from that baseline is
measured:

| Quantity | Definition |
|---|---|
| amplitude *A* | the deviation with the largest absolute value (signed) |
| rise time, return time | time from 10 % to 90 % of \|*A*\| at the start, and from 90 % back to 10 % at the end |
| P<sub>n</sub> | share of the episode's spectral power above 0.4 Hz (Hossain et al., 2022) |
| level shift | median after minus median before |
| duration | number of windows × 5 s |

Each episode is assigned to one class by fixed rules, applied in this order
**[decision]**:

1. **C5, dropout or saturation**: a value ≤ 0.01 µS or ≥ 60 µS (60 µS is the upper
   range limit of Kleckner et al., 2018; 0.01 µS is a decision), or a flat run of at
   least 1 s near zero;
2. unclassified if \|*A*\| < 0.05 µS;
3. **C3, level shift**: \|shift\| > 3 σ<sub>eff</sub> and ≥ 0.05 µS;
4. **C4, walking oscillation**: P<sub>n</sub> ≥ 0.5 during a walking stage;
5. **C2, drop** if *A* < 0, otherwise **C1, peak**.

The measured rise, return, amplitude, shift and duration become the parameters of
the artifact templates of §4.3. For linear edges the 10–90 % time is 0.8 of the full
edge time, so edge times are rise time / 0.8.

These rules were written before looking at the UTD episodes. On the AWW development
data they put 68 % of episodes into C3, while the AWW traces mostly show short drops,
so the rules must be inspected on UTD before their class shares are used
([STATUS.md](STATUS.md)).

### 3.7 Movement events

Sudden wrist movements are found from the jerk, the rate of change of acceleration:

$$
j[n] = f_s\,\lVert \mathbf a[n] - \mathbf a[n-1] \rVert \quad (\text{g/s}),
$$

with *f<sub>s</sub>* = 8 Hz. A **movement event** is a run of at least 4 samples
(0.5 s) with *j* above θ, the participant's 99th percentile of *j* during rest
stages; runs less than 1 s apart are merged. The event's **intensity** *I* is its
largest jerk. All thresholds are **[decision]**. Using each participant's own
percentile makes the definition independent of how much a person moves at rest.

### 3.8 Walking gain for class C4

For the walking-oscillation class, the fast part of the EDA (signal minus its 2-s
moving average) is regressed on the fast part of the acceleration magnitude, through
the origin, over all samples of the walking stages. The slope *G* (µS per g) is
pooled across participants by the median **[measured]**. On AWW, *G* ≈ 0.0003 µS/g,
too small to produce visible artifacts.

### 3.9 How the expert labels relate to movement

Two descriptive numbers are stored: the share of artifact episodes that start between
1 s before and 5 s after a movement onset, and the share of movement events followed
by an episode in that window **[decision on the window]**. They describe the expert
labels, which partly depend on movement; they are not estimates of how often
movement really disturbs the electrode.

---

## 4. Generating semi-synthetic participants

Module `edasynth.generator`.

### 4.1 Donors

Synthetic participant *j* combines the real acceleration record, protocol stages and
movement events of a **movement donor** *u*(*j*) with the measured physiology
(τ<sub>0</sub>, response rates and sizes, tonic fits, σ, *q*) of a different
**physiology donor** *d*(*j*) ≠ *u*(*j*). The pairing is a random derangement, a
permutation without fixed points. Because the two donors are different people, any
relation between movement and EDA in the synthetic data was put there by the
generator, and the model cannot recognise a person from two matching signals. One
**replicate** contains as many synthetic participants as there are movement donors
(20 for UTD), each movement record used once.

### 4.2 Clean signal

The clean signal *y*<sub>0</sub> is built on a 64 Hz grid and then averaged in blocks
of 8 samples to 8 Hz **[decision]**. The finer grid lets response onsets and artifact
edges fall between the 8 Hz samples instead of only on them.

- **Tonic level.** For each protocol stage of the movement donor, one fitted tonic
  (§3.5) of the same stage type is drawn from the physiology donor with a random
  phase φ. Stages are joined without jumps: each starts where the previous ended.
- **Spontaneous responses.** In each stage, the number of responses is drawn from a
  Poisson distribution with the donor's rate for that stage type, onset times are
  uniform within the stage, and sizes are drawn from the donor's measured sizes.
- **Movement-evoked responses.** Each movement event is followed by a genuine
  response with probability π<sub>r</sub>, after a latency drawn uniformly from 1 to
  5 s. The latency range is the usual one for stimulus-evoked responses
  [source: Boucsein et al., 2012, as cited in Posada-Quintero & Chon, 2020]; applying it to
  movement is an assumption.
- **Noise and quantisation.** White Gaussian noise with the donor's σ is added; after
  the artifacts the signal is rounded to the donor's step *q*.

### 4.3 Artifact classes

The artifact shapes are the simplest shapes consistent with the published
descriptions **[decision]**; their sizes and durations are drawn from the measured
episodes of §3.6, pooled over all participants.

| Class | Shape | Origin of the idea |
|---|---|---|
| C1 peak | linear rise, linear fall (triangle) | peaks without exponential decay (Taylor et al., 2015; Zhang et al., 2017) |
| C2 drop | downward trapezoid (linear fall, plateau, linear return) | sudden decrease (Zhang et al., 2017) |
| C3 level shift | linear ramp to a lasting offset | persistent contact change (assumption) |
| C4 walking oscillation | *G* × (\|**a**\| − its 2-s moving mean) inside a movement event | walking artifacts (Hossain et al., 2022, foot EDA) |
| C5 dropout / saturation | values replaced by 0 µS or the upper limit | range rules of Kleckner et al. (2018) |

C4 is built directly from the acceleration, so the accelerometer must help with it.
It serves as a positive control (§7).

### 4.4 The movement–artifact link

Every artifact of classes C1, C2, C3 and C5 is **caused by movement** with
probability *c*. A caused artifact is attached to a random movement event *e* of
the movement donor: it starts at *t<sub>e</sub>* + *U*(0, δ<sub>max</sub>) and its
size is multiplied by

$$
\left(\frac{I_e}{\operatorname{median}(I)}\right)^{\kappa},
$$

so with κ = 0 the size does not depend on how vigorous the movement was, and with
κ = 1 it is proportional to its intensity. An artifact that is not caused by movement
starts at a uniformly random time with its measured size. C4 artifacts are always
attached to a movement event and scaled in the same way.

In a single-class cell every artifact has that class. In the mixture the class of each
new artifact is drawn with the class shares measured on the real episodes (§3.6),
C4 excluded.

### 4.5 How artifacts change the signal

- **V1, additive (main variant).** *y* = *y*<sub>0</sub> + Σ<sub>m</sub>
  *A<sub>m</sub>*(*t*): each artifact is added to the clean signal. C5 replaces
  values instead of adding to them.
- **V2, multiplicative (C2 and C3 only).** *y* = *y*<sub>0</sub>(1 − ρ *g*(*t*)),
  where *g* is the template shape scaled to 0…1 and ρ, the lost share of contact,
  is the measured size divided by the local level (median of the 5 s before) and
  capped at 0.95. The signal is always lowered, also for level shifts that were
  measured as upward. This variant represents the idea that losing contact removes a
  fraction of the signal; it is an assumption to be tested in the laboratory.

### 4.6 Prevalence and true labels

Artifacts are added one by one until 3.8 % of windows are positive, the UTD
proportion (Zhang et al., 2017) **[source]**, with at most 400 artifacts per
participant. A variant uses 15 %, the share in the AWW resting data.

A sample is **corrupted** when the artifact moves it by at least 0.05 µS, the minimum
visible response **[decision]**. For level shifts (C3) only the ramp and the first
sample at the new level count; the lasting offset afterwards is not labelled
**[decision]**. After the change the signal looks clean again at its new level, and
the published labelling criteria describe sudden changes, not shifted levels. A window is
**positive** when any of its samples is corrupted (main rule) or, as a variant, when
more than half are.

---

## 5. Features

Module `edasynth.features`. All features are computed for non-overlapping 5-second
windows (40 samples), the unit used by Taylor et al. (2015) and Zhang et al. (2017).
The same code is applied to real and synthetic recordings.

**F0, 24 EDA features** (Zhang et al., 2017). For six signals of each window, the
mean, standard deviation, maximum and minimum are computed. The six signals are the
EDA itself, its first difference, its second difference and its Haar wavelet detail
coefficients at levels 3, 2 and 1 (which Zhang et al. describe as 1, 2 and 4 Hz).
The first difference is the central difference (*x*<sub>n+1</sub> −
*x*<sub>n−1</sub>)/2, as in the published script. The stored UTD feature matrices
were found to use the forward difference instead (reconstruction from the raw
files); `features.derivative: forward` reproduces them.

**FA, 96 acceleration features** (Zhang et al., 2017): the same 24 statistics for
the acceleration magnitude and for each of the three axes. **FA4** is its first four
columns: mean, SD, maximum and minimum of the acceleration magnitude.

**FB, 4 coupling features** **[decision]**. These describe whether the EDA and the
movement change together, which FA and FA4 do not. With *v*[*n*] =
*f<sub>s</sub>*\|*y*[*n*] − *y*[*n*−1]\| the speed of EDA change (µS/s) and *j* the
jerk (§3.7):

| Feature | Definition | Meaning |
|---|---|---|
| B1 | max over lags ℓ = 0 … 5 s of corr(*j*[*n* − ℓ], *v*[*n*]) within the window | how closely EDA changes follow jerks |
| B2 | the lag ℓ of B1, in seconds | how long after the movement the EDA changes |
| B3 | [max *v* / (median₆₀ *v* + *qf<sub>s</sub>*)] × [max *j* / (median₆₀ *j* + *q*<sub>acc</sub>*f<sub>s</sub>*)] | unusually fast EDA change in the window together with an unusually strong jerk in the window or the 5 s before; medians over the preceding 60 s |
| B4 | mean magnitude-squared coherence of \|**a**\| and *y* in 0.5–3 Hz | shared rhythmic content, e.g. during walking |

All four use only the window and the 5 s before it, so they could run in real time.
The terms *qf<sub>s</sub>* keep the ratios finite on flat, quantised signals; B4 is 0
when either signal is constant.

---

## 6. Models and performance measures

Module `edasynth.models` and `edasynth.experiment`.

### 6.1 Models

The primary model is **logistic regression**: it estimates the probability that a
window is an artifact as a logistic function of a weighted sum of its features. The
settings are L2 regularisation, standardised features, class weights inversely
proportional to class frequencies (so the rare artifacts are not ignored) and the
`newton-cholesky` solver. The regularisation strength *C* ∈ {0.01, 0.1, 1, 10} is
chosen inside each training set by grouped 5-fold cross-validation, maximising
AUPRC **[decision]**; leave-one-participant-out inside each training set would be
slower. A **random
forest** of 500 trees is used as a robustness check.

**Leave-one-participant-out (LOSO).** The model is trained on all participants but
one and tested on the one left out; this is repeated for every participant. The test
participant never contributes to standardisation, the choice of *C* or the fit. The
predictions of all folds are pooled before the measures below are computed.

### 6.2 Measures

- **AUPRC, the main measure** (area under the precision–recall curve, computed as
  average precision). For every threshold on the predicted probability, *precision*
  is the share of flagged windows that are real artifacts and *recall* is the share
  of real artifacts that are flagged. AUPRC summarises this trade-off from 0 to 1.
  A detector that guesses scores about the artifact share (0.038 here), and a perfect
  detector scores 1. It is preferred over the ROC area when the positive class is
  rare, because it is sensitive to false alarms among the many clean windows.
- **Gain Δ.** Δ<sub>X</sub> = AUPRC(F0 + F<sub>X</sub>) − AUPRC(F0) for X = A, A4,
  B. A positive Δ means that adding the acceleration features improved detection.
  Δ<sub>B</sub> − Δ<sub>A4</sub> compares the kind of feature at equal count.
- **AUROC** (area under the ROC curve): the probability that a random artifact
  window gets a higher score than a random clean window. Reported for comparison
  with Zhang et al. (2017).
- **Sensitivity at 95 % specificity**: the share of artifact windows flagged when
  the threshold is set so that 95 % of clean windows are not flagged.
- **Lost responses (cost)**: at the same threshold, the share of genuine responses
  (height ≥ 0.05 µS) whose peak falls in a flagged window, separately for
  spontaneous and movement-evoked responses. It measures how much physiology a
  detector would discard.
- **Detection by class**: at the same threshold, the share of windows of each
  artifact class that are flagged.

### 6.3 Uncertainty

- **Participant bootstrap** (real data and power analysis): participants are drawn
  with replacement 1 000 times, Δ is recomputed on each sample, and the 2.5th and
  97.5th percentiles form the 95 % confidence interval. Resampling whole participants
  respects the fact that windows of one person are not independent.
- **Per-participant test** (real data): Δ is computed for each participant with at
  least one artifact, and a Wilcoxon signed-rank test asks whether the median is
  different from zero.
- **Across replicates** (synthetic data): each replicate is a new synthetic data
  set; the mean Δ over *R* replicates gets a 95 % *t*-interval.

---

## 7. Validation checks (G1–G7)

Module `edasynth.validation`, command `validate`. Thresholds are in
`configs/default.yaml` (section `validation`) and may not be changed to make a check
pass ([PREREGISTRATION.md](PREREGISTRATION.md)).

| Check | What it tests | Pass condition |
|---|---|---|
| G1 | the code: unit tests of readers, features, generator, models, estimation | all tests pass |
| G2 | parameter recovery: 20 artifact-free synthetic participants with known physiology are re-estimated with §3 | median \|τ<sub>0</sub> error\| ≤ 0.3 s; median \|relative rate error\| ≤ 0.2; two-sample Kolmogorov–Smirnov test of response sizes *p* > 0.05 |
| G3 | realism of clean windows: can a classifier tell real from synthetic? | logistic regression on F0, grouped 5-fold cross-validation, AUC ≤ 0.70 |
| G4 | realism of artifact windows, same test | AUC ≤ 0.70 (the real labels partly depend on movement, which limits this check) |
| G5 | negative controls: no movement–artifact link (*c* = 0), and a full link (*c* = 1) with the acceleration shifted in time by at least 60 s after generation | for Δ<sub>A</sub>, Δ<sub>A4</sub> and Δ<sub>B</sub> over 10 replicates: the 95 % interval contains 0 and \|mean\| < 0.01 |
| G6 | positive controls: C4 only, and the mixture with *c* = 1, κ = 1, δ<sub>max</sub> = 1 s | in each control, Δ<sub>A</sub> or Δ<sub>B</sub> has mean > 0 and lower 95 % bound > 0 |
| G7 | our features on real UTD reproduce the published-feature baseline | same sign of Δ, both AUPRC values within 0.10 of 0.449 and 0.475 |

A classifier two-sample test (G3, G4) is a direct check of realism: if a simple
classifier cannot separate real from synthetic windows better than AUC 0.70, the
synthetic windows resemble the real ones in the features the detectors use. The
threshold 0.70 is a **[decision]**. A linear classifier is a lenient test; passing
it does not prove that the data are realistic in every respect.

The negative controls test whether the pipeline invents a gain: when movement and
artifacts are unrelated, adding movement features must not help. The positive
controls test whether it can detect a gain that exists.

---

## 8. Experiment and analysis

Module `edasynth.experiment`, commands `pilot`, `grid`, `analyze`, `power`.

**Grid.** *c* ∈ {0, 0.25, 0.5, 0.75, 1} × κ ∈ {0, 1} × δ<sub>max</sub> ∈ {1, 5} s ×
π<sub>r</sub> ∈ {0, 0.5} × artifact composition ∈ {C1, C2, C3, C5, mixture in the
measured shares} = 200 cells **[decision]**. Variants on part of the grid: V2 for C2
and C3, 15 % prevalence, the majority labelling rule and C4 alone.

**Number of replicates.** A pilot runs 10 replicates in four cells (*c* ∈ {0, 1} ×
κ ∈ {0, 1}, mixture). With *s* the largest standard deviation of Δ<sub>A4</sub> or
Δ<sub>B</sub> among them, the number of replicates needed for a 95 % interval of
half-width *h* = 0.005 is

$$
R = \left\lceil \left(\frac{1.96\, s}{h}\right)^2 \right\rceil ,
$$

limited to 2 ≤ *R* ≤ 50 **[decision]**. Each replicate uses its own seed derived from
the cell and the replicate number, so every run is reproducible.

**Decision rule.** In a cell, the accelerometer counts as *useful* when the mean Δ
is at least 0.01 and the lower 95 % bound is above 0 **[decision]**. The whole map is
reported; no cell is singled out after the fact.

**Goal 5, compatibility with the real data.** A cell is *compatible* with UTD when
the real gain (Δ = 0.027, §1 of [RESULTS.md](RESULTS.md)) lies within the central
95 % of the cell's replicate values of Δ<sub>A4</sub>. A variant flips synthetic
labels at random until Cohen's κ between true and noisy labels is about 0.62, the
agreement among the UTD experts, and trains on the noisy labels; using the expert
agreement for this purpose is an **assumption**.

**Goal 6, what to measure first.** Starting from the reference cell (*c* = 0.5,
κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0, mixture), each parameter is varied
alone over its grid values, and the range of the mean Δ shows which parameter
matters most. For the first three compatible cells, the power for 5, 10, 15 and 20
participants is the share of 50 replicates whose participant-bootstrap interval of
Δ lies above zero.

---

## 9. Untrained reference detectors

Module `edasynth.detectors`, command `detectors`. These detectors use EDA only and
are not trained. They show which artifact classes are caught without movement
information (goal 2); they are not part of the map.

| Detector | Rule | Source |
|---|---|---|
| slope | largest \|d*y*/d*t*\| in the window (score) | common baseline |
| Kleckner | rules 1, 2 and 4 on a low-pass filtered signal (FIR, 0.35 Hz cut-off, 64 s long); rule 3 needs temperature and is omitted | Kleckner et al. (2018) |
| Kocielnik | relative change faster than +20 % or −10 % per second | Kocielnik et al., as cited in Taylor et al. (2015) |
| P<sub>n</sub> | share of power above 0.4 Hz (score) | Hossain et al. (2022) |
| SWT | stationary Haar wavelet transform, 8 levels, Gaussian-mixture thresholds per 400-s segment; a window is flagged if any coefficient in it exceeds the thresholds | after Chen et al. (2015); the window rule is a **[decision]** |

The SWT detector follows the description of Chen et al. (2015) as implemented
here; it has not been compared with the authors' own code.

---

## 10. Assumptions and limitations of the method

- **Conditional results.** Every synthetic result holds only under the assumed
  movement–artifact link and the assumed artifact shapes. It does not tell how often
  movement disturbs the electrode in people.
- **One device, one protocol.** All parameters come from 20 participants recorded
  with one wristband in one laboratory protocol.
- **Labels depend on movement.** The real labels used for estimation, realism checks
  and the real-data comparison were made with a criterion that refers to movement.
- **Model of the clean signal.** The tonic and response models are simplifications.
  Different decomposition methods give different tonic and phasic parts
  (Posada-Quintero & Chon, 2020); parameters measured with cvxEDA inherit its view.
- **Artifact shapes and classes** are the simplest shapes consistent with published
  descriptions, and the class rules have not yet been checked on UTD episodes.
- **Feature scaling differs between real and synthetic data.** The published UTD
  features are z-scored within each participant; the synthetic pipeline standardises
  over the training set only. This matters for the comparison in goal 5.
- **Label noise model.** Symmetric random flips are the simplest model of expert
  disagreement; real disagreements concentrate at artifact edges.

## 11. Parameter table

| Parameter | Value | Origin |
|---|---|---|
| window length | 5 s = 40 samples | Taylor et al. (2015); Zhang et al. (2017) |
| artifact window (real data) | ≥ 2 of 3 experts | decision |
| alignment correlation | ≥ 0.999 | decision |
| clean window | unanimously clean, neighbours too | decision, after Kleckner et al. (2018) rule 4 |
| noise σ | median\|*d*<sub>1</sub>\|/0.6745 | method: Donoho & Johnstone (1994); value measured |
| τ<sub>1</sub>, τ<sub>0</sub> range | 0.7 s; 2–4 s | Greco et al. (2016) |
| τ<sub>0</sub> step and selection | 0.1 s; smallest objective | decision (departs from Greco et al.) |
| α, γ, knot spacing | 8·10⁻⁴, 10⁻², 10 s | Greco et al. (2016) |
| z-scoring before cvxEDA | yes | decision |
| minimum response | 0.05 µS | Kleckner et al. (2018); Boucsein et al. (2012) via Posada-Quintero & Chon (2020) |
| tonic form, period range | line + sinusoid, 45–90 s | Greco et al. (2016), synthetic tonic |
| range limits for C5 | 0.01 µS (dropout), 60 µS | 60 µS: Kleckner et al. (2018); 0.01 µS: decision |
| P<sub>n</sub> band | > 0.4 Hz | Hossain et al. (2022) |
| movement event | jerk > 99th rest percentile, ≥ 0.5 s, merge < 1 s | decision |
| evoked-response latency | 1–5 s | Boucsein et al. (2012) via Posada-Quintero & Chon (2020); use for movement: assumption |
| generator grid | 64 Hz | decision |
| label threshold | 0.05 µS | decision |
| prevalence | 3.8 % (15 % variant) | Zhang et al. (2017) |
| V2 maximum lost contact | 0.95 | decision |
| FB lags, context, medians, band | 0–5 s, 5 s, 60 s, 0.5–3 Hz | decision |
| logistic regression *C* grid | 0.01, 0.1, 1, 10 | decision |
| bootstrap samples | 1 000 | decision |
| decision rule | mean Δ ≥ 0.01 and lower bound > 0 | decision |
| realism threshold | AUC ≤ 0.70 | decision |
| replicate half-width | 0.005, *R* ≤ 50 | decision |
| label-noise target | Cohen κ ≈ 0.62 | measured Fleiss κ of UTD experts; use as target: assumption |
