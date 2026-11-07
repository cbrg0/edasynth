# Method

Every step, formula and number of the study. Each number is labelled **[source]**
(from the cited publication), **[measured]** (estimated from the recordings by this
code) or **[decision]** (chosen here and fixed in `settings/default.yaml` before the
main results; see [PREREGISTRATION.md](PREREGISTRATION.md)).

## 1. Design

**Question.** Does a synchronous wrist accelerometer improve the detection of motion
artifacts in wrist EDA compared with a detector that uses EDA alone?

**Why a conditional answer.** The benefit depends on how often and how movement
disturbs the electrode contact, and no data set measures that. The link is therefore
treated as a set of parameters; data are generated for many values and the benefit is
measured for each. Zhang et al. (2017) suggested two reasons for the small benefit they
saw: unsuitable acceleration features, and movements that do not affect the contact.
These are the two axes of the design: four feature sets (§5) and five link parameters
(§4.4).

## 2. Data

**UTD** (module `data`, reader `wfdbio`). PhysioNet *Non-EEG Dataset for Assessment of
Neurological Status* v1.0.0 (Birjandtalab et al., 2016, 2017): 20 healthy participants,
wrist acceleration (3 axes, g), skin temperature and EDA at 8 Hz [source]. Seven
protocol stages: relaxation, physical stress (standing, treadmill), relaxation,
cognitive stress, relaxation, emotional stress, relaxation. The header gives the EDA
unit as "NU"; values are treated as µS (medians 0.03–7.8 per participant, plausible).
Every channel is checked against the checksum in its header.

**Labels.** Zhang et al. (2017) published labels for UTD (IdeasLabUT, commit
`7c4bd46`): three experts labelled each non-overlapping 5-s window. A window is an
artifact when at least two experts say so **[decision]**: 349 of 9 211 windows (3.79 %)
[measured]; Fleiss' κ between experts 0.62 [measured]. One labelling criterion is
"a sudden change correlated with motion", so the labels partly depend on movement.

**Linking labels to signals** (module `alignment`, command `align`). Window *i* is
taken to cover samples 40*i* … 40*i*+39. This is checked per participant by
correlating the window mean, SD, maximum and minimum recomputed from the raw signal
with the published feature matrix; accepted when all four correlations are ≥ 0.999
**[decision]**. All 20 participants were accepted with no offset [measured].

**Independent data** (used only for the checks of RESULTS.md §5). EDABE
(Llanes-Jurado et al., 2023): 43 Shimmer3 GSR+ recordings at 128 Hz with an artifact
label per sample, no acceleration. WESAD (Schmidt et al., 2018): 15 participants,
Empatica E4 on the wrist, EDA at 4 Hz and acceleration at 32 Hz, protocol stages, no
artifact labels. Preparation: `data-download-and-external-checks/edabe_pack.py`, `data-download-and-external-checks/wesad_pack.py`.

## 3. Parameters measured on the recordings

Module `estimation`, command `estimate`, output `result-tables/bank-utd.json/.npz` (the
"bank"). The generator later resamples measured values instead of assuming
distributions.

**Clean windows.** All three experts say clean, and so do both neighbouring windows
(one-window margin after rule 4 of Kleckner et al., 2018) **[decision]**.

**Noise.** The Haar-wavelet estimate σ = median|d₁|/0.6745 (Donoho & Johnstone, 1994)
is 0 for all UTD participants because the noise is below the quantisation step
(*q* ≈ 0.001 µS; 0.006 µS for S10) [measured]. The generator therefore uses the
**decomposition residual** *e* = *x* − *t* − *r* (below), cut into 5-s blocks inside
clean rest windows, centred, and resampled with replacement **[decision]**.

**Tonic and phasic parts: cvxEDA** (Greco et al., 2016). A clean signal is modelled as

$$
y(t) = s(t) + \sum_k a_k\,h(t-t_k) + \varepsilon(t),\qquad
h(t) = \left(e^{-t/\tau_0} - e^{-t/\tau_1}\right)[t\ge 0],
$$

with tonic level *s*, responses of size *a<sub>k</sub>* starting at *t<sub>k</sub>*
and the Bateman shape *h*. The parts are found by the convex problem

$$
\min_{p\ge 0,\ \ell,\ d}\ \tfrac12\lVert Kp + B\ell + Cd - y\rVert^2 + \alpha\sum_n p_n + \tfrac{\gamma}{2}\lVert\ell\rVert^2
$$

(driver *p*, Bateman matrix *K*, cubic spline *Bℓ* with 10-s knots, offset and trend
*Cd*). τ₁ = 0.7 s, τ₀ ∈ [2, 4] s, α = 8·10⁻⁴, γ = 10⁻² [source]. The signal is
z-scored before solving and rescaled after **[decision]**. τ₀ is searched in steps of
0.1 s and chosen by the **smallest objective**, departing from Greco et al. (residual
norm), because on synthetic records only the objective recovered the true τ₀ (test
`testTau0SelectionRecoversKnownValue`) **[decision]**. A built-in solver
(`decompose.numpySolve`) is used when the authors' `cvxopt` code is unavailable.

**Responses.** Peaks of the phasic part with prominence ≥ 0.05 µS, the accepted
minimum amplitude (Kleckner et al., 2018; Boucsein et al., 2012, via Posada-Quintero
& Chon, 2020). Sizes only from responses whose onset and peak lie in clean windows
**[decision]**. Rates λ = *K*/*T* separately for rest and load stages [measured].

**Tonic level.** Each stage longer than 180 s is fitted by
*s*(τ) = *s*₀ + *b*τ + *A<sub>s</sub>* sin(2πτ/*T<sub>t</sub>* + φ), *T<sub>t</sub>* =
45…90 s, the synthetic tonic of Greco et al. (2016) [source]; using it to describe real
stages is a **[decision]**. The fit is made to the signal minus the detected responses.

**Artifact episodes.** An episode is a run of consecutive artifact windows. Relative
to a straight baseline between the medians of the 2 s before and after, its amplitude,
edge times, level shift, high-frequency power share P<sub>n</sub> (> 0.4 Hz; Hossain et
al., 2022) and its **course** sample by sample are stored. Sizes are stored relative
to the level before the episode, max(level, 0.05 µS) **[decision]**. Classes, by fixed
rules applied in order **[decision]**:

1. **C5 dropout/saturation**: a value ≤ 0.01 µS or ≥ 60 µS (60 µS: Kleckner et al.,
   2018), or a flat run ≥ 1 s near zero;
2. unclassified if |amplitude| < 0.05 µS;
3. **C3 level shift**: |shift| > 3 σ<sub>eff</sub> and ≥ 0.05 µS;
4. **C4 walking oscillation**: P<sub>n</sub> ≥ 0.5 during walking;
5. **C2 drop** if the amplitude is negative, otherwise **C1 peak**.

UTD: 208 episodes, C3 147, C2 43, C1 14, C5 2, C4 2 [measured].

**Movement events.** Jerk *j*[*n*] = *f<sub>s</sub>*‖**a**[*n*] − **a**[*n*−1]‖ (g/s).
An event is a run ≥ 0.5 s with *j* above the participant's 99th percentile at rest;
runs < 1 s apart are merged; intensity = the largest jerk **[decision]**.

**Walking gain for C4.** The measured gain of fast EDA on fast acceleration is
−0.0009 µS/g on UTD, too small to produce artifacts. The C4 class therefore uses a
**control gain** that gives a median artifact of 0.1 µS **[decision]**; C4 is a
positive control, not a claim about real walking artifacts.

## 4. Generator

Module `generator`.

**4.1 Donors.** A synthetic participant combines the real acceleration and stages of a
**movement donor** with the measured physiology of a **different** physiology donor.
Any link between movement and EDA in the synthetic data is therefore put there by the
generator, at a known strength. A recording with more than half of its EDA below
0.05 µS is not a physiology donor (S18 on UTD) **[decision]**; 0.05 µS is rule 1 of
Kleckner et al. (2018). One replicate has 20 synthetic participants.

**4.2 Clean signal**, built at 64 Hz and averaged to 8 Hz **[decision]**:

- tonic: the physiology donor's own stage fits joined without jumps, kept ≥ 0.05 µS;
- spontaneous responses: Poisson counts with the donor's rates, sizes drawn from the
  donor's measured sizes;
- movement-evoked responses: after a movement event with probability π<sub>r</sub>,
  latency 1–5 s (the usual stimulus latency, Boucsein et al., 2012, via Posada-Quintero
  & Chon, 2020; applying it to movement is an assumption);
- noise: the donor's residual blocks; rounding to the donor's step *q*.

**4.3 Artifacts.** C1–C3 add the measured course of a real episode of that class,
scaled by the local level of the clean synthetic signal (median of the 5 s before,
floor 0.05 µS) **[decision]**. A level shift (C3) sets a contact state: the next shift
replaces the previous offset instead of adding to it. Downward changes are limited to
95 % of the signal, and the result is kept between 5 % of the clean signal and 60 µS.
C4 = control gain × (|**a**| minus its 2-s moving mean) inside a movement event. C5
replaces values with 0 µS or the upper limit.

**4.4 The movement–artifact link** (the grid axes):

| Parameter | Meaning | Grid |
|---|---|---|
| *c* | probability that an artifact is caused by movement | 0, 0.25, 0.5, 0.75, 1 |
| κ | size scales as (*I*/median *I*)<sup>κ</sup> with movement intensity *I* | 0, 1 |
| δ<sub>max</sub> | a caused artifact starts *U*(0, δ<sub>max</sub>) after the movement | 1, 5, 10 s |
| π<sub>r</sub> | share of movements followed by a genuine response | 0, 0.5 |
| composition | C1, C2, C3, C5, or their mixture in the measured shares | 5 options |

Artifacts not caused by movement start at random times. **4.5 Prevalence and labels.**
Artifacts are added until 3.8 % of windows are positive (the UTD share; Zhang et al.,
2017) [source]. A sample is corrupted when the artifact moves it by ≥ 0.05 µS; for C3
only the change itself, not the lasting offset **[decision]**. A window is positive
when any sample is corrupted (variant: more than half).

**Variants:** V2 multiplicative contact loss, *y* = *y*₀(1 + course) (C2, C3); 15 %
prevalence; majority labelling rule; C4 only.

## 5. Features

Module `features`; non-overlapping 5-s windows (40 samples), identical code for real
and synthetic data.

| Set | Features | Source |
|---|---|---|
| F0 | 24 EDA features: mean, SD, max, min of the EDA, its 1st and 2nd differences and its Haar details at three levels | Zhang et al. (2017) |
| FA | the same 24 statistics for the acceleration magnitude and each axis (96) | Zhang et al. (2017) |
| FA4 | mean, SD, max, min of the acceleration magnitude | first 4 of FA |
| FB | 4 coupling features (below) | **[decision]** |

With *v* = *f<sub>s</sub>*|Δ*y*| the speed of EDA change and *j* the jerk:
**B1** = max over lags 0–5 s of corr(*j*[*n*−ℓ], *v*[*n*]); **B2** = that lag;
**B3** = unusually fast EDA change times unusually strong jerk (each relative to its
median over the preceding 60 s); **B4** = mean coherence of |**a**| and *y* in
0.5–3 Hz. All use only the window and the 5 s before it, so they could run in real time.

## 6. Models and measures

Modules `models`, `experiment`.

- **Model:** L2 logistic regression, standardised features, balanced class weights,
  *C* = 1 (as in the real-data comparison) **[decision]**; random forest (500 trees)
  as a robustness check on real data.
- **Leave one participant out:** train on all but one, test on that one, pool the
  predictions.
- **AUPRC** (average precision), the main measure: guessing scores the artifact share
  (0.038), a perfect detector 1. Preferred to AUROC because artifacts are rare.
- **Gain** Δ<sub>X</sub> = AUPRC(F0 + F<sub>X</sub>) − AUPRC(F0), X = A, A4, B.
- **Cost:** at 95 % specificity, the share of genuine responses whose peak lies in a
  flagged window (this also counts windows that contain a real artifact).
- **Uncertainty:** participant bootstrap (1 000 resamples) on real data; 95 %
  *t*-interval over replicates on synthetic data.

## 7. Validation checks

Module `validation`, command `validate`. Thresholds may not be changed to make a check
pass.

| Check | Tests | Pass condition |
|---|---|---|
| G1 | unit tests | all pass |
| G2 | 20 artifact-free synthetic participants with known physiology re-estimated; also with varied response shapes (τ₀ ± 25 %, τ₁ 0.35–1.05 s) | median \|τ₀ error\| ≤ 0.3 s, median \|rate error\| ≤ 0.2, KS test of sizes *p* > 0.05 |
| G3 | can a classifier tell real from synthetic **clean** windows? The generator never sees the compared participant (built from the other 19) | logistic regression on F0, grouped 5-fold CV, median AUC over 10 repetitions ≤ 0.70 |
| G4 | the same for **artifact** windows | as G3 |
| G5 | no link (*c* = 0), and full link with the acceleration shifted by ≥ 60 s | 95 % interval contains 0 and \|mean\| < 0.01 |
| G6 | C4 only; mixture with *c* = 1 | Δ<sub>A</sub> or Δ<sub>B</sub> with lower bound > 0 |
| G7 | our features on raw UTD reproduce the published-feature baseline | same sign of Δ, AUPRC within 0.10 |

G5 also fails on a significant **loss**. Feature sets lose slightly without a link
(accepted exception, PREREGISTRATION.md), so their gains in the map are also read
against that loss.

## 8. Experiment and analysis

Commands `pilot`, `grid`, `analyze`, `power`.

- **Grid:** 5 × 2 × 3 × 2 × 5 = 300 cells. Cells where a parameter cannot act (κ and
  δ<sub>max</sub> at *c* = 0; κ for C5) are computed once and copied: 226 computed.
- **Replicates:** a pilot gives the largest SD *s* of Δ;
  *R* = ⌈(1.96 *s* / 0.005)²⌉, capped at 50 **[decision]**. *R* = 50.
- **Decision rule:** a cell is *useful* when the mean Δ ≥ 0.01 and its lower 95 %
  bound > 0 **[decision]**. All cells are reported.
- **Compatibility with the real data:** a mixture cell is compatible when the real
  Δ = 0.027 lies within the central 95 % of its replicate values of Δ<sub>A4</sub>.
  Repeated with training labels flipped at random to Cohen's κ ≈ 0.62 (the expert
  agreement; using it as a target is an assumption).
- **Sensitivity:** from the reference cell (*c* = 0.5, κ = 1, δ<sub>max</sub> = 1 s,
  π<sub>r</sub> = 0, mixture) each parameter is varied alone.
- **Power:** for three compatible mixture cells at the quartiles of their mean
  Δ<sub>B</sub>, the share of 50 simulated studies with 5, 10, 15 or 20 participants
  whose participant-bootstrap interval of Δ lies above zero.

## 9. Limitations of the method

- Synthetic results are conditional on the assumed link and artifact shapes; they do
  not say how often movement disturbs the electrode in people.
- The generator is built from one device, one protocol and 20 people; checks on two
  other data sets (RESULTS.md §5) support, but do not prove, transfer.
- Expert labels partly depend on movement and are not an independent ground truth.
- Relative artifact sizes assume that artifacts scale with the signal level.
- Residual-block noise keeps correlations only up to 5 s.
- Published UTD features are z-scored within participants; the synthetic pipeline
  standardises over the training set. This matters for the compatibility comparison.
