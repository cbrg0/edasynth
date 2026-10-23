# Glossary

Terms are in alphabetical order. Section references point to [METHOD.md](METHOD.md)
unless stated otherwise.

**Accelerometer (ACC).** A sensor that measures acceleration along three
perpendicular axes, in units of g (1 g = 9.81 m/s², the acceleration due to gravity).
At rest the magnitude of the three-axis vector is about 1 g because gravity is always
measured. Most EDA wristbands contain one.

**Artifact (motion artifact).** A part of the recording caused by the measurement,
here by movement of the electrode against the skin, and not by the body. It can look
like a peak, a drop, a lasting step, a fast oscillation or a loss of signal.

**Artifact class.** The five shapes used by the generator (§4.3): C1 peak, C2 drop,
C3 level shift, C4 walking oscillation, C5 dropout or saturation.

**AUPRC.** Area under the precision–recall curve, the main performance measure
(§6.2). It ranges from 0 to 1; guessing gives about the share of artifact windows
(0.038 here) and a perfect detector gives 1.

**AUROC.** Area under the ROC curve: the probability that a randomly chosen artifact
window receives a higher score than a randomly chosen clean window. 0.5 means
guessing.

**AWW.** *Alan Walks Wales*: wrist recordings of one person on 11 days, labelled by
Zhang et al. (2017). Used here only to develop and test the code.

**Bank.** The file of parameters measured on real recordings
(`results/bank-<source>.json` and `.npz`), from which the generator draws.

**Bateman function.** The shape used for one skin-conductance response,
*h*(*t*) = e<sup>−*t*/τ0</sup> − e<sup>−*t*/τ1</sup>: a rise with time constant τ<sub>1</sub>
(0.7 s) followed by a slower decay with time constant τ<sub>0</sub> (2–4 s) (§3.3).

**Bootstrap (participant).** A way to estimate uncertainty: whole participants are
drawn at random with replacement many times, the measure is recomputed each time,
and the spread of the results gives a confidence interval (§6.3).

**Causal feature.** A feature that uses only the current window and earlier data, so
it could be computed in real time.

**Classifier two-sample test.** A realism check: a classifier is trained to tell real
from synthetic windows. If it cannot do much better than chance (here AUC ≤ 0.70), the
two kinds of windows are similar in the features used (§7).

**Clean window.** In the real data, a window that all three experts labelled clean and
whose neighbours were also labelled clean (§3.1).

**Coherence (magnitude-squared).** A number between 0 and 1 that measures, for each
frequency, how consistently two signals oscillate together. Used in feature B4.

**Confidence interval (95 %).** A range of values that are compatible with the data at
the 95 % level. If it contains zero, the data do not show a clear effect in either
direction.

**cvxEDA.** A method by Greco et al. (2016) that splits an EDA recording into tonic
level, phasic responses and noise by solving a convex optimisation problem (§3.3).
Used here only to measure properties of real recordings.

**Derangement.** A random reordering in which no element stays in its place. Used to
pair each movement donor with a different physiology donor.

**Donor.** A real participant who supplies part of a synthetic participant: the
movement donor supplies the acceleration record and movement events, a different
physiology donor supplies the EDA parameters (§4.1).

**Driver (sudomotor driver).** In cvxEDA, the sparse non-negative input signal that
produces the phasic responses when passed through the Bateman shape.

**EDA (electrodermal activity, skin conductance).** How easily a very small current
passes through the skin, in microsiemens (µS). It rises with sweat-gland activity,
which is controlled by the sympathetic nervous system ([BACKGROUND.md](BACKGROUND.md)).

**Episode.** A run of consecutive artifact windows in the real data, measured and
classified in §3.6.

**Feature.** A number computed from one window, such as the maximum of the EDA or the
standard deviation of the acceleration magnitude. **Feature sets** (§5): F0, 24 EDA
features; FA, 96 acceleration features; FA4, 4 acceleration-magnitude features; FB,
4 coupling features describing whether EDA and movement change together.

**Gain (Δ).** AUPRC with the acceleration features minus AUPRC with the EDA features
only. Positive means that adding acceleration helped.

**Ground truth.** The correct answer for each window. In the synthetic data it is known
exactly because the generator placed the artifacts; in the real data only expert
labels exist.

**Haar wavelet details.** Differences between neighbouring samples (level 1) or
between neighbouring averages (levels 2, 3), scaled by 1/√2. They capture fast
changes at different time scales and are used in F0, FA and the noise estimate.

**Jerk.** The rate of change of acceleration, *j*[*n*] = *f<sub>s</sub>*‖**a**[*n*] −
**a**[*n*−1]‖, in g/s. It is large during sudden movements (§3.7).

**KKT conditions.** The optimality conditions of a constrained optimisation problem.
The built-in cvxEDA solver reports how far its solution is from meeting them.

**Label.** The answer "artifact" or "clean" for one window. Real labels come from
experts; synthetic labels come from the generator.

**Leave-one-participant-out (LOSO).** Train on all participants but one, test on the
one left out, repeat for everyone. It measures performance on people the model has
never seen.

**Logistic regression.** A model that estimates the probability of an artifact as a
logistic function of a weighted sum of the features; the main model here (§6.1).

**Movement event.** A period of at least 0.5 s in which the jerk exceeds the
participant's 99th percentile at rest (§3.7). Its intensity is its largest jerk.

**Movement–artifact link.** The assumed relation between movement and artifacts,
described by five parameters (§4.4): *c*, the share of artifacts caused by movement;
κ, how strongly artifact size grows with movement intensity; δ<sub>max</sub>, the
longest delay after a movement; π<sub>r</sub>, the share of movements followed by a
genuine response; and the artifact class.

**Negative and positive control.** Runs in which the right answer is known: no link
between movement and artifacts (the accelerometer must not help) or artifacts made
from movement (it must help). The pipeline has to pass both (§7).

**Phasic response (skin conductance response, SCR).** A short rise of EDA starting
1–5 s after a stimulus, followed by a slower decay. A rise of at least 0.05 µS is
conventionally counted.

**Precision and recall.** Among windows flagged as artifacts, precision is the share
that are artifacts; among all artifact windows, recall is the share that are flagged.

**Preregistration.** Writing down the analysis rules, with a public date, before the
results exist ([PREREGISTRATION.md](PREREGISTRATION.md)).

**Prevalence.** The share of windows that are artifacts: 3.8 % in UTD.

**Quantisation step.** The smallest change a digital sensor can record. Values are
rounded to multiples of it, which appears as small steps in the signal.

**Random forest.** A model that averages many decision trees; used here as a
robustness check.

**Replicate.** One complete synthetic data set for one grid cell, generated with its
own seed. Each cell is repeated *R* times.

**Semi-synthetic data.** Data in which some parts are real (here the acceleration)
and others are generated (here the EDA and its artifacts).

**Sensitivity and specificity.** Sensitivity is the share of artifact windows that are
flagged; specificity is the share of clean windows that are not flagged.

**Stage.** A part of the UTD protocol: relaxation, physical stress, cognitive stress
or emotional stress. Relaxation stages are treated as *rest*, all others as *load*.

**Tonic level.** The slowly changing baseline of skin conductance, on a scale of tens
of seconds to minutes.

**Two-hand protocol.** A planned laboratory study in which EDA is recorded on both
wrists while only one hand moves. The still hand serves as a reference for judging
changes on the moving hand. It is a reference with its own errors, not a perfect
ground truth.

**UTD.** The PhysioNet *Non-EEG Dataset for Assessment of Neurological Status*
(Birjandtalab et al., 2016, 2017): wrist recordings of 20 participants with expert
artifact labels from Zhang et al. (2017).

**Validation checks G1–G7.** The tests of correctness, realism and controls that the
pipeline must pass before any synthetic result is interpreted (§7).

**Window.** A non-overlapping 5-second piece of a recording (40 samples at 8 Hz). Every
label, feature and prediction refers to one window.

**z-score.** A value expressed as the number of standard deviations from the mean,
(*x* − mean)/SD.
