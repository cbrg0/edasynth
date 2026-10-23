# Glossary (for readers outside signal processing)

**Electrodermal activity (EDA), skin conductance.** How easily a tiny current
passes through the skin, measured in microsiemens (uS). Sweat glands controlled
by the sympathetic nervous system change it, so EDA tracks arousal.

**Tonic level.** The slowly drifting baseline of skin conductance (minutes).

**Phasic response (skin conductance response, SCR).** A short rise of 1-5 s
followed by a slower return, triggered by a burst of sympathetic activity. Its
smallest visible size is usually taken as 0.05 uS.

**Bateman function.** The difference of two decaying exponentials,
h(t) = exp(-t/tau0) - exp(-t/tau1), used as the shape of one response
(tau1 = 0.7 s rise, tau0 = 2-4 s decay).

**cvxEDA.** A method (Greco et al., 2016) that splits a recording into tonic level,
phasic responses and noise by convex optimisation. We use it to *measure* real
recordings, not to clean them.

**Motion artifact.** A change in the recorded signal caused by the electrode
moving, being pressed or losing contact, not by the body. It can look like a
response, a drop, a step or a burst of fast oscillation.

**Accelerometer (ACC).** A sensor in the same wristband that records
acceleration along three axes (in g, 1 g = gravity). **Jerk** is how fast the
acceleration itself changes; it is large during sudden wrist movements.

**Window.** A 5-s piece of the recording (40 samples at 8 Hz). Every label,
feature and prediction refers to one window.

**Label / ground truth.** The answer "artifact or clean" for a window. In real
data experts give it; in synthetic data we know it exactly because we added the
artifacts ourselves.

**Synthetic / semi-synthetic data.** Data built by a program. Here the
acceleration is real (recorded from people) and the EDA is synthetic, so the true
answer is known.

**Donor.** A real participant whose recording supplies part of a synthetic
participant: one donor gives the movement, a *different* donor gives the
physiological parameters.

**Feature.** A number computed from a window (e.g. the maximum slope).
F0 = 24 EDA features, FA = 96 acceleration features, FA4 = 4 acceleration
features, FB = 4 features describing how movement and EDA change *together*.

**Logistic regression, random forest.** Two standard ways of turning features
into a probability that a window is an artifact.

**Leave-one-participant-out (LOSO).** Train on all participants but one, test on
the one left out, repeat for everyone. It tells how well the method works on a
person it has never seen.

**AUPRC (area under the precision-recall curve).** Main quality score. 1 is
perfect; a method that guesses scores about the share of artifacts (here ~0.04).
It is preferred over ROC AUC when artifacts are rare.

**Delta (Δ).** AUPRC with acceleration minus AUPRC without it. Positive means the
accelerometer helped.

**Confidence interval (95 %).** The range of values compatible with the data; if
it contains 0, the data do not show a clear gain.

**Bootstrap by participant.** Re-drawing whole participants at random many times
to see how much a result would vary with different people.

**Negative / positive control.** A run where we know there should be no effect
(movement unrelated to artifacts) or a clear effect (artifacts made from
movement). The method must pass both before its other results are trusted.

**Preregistration.** Writing down the analysis plan and decision rules, with a
public date, *before* seeing the results.

**Two-hand protocol.** The planned laboratory study: EDA on both wrists while
only one hand moves, so that the still hand helps judge whether a change on the
moving hand is physiological or mechanical. It does not give a perfect ground
truth (see the method).
