# Background

This page explains the measurement, the problem and the research gap for readers
who are new to the field. The method itself is described in [METHOD.md](METHOD.md).

## 1. What electrodermal activity measures

Electrodermal activity (EDA) is the electrical conductance of the skin. Two
electrodes touch the skin, the device applies a very small voltage or current, and
it records how easily current flows between them. Conductance is the inverse of resistance and is
measured in microsiemens (µS).

Skin conductance rises when eccrine sweat glands fill with sweat, because sweat is a
good conductor. These glands are controlled by the sympathetic branch of the
autonomic nervous system, the part that prepares the body for effort or threat.
EDA is therefore used as an indirect measure of sympathetic activity, for example
in stress research, emotion research and some clinical monitoring
(Posada-Quintero & Chon, 2020).

A recording is usually described as the sum of two parts:

- the **tonic level**, a slowly drifting baseline that changes over tens of seconds
  to minutes;
- **phasic responses** (skin conductance responses, SCRs), short rises that start
  about 1 to 5 s after a stimulus, peak after a few seconds and then decay more
  slowly. A rise of at least 0.05 µS is conventionally counted as a response
  (Kleckner et al., 2018; Posada-Quintero & Chon, 2020, citing Boucsein et al., 2012).

Separating the two parts is not straightforward, and different decomposition
methods give noticeably different tonic and phasic components for the same signal
(Posada-Quintero & Chon, 2020). This project uses one of them, cvxEDA
(Greco et al., 2016), only to *measure* properties of real recordings; see
[METHOD.md §3](METHOD.md#3-measuring-parameters-on-real-recordings).

## 2. Why motion disturbs wrist recordings

Laboratory EDA is often recorded from the fingers or palm while the person sits
still. Wristbands record from the wrist during everyday life, where the arm moves
all the time. Movement can change how the electrodes touch the skin:

- pressure on the electrode can increase the contact area and raise the
  conductance for a moment;
- a sliding or lifted electrode can reduce the contact and make the signal drop;
- a contact change that persists can shift the whole level up or down;
- complete loss of contact can produce values at or near zero.

These four mechanisms are the ones represented by the artifact classes of this
project ([METHOD.md §4.3](METHOD.md#43-artifact-classes)). How often each of them
occurs in real wrist recordings is not known.

These changes are **motion artifacts**: features of the recording that come from
the measurement, not from the body. Some of them have shapes that resemble genuine
responses. Expert labellers in two published studies used criteria such as
"a peak which does not show exponential decay", "a sudden change in EDA correlated
with motion" and, for the data used here, a sudden drop of more than 0.1 µS,
because skin conductance has no physiological mechanism for sudden decreases
(Taylor et al., 2015; Zhang et al., 2017).

Movement also has a second, opposite effect: effort and sudden actions can trigger
*real* sympathetic responses. Taylor et al. (2015) gave this as a reason for not
using acceleration in their detector: in applications such as seizure detection a
strong accelerometer signal occurs together with a genuine rise in EDA. A detector
that treats every change during movement as an artifact will therefore also discard
real physiology.

## 3. What an accelerometer adds

An accelerometer records the acceleration of the device along three perpendicular
axes, in units of g (1 g is the acceleration due to gravity, 9.81 m/s²). When the
wrist is still, the magnitude of the three-axis vector is close to 1 g because
gravity is always measured; movement adds to or subtracts from it. Many wristbands
that record EDA, including the one used for the data here, also record acceleration.

The appeal is clear: if the wrist moved at the moment the EDA changed, the change is
suspicious. The difficulty is that the accelerometer measures movement of the
device, not the state of the electrode–skin contact. A movement may leave the
contact unchanged, and a contact change could also happen without strong
acceleration, for example through slow pressure on the band.

## 4. What is known

Several groups have built artifact detectors or cleaning methods:

- Taylor et al. (2015) trained classifiers on features of 5-second EDA windows
  labelled by experts. They deliberately did not give the classifier any
  acceleration data.
- Chen et al. (2015) proposed a wavelet-based method for removing motion artifacts
  from EDA.
- Kleckner et al. (2018) proposed simple, transparent rules: EDA outside 0.05–60 µS,
  changes faster than ±10 µS/s, skin temperature outside 30–40 °C, and the 5 s
  around any of these are marked invalid.
- Hossain et al. (2022) compared EDA from several body sites and defined an
  artifact index, P<sub>n</sub>, the share of signal power above 0.4 Hz. They observed that
  walking artifacts in foot EDA were largely high-frequency and could be filtered
  out.
- Amin et al. (2024) used the accelerometer as a noise reference for reducing
  artifacts in skin conductance responses.

The study closest to this project is Zhang et al. (2017). They collected expert
labels for two public wrist data sets recorded with the Affectiva Q Sensor: UTD
(20 people in a laboratory stress protocol, 3.8 % of windows labelled as artifact)
and AWW (one person on several days, 33 % of walking and 15 % of resting windows
labelled as artifact). They computed 24 EDA features and 96 acceleration features
per 5-second window and compared detectors with and without acceleration. For
supervised detectors, adding acceleration raised the area under the ROC curve by
about 0.4 % on average; for unsupervised detectors it lowered it. They offered two
explanations: many movements seen by the experts did not change the EDA, so they may
not have affected the electrode contact, and the acceleration features may not suit
the task.

## 5. The gap this project addresses

Whether a synchronous accelerometer improves artifact detection on wrist EDA is
therefore unsettled, and the available data cannot settle it:

1. **No independent ground truth.** In UTD and AWW the artifact labels come from
   experts. One of their criteria explicitly uses movement, so a detector
   that uses movement is partly rewarded for agreeing with how the labels were
   made, not necessarily with what happened at the electrode.
2. **The decisive quantity is unmeasured.** The benefit of the accelerometer depends
   on how often, how strongly and with what delay movement disturbs the contact.
   None of the data sets used here measures this, and we know of no data set that
   does.
3. **Small data.** UTD has 349 artifact windows from 20 people. Our re-analysis
   gives a confidence interval for the gain that runs from slightly below zero to
   about +0.06 in AUPRC ([RESULTS.md](RESULTS.md)).

The project responds in two ways. It quantifies what the real data do and do not
show, and it builds semi-synthetic data in which the movement–artifact link is a
parameter. Varying that parameter replaces a question the current data cannot answer ("does
the accelerometer help?") with questions they can help to answer: under which links it helps, by how much,
at what cost, which links are consistent with the real data, and what a laboratory
study would need to measure to decide between them.

## 6. Why semi-synthetic data, and what they cannot do

Purely synthetic data would need a model of human movement, which is hard to make
realistic. Purely real data lack a ground truth. The semi-synthetic design keeps
the real acceleration recordings unchanged and generates only the EDA, using
properties measured on real EDA from the same data set. Because the generator
places every artifact, the correct label of every window is known.

The cost is that every synthetic result is conditional on the generator. If the
synthetic signals differ from real ones in ways that matter to the detectors,
the map may not describe reality. This is why the pipeline includes realism checks
(G3, G4) and control experiments (G5, G6) that must pass before any synthetic result
is interpreted ([METHOD.md §7](METHOD.md#7-validation-checks-g1g7)). It is also why
the conditional map is meant to *plan* a laboratory study, not to replace one.

## 7. Why the question matters

Wearable EDA is used outside the laboratory, where movement cannot be avoided. The
U.S. Food and Drug Administration (2018) cleared the Embrace wristband, which detects
generalized tonic–clonic seizures from EDA and motion data, as an adjunct to seizure
monitoring "during periods of rest". Knowing when acceleration genuinely helps to
separate artifacts from physiology, and when it would remove real responses, is
part of what is needed before such devices can be trusted during movement.
