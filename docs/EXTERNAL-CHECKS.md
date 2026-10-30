# Checks on independent data (EDABE and WESAD)

This plan was written on 30 October 2026, before any synthetic data were compared
with these data sets and before any detector was applied to them. At that point only
the following had been seen: the column names of EDABE, the number of samples, the
share of artifact samples and the range of the raw signal per EDABE file
(`edabe-8hz-info.txt`), and the file list, one questionnaire and one `info.txt` of
WESAD. The code and the generator of version 1.4.0 are used unchanged. Results are
reported in [RESULTS.md §5](RESULTS.md#5-checks-on-independent-data) whatever they are.

## Why

All checks so far used UTD only: the generator was built from UTD, and its realism
was tested against UTD. A generator that only reproduces the recordings it was built
from could pass those checks and still be unrealistic for other devices. Two public
data sets allow independent checks:

| Data set | Device, site | People | What it has | What it lacks |
|---|---|---|---|---|
| EDABE (Llanes-Jurado et al., 2023, doi:10.1016/j.eswa.2023.120581) | Shimmer3 GSR+, 128 Hz | 43 recordings, about 74 h | raw EDA, EDA corrected by an expert, artifact label per sample (one expert per file) | acceleration |
| WESAD (Schmidt et al., 2018, doi:10.1145/3242969.3242985) | Empatica E4 on the wrist, EDA 4 Hz, acceleration 32 Hz | 15 | wrist EDA and acceleration, protocol stages (baseline, stress, amusement, meditation) | artifact labels |

Without acceleration, EDABE cannot test whether the accelerometer helps. It can test
whether the synthetic artifacts look like artifacts marked by other experts on another
device. Without artifact labels, WESAD cannot test detection. It can show whether the
map changes when the movement comes from another device and another protocol.

## Preparation (done on macOS by the user, scripts in `scripts/`)

- **EDABE**: each file averaged over blocks of 16 samples (128 Hz → 8 Hz, the rate of
  the UTD data). A block counts as artifact when any of its samples is labelled
  (`binarytarget`). Columns kept: `rawdata`, `cleandata`, `binarytarget`. The unit is
  taken to be µS (**assumption**; to be checked against the EDABE paper).
- **EDABE start and end**: the first rows of a file show `rawdata` = 0 while the expert
  signal is not zero, which looks like the device before recording started. Leading and
  trailing samples with `rawdata` = 0 are therefore removed **[decision, added after
  seeing the first three rows of one file and before any comparison]**; zeros inside a
  recording are kept (they are signal loss, class C5).
- **WESAD**: from `S*.pkl`, the wrist acceleration (32 Hz, units of 1/64 g as stored by
  the E4, **assumption** to be checked against the WESAD readme), the wrist EDA (4 Hz)
  and the protocol label, taken every 175th sample of the 700 Hz label (4 Hz). The
  `.pkl` signals are synchronised by the WESAD authors; the separate E4 files are not
  synchronised with the protocol times.

Windows are 5 s long and do not overlap, as for UTD. A window is positive when any of
its samples is an artifact (the main rule of the grid).

## E1, E2: realism against EDABE

The same classifier two-sample test as G3 and G4 (logistic regression, 5-fold
cross-validation grouped by recording, AUC), comparing EDABE windows with synthetic
windows:

- **E1**: EDABE artifact windows against synthetic artifact windows.
- **E2**: EDABE clean windows against synthetic clean windows.

Synthetic windows: 10 repetitions, each with 20 synthetic participants from the full
UTD bank, mixture C1+C2+C3+C5, *c* = 0.5, κ = 1, δ<sub>max</sub> = 1 s,
π<sub>r</sub> = 0, prevalence 3.8 % (the setting of G3 and G4). Features: the 24 EDA
features F0.

The two devices and populations differ in level, so the **main version z-scores each
feature within each recording** (real and synthetic alike), as the published UTD
features are; the unscaled version is also reported. **Criterion [decision]: median
AUC over the 10 repetitions at most 0.70 in the main version**, the limit of G3 and G4.

## E3: does a detector trained on synthetic data find EDABE artifacts?

Three detectors with the 24 EDA features (z-scored within recording) and logistic
regression with *C* = 1, each scored on all 43 EDABE recordings by AUPRC:

- **(a) synthetic**: trained on one synthetic replicate of the E1 setting; 10 replicates,
  median reported;
- **(b) real UTD**: trained on the 20 raw UTD recordings with the majority expert label;
- **(c) EDABE itself**: leave one recording out (reference for what the features allow).

**Criterion [decision]: AUPRC of (a) at least AUPRC of (b) minus 0.05.** That is,
training on synthetic data should not be clearly worse than training on the real UTD
data from which the generator was built. The share of artifact windows (the AUPRC of
guessing) is reported with them.

## W1: the map with WESAD movement

The movement donors of UTD are replaced by the 15 WESAD wrist recordings; the
physiology donors stay those of UTD. Preparation of each WESAD recording:

- acceleration averaged over blocks of 4 samples (32 → 8 Hz) and divided by 64 to give g;
- cropped from the first to the last sample with a protocol label 1 to 4;
- stage "rest" for baseline (1) and meditation (4), "load" for stress (2), amusement (3)
  and every other sample in between;
- movement events found by the unchanged `motionEvents`, with the threshold from the
  "rest" samples.

Cells: mixture, κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0, *c* ∈ {0, 0.5, 1};
50 replicates each, 15 synthetic participants per replicate (one per movement donor).
They are compared with the same cells of the UTD grid.

**Criterion [decision]: the usefulness decision for FB (mean Δ ≥ 0.01 and lower 95 %
bound > 0) is the same as on UTD in all three cells.** The differences of Δ<sub>B</sub>
and Δ<sub>A4</sub> are reported. With 15 instead of 20 participants the intervals are
wider; this is a known difference, not a criterion.

## What the checks cannot show

- EDABE labels are corrections by experts who saw the EDA only; like UTD labels, they
  are not an independent ground truth of contact disturbances.
- W1 changes the source of movement but keeps the assumed link between movement and
  artifacts; it cannot show what *c* is on the E4.
- If E1 or E3 fails, the map is not withdrawn; the failure is reported with the map and
  limits how far its numbers transfer to other devices.
