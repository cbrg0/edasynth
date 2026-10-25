# Preregistration of the synthetic experiment

## Why this file exists

The main experiment produces 200 cells, several measures and several variants. With
that many numbers, some will look interesting by chance. To prevent the analysis from
being shaped by the results, the rules for running and judging the experiment are
fixed **before** the main grid is run, and the moment of fixing is made public with a
git tag pushed to GitHub. The tag's date shows that the rules existed before the
results.

**Current state: complete for version 1.3.0.** The design questions are settled
(below and in [STATUS.md](STATUS.md#open-design-decisions-to-settle-before-the-tag-prereg-synth-v1)),
the checks have been run on UTD, and the pilot has fixed *R* = 40. The tag
`prereg-synth-v1` has not been created yet.

## Order of work

1. Before the tag: data checks (`inspect`, `align`), parameter estimation
   (`estimate`), the validation checks G1–G7 (`validate`), the untrained detectors and
   the pilot that sets the number of replicates (`pilot`). This is `make before-tag`.
2. The tag `prereg-synth-v1` is created and pushed.
3. After the tag: the main grid, the variants, the label-noise runs, the analysis, the
   power analysis and the figures. This is `make after-tag`.

The validation checks and the pilot are run before the tag because they test the
generator and determine the number of replicates; they do not use the outcome of the
main grid.

## What the tag fixes

1. `configs/default.yaml` in full: grid, thresholds of the checks, seeds, model
   settings, feature parameters, decision rule, compatibility rule, the reference cell
   of the sensitivity analysis and the power settings.
2. The code in `src/edasynth/` at the tagged commit.
3. The number of replicates *R* written by the pilot to `results/pilot-utd-R.json`.

## Primary outcome

For each of the 200 cells: the mean gains Δ<sub>A4</sub> and Δ<sub>B</sub> (AUPRC with
the feature set minus AUPRC with EDA features only; logistic regression;
leave-one-participant-out) over *R* replicates, with a 95 % *t*-interval. A cell counts
as **accelerometer useful** when the mean gain is at least 0.01 and the lower bound of
the interval is above 0. All cells are reported; no single cell is presented as the
finding after the results are known.

## Secondary outcomes

All of these are specified in advance and reported whatever their values:

- Δ<sub>A</sub> (the 96 acceleration features of Zhang et al., 2017) and
  Δ<sub>B</sub> − Δ<sub>A4</sub> (kind of acceleration information at equal number of
  features);
- the random-forest gains in the same cells, as a robustness check;
- differences in the area under the ROC curve, for qualitative comparison with Zhang
  et al. (2017);
- sensitivity at 95 % specificity; the share of genuine responses flagged as artifacts,
  separately for spontaneous and movement-evoked responses; the detection rate per
  artifact class;
- the variants: multiplicative contact loss (V2), 15 % prevalence, majority labelling
  rule, C4 only;
- goal 5: whether each cell is compatible with the real gain Δ = 0.027, with and
  without the label-noise variant;
- goal 6: the one-factor sensitivity ranges and the power for 5, 10, 15 and 20
  participants.

Definitions of all measures are in [METHOD.md §6](METHOD.md#6-models-and-performance-measures).

## Checks that must pass before the main grid

The checks G1–G7 of [METHOD.md §7](METHOD.md#7-validation-checks-g1g7), with the
thresholds in `configs/default.yaml`. If a check fails, the generator or the estimation
is corrected, the correction is recorded in `CHANGELOG.md`, and all checks are run
again. **Thresholds are not changed to make a check pass.**

## Decisions taken after the first checks on UTD, before the tag

The checks were run on UTD before the tag, as planned, and several failed. The
corrections below were made before any main-grid result existed; each is recorded in
`CHANGELOG.md` (versions 1.2.0 to 1.3.0) with the result that led to it. They were
chosen after seeing the checks, so they are fitted to the checks in the sense that a
check that failed was the reason for each change; no threshold of a check was changed.

| Version | Change | Check that led to it |
|---|---|---|
| 1.2.0 | noise from blocks of the decomposition residual | G3 (σ = 0 on UTD) |
| 1.2.0 | artifact sizes relative to the level | sizes measured at one level added to another |
| 1.2.0 | physiology donors with EDA mostly below 0.05 µS excluded (S18) | inspection |
| 1.2.0 | response sizes only from clean windows | heights larger than the recording range |
| 1.3.0 | C1–C3 artifacts use the measured course of real episodes | G4 |
| 1.3.0 | tonic follows the donor's own level from stage to stage, floored at 0.05 µS | negative clean signals found while checking G4 |
| 1.3.0 | control gain for C4 | G6 (C4 produced no artifacts with the measured gain) |
| 1.3.0 | equivalent grid cells computed once; power cells chosen by quartiles of Δ<sub>B</sub>; label-noise compatibility analysed | computing cost; known problems 4 and 7 |

**G5 is an accepted exception.** With 1.3.0 every other check passes on UTD. G5 fails
because the feature sets lower AUPRC slightly when movement and artifacts are
unrelated: FA by 0.032, FA4 by 0.002 and FB by 0.002 (no link); no spurious gain was
found. The rule of G5 was kept unchanged (25 October 2026), and the main grid is run
although G5 fails in this way. The consequence is fixed here: in the map, a gain of a
feature set is also compared with the loss of the same set in the no-link control of
G5, and this comparison is reported next to the decision rule. For FB and FA4 the
loss is five times smaller than the 0.01 that the decision rule requires.

**Kept unchanged on purpose:** the τ<sub>0</sub> range of 2–4 s (12 of 20 participants
at an edge; G2 passes, and a wider range has no source in this project), and the
episode class rules (71 % of UTD episodes are C3; changing the rules after seeing
the shares would be fitted to the data, and the map is also reported per class).

## Changes after the tag

Any change after the tag is made openly: a new tag (`prereg-synth-v2`, and so on) with
the reason in `CHANGELOG.md`. Any analysis designed after seeing results is labelled
*post hoc* wherever it is reported.
