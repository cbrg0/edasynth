# Preregistration of the synthetic experiment

## Why this file exists

The main experiment produces 200 cells, several measures and several variants. With
that many numbers, some will look interesting by chance. To prevent the analysis from
being shaped by the results, the rules for running and judging the experiment are
fixed **before** the main grid is run, and the moment of fixing is made public with a
git tag pushed to GitHub. The tag's date shows that the rules existed before the
results.

**Current state: draft.** The tag `prereg-synth-v1` does not exist yet. Several design
questions are still open ([STATUS.md](STATUS.md#open-design-decisions-to-settle-before-the-tag-prereg-synth-v1));
they have to be settled, and this file updated, before the tag is created.

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

## Changes after the tag

Any change after the tag is made openly: a new tag (`prereg-synth-v2`, and so on) with
the reason in `CHANGELOG.md`. Any analysis designed after seeing results is labelled
*post hoc* wherever it is reported.
