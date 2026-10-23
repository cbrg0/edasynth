# Preregistration of the synthetic experiment

Status: **draft until the git tag `prereg-synth-v1` exists on GitHub.** The tag's
date is the public proof that everything below was fixed before the main grid
(`make grid`) was run. Validation (`make validate`) and the pilot (`make pilot`)
are run before tagging; the grid, variants, analysis and power runs after it.

## Fixed by the tag

1. `configs/default.yaml` in full: grid, thresholds, seeds, model settings,
   feature parameters (including the four coupling features F_B), decision rule,
   compatibility rule, sensitivity reference, power settings.
2. The code in `src/edasynth/` at the tagged commit.
3. The number of replicates R written by the pilot to `results/pilot_utd_R.json`.

## Primary outcome

For each of the 200 cells: mean Δ_A4 and Δ_B (AUPRC difference against F0,
logistic regression, leave-one-participant-out) over R replicates with a 95 %
t-interval. A cell counts as "accelerometer useful" when mean Δ ≥ 0.01 and the
lower bound > 0. All cells are reported; no cell is singled out post hoc as a
finding.

## Secondary outcomes (all pre-specified)

* Δ_A (96 features, Zhang's setting) and Δ_B − Δ_A4 (kind of features at equal count).
* Random forest Δ on the same cells (robustness).
* AUROC differences (comparison with Zhang et al., qualitative only).
* Sensitivity at specificity 0.95; share of real responses flagged as artifacts
  (spontaneous and motion-evoked separately); share of each artifact class detected.
* Variants: V2 multiplicative contact loss, 15 % prevalence, majority labels, C4 only.
* Goal 5: compatibility of each cell with the real Δ = 0.027; label-noise variant.
* Goal 6: one-factor sensitivity ranges; power for 5, 10, 15, 20 participants.

## Gates that must pass before the grid

G1-G7 as defined in `docs/METHOD.md` §8 with the thresholds in the config. If a
gate fails, the generator or estimation is corrected, the correction is recorded
in `CHANGELOG.md`, and all gates are rerun. Thresholds are not changed to make a
gate pass.

## Changes after the tag

Made openly: a new tag `prereg-synth-v2` with the reason in `CHANGELOG.md`. Any
analysis designed after seeing results is labelled *post hoc* wherever reported.
