# Preregistration

With 300 cells, several measures and variants, some numbers will look interesting by
chance. So the rules for running and judging the experiment were fixed **before** the
main results, and the moment of fixing was made public with git tags. Every later
change is listed here with its reason. The full original texts are in the repository
at each tag.

| Tag | Commit date | What it fixed |
|---|---|---|
| `prereg-synth-v1` | 25 Oct 2026 | configuration, code and number of replicates before the main grid (version 1.3.0) |
| `prereg-synth-v2` | 27 Oct 2026 | label-noise grid on the mixture cells only (1.3.1) |
| `prereg-synth-v3` | 28 Oct 2026 | the corrected generator and checks (1.4.0), after the 1.3.0 grid was withdrawn |

## What is fixed

1. `settings/default.yaml` in full: grid, thresholds, seeds, model settings, features,
   decision and compatibility rules, sensitivity reference cell, power settings.
2. The code in `program-code/edasynth/` at the tagged commit.
3. The number of replicates *R* from the pilot (`result-tables/pilot-utd-R.json`).

**Order of work.** Before the tag: data checks, estimation, checks G1–G7, pilot
(`make before-tag`). After the tag: grid, variants, label-noise runs, analysis, power,
figures (`make after-tag`). The checks and the pilot test the generator; they do not
use the outcome of the grid.

**Primary outcome.** For each cell, the mean Δ<sub>A4</sub> and Δ<sub>B</sub> over *R*
replicates with a 95 % *t*-interval. A cell is *useful* when the mean is ≥ 0.01 and the
lower bound > 0. All cells are reported; no single cell is singled out afterwards.

**Secondary outcomes,** reported whatever their values: Δ<sub>A</sub> and
Δ<sub>B</sub> − Δ<sub>A4</sub>; sensitivity at 95 % specificity; lost genuine responses;
the variants (V2, 15 % prevalence, majority rule, C4); compatibility with the real
Δ = 0.027 with and without label noise; sensitivity ranges and power for 5–20
participants.

**Checks before the grid.** G1–G7 with the thresholds in the configuration. If a check
fails, the generator is corrected, the correction is recorded, and all checks are run
again. **Thresholds are never changed to make a check pass.**

## Changes made after seeing check results, before the first tag

Each was triggered by a failed check (see RESULTS.md §2); no threshold changed.

| Version | Change | Triggered by |
|---|---|---|
| 1.2.0 | residual-block noise; relative artifact sizes; S18 excluded as physiology donor; response sizes from clean windows | G3; level mismatch; inspection |
| 1.3.0 | measured artifact courses; tonic from the donor's own stages; C4 control gain; equivalent cells computed once; power cells at the quartiles | G4; negative signals; G6; computing cost |

**G5 is an accepted exception** (decided 25 October 2026): without a link every feature
set loses slightly (FA −0.032, FA4 −0.002, FB −0.002) and none gains. The rule was kept;
in the map a gain is also compared with that no-link loss.

**Kept on purpose:** the τ₀ range of 2–4 s and the episode class rules, although 71 %
of UTD episodes fall into C3; changing them after seeing the shares would fit the data.

## Changes after the tags

**v2 (27 October).** The label-noise grid is computed for the mixture cells only,
because only they enter its analysis. Made while the main grid was running and before
any of its results were seen.

**v3 (28 October), made after the summary of the 1.3.0 grid had been seen.** An odd
result in the 15 % variant revealed that overlapping artifacts multiplied; the 1.3.0
grid was withdrawn. While correcting this, three circular checks were also fixed:
synthetic artifacts had been judged against the episodes they were built from (G4),
the generator only produced delays inside the 5-s window of the FB features, and G2
used data with exactly the response shape the estimation assumes. To make the rerun
feasible, *C* = 1 is fixed as in the real-data comparison.

Changes: clean reference level and contact state for artifacts, range limits;
δ<sub>max</sub> = 10 s added (300 cells); fixed *C*; G2 also with varied shapes; G3/G4
with a generator that never sees the compared participant, median over 10 repetitions
(chosen after seeing a split-half result and the variation between seeds).

Unchanged: primary outcome, decision rule, every threshold, compatibility rule, the
reading of gains against the G5 loss, the rule for *R*.

Because the 1.3.0 summary was seen first, a reader cannot rule out that it influenced
these changes. None of them targets a particular cell or direction of the map.

**Folder names (7 November 2026).** After all results were obtained, the folders were
renamed so that their content is clear from the name, for example `src` became
`program-code`, `configs` became `settings` and `results` became `result-tables`. The
tags keep the old names. No number changed: a replicate recomputed after the renaming
equals the reported value exactly.

## Plan of the checks on independent data (30 October 2026)

Written before any synthetic data were compared with EDABE or WESAD and before any
detector was applied to them; at that point only the EDABE column names, sample counts,
artifact shares and value ranges, and the WESAD file list had been seen. Version 1.4.0
unchanged.

- **E1, E2:** the realism test of G3/G4 between EDABE and synthetic artifact (E1) and
  clean (E2) windows; 10 repetitions of 20 synthetic participants (mixture, *c* = 0.5,
  κ = 1, δ<sub>max</sub> = 1 s, π<sub>r</sub> = 0). Features z-scored within each
  recording (main version) because devices differ in level; unscaled also reported.
  **Criterion: median AUC ≤ 0.70.**
- **E3:** AUPRC on EDABE of a detector trained on synthetic data (a), on real UTD (b)
  and on EDABE itself (c). **Criterion: (a) ≥ (b) − 0.05.**
- **W1:** the map with the 15 WESAD wrist recordings as movement donors, *c* ∈ {0, 0.5,
  1}, 50 replicates. **Criterion: the usefulness decision for FB is the same as with
  UTD movement in all three cells.**
- Preparation: EDABE averaged to 8 Hz (a block is an artifact if any sample is
  labelled); leading and trailing zeros removed (added after seeing the first three rows
  of one file, before any comparison); WESAD acceleration in 1/64 g averaged to 8 Hz,
  stages rest (baseline, meditation) and load (others).

**Limit of this evidence.** These plan commits are dated 30 October 2026, 21:34 and
21:35, hours before the runs, but they reached GitHub together with the results on
31 October. Their timing therefore rests on the commit dates, which a reader cannot
verify independently, unlike the tags above.
