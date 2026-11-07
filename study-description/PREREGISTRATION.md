# Preregistration

The study computes 300 cells, several measures and several variants. Among so many
numbers some will look interesting by chance, and an author who chooses the rules
after seeing the numbers can find almost any result. To prevent this, the rules for
running and judging the experiment were fixed **before** the main results existed. The
moment of fixing was made public with git tags: a tag marks one exact state of the
repository with its date, and that state can no longer be changed without the change
being visible. Every later change is listed in this file with its reason. The full
original texts can be read in the repository at each tag.

| Tag | Date of the commit | What it fixed |
|---|---|---|
| `prereg-synth-v1` | 25 October 2026 | the settings, the code and the number of replicates before the main grid, in version 1.3.0 |
| `prereg-synth-v2` | 27 October 2026 | that the grid with noisy labels is computed for the mixture cells only, in version 1.3.1 |
| `prereg-synth-v3` | 28 October 2026 | the corrected generator and checks of version 1.4.0, after the grid of version 1.3.0 was withdrawn |

## Preregistered elements

1. The whole settings file, now `settings/default.yaml`, which contains the grid, the
   thresholds of all checks, the random seeds and the decision rules.
2. The code of the program at the tagged commit, now in `program-code/edasynth/`.
3. The number of replicates *R*, which the pilot run determined and wrote to
   `result-tables/pilot-utd-R.json`.

**Primary outcome.** For each cell, the mean gains Δ<sub>A4</sub> and
Δ<sub>B</sub> over the *R* replicates are reported with their 95 % *t*-intervals. A
cell is useful when the mean gain is at least 0.01 and the lower bound of the interval
is above zero. All cells are reported, and no single cell is picked out afterwards.

**Secondary outcomes** are reported whatever their values turn out to be. They are
the gain Δ<sub>A</sub>, the difference Δ<sub>B</sub> − Δ<sub>A4</sub>, the sensitivity at
95 % specificity, the share of lost genuine responses, the variants, the agreement with
the real gain of 0.027 with and without noisy labels, the effect of each parameter
alone, and the power for 5 to 20 participants.

**Validation before the grid.** The checks G1 to G7 use the thresholds written in the
settings file. If a check fails, the generator is corrected, the correction is
recorded, and all checks are run again. **A threshold is never changed to make a check
pass.**

## Changes before the first tag

Each of these changes was caused by a failed check.
No threshold was changed.

| Version | What was changed | What caused it |
|---|---|---|
| 1.2.0 | the noise is taken from blocks of the remainder of the decomposition, artifact sizes are stored relative to the signal level, participant S18 is no longer a physiology donor, and response sizes are measured in clean windows only | check G3 failed, the artifact sizes did not match the signal level, and an inspection of the data showed the other two problems |
| 1.3.0 | artifacts copy the measured course of real episodes, the tonic level follows the donor's own stages, C4 uses a control gain, equivalent cells are computed once, and the power is computed for cells at the quartiles | check G4 failed, some clean signals became negative, check G6 failed, and the full grid would have taken too long |

**Exception for G5.** This was decided on 25 October 2026. Without a
link between movement and artifacts, no feature set gains, but every set loses
slightly: FA by 0.032, FA4 by 0.002 and FB by 0.002. The rule of G5 was kept as it was,
and in the map every gain is also compared with this loss.

## Changes after the tags

**Tag v2, 27 October 2026.** The grid with noisy labels is computed for the mixture
cells only, because only these cells enter its analysis. The change was made while the
main grid was running and before any of its results had been seen.

**Tag v3, 28 October 2026, made after the summary of the grid of version 1.3.0 had been
seen.** An odd result in the variant with 15 % artifacts revealed that overlapping
artifacts multiplied each other, and the grid of version 1.3.0 was withdrawn. While
this error was being corrected, three checks that had been circular were fixed as
well. First, synthetic artifacts had been compared with the very episodes they were
built from (check G4). Second, the generator produced delays only within the five
seconds that the coupling features look at. Third, check G2 used data with exactly the
response shape that the measurement assumes. To make the new run possible in the time
available, the regularisation constant was fixed at *C* = 1, as in the comparison on
real data, and checks G3 and G4 now compare each real participant with a generator that
never saw this participant.

The following stayed unchanged: the primary outcome, the decision rule, every
threshold, the rule for agreement with the real data, the reading of gains against
the loss of G5, and the rule for *R*.

Because the summary of version 1.3.0 had been seen first, a reader cannot rule out
that it influenced these changes. None of them, however, targets a particular cell or
a particular direction of the map.

**New folder names, 7 November 2026.** After all results had been obtained, the
folders were renamed so that their content is clear from the name; for example, `src`
became `program-code`, `configs` became `settings` and `results` became
`result-tables`. The tags keep the old names. No number changed: one replicate
recomputed after the renaming equals the reported value exactly.

## Plan of the external validation

This plan was written before any synthetic data had been compared with EDABE or WESAD
and before any detector had been applied to them. At that point only the column names
of EDABE, its numbers of samples, its shares of artifacts and its value ranges had been
seen, as well as the list of files in WESAD. Version 1.4.0 stayed unchanged.

- **E1 and E2.** The realism test of checks G3 and G4 is repeated between EDABE and
  synthetic data, for artifact windows in E1 and for clean windows in E2. It uses 10
  repetitions of 20 synthetic participants in the cell with the mixture, *c* = 0.5,
  κ = 1, δ<sub>max</sub> = 1 s and π<sub>r</sub> = 0. Because devices differ in their
  signal level, the features are scaled within each recording in the main version, and
  the unscaled result is reported as well. **The check passes when the median AUC is
  at most 0.70.**
- **E3.** Three detectors are tested on EDABE by their AUPRC: (a) one trained on
  synthetic data, (b) one trained on the real UTD recordings, and (c) one trained on
  EDABE itself. **The check passes when (a) is at least (b) minus 0.05.**
- **W1.** The map is computed with the 15 WESAD wrist recordings as movement donors,
  for *c* equal to 0, 0.5 and 1, with 50 replicates. **The check passes when the
  decision about the usefulness of FB is the same as with the UTD movements in all
  three cells.**
- **Preparation of the data.** Both data sets are averaged down to 8 samples per
  second. Zeros at the start and the end of each EDABE recording are removed; this rule
  was added after the first three rows of one file had been seen, before any
  comparison.

**Limitation.** The commits with this plan are dated 30 October 2026
at 21:34 and 21:35, hours before the runs. However, they reached GitHub together with
the results on 31 October. Their timing therefore rests on the dates written in the
commits, which a reader cannot verify independently, unlike the dates of the tags.
