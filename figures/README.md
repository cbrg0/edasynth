# Figures folder

Figures are drawn by `python -m edasynth figures` from the files in `results/`.
F1 and F2 generate a fresh example from the parameter bank with a fixed seed; all
other figures plot stored numbers only. Each figure is written as
PDF and PNG. None is committed to the repository.

| File | Shows | Goal |
|---|---|---|
| `F1-real-vs-synthetic` | a real recording next to a synthetic participant whose physiology comes from that person: EDA, acceleration magnitude and labels | realism, by eye |
| `F2-artifact-classes` | each artifact class added to the same clean segment | the generator |
| `F3-realism` | for each EDA feature, how differently it is distributed in real and synthetic clean windows (Kolmogorov–Smirnov D), with the classifier AUC of check G3 | realism (G3) |
| `F4-map-<composition>` | mean gain of F0+FA4 and F0+FB over *c* and κ, one panel row per δ<sub>max</sub> and π<sub>r</sub>; an asterisk marks cells where the accelerometer counts as useful | the map (goal 3) |
| `F5-by-class` | gain against *c* for each single artifact class | goal 3 |
| `F6-cost` | share of genuine movement-evoked responses flagged as artifacts, by feature set | cost (goal 4) |
| `F7-detectors-by-class` | share of artifact windows flagged by each untrained detector, by class | blind spots (goal 2) |
| `alignment/step2-S01` … | one plot per UTD participant: EDA, acceleration, the three experts' labels and stage boundaries | label alignment |

`figures/smoke/` holds the figures of the development run on AWW, which are not
results.
