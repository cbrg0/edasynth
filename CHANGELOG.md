# Changelog

One line per version; the full entries are in the git history. Changes to the method
after a preregistration tag are released under a new tag
([docs/PREREGISTRATION.md](docs/PREREGISTRATION.md)).

| Version | Date (2026) | Change |
|---|---|---|
| – | 7 Nov | Colab notebook `notebooks/quickstart.ipynb`: the map and all checks recomputed from the committed tables, the real-data result from the public data, the unit tests. |
| 1.4.1 | 7 Nov | Repository shortened: documentation condensed into README, METHOD, RESULTS and PREREGISTRATION; unused code removed (untrained reference detectors, `inspect`, the one-person development run, a one-off power script); reported results and figures committed; summary figure F0. No change to any computation: a grid replicate recomputed after the cut matches the reported value exactly. |
| – | 31 Oct | Checks on EDABE and WESAD, planned before running (all pass). |
| – | 30 Oct | Results of the 1.4.0 grid, variants, label-noise grid, power and figures. |
| 1.4.0 | 28 Oct | Overlapping artifacts no longer compound (contact state, clean reference level, range limits); realism checks without circularity; G2 with varied response shapes; δ<sub>max</sub> up to 10 s; *C* = 1. Tag `prereg-synth-v3`; the 1.3.0 grid is withdrawn. |
| 1.3.2 | 28 Oct | `power` uses all requested cores. |
| 1.3.1 | 27 Oct | Label-noise grid on the mixture cells only. Tag `prereg-synth-v2`. |
| 1.3.0 | 25 Oct | Artifacts follow the measured course of real episodes; tonic from the donor's own stages; C4 control gain; equivalent cells computed once. Tag `prereg-synth-v1`. |
| 1.2.1 | 25 Oct | A results file made with another parameter bank is refused. |
| 1.2.0 | 25 Oct | Residual-block noise; relative artifact sizes; low-EDA donors excluded; response sizes from clean windows. |
| 1.1.x | 23 Oct | Alignment of labels to raw signals; reader checks; fixes found on the development run. |
| 1.0.0 | 23 Oct | First complete pipeline. |
