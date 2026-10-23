# Method, step by step, with the code that implements it

This file is the technical counterpart of the methods document. Each step names
the module and command that implement it and marks every number as either
**[source]** (taken from the cited work), **[measured]** (estimated from real data
by this code) or **[decision]** (chosen by this project, fixed in
`configs/default.yaml` before the main run).

## 0. Question and what synthetic data can answer

*Question.* Does a synchronous wrist accelerometer improve the detection of motion
artifacts in wrist EDA compared with a model that uses EDA only?

Synthetic data cannot answer this for real people: the benefit depends on how
often wrist movement actually disturbs the electrode contact, and no available
data set measures that. The study therefore builds a **conditional map**: the
movement-contact link is varied over a grid, and for each cell the gain in AUPRC
with its confidence interval is computed. The map's axes follow the two
explanations Zhang et al. (2017) gave for the small gain they observed:
(1) the kind of acceleration features, (2) how strongly movement is tied to
contact changes. Six goals follow (README, "Goals"); goal 1 (validation) gates all
others.

## 1. Data — `edasynth.data`, `edasynth.io_wfdb`, `make data`

* UTD raw signals: PhysioNet `noneeg` 1.0.0, 20 participants, 8 Hz, ax/ay/az/temp/EDA
  [source: Birjandtalab et al., 2016, 2017]. Read by our format-16 reader, which
  checks every channel against the header's initial value and checksum.
* Expert labels: Zhang et al. repository, majority of 3 experts per 5-s window
  [source]. 9 211 windows, 349 artifacts, Fleiss κ = 0.619 [measured from the files].
* AWW: one person, used for checking feature code and for the smoke run only.

## 2. Linking labels to signals — `edasynth.alignment`, `make align`

Window *i* = samples 40*i* ... 40*i*+39 [source: Zhang's script]. Group *g* ↔
`Subject(g+1)` is a working hypothesis checked per participant: equal window
counts and r ≥ 0.999 [decision] between recomputed and published amplitude
features (z-scored within participant). Otherwise shifts of ±1-10 windows and
0-39 samples are searched; failures are excluded with the reason recorded.

## 3. Parameters from real data — `edasynth.estimation`, `make estimate`

| Quantity | How | Status |
|---|---|---|
| clean windows | all 3 experts "clean" and both neighbours too (5-s transition zone of Kleckner et al., 2018) | [decision] |
| noise σ | median(\|d1\|)/0.6745 of level-1 Haar details in clean rest windows (Donoho & Johnstone, 1994) | [measured] |
| quantisation q | most frequent non-zero \|Δy\| | [measured] |
| decomposition | z-score record; cvxEDA with τ1 = 0.7 s, α = 8e-4, γ = 1e-2, 10-s knots; τ0 ∈ {2.0, 2.1, ..., 4.0} s, chosen by the smallest cvxEDA objective value | model and settings [source: Greco et al., 2016]; τ0 grid and selection rule [decision] |
| responses | peaks of the phasic component r with trough-to-peak height A ≥ 0.05 uS; area for the generator a = A/h_max(τ0) (runs of non-zero driver p were rejected: the driver is dense on noisy data, see CHANGELOG) | [source] threshold; detection rule [decision]; values [measured] |
| response rate | λ = K/T per stage type (rest, load) in clean windows (Poisson MLE) | [measured] |
| plausibility | warning (not a gate) when a response rate exceeds 10/min, τ0 sits on the grid edge or the solver's KKT violation exceeds 1e-3 | [decision] |
| tonic shape | s0 + b·t + A_s·sin(2πt/T_t + φ), T_t ∈ [45, 90] s, per protocol stage, fitted to the signal minus the detected responses (t + r − r_events) so that the generator's level matches the real one | form and target [decision; a published source for the sinusoidal term is still to be checked]; values [measured] |
| artifact episodes | runs of majority-artifact windows; baseline = line between 2-s medians before/after; amplitude, slopes, 10-90 % rise and return, power share > 0.4 Hz, residual shift | [measured]; class rules [decision] |
| motion events | jerk j = fs‖a[n]−a[n−1]‖ above the participant's 99th rest percentile for ≥ 0.5 s, merged if < 1 s apart | [decision] |
| UTD reference point | share of episodes starting −1...+5 s after a motion onset; share of motion events followed by an artifact | [measured]; describes the labels, not reality |

The numpy decomposition backend solves the same convex problem as the reference
cvxEDA code (Greco et al., 2016); a test compares their objective values when
that code is on the Python path. Until then the agreement is not verified.

## 4. Generator — `edasynth.generator`

Synthetic participant j = real acceleration, stages and motion events of motion
donor u(j) + physiology of donor d(j) ≠ u(j) (derangement).
Signal built at 64 Hz and block-averaged to 8 Hz [decision; the ARMA form at 8 Hz
deviates 6-7 % from the exact Bateman shape].

y0 = s + r + ε; r = Σ a_k h(t − t_k) + Σ_{e∈E_r} a_e h(t − t_e − L_e);
spontaneous t_k ~ Poisson(λ of the donor and stage type); areas resampled from the
donor; E_r: each motion event with probability π_r ∈ {0, 0.5} [decision];
L_e ~ U(1, 5) s [source: Boucsein, 2012, as cited in Posada-Quintero & Chon, 2020;
applying it to movement is an assumption]; ε ~ N(0, σ²); output rounded to q.

Artifacts (V1 additive, y = y0 + ΣA_m):

| Class | Shape | Basis |
|---|---|---|
| C1 pressure peak | triangle (linear rise, linear fall) | Chen et al. (2015); Taylor et al. (2015) |
| C2 contact drop | downward trapezoid, edges w ∈ {1/8, 1} s | Taylor et al. (2015); Zhang et al. (2017) |
| C3 level shift | linear ramp to a lasting offset | contact change |
| C4 walking oscillation | G·(‖a‖ − 2-s mean of ‖a‖) inside a motion event; positive control | Hossain et al. (2022) |
| C5 dropout / saturation | value replaced by 0 or the upper limit | Taylor et al. (2015); Kleckner et al. (2018) |

Shapes are the simplest consistent with the descriptions [decision]; their
numbers come from UTD episodes [measured]. V2 (C2, C3 only): y = y0(1 − ρ), ρ the
lost contact fraction (assumption, to be tested on the bench).
Prevalence: artifacts are added until 3.8 % of windows are positive (15 % variant).
True label: sample corrupted if |y − y0| ≥ 0.05 uS; window positive if any sample
is (majority rule as a variant).

## 5-6. Feature sets and the movement-contact link

Axis 1 (feature kind): F0 (24 EDA), F0+FA (96 acc, Zhang's setting), F0+FA4
(4 acc, earlier real-data analysis), F0+FB (4 coupling features: lagged
correlation of jerk and EDA speed, its lag, co-occurrence of unusual EDA speed and
jerk, EDA-|a| coherence 0.5-3 Hz; all causal over 5 s of past + the window).

Axis 2 (link): each artifact of C1/C2/C3/C5 is caused by a random motion event
with probability c; its onset is t_e + U(0, δ_max); its size is multiplied by
(I_e / median I)^κ.

Grid: c ∈ {0, .25, .5, .75, 1} × κ ∈ {0, 1} × δ_max ∈ {1, 5} s × π_r ∈ {0, .5}
× composition ∈ {C1, C2, C3, C5, mixture} = 200 cells [decision]. Variants on part
of the grid: C4 only, V2, 15 % prevalence, majority labels.

## 7. Models — `edasynth.models`

Logistic regression (L2, balanced class weights, standardisation inside the
pipeline, C ∈ {0.01, 0.1, 1, 10} chosen by grouped 5-fold CV inside each
training set [decision: inner LOSO would multiply run time by ~4]) is primary;
random forest (500 trees) checks robustness. Outer evaluation is always LOSO.
Δ_X = AUPRC(F0 + F_X) − AUPRC(F0) on pooled out-of-fold predictions.
Untrained EDA-only detectors (goal 2): slope rule, Kleckner rules 1, 2, 4,
Kocielnik relative-rate rule, P_n (Hossain et al., 2022), SWT with Gaussian-mixture
thresholds (Chen et al., 2015).

## 8. Validation gates — `edasynth.validation`, `make validate`

G1 unit tests; G2 parameter recovery; G3/G4 classifier two-sample tests
(AUC ≤ 0.70 [decision]); G5 negative controls (c = 0; acceleration circularly
shifted ≥ 60 s); G6 positive controls (C4 only; c = 1, κ = 1, δ_max = 1 s);
G7 our features on real UTD reproduce the published-feature baseline in direction
and size (|ΔAUPRC| ≤ 0.10 [decision]).

## 9. Experiment and analysis — `edasynth.experiment`, `make pilot grid analyze power`

One replicate = 20 synthetic participants (each UTD acceleration once). Pilot:
10 replicates in 4 cells → R = ⌈(1.96 s_Δ/0.005)²⌉ ≤ 50. Per cell: mean Δ and
95 % t-interval over R replicates; decision rule "useful" = mean ≥ 0.01 and lower
bound > 0 [decision]; the whole map is reported. Goal 5: a cell is compatible with
UTD if the real Δ = 0.027 lies in the central 95 % of its replicates (the main
replicates already have UTD's size and prevalence); a label-noise variant flips
labels to Cohen κ ≈ 0.62 [assumption]. Goal 6: one-factor ranges of mean Δ around
(c = .5, κ = 1, δ_max = 1, π_r = 0) and power (share of replicates whose
participant-bootstrap CI excludes 0) for 5, 10, 15, 20 participants.
