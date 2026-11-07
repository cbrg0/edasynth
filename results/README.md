# Results

The files reported in [docs/RESULTS.md](../docs/RESULTS.md). The parameter bank
(`bank-utd.*`, written by `estimate`) is not committed because of its size; it is
rebuilt by `make estimate`.

| File | Command | Contents |
|---|---|---|
| `alignment.csv` | `align` | label-to-signal alignment per participant |
| `realdata.json` | `realdata --stored --raw` | real-data comparison (RESULTS §1) |
| `validation-utd.json` | `validate --pytest` | checks G1–G7 (RESULTS §2) |
| `pilot-utd-R.json` | `pilot` | number of replicates *R* |
| `grid-metrics-utd.csv` | `grid` | every cell and replicate, all measures (11 300 rows) |
| `grid-labelnoise-metrics-utd.csv` | `grid --label-noise` | mixture cells with noisy training labels |
| `variants-metrics-utd.csv` | `grid --variants` | the variants |
| `grid-summary-utd.csv` | `analyze` | per cell: mean gain, 95 % interval, decision |
| `compatibility*.csv` | `analyze` | which cells agree with the real gain |
| `sensitivity-utd.csv` | `analyze` | effect of each link parameter alone |
| `protocol-planning-utd.csv` | `power` | power by number of participants |
| `external-edabe.json` | `scripts/external_checks.py edabe` | checks E1–E3 (RESULTS §4) |
| `external-wesad*.json/csv` | `scripts/external_checks.py wesad` | check W1 |

All commands are `python -m edasynth <command> --source utd` unless noted.
