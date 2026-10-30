"""Finish a power run that stopped before its end (STATUS.md, known problem 13).

Same code path as `python -m edasynth power --source utd`; the (condition, n) points
already printed in the log of the stopped run are read back instead of recomputed.

Usage: python scripts/power_resume.py <log of the stopped run> [--jobs N]
"""

import argparse
import ast
import re
from pathlib import Path

import pandas as pd

from edasynth import cli
from edasynth.config import key, loadConfig
from edasynth.experiment import bootstrapPower, loadBank
from edasynth.generator import Condition


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    done = {}
    for line in Path(args.log).read_text().splitlines():
        m = re.match(r"\[[^\]]+\] (.+) n=(\d+): (\{.*\})$", line)
        if m:
            done[(m.group(1), int(m.group(2)))] = ast.literal_eval(m.group(3))
    cfg = loadConfig(None)
    cfg["repoRoot"] = str(Path(cli.__file__).resolve().parents[2])
    res = cli.resultsDir(cfg)
    bank = loadBank(res / "bank-utd", cfg)
    mix = tuple(cfg["grid"]["compositions"][-1])
    gm = cli.readGrid(res / "grid-metrics-utd.csv", cfg)
    gm = gm[gm.classes == "+".join(mix)]
    by = ["c", "kappa", "deltaMax", "piR"]
    comp = pd.read_csv(res / "compatibility-utd.csv")
    comp = comp[comp.compatible].merge(gm.groupby(by, as_index=False)["dB"].mean(), on=by)
    comp = cli.selectPowerConditions(comp, cfg["power"]["nConditions"])
    rows = []
    for _, c in comp.iterrows():
        cond = Condition(c=c.c, kappa=c.kappa, deltaMax=c.deltaMax, piR=c.piR, classes=mix)
        for n in cfg["power"]["participants"]:
            pw = done.get((cond.key(), n))
            if pw is None:
                pw = bootstrapPower(bank, cond, cfg, n, cfg["power"]["reps"], nJobs=args.jobs)
                cli.log(f"{cond.key()} n={n}: {pw}")
            rows.append({"condition": cond.key(), "nParticipants": n, **{key("power", k): v for k, v in pw.items()}})
            pd.DataFrame(rows).to_csv(res / "protocol-planning-utd.csv", index=False)


if __name__ == "__main__":
    main()
