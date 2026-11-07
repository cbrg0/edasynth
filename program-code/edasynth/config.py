"""This file reads the settings of the study and is used by every other step.

All numbers that can change a result, such as thresholds, random seeds and the values of
the grid, are kept in the file settings/default.yaml. This file loads them so that the
rest of the program takes every number from one place and never from its own code.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

repoRoot = Path(__file__).resolve().parents[2]


def loadConfig(path: str | Path | None = None, overrides: dict | None = None) -> dict:
    path = Path(path) if path else repoRoot / "settings" / "default.yaml"
    with path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if overrides:
        cfg = merge(cfg, overrides)
    cfg["configPath"] = str(path)
    return cfg


def merge(base: dict, extra: dict) -> dict:
    out = copy.deepcopy(base)
    for key, val in extra.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


def key(*parts) -> str:
    words = [str(p) for p in parts if str(p)]
    return words[0] + "".join(w[:1].upper() + w[1:] for w in words[1:]) if words else ""


def resolve(cfg: dict, key: str) -> Path:
    p = Path(cfg["paths"][key])
    return p if p.is_absolute() else repoRoot / p


def get(cfg: dict, dotted: str, default: Any = None) -> Any:
    node: Any = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node
