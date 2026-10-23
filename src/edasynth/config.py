"""Configuration loading.

Every number that changes a result lives in a YAML file under ``configs/``.
Code never hard-codes thresholds; it reads them from the dictionary returned here.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

repoRoot = Path(__file__).resolve().parents[2]


def loadConfig(path: str | Path | None = None, overrides: dict | None = None) -> dict:
    """Load a YAML config and apply optional nested overrides.

    Args:
        path: config file; defaults to ``configs/default.yaml`` in the repository.
        overrides: nested dict merged on top of the file (used by tests and smoke runs).

    Returns:
        Plain nested dict.
    """
    path = Path(path) if path else repoRoot / "configs" / "default.yaml"
    with path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if overrides:
        cfg = merge(cfg, overrides)
    cfg["configPath"] = str(path)
    return cfg


def merge(base: dict, extra: dict) -> dict:
    """Return a deep copy of ``base`` with ``extra`` merged in recursively."""
    out = copy.deepcopy(base)
    for key, val in extra.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


def key(*parts) -> str:
    """Join name parts without separators, capitalising every part after the first.

    ``key("auprc", "F0")`` gives ``"auprcF0"``; used for result columns and keys.
    """
    words = [str(p) for p in parts if str(p)]
    return words[0] + "".join(w[:1].upper() + w[1:] for w in words[1:]) if words else ""


def resolve(cfg: dict, key: str) -> Path:
    """Resolve a path entry from ``cfg['paths']`` relative to the repository root."""
    p = Path(cfg["paths"][key])
    return p if p.is_absolute() else repoRoot / p


def get(cfg: dict, dotted: str, default: Any = None) -> Any:
    """Read ``cfg['a']['b']`` as ``get(cfg, 'a.b')``."""
    node: Any = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node
