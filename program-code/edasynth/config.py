"""This module reads the settings of the study.

Every number that can change a result is stored in a YAML file in the folder settings.
The code contains no thresholds of its own and reads them from the dictionary that this
module returns.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

repoRoot = Path(__file__).resolve().parents[2]


def loadConfig(path: str | Path | None = None, overrides: dict | None = None) -> dict:
    """Read a YAML settings file and apply optional nested changes on top of it.

    The argument path names the settings file; without it, settings/default.yaml in the
    repository is used. The argument overrides is a nested dictionary that is merged
    over the file; the tests use it. The function returns an ordinary nested dictionary.
    """
    path = Path(path) if path else repoRoot / "settings" / "default.yaml"
    with path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if overrides:
        cfg = merge(cfg, overrides)
    cfg["configPath"] = str(path)
    return cfg


def merge(base: dict, extra: dict) -> dict:
    """Return a deep copy of base into which extra has been merged at every level."""
    out = copy.deepcopy(base)
    for key, val in extra.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


def key(*parts) -> str:
    """Join parts of a name without separators and capitalise every part after the first.

    For example, key("auprc", "F0") gives "auprcF0". The names of the result columns are
    built this way.
    """
    words = [str(p) for p in parts if str(p)]
    return words[0] + "".join(w[:1].upper() + w[1:] for w in words[1:]) if words else ""


def resolve(cfg: dict, key: str) -> Path:
    """Turn a path from the section paths of the settings into a full path inside the
    repository.
    """
    p = Path(cfg["paths"][key])
    return p if p.is_absolute() else repoRoot / p


def get(cfg: dict, dotted: str, default: Any = None) -> Any:
    """Read a nested setting by a dotted name, so that get(cfg, 'a.b') returns
    cfg['a']['b'].
    """
    node: Any = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node
