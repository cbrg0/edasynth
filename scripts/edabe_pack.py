"""Reduce the EDABE CSV files to 8 Hz (docs/EXTERNAL-CHECKS.md, preparation).

Run in the folder that holds the unpacked EDABE data set; writes edabe-8hz.npz and
edabe-8hz-info.txt. A block of 16 samples (128 Hz -> 8 Hz) is averaged; it counts as
artifact when any of its samples is labelled.
"""

import glob
import os

import numpy as np
import pandas as pd

root = "Electrodermal Activity artifact correction BEnchmark (EDABE)-EDABE dataset"
out, info = {}, []
for path in sorted(glob.glob(os.path.join(root, "*", "*.csv"))):
    split = os.path.basename(os.path.dirname(path))
    name = os.path.basename(path)[:-4]
    df = pd.read_csv(path, sep=";")
    n = len(df) // 16 * 16
    raw = df["rawdata"].to_numpy(float)[:n].reshape(-1, 16).mean(1)
    cln = df["cleandata"].to_numpy(float)[:n].reshape(-1, 16).mean(1)
    lab = df["binarytarget"].to_numpy(float)[:n].reshape(-1, 16).max(1)
    key = f"{split}__{name}"
    out[key + "__raw"] = raw.astype(np.float32)
    out[key + "__clean"] = cln.astype(np.float32)
    out[key + "__label"] = lab.astype(np.int8)
    rng = f"{np.nanmin(raw):.3f}..{np.nanmax(raw):.3f}"
    info.append(f"{key} rows={len(df)} artifactShare={lab.mean():.3f} rawRange={rng}")
    print(info[-1], flush=True)
np.savez_compressed("edabe-8hz.npz", **out)
open("edabe-8hz-info.txt", "w").write("\n".join(info))
print("done", len(info), "files")
