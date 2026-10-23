#!/usr/bin/env bash
# Download the two public data sources at fixed versions and verify them.
#   1. PhysioNet "Non-EEG Dataset for Assessment of Neurological Status" v1.0.0 (ODC-By 1.0)
#   2. Zhang et al. (2017) repository IdeasLabUT/EDA-Artifact-Detection at commit 7c4bd46 (BSD)
# Nothing in data/raw or data/external is ever edited by hand or by code.
# Requires the package to be installed (pip install -e . or the conda environment).
set -euo pipefail
cd "$(dirname "$0")/.."

noneegDir=data/raw/noneeg-1.0.0
zhangDir=data/external/eda-artifact-detection
zhangCommit=7c4bd467bec9099669ddf3581b00a9a843a2f136

mkdir -p data/raw data/external

if [ ! -f "$noneegDir/SHA256SUMS.txt" ]; then
  echo ">> downloading PhysioNet noneeg 1.0.0 (3.7 MB)"
  tmp=$(mktemp -d)
  curl -L --fail -o "$tmp/noneeg.zip" "https://physionet.org/content/noneeg/get-zip/1.0.0/"
  unzip -q "$tmp/noneeg.zip" -d "$tmp/unz"
  src=$(dirname "$(find "$tmp/unz" -name SHA256SUMS.txt | head -n1)")
  mkdir -p "$noneegDir"
  cp -R "$src"/. "$noneegDir"/
  rm -rf "$tmp"
fi

if [ ! -d "$zhangDir/.git" ]; then
  echo ">> cloning Zhang et al. repository"
  git clone --quiet https://github.com/IdeasLabUT/EDA-Artifact-Detection "$zhangDir"
fi
git -C "$zhangDir" checkout --quiet "$zhangCommit"

# Both sources are verified by the package itself. PhysioNet's SHA256SUMS.txt
# separates hash and file name by a single space, which `shasum -c` rejects, so
# shasum is not used. A failed check stops the script before PROVENANCE.txt is written.
echo ">> verifying checksums of both sources"
if ! "${PYTHON:-python}" -m edasynth check-data; then
  echo "checksum verification failed: see the MISSING or MISMATCH lines above" >&2
  exit 1
fi

cat > data/PROVENANCE.txt <<EOF
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
noneeg: https://physionet.org/content/noneeg/1.0.0/ (zip), checksums verified against SHA256SUMS.txt
zhang:  https://github.com/IdeasLabUT/EDA-Artifact-Detection @ $zhangCommit, checksums verified
EOF
echo ">> done; provenance written to data/PROVENANCE.txt"
