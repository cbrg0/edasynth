#!/usr/bin/env bash
# Download the two public data sources at fixed versions and verify them.
#   1. PhysioNet "Non-EEG Dataset for Assessment of Neurological Status" v1.0.0 (ODC-By 1.0)
#   2. Zhang et al. (2017) repository IdeasLabUT/EDA-Artifact-Detection at commit 7c4bd46 (BSD)
# Nothing in data/raw or data/external is ever edited by hand or by code.
set -euo pipefail
cd "$(dirname "$0")/.."

NONEEG_DIR=data/raw/noneeg-1.0.0
ZHANG_DIR=data/external/eda-artifact-detection
ZHANG_COMMIT=7c4bd467bec9099669ddf3581b00a9a843a2f136

mkdir -p data/raw data/external

if [ ! -f "$NONEEG_DIR/SHA256SUMS.txt" ]; then
  echo ">> downloading PhysioNet noneeg 1.0.0 (3.7 MB)"
  tmp=$(mktemp -d)
  curl -L --fail -o "$tmp/noneeg.zip" "https://physionet.org/content/noneeg/get-zip/1.0.0/"
  unzip -q "$tmp/noneeg.zip" -d "$tmp/unz"
  src=$(dirname "$(find "$tmp/unz" -name SHA256SUMS.txt | head -n1)")
  mkdir -p "$NONEEG_DIR"
  cp -R "$src"/. "$NONEEG_DIR"/
  rm -rf "$tmp"
fi
echo ">> verifying PhysioNet checksums"
( cd "$NONEEG_DIR" && shasum -a 256 -c SHA256SUMS.txt --quiet ) && echo "   all PhysioNet files OK"

if [ ! -d "$ZHANG_DIR/.git" ]; then
  echo ">> cloning Zhang et al. repository"
  git clone --quiet https://github.com/IdeasLabUT/EDA-Artifact-Detection "$ZHANG_DIR"
fi
git -C "$ZHANG_DIR" checkout --quiet "$ZHANG_COMMIT"
echo ">> verifying Zhang et al. files"
( cd "$ZHANG_DIR" && shasum -a 256 -c ../../checksums/zhang-eda-artifact-detection-7c4bd46.sha256 --quiet ) && echo "   all Zhang files OK"

cat > data/PROVENANCE.txt <<EOF
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
noneeg: https://physionet.org/content/noneeg/1.0.0/ (zip), checksums verified against SHA256SUMS.txt
zhang:  https://github.com/IdeasLabUT/EDA-Artifact-Detection @ $ZHANG_COMMIT, checksums verified
EOF
echo ">> done; provenance written to data/PROVENANCE.txt"
