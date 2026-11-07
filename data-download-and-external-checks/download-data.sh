#!/usr/bin/env bash
# This script downloads the two public data sources at fixed versions and checks them.
# The first source is the PhysioNet "Non-EEG Dataset for Assessment of Neurological Status",
# version 1.0.0, under the Open Data Commons Attribution License 1.0. The second is the
# repository IdeasLabUT/EDA-Artifact-Detection of Zhang et al. (2017) at commit 7c4bd46,
# under the BSD licence. Nothing in downloaded-data/raw or downloaded-data/external is ever
# edited, by hand or by code. The program must be installed first, with pip install -e .
set -euo pipefail
cd "$(dirname "$0")/.."

noneegDir=downloaded-data/raw/noneeg-1.0.0
zhangDir=downloaded-data/external/eda-artifact-detection
zhangCommit=7c4bd467bec9099669ddf3581b00a9a843a2f136

mkdir -p downloaded-data/raw downloaded-data/external

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

# The program itself checks both sources. The SHA256SUMS.txt file of PhysioNet separates the
# checksum and the file name by a single space, which the command shasum -c rejects, so shasum
# is not used. If a check fails, the script stops before it writes PROVENANCE.txt.
echo ">> verifying checksums of both sources"
if ! "${PYTHON:-python}" -m edasynth check-data; then
  echo "checksum verification failed: see the MISSING or MISMATCH lines above" >&2
  exit 1
fi

cat > downloaded-data/PROVENANCE.txt <<EOF
downloaded: $(date -u +%Y-%m-%dT%H:%M:%SZ)
noneeg: https://physionet.org/content/noneeg/1.0.0/ (zip), checksums verified against SHA256SUMS.txt
zhang:  https://github.com/IdeasLabUT/EDA-Artifact-Detection @ $zhangCommit, checksums verified
EOF
echo ">> done; provenance written to downloaded-data/PROVENANCE.txt"
