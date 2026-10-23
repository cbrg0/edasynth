# One command per step of the method. `make all` runs the whole study in order.
# JOBS = number of CPU cores used by the long steps (grid, variants).
PY   ?= python
JOBS ?= 4
RUN  = $(PY) -m edasynth
SRC  = --source utd

.PHONY: help install data check test test-all inspect align realdata estimate validate detectors \
        pilot grid variants labelnoise analyze power figures before-tag after-tag all smoke clean-results

help:               ## list the targets with a one-line description
	@grep -E '^[a-zA-Z-]+:.*## ' Makefile | sed -E 's/:.*## /\t/'

install:            ## editable install with test and reference extras
	$(PY) -m pip install -e ".[dev,reference]"

data:               ## download and verify both data sources
	PYTHON=$(PY) bash scripts/download-data.sh

check:              ## re-verify checksums of downloaded data
	$(RUN) check-data

test:               ## G1: unit tests (fast)
	$(PY) -m pytest -m "not slow"

test-all:           ## all tests, including those that need the downloaded data
	$(PY) -m pytest

inspect:            ## METHOD 2.1: units, ranges and annotation codes of the raw files
	$(RUN) inspect

align:              ## METHOD 2.4: expert labels vs raw signal, one plot per participant
	$(RUN) align --plots

realdata:           ## real-data baseline: published features and recomputed features (G7)
	$(RUN) realdata --stored --raw

estimate:           ## METHOD 3: measure generator parameters on UTD
	$(RUN) estimate $(SRC) --jobs $(JOBS)

validate:           ## METHOD 7: validation checks G1-G7
	$(RUN) validate $(SRC) --pytest

detectors:          ## METHOD 9 / goal 2: untrained EDA-only detectors by artifact class
	$(RUN) detectors $(SRC) --reps 5

pilot:              ## METHOD 8: pilot and number of replicates R
	$(RUN) pilot $(SRC) --jobs $(JOBS)

grid:               ## METHOD 8: main grid (resumable)
	$(RUN) grid $(SRC) --jobs $(JOBS)

variants:           ## METHOD 8: variants (V2, 15 % prevalence, majority rule, C4)
	$(RUN) grid $(SRC) --variants --jobs $(JOBS)

labelnoise:         ## METHOD 8 / goal 5: main grid with noisy training labels
	$(RUN) grid $(SRC) --label-noise --jobs $(JOBS)

analyze:            ## METHOD 8: map summary, compatibility (goal 5), sensitivity (goal 6)
	$(RUN) analyze $(SRC)

power:              ## METHOD 8 / goal 6: power by number of participants
	$(RUN) power $(SRC) --jobs $(JOBS)

figures:            ## F1-F7
	$(RUN) figures $(SRC)

# The preregistration tag must be pushed between the two halves:
#   make before-tag  ->  settle docs/STATUS.md open questions, git tag prereg-synth-v1, git push --tags  ->  make after-tag
before-tag: data test inspect align realdata estimate validate detectors pilot
after-tag: grid variants labelnoise analyze power figures
all: before-tag after-tag

smoke:              ## ~50-min end-to-end run on AWW (development check, NOT results)
	$(RUN) --config configs/smoke.yaml estimate --source aww --jobs $(JOBS)
	$(RUN) --config configs/smoke.yaml validate --source aww --gates G2 G3 G4 G5 G6
	$(RUN) --config configs/smoke.yaml detectors --source aww --reps 1
	$(RUN) --config configs/smoke.yaml pilot --source aww --jobs $(JOBS)
	$(RUN) --config configs/smoke.yaml grid --source aww --jobs $(JOBS)
	$(RUN) --config configs/smoke.yaml analyze --source aww
	$(RUN) --config configs/smoke.yaml power --source aww --jobs $(JOBS)
	$(RUN) --config configs/smoke.yaml figures --source aww

clean-results:      ## delete generated outputs, keep the README files
	find results figures -mindepth 1 -maxdepth 1 ! -name README.md -exec rm -rf {} +
