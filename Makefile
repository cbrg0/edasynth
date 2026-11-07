# One command per step of the method. `make all` runs the whole study in order.
# JOBS = number of CPU cores used by the long steps (grid, variants).
PY   ?= python
JOBS ?= 4
RUN  = $(PY) -m edasynth
SRC  = --source utd

.PHONY: help install data check test test-all align realdata estimate validate \
        pilot grid variants labelnoise analyze power figures before-tag after-tag all clean-results

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

align:              ## METHOD 2.4: expert labels vs raw signal, one plot per participant
	$(RUN) align --plots

realdata:           ## real-data baseline: published features and recomputed features (G7)
	$(RUN) realdata --stored --raw

estimate:           ## METHOD 3: measure generator parameters on UTD
	$(RUN) estimate $(SRC) --jobs $(JOBS)

validate:           ## METHOD 7: validation checks G1-G7
	$(RUN) validate $(SRC) --pytest

pilot:              ## METHOD 8: pilot and number of replicates R
	$(RUN) pilot $(SRC) --jobs $(JOBS)

grid:               ## METHOD 8: main grid (resumable)
	$(RUN) grid $(SRC) --jobs $(JOBS)

variants:           ## METHOD 8: variants (V2, 15 % prevalence, majority rule, C4)
	$(RUN) grid $(SRC) --variants --jobs $(JOBS)

labelnoise:         ## METHOD 8 / goal 5: mixture cells of the grid with noisy training labels
	$(RUN) grid $(SRC) --label-noise --jobs $(JOBS)

analyze:            ## METHOD 8: map summary, compatibility (goal 5), sensitivity (goal 6)
	$(RUN) analyze $(SRC)

power:              ## METHOD 8 / goal 6: power by number of participants
	$(RUN) power $(SRC) --jobs $(JOBS)

figures:            ## F1-F6
	$(RUN) figures $(SRC)

# The preregistration tag must be pushed between the two halves:
#   make before-tag  ->  git tag prereg-synth-vN, git push --tags  ->  make after-tag
before-tag: data test align realdata estimate validate pilot
after-tag: grid variants labelnoise analyze power figures
all: before-tag after-tag


clean-results:      ## delete generated outputs, keep the README files
	find results figures -mindepth 1 -maxdepth 1 ! -name README.md -exec rm -rf {} +
