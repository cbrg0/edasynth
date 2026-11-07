# One command per step of the method. `make all` runs the whole study in order.
# JOBS = number of CPU cores used by the long steps (grid, variants).
PY   ?= python
JOBS ?= 4
RUN  = $(PY) -m edasynth
SRC  = --source utd

.PHONY: help install data check test test-all align realdata estimate validate \
        pilot grid variants labelnoise analyze power figures before-tag after-tag all clean-results

help:               ## shows this list of steps with a short description of each
	@grep -E '^[a-zA-Z-]+:.*## ' Makefile | sed -E 's/:.*## /\t/'

install:            ## installs the program together with the packages for tests and comparisons
	$(PY) -m pip install -e ".[dev,reference]"

data:               ## downloads both public data sources and checks every file by its checksum
	PYTHON=$(PY) bash data-download-and-external-checks/download-data.sh

check:              ## checks the checksums of the downloaded files again
	$(RUN) check-data

test:               ## runs the fast automatic tests of the program, which form check G1
	$(PY) -m pytest -m "not slow"

test-all:           ## runs all automatic tests, including those that need the downloaded data
	$(PY) -m pytest

align:              ## checks that the expert labels match the raw signal and draws one plot per participant (METHOD section 2)
	$(RUN) align --plots

realdata:           ## compares the detectors on real data, with the published and with the recomputed features (RESULTS section 1)
	$(RUN) realdata --stored --raw

estimate:           ## measures the parameters of the generator on the UTD recordings (METHOD section 3)
	$(RUN) estimate $(SRC) --jobs $(JOBS)

validate:           ## runs the checks G1 to G7 of the generator (METHOD section 7)
	$(RUN) validate $(SRC) --pytest

pilot:              ## runs the pilot that determines the number of replicates (METHOD section 8)
	$(RUN) pilot $(SRC) --jobs $(JOBS)

grid:               ## computes the main grid of 300 scenarios; an interrupted run continues where it stopped
	$(RUN) grid $(SRC) --jobs $(JOBS)

variants:           ## computes the variants of the grid (METHOD section 4.5)
	$(RUN) grid $(SRC) --variants --jobs $(JOBS)

labelnoise:         ## computes the mixture cells again with noisy training labels (METHOD section 8)
	$(RUN) grid $(SRC) --label-noise --jobs $(JOBS)

analyze:            ## summarises the grid, its agreement with the real data and the effect of each parameter
	$(RUN) analyze $(SRC)

power:              ## estimates how many participants a laboratory study needs (METHOD section 8)
	$(RUN) power $(SRC) --jobs $(JOBS)

figures:            ## draws the figures F0 to F6 from the result tables
	$(RUN) figures $(SRC)

# The preregistration tag must be pushed between the two halves:
#   make before-tag  ->  git tag prereg-synth-vN, git push --tags  ->  make after-tag
before-tag: data test align realdata estimate validate pilot ## runs every step that must be finished before the preregistration tag
after-tag: grid variants labelnoise analyze power figures ## runs every step that follows the preregistration tag; this takes one to two days
all: before-tag after-tag


clean-results:      ## deletes all generated outputs
	find result-tables result-figures -mindepth 1 -maxdepth 1 -exec rm -rf {} +
