# Build the ISO 532-1 reference oracle used to generate the golden test data.
#
# The ISO 532-1:2017 Annex A.4 reference sources are NOT part of this repository:
# they are distributed with the standard and are not ours to redistribute. Point
# ISO_SRC at your own copy of the ISO_532-1_LIB folder, e.g.
#
#   make ISO_SRC="$$HOME/ISO 532-1 - Program etc/Annex A.4/ISO_532-1_LIB"
#
# The generated CSVs under test/golden/ are committed, so running the MATLAB test
# suite needs neither this build nor the ISO sources. You only need them to
# regenerate the golden data.

ISO_SRC ?= $(HOME)/Documents/git/ISO 532-1 - Program etc/Annex A.4/ISO_532-1_LIB

CC     ?= cc
CFLAGS ?= -O2 -std=c99 -Wall
BIN     = iso532_oracle

.PHONY: all golden clean check

all: $(BIN)

$(BIN): iso532_oracle.c check
	$(CC) $(CFLAGS) -I"$(ISO_SRC)/include" -o $@ iso532_oracle.c "$(ISO_SRC)/src/ISO_532-1.c" -lm

check:
	@test -f "$(ISO_SRC)/src/ISO_532-1.c" || { \
		echo "ERROR: ISO reference sources not found at:"; \
		echo "  $(ISO_SRC)"; \
		echo "Set ISO_SRC to your copy of the Annex A.4 ISO_532-1_LIB folder."; \
		exit 1; }

golden: $(BIN)
	./gen_golden.sh

clean:
	rm -f $(BIN)
