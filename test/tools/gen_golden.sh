#!/bin/sh
# Regenerate the golden reference data in test/golden/ from the ISO 532-1:2017
# Annex A.4 reference implementation.
#
#   make golden                     # uses the default paths below
#   SOUNDS=... ISO_DIR=... make golden
#
# SOUNDS  : folder holding the ISO test signals (the Zenodo validation set,
#           sound_files/validation_SQAT_v1_0/Loudness_ISO532_1 by default)
# ISO_DIR : root of the "ISO 532-1 - Program etc" folder shipped with the
#           standard, used for the Annex B.2 level vector and the Annex C
#           1 kHz 40 dB anchor.
#
# The resulting CSVs are committed, so the MATLAB test suite does not need
# either this script or the ISO sources.

set -e

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
OUT="$REPO/test/golden"

SOUNDS=${SOUNDS:-"$REPO/sound_files/validation_SQAT_v1_0/Loudness_ISO532_1"}
ISO_DIR=${ISO_DIR:-"$HOME/Documents/git/ISO 532-1 - Program etc"}
CAL="$SOUNDS/calibration signal sine 1kHz 60dB.wav"
ORACLE="$HERE/iso532_oracle"

[ -x "$ORACLE" ] || { echo "build the oracle first: make"; exit 1; }
[ -f "$CAL" ]    || { echo "calibration signal not found: $CAL"; exit 1; }

mkdir -p "$OUT"

# --- signal 1: stationary, from 28 third-octave levels (Annex B.2) ----------
LEVELS="$ISO_DIR/Annex B.2/Test signal 1.txt"
if [ -f "$LEVELS" ]; then
    "$ORACLE" lv "$LEVELS" "$OUT/sig01_stationary_levels.csv"
else
    echo "skip: $LEVELS not found"
fi

# --- signals 2-5: stationary, from audio (Annex B.3), time_skip = 0 --------
for n in 2 3 4 5; do
    f=$(ls "$SOUNDS"/"Test signal $n ("*.wav 2>/dev/null | head -1) || true
    [ -n "$f" ] || { echo "skip: signal $n not found"; continue; }
    "$ORACLE" st "$f" "$CAL" 60 0 "$OUT/$(printf 'sig%02d' $n)_stationary.csv"
done

# --- signals 6-25: time varying (Annex B.4 synthetic, B.5 technical) -------
for n in 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25; do
    f=$(ls "$SOUNDS"/"Test signal $n ("*.wav 2>/dev/null | head -1) || true
    [ -n "$f" ] || { echo "skip: signal $n not found"; continue; }
    tag=$(printf 'sig%02d' $n)
    "$ORACLE" tv "$f" "$CAL" 60 "$OUT/${tag}_timevarying.csv"
done

# --- anchor: 1 kHz at 40 dB SPL, must yield 1 sone (Annex C) ---------------
ANCHOR="$ISO_DIR/Annex C/sine 1kHz 40dB 16bit.wav"
if [ -f "$ANCHOR" ]; then
    "$ORACLE" tv "$ANCHOR" "$CAL" 60 "$OUT/anchor_1kHz_40dB_timevarying.csv"
    "$ORACLE" st "$ANCHOR" "$CAL" 60 0 "$OUT/anchor_1kHz_40dB_stationary.csv"
else
    echo "skip: $ANCHOR not found"
fi

echo "golden data written to $OUT"
