# SQAT test suite

Automated tests for the SQAT psychoacoustic metrics, built on `matlab.unittest`.

```matlab
cd test
run_all_tests                      % everything
run_all_tests('Name', '*anchor*')  % one class
```

or headlessly:

```sh
matlab -batch "cd test; run_all_tests"
```

`run_all_tests` bootstraps `startup_SQAT` itself, so it works from a clean session.

## What is here

| file | purpose |
| --- | --- |
| `tLoudness_ISO532_1_anchor.m` | Self-contained sanity tests. Signals are synthesised or hard-coded, so these need no external data and run in seconds. |
| `tLoudness_ISO532_1_reference.m` | Conformance against the ISO 532-1:2017 Annex A.4 reference implementation, sample by sample. |
| `golden/` | Reference results produced by the ISO reference C. Committed. |
| `tools/` | The oracle used to regenerate `golden/`. |

## Why compare against the reference implementation

The obvious way to check a loudness model is against the single-value statistics the
standard tabulates (`Nmax`, `N5`). That turns out to be far too weak: it compresses a
whole time series into two numbers, and it is entirely possible for both to look right
while the loudness-versus-time curve is wrong.

That is exactly what happened here. Before the ISO 532-1 conformance fixes, the
`Nmax`/`N5` comparison reported agreement for seven of the eight synthetic time-varying
signals, while **every one of the twenty time-varying signals** in fact deviated from the
reference implementation by 3–18 % of peak loudness. Only signal 10 - the 10 ms tone
pulse, where the errors have the least room to average out - was visibly out of tolerance.

So the tests here compare the full `N(t)` series against the reference implementation
point by point. It is a stricter oracle than the standard's own tables, it localises a
regression in time instead of hiding it in an average, and it is deterministic.

The anchor tests remain valuable for a different reason: they are the human-readable
statement of what the model means (a 1 kHz tone at 40 dB SPL is 1 sone), they need no
external data, and they run fast enough to use while iterating.

## Test data

`tLoudness_ISO532_1_reference` needs the ISO test signals, which are not redistributable
here. Download the dataset from <https://doi.org/10.5281/zenodo.7933206> and place the
`validation_SQAT_v1_0` folder inside `sound_files/`. Without it those tests report as
**skipped**, not failed; the anchor tests still run.

## Regenerating the golden data

Only needed if the oracle or the signal set changes - the CSVs are committed.

```sh
cd test/tools
make ISO_SRC="/path/to/ISO 532-1 - Program etc/Annex A.4/ISO_532-1_LIB"
ISO_DIR="/path/to/ISO 532-1 - Program etc" ./gen_golden.sh
```

The ISO 532-1:2017 Annex A.4 reference sources are **not** vendored: they ship with the
standard and are not ours to redistribute. `tools/iso532_oracle.c` is a thin driver around
them - it reads a WAV, applies the same calibration as `utilities/calibrate.m`, calls the
unmodified reference library, and writes CSV. Time-varying results are emitted on the
500 Hz grid the reference program itself reports on, matching `OUT.time`.

## Adding tests for other metrics

Name new files `t<Something>.m` and they are picked up automatically. Prefer an oracle
that is independent of the implementation under test - a reference implementation, a
closed-form result, or a documented invariant - over a snapshot of current behaviour,
which only ever proves the code has not changed.
