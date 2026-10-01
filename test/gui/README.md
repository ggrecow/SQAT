# Tests of the SQAT GUI

The tests of the graphical interface in `gui/` are split in three files, from
the fastest to the slowest. They use MATLAB's own unit test framework
(`functiontests`) and need nothing besides SQAT: the test signals are
synthetic, written to a temporary folder when a file needs WAV files, and the
audio player runs with a silent buffer, so the tests make no sound.

| File | Tests | Time | What it covers |
|---|---:|---:|---|
| `SQAT_GUI_unit_test.m` | 29 | 13 s | the functions of `gui/` alone: catalogue, loading, share plan, windows, spectrogram, weighting, enhanced STFT. No window opens and no SQAT metric runs. |
| `SQAT_GUI_tool_calling_test.m` | 8 | 70 s | the functions that call the SQAT metrics: the output of the GUI must equal a direct call of each metric, and the analyses and single values are read from it. No window opens. |
| `SQAT_GUI_integration_test.m` | 85 | 349 s | the interface itself: main window, runs, results, graphs windows, saving, waveform window, player, enhanced STFT with the background pool. |

Times measured on a Mac with 12 cores and MATLAB R2026a, one run each (432 s
for the three). Run the unit tests after any change, and the three files
before a commit; together they must stay under 5 minutes.

The other files in this folder are local scripts and stay out of git.

## Running

From the SQAT root, in a terminal:

```
matlab -batch "startup_SQAT; cd test/gui; r = runtests({'SQAT_GUI_unit_test','SQAT_GUI_tool_calling_test','SQAT_GUI_integration_test'}); disp(table(r)); assert(all([r.Passed]))"
```

or at the MATLAB prompt, with the SQAT root as the current folder:

```matlab
startup_SQAT
r = runtests('test/gui/SQAT_GUI_unit_test.m');     % one file
table(r)
```

To run part of a file, filter by name:

```matlab
runtests('test/gui/SQAT_GUI_integration_test.m', 'Name', 'SQAT_GUI_integration_test/test_waveform_*')       % waveform window
runtests('test/gui/SQAT_GUI_integration_test.m', 'Name', 'SQAT_GUI_integration_test/test_enhanced_stft_*')  % enhanced STFT
runtests('test/gui/SQAT_GUI_integration_test.m', 'Name', 'SQAT_GUI_integration_test/test_gui_*')            % main and graphs windows
```

The windows of the integration tests open hidden, so the screen stays free.
The one exception is `test_gui_run_shows_a_progress_dialog`: the progress
dialog needs a visible window, which flashes on the screen for a few seconds.

The integration tests keep the background pool off (application data
`sqat_no_background` on `groot`), so that no enhanced map is computed behind
them; the four tests of the pool turn it on for themselves (`il_use_pool`).

## MATLAB versions

The tests pass on R2026a (macOS). Two parts of the GUI depend on the version,
and the tests follow them:

- the dark theme needs `theme` (R2025a or newer); before that the GUI opens
  in the light look with the theme button disabled;
- the enhanced maps go to `backgroundPool` (R2021b or newer).
