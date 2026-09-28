# Tests of the SQAT GUI

`tSQAT_GUI.m` holds 103 tests of the graphical interface in `gui/`. They use
MATLAB's own unit test framework (`functiontests`) and need nothing besides
SQAT: the test signals are synthetic WAV files written to a temporary folder
at the start, and the audio player runs with a silent buffer, so the tests
make no sound.

The other files in this folder are local scripts and stay out of git.

## Running

From the SQAT root, in a terminal:

```
matlab -batch "startup_SQAT; cd test; r = runtests('tSQAT_GUI'); disp(table(r)); assert(all([r.Passed]))"
```

or at the MATLAB prompt, with the SQAT root as the current folder:

```matlab
startup_SQAT
r = runtests('test/tSQAT_GUI.m');
table(r)
```

Windows open hidden, so the screen stays free while the suite runs. The one
exception is `test_gui_run_shows_a_progress_dialog`: the progress dialog needs
a visible window, which flashes on the screen for a few seconds.

To run part of the suite, filter by name:

```matlab
runtests('test/tSQAT_GUI.m', 'Name', 'tSQAT_GUI/test_waveform_*')       % waveform window
runtests('test/tSQAT_GUI.m', 'Name', 'tSQAT_GUI/test_enhanced_stft_*')  % enhanced STFT
runtests('test/tSQAT_GUI.m', 'Name', 'tSQAT_GUI/test_gui_*')            % main and graphs windows
runtests('test/tSQAT_GUI.m', 'Name', 'tSQAT_GUI/test_gui_run_shows_a_progress_dialog')   % one test
```

## How long it takes

The whole suite takes 4 to 6 minutes. Measured on a Mac with 12 cores and
MATLAB R2026a: 257 s in one run and 341 s in another; in low power mode it
took up to 576 s. The time of each section in the 341 s run, which had 101
tests (the channel test of the graphs windows and the tab test of the waveform
window came later and take a few seconds each):

| Section | Tests | Time (s) |
|---|---:|---:|
| Catalogue | 3 | 1 |
| Running a metric | 3 | 32 |
| Loading audio | 4 | 0.3 |
| Extracting results | 5 | 27 |
| Sharing results between metrics | 2 | 7 |
| Spectrogram, windows, filters and weighting | 11 | 6 |
| Main window | 13 | 39 |
| Running the analyses | 15 | 69 |
| Graphs windows | 20 | 71 |
| Waveform window | 17 | 44 |
| Enhanced STFT in the waveform window | 10 | 42 |

The slowest single test is `test_run_equals_the_direct_call` (about 28 s),
which calls every metric of SQAT once.

## Sections

The file is split into sections (`%%` headers) that follow the parts of the
GUI:

- **Catalogue, Running a metric, Loading audio, Extracting results, Sharing
  results between metrics**: the helper functions behind the analyses. The
  values shown by the GUI must equal a direct call of the SQAT function.
- **Spectrogram, windows, filters and weighting**: the signal processing
  helpers of the waveform window, tested as functions.
- **Main window**: controls, signal list, channels, theme, status bar.
- **Running the analyses**: runs, results tabs, parameters, export, errors.
- **Graphs windows**: plots, SQAT figures, overlays, statistics, pinned
  windows.
- **Waveform window**: playback, seeking, loop, filter box, weighting,
  spectrogram options, colour scale.
- **Enhanced STFT in the waveform window**: zoom, cache, colour floor, the
  background pool and the preview of long signals.
- **Helpers**: small functions shared by the tests. They call the callbacks
  of the controls directly, so no real mouse click is involved.

## MATLAB versions

The suite passes on R2026a (macOS). Two parts of the GUI depend on the
version, and the tests follow them:

- the dark theme needs `theme` (R2025a or newer); before that the GUI opens
  in the light look with the theme button disabled;
- the enhanced maps go to `backgroundPool` (R2021b or newer). Several tests
  set the application data `sqat_no_background` on the main window to test
  the path that computes the map at once.
