# Tests of the SQAT GUI

The tests of the graphical interface in `gui/` are split in three files, from
the fastest to the slowest. They use MATLAB's own unit test framework
(`functiontests`) and need nothing besides SQAT: the test signals are
synthetic, written to a temporary folder when a file needs WAV files, and the
audio player runs with a silent buffer, so the tests make no sound.

| File | Tests | Time | What it covers |
|---|---:|---:|---|
| `SQAT_GUI_unit_test.m` | 20 | 10 s | the functions of `gui/` alone: catalogue, loading and calibration, export to Excel, spectrogram, spectral filter, weighting, enhanced STFT, the image PDF. No window opens and no SQAT metric runs. |
| `SQAT_GUI_tool_calling_test.m` | 7 | 54 s | the functions that call the SQAT metrics: the output of the GUI must equal a direct call of each metric, and the analyses and single values are read from it. No window opens. |
| `SQAT_GUI_integration_test.m` | 23 | 143 s | the interface itself, by what it does: opens, loads, runs, exports, saves, calibrates, the graphs windows and the waveform window. |

Times measured on a Mac with 12 cores and MATLAB R2026b, one run each (207 s
for the three). Run the unit tests after any change, and the three files
before a commit; together they must stay under 5 minutes.

## What is a test here

A test checks what the GUI does: a value with a source (a standard, the help
of a function, a signal of known level), a file that comes out with the right
content, or the main path of use. A test that breaks because a control moved,
a row changed height or a list got one more item is no test of the GUI, and
is left out.

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
They keep the background pool off (application data `sqat_no_background` on
`groot`), so that no enhanced map is computed behind them, and the player
runs with a silent buffer (`sqat_gui_mute`), so nothing is heard. Two tests
press Play and need an audio output on the machine; without one they are
skipped.

The integration tests read the state of the GUI, never its look. The rules,
so that a change of layout, wording, font or colour leaves the tests
standing:

- a control is found by its `Tag`, never by its text or title;
- the state of a run comes from the GUI itself through `il_gui_state` (the
  `results`, `matrix`, `stale`, `loaded`, `analyses` and `active_idx`
  variables that `SQAT_GUI` shares with its callbacks), the state of the
  player through `il_playing` (the `sqat_play` application data);
- a menu, tree or tab is read through the data behind its entries
  (`ItemsData`, `NodeData`, `UserData`), never through the entries' text;
- where a label shows a value (a calibration, a time skip, a count, a unit),
  the test looks for the value inside the text and leaves the wording alone;
- the list of controls, the theme the window opens with, the wording of
  titles and captions, fonts, colours and the sizes of the layout are not
  tested.

## MATLAB versions

The tests pass on R2026a (macOS). Two parts of the GUI depend on the version,
and the tests follow them:

- the dark theme needs `theme` (R2025a or newer); before that the GUI opens
  in the light look with the theme button disabled;
- the enhanced maps go to `backgroundPool` (R2021b or newer).
