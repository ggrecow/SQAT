function tests = SQAT_GUI_integration_test
% Integration tests of the SQAT graphical interface: the windows open (hidden),
% the callbacks run, the metrics run, and the player and the background pool work.
% How to run the three test files and how long they take: test/gui/README.md.
tests = functiontests(localfunctions);
end

%% Fixtures ----------------------------------------------------------------

function setupOnce(tc)
setappdata(groot, 'sqat_gui_mute', true);           % the player runs, with a silent buffer
setappdata(groot, 'sqat_no_background', true);      % no enhanced maps computed behind the tests (il_use_pool)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and gui/, also on a clean MATLAB path
fs = 48000;
t = (0:1/fs:3-1/fs)';
% 1 kHz tone, 60 dB SPL, amplitude-modulated at 4 Hz (roughness and FS > 0)
x_mono = sqrt(2)*2e-5*10^(60/20) * (1 + 0.5*sin(2*pi*4*t)) .* sin(2*pi*1000*t);
% stereo: channel 2 is channel 1 attenuated by 10 dB
x_stereo = [x_mono, x_mono*10^(-10/20)];
% flyover-like noise burst for EPNL
rng(1);
x_burst = sqrt(2)*2e-5*10^(80/20) * sin(pi*t/t(end)).^2 .* randn(size(t))/3;
dir_tmp = tempname; mkdir(dir_tmp);
tc.TestData.dir_tmp = dir_tmp;
tc.TestData.fs = fs;
tc.TestData.x_mono = x_mono;
tc.TestData.x_burst = x_burst;
tc.TestData.wav_mono = fullfile(dir_tmp, 'tone_mono.wav');
tc.TestData.wav_stereo = fullfile(dir_tmp, 'tone_stereo.wav');
tc.TestData.wav_burst = fullfile(dir_tmp, 'burst.wav');
% the first second of the mono and the stereo tone: enough for the tests that
% run the ECMA-418-2 metrics twice (through the GUI and directly), which take
% most of the time of the suite on 3 s
tc.TestData.wav_mono_1s = fullfile(dir_tmp, 'tone_mono_1s.wav');
tc.TestData.wav_stereo_1s = fullfile(dir_tmp, 'tone_stereo_1s.wav');
audiowrite(tc.TestData.wav_mono, x_mono, fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_stereo, x_stereo, fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_burst, x_burst, fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_mono_1s, x_mono(1:fs), fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_stereo_1s, x_stereo(1:fs, :), fs, 'BitsPerSample', 32);
% 1 kHz, 60 dB SPL, 100 % modulated at 70 Hz: nearly constant roughness
x_rough = sqrt(2)*2e-5*10^(60/20) * (1 + sin(2*pi*70*t)) .* sin(2*pi*1000*t) / sqrt(1.5);
tc.TestData.wav_rough = fullfile(dir_tmp, 'am70.wav');
audiowrite(tc.TestData.wav_rough, x_rough, fs, 'BitsPerSample', 32);
% two steady tones, 500 Hz and 2 kHz, for the filter of the waveform window
tc.TestData.wav_two = fullfile(dir_tmp, 'two_tones.wav');
audiowrite(tc.TestData.wav_two, 0.1*sin(2*pi*500*t) + 0.1*sin(2*pi*2000*t), fs, 'BitsPerSample', 32);
% a short tone, to reach the end of the file while a test waits
tc.TestData.wav_short = fullfile(dir_tmp, 'short.wav');
audiowrite(tc.TestData.wav_short, 0.1*sin(2*pi*440*t(1:round(0.3*fs))), fs, 'BitsPerSample', 32);
% steady 1 kHz tone at 60 dB SPL, for the level of the spectrogram
tc.TestData.wav_tone = fullfile(dir_tmp, 'tone_1k_60dB.wav');
audiowrite(tc.TestData.wav_tone, sqrt(2)*2e-5*10^(60/20)*sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
end

function teardownOnce(tc)
for name = {'sqat_gui_mute', 'sqat_no_background'}
    if isappdata(groot, name{1})
        rmappdata(groot, name{1});
    end
end
rmdir(tc.TestData.dir_tmp, 's');
end

function teardown(~)
delete(findall(groot, 'Type', 'figure'));
end

%% Main window -------------------------------------------------------------

function test_gui_opens_with_the_expected_controls(tc)
% The main window opens with every control of the layout: logo, file count,
% Load files, signal list, theme, analyses, Run, the two windows, export,
% the results matrix with its plot, the table, the log, status and progress;
% the empty signal list says how to start.
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyClass(fig, 'matlab.ui.Figure');
for tag = {'logo','file_count','load_files','signals_list','theme', ...
           'add_analysis','analysis_list','run','open_graphs','export', ...
           'console','results_table','results_matrix','waveform_dock','status','progress'}
    tc.verifyNotEmpty(findobj(fig, 'Tag', tag{1}), ['missing control: ' tag{1}]);
end
tc.verifyNotEmpty(findobj(fig, 'Tag', 'signals_hint'));        % the empty list says how to start
tc.verifyEqual(findobj(fig, 'Tag', 'run').Text, 'Run Analysis');
for tag = {'show_plots','save_figures','split_figures','stop_run'}   % Stop is in the progress dialog
    tc.verifyEmpty(findobj(fig, 'Tag', tag{1}), ['control still there: ' tag{1}]);
end
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, 'No files loaded');
tc.verifyEqual(findobj(fig, 'Tag', 'status').Text, 'Ready');
tc.verifyEmpty(findobj(fig, 'Tag', 'sqat_figures'), 'the graphs window shows the SQAT figures');
dd = findobj(fig, 'Tag', 'analysis_metric_1');
m = SQAT_GUI_metrics;
tc.verifyEqual(dd.ItemsData, {m.id});
tc.verifyEqual(dd.Value, 'Loudness_ISO532_1');                 % one analysis at start
tc.verifyEmpty(findobj(fig, 'Tag', 'analysis_metric_2'));
tc.verifyEmpty(findobj(fig, 'Tag', 'channel'), 'the channel is set per signal');
tc.verifyEmpty(findobj(fig, 'Tag', 'dbfs'), 'the dBFS is set per signal');
end

function test_gui_called_without_output_keeps_working(tc)
% The usual call, SQAT_GUI at the prompt, returns nothing: the callbacks must still work.
evalc('SQAT_GUI({tc.TestData.wav_mono}, ''Visible'', ''off'')');
fig = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI');
tc.assertNumElements(fig, 1);
tc.addTeardown(@() delete(fig));
if il_has_theme()
    il_press(fig, 'theme');
    tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'light');
end
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
il_press(fig, 'open_graphs');
tc.verifyNumElements(il_window('SQAT_GUI_graphs'), 1);
end

function test_gui_windows_skip_a_user_default_CreateFcn(tc)
% A startup.m may set a default CreateFcn that fails on uifigures (e.g.
% addToolbarExplorationButtons in R2024b); the GUI windows must not run it.
old = get(groot, 'defaultFigureCreateFcn');
tc.addTeardown(@() set(groot, 'defaultFigureCreateFcn', old));
setappdata(groot, 'sqat_gui_created', {});
set(groot, 'defaultFigureCreateFcn', @(f, ~) setappdata(groot, 'sqat_gui_created', ...
    [getappdata(groot, 'sqat_gui_created'), {f}]));
tc.addTeardown(@() rmappdata(groot, 'sqat_gui_created'));
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
wins = [fig, il_window('SQAT_GUI_graphs')];
tc.assertNumElements(wins, 2);
created = getappdata(groot, 'sqat_gui_created');
for k = 1:numel(wins)
    tc.verifyFalse(any(cellfun(@(f) isequal(f, wins(k)), created)), ...
        ['the default CreateFcn ran on ' wins(k).Tag]);
end
end

function test_gui_play_and_stop_never_break_the_interface(tc)
% Audio output may be missing (batch mode); the buttons must then only log it.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_press(w, 'play');
il_press(w, 'stop');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyTrue(contains(log, 'Playing') || contains(log, 'Audio output unavailable'));
end

function test_gui_loads_files_and_draws_the_waveform(tc)
% Two files given at the start appear in the list in order, ticked, each with
% its bin; a mono file starts on channel 1 and a stereo one on All; the
% calibration shows 94 dBFS in italics until it is set, even to the same
% value.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(il_signal_names(fig), {'tone_mono.wav', 'tone_stereo.wav'});
tc.verifyTrue(findobj(fig, 'Tag', 'signal_tick_1').Value);   % new signals come ticked
tc.verifyTrue(findobj(fig, 'Tag', 'signal_tick_2').Value);
tc.verifyNotEmpty(findobj(fig, 'Tag', 'signal_remove_2'));
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, '2 files loaded');
% each signal has its channel (All by default for a stereo file) and its calibration
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Value, '1');
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_2').Value, 'All');
c = findobj(fig, 'Tag', 'signal_cal_1');
tc.verifyEqual(c.Text, '94 dBFS');
tc.verifyEqual(c.FontAngle, 'italic');                      % the default is marked
il_signal_dbfs(fig, 1, 94);
tc.verifyEqual(findobj(fig, 'Tag', 'signal_cal_1').FontAngle, 'normal');   % set, even to the same value
end

function test_gui_each_signal_offers_only_its_own_channels(tc)
% The channel menu of each signal lists only the channels of its file: 1 for
% a mono file; 1, 2 and All for a stereo one.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Items, {'1'});
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_2').Items, {'1', '2', 'All'});
end

function test_gui_signal_list_keeps_the_x_in_view_for_long_names(tc)
% The name takes the space that is left, so a long name cannot push the bin out.
long = fullfile(tc.TestData.dir_tmp, [repmat('a_very_long_signal_name_', 1, 5) '.wav']);
copyfile(tc.TestData.wav_mono, long);
fig = SQAT_GUI({long}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
g = findobj(fig, 'Tag', 'signals_list');
tc.verifyEqual(g.ColumnWidth{3}, '1x');
tc.verifyTrue(all(cellfun(@isnumeric, g.ColumnWidth([1 2 4 5 6]))));
tc.verifyLessThanOrEqual(sum([g.ColumnWidth{[1 2 4 5 6]}]), 250);   % the name keeps most of the 600 px
tc.verifyEqual(g.Parent.Parent.Parent.Parent.ColumnWidth{1}, 600);   % list, box, panel, left column, main grid
end

function test_gui_lists_remove_with_a_bin_and_show_the_parameters_in_full(tc)
% The remove buttons carry the bin icon; the short parameters of an analysis
% show in full in their tooltip.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
for tag = {'signal_remove_1', 'analysis_remove_1'}
    b = findobj(fig, 'Tag', tag{1});
    tc.verifyEmpty(b.Text, tag{1});
    tc.verifySubstring(b.Icon, 'trash.png', tag{1});
end
s = findobj(fig, 'Tag', 'analysis_summary_1');
tc.verifySubstring(s.Tooltip, 'Time skip');
pw = il_open_params(fig, 1);
c = findobj(pw, 'Tag', 'param_time_skip'); c.Value = 0.7; c.ValueChangedFcn(c, []);
tc.verifySubstring(s.Tooltip, '0.7');
end

function test_gui_channels_offer_all_only_with_a_stereo_file(tc)
% A mono file offers channel 1 only, without All, and its bin empties the
% list.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Items, {'1'});
il_remove_signal(fig, 1);
tc.verifyEmpty(il_signal_names(fig));
end

function test_gui_signal_list_ticks_and_removes_signals(tc)
% Only the ticked signals are analysed; with none ticked the run warns and
% does nothing; removing a signal takes its results with it and updates the
% file count down to "No files loaded".
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_mark_signal(fig, 2, false);                       % only the first is analysed
tc.verifyEqual(findobj(fig, 'Tag', 'run').Text, ['Run 1 signal ' char(215) ' 1 analysis']);
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_mono.wav'});
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, '1 file(s)');
il_mark_signal(fig, 2, true);
il_mark_signal(fig, 1, false);
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'});
% nothing ticked: a warning and no run
il_mark_signal(fig, 2, false);
il_press(fig, 'run');
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No signals ticked');
% removing a signal takes its results away
il_mark_signal(fig, 1, true); il_mark_signal(fig, 2, true);
il_press(fig, 'run');
il_remove_signal(fig, 1);
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'});
tc.verifyEqual(il_signal_names(fig), {'tone_1k_60dB.wav'});
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, '1 file loaded');
il_remove_signal(fig, 1);
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, 'No files loaded');
tc.verifyEmpty(findobj(fig, 'Tag', 'results_table').Data);
il_press(fig, 'run');                                % the interface goes on
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No files loaded');
end

function test_gui_signals_carry_a_number_and_the_plots_use_it_as_a_tag(tc)
% Signals are numbered #1, #2 in the list, and the legends and the statistics
% table of the graphs window name them "Signal #1, ch1". After a removal the
% remaining signal keeps its number.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(il_signal_numbers(fig), {'#1', '#2'});
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'loudness');
lg = findobj(g, 'Type', 'legend');
tc.assertNumElements(lg, 1);
tc.verifyEqual(lg.String, {'Signal #1, ch1', 'Signal #2, ch1'});
il_set(g, 'graph_analysis', 'stats');
tc.verifyEqual(findobj(g, 'Tag', 'graph_stats').ColumnName(:)', ...
    {'Quantity', 'Signal #1, ch1', 'Signal #2, ch1'});
% a number is not given twice: the signal that is left keeps its own
il_remove_signal(fig, 1);
tc.verifyEqual(il_signal_numbers(fig), {'#2'});
end

function test_gui_theme_switches_between_dark_and_light(tc)
% The window starts dark, as pySQAT, with the white logo and a sun on the
% theme button; the button switches to light, with the dark logo, a moon and
% white in place of the light grey of MATLAB, and back to dark.
% Before R2025a the window stays light and the button is disabled, with a
% tooltip naming R2025a.
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
logo = findobj(fig, 'Tag', 'logo');
b = findobj(fig, 'Tag', 'theme');
if ~il_has_theme()                                  % before R2025a: the light look, and no toggle
    tc.verifyTrue(endsWith(logo.ImageSource, 'logo.png'));
    tc.verifyTrue(isfile(logo.ImageSource));
    tc.verifyEqual(char(b.Enable), 'off');
    tc.verifySubstring(b.Tooltip, 'R2025a');
    return
end
tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'dark');         % starts dark, like pySQAT
tc.verifyTrue(endsWith(logo.ImageSource, 'logo_white.png'));
tc.verifyTrue(isfile(logo.ImageSource));
tc.verifyEqual(b.Text, char(9788));                 % a sun
tc.verifySubstring(b.Tooltip, 'Light');
il_press(fig, 'theme');
tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'light');
tc.verifyTrue(endsWith(logo.ImageSource, 'logo.png'));
tc.verifyTrue(isfile(logo.ImageSource));
tc.verifyEqual(b.Text, char(9790));                 % a moon
tc.verifySubstring(b.Tooltip, 'Dark');
dock = findobj(fig, 'Tag', 'waveform_dock');         % the light theme in white, not the grey of MATLAB
tc.verifyEqual(fig.Color, [1 1 1]);
tc.verifyEqual(dock.BackgroundColor, [1 1 1]);
tc.verifyEqual(findobj(fig, 'Tag', 'signals_list').BackgroundColor, [1 1 1]);
il_press(fig, 'theme');                              % and dark again: the colours of the theme
tc.verifyLessThan(fig.Color, 0.5);
tc.verifyLessThan(dock.BackgroundColor, 0.5);
end

function test_gui_status_bar_follows_the_run(tc)
% The progress bar starts at 0 and ends at 100 after a run, with "Done" in
% the status line.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(findobj(fig, 'Tag', 'progress').Value, 0);
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
tc.verifyEqual(findobj(fig, 'Tag', 'progress').Value, 100);
end

%% Running the analyses ----------------------------------------------------

function test_gui_run_shows_a_progress_dialog(tc)
% the dialog needs a visible window, so this one is shown
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'on');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Roughness_Daniel1997'});
% the run gives the timer a turn between steps: record what the dialog shows then
tm = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.05, 'TimerFcn', @(~, ~) record());
tc.addTeardown(@() delete(tm));
setappdata(groot, 'sqat_gui_dialogs', {});
tc.addTeardown(@() rmappdata(groot, 'sqat_gui_dialogs'));
    function record()
        if isappdata(fig, 'sqat_progress')
            d = getappdata(fig, 'sqat_progress');
            setappdata(groot, 'sqat_gui_dialogs', [getappdata(groot, 'sqat_gui_dialogs'), ...
                {struct('value', d.Value, 'message', d.Message, 'cancelable', d.Cancelable)}]);
        end
    end
start(tm);
il_press(fig, 'run');
stop(tm);
seen = getappdata(groot, 'sqat_gui_dialogs');
tc.assertNotEmpty(seen, 'no progress dialog during the run');
tc.verifyEqual(seen{1}.cancelable, matlab.lang.OnOffSwitchState.on);
tc.verifyTrue(any(cellfun(@(r) contains(r.message, 'Roughness_Daniel1997') || contains(r.message, 'Loudness'), seen)));
tc.verifyGreaterThan(max(cellfun(@(r) r.value, seen)), 0);
tc.verifyFalse(isappdata(fig, 'sqat_progress'), 'the dialog stays after the run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
end

function test_gui_results_matrix_sets_the_signals_side_by_side(tc)
% After a run the Results tab holds a matrix: one row per analysis and
% quantity, one column per signal and channel, each cell the value of the full
% table; the Waveform tab stays on screen; a removed signal takes its column
% away.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifyEqual(findobj(fig, 'Type', 'uitab', 'Title', 'Results').Parent.SelectedTab.Title, 'Waveform');   % the player stays on screen
M = findobj(fig, 'Tag', 'results_matrix');
tc.verifyEqual(M.ColumnName(:)', {'Analysis', 'Quantity', 'Signal #1, ch1', 'Signal #2, ch1'});
T = findobj(fig, 'Tag', 'results_table').Data;
row = find(strcmp(M.Data(:, 2), 'N5 (sone)'));
tc.assertNumElements(row, 1);
tc.verifyEqual(M.Data{row, 3}, T.Value(strcmp(T.File, 'tone_mono.wav') & strcmp(T.Quantity, 'N5')));
tc.verifyEqual(M.Data{row, 4}, T.Value(strcmp(T.File, 'tone_1k_60dB.wav') & strcmp(T.Quantity, 'N5')));
tc.verifyEmpty(find(strcmp(M.Data(:, 2), 'N10 (sone)'), 1));     % percentiles: only 5 and 90 %
tc.verifyNotEmpty(find(strcmp(M.Data(:, 2), 'N90 (sone)'), 1));
il_remove_signal(fig, 1);
tc.verifyEqual(M.ColumnName(:)', {'Analysis', 'Quantity', 'Signal #2, ch1'});
end

function test_gui_results_put_the_signals_and_channels_side_by_side(tc)
% In the results table each quantity is listed for every signal and channel
% before the next quantity, the analyses come in order and with their units,
% and the matrix gives each channel of a stereo signal its own column.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Roughness_Daniel1997'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
% the first quantity of #1: the mono file, then channels 1 and 2 of the stereo one
tc.verifyEqual(T.File(1:3), {'tone_mono.wav'; 'tone_stereo.wav'; 'tone_stereo.wav'});
tc.verifyEqual(T.Channel(1:3), {'1'; '1'; '2'});
tc.verifyEqual(T.Quantity(2:3), repmat(T.Quantity(1), 2, 1));
tc.verifyEqual(T.Signal(1:3), {'#1'; '#2'; '#2'});
tc.verifyEqual(unique(T.Unit(strcmp(T.Quantity, 'N5'))), {'sone'});
tc.verifyEqual(unique(T.Unit(strcmp(T.Metric, 'Roughness_Daniel1997'))), {'asper'});
tc.verifyNotEqual(T.Quantity{4}, T.Quantity{1});
% every analysis comes whole, #1 before #2
a = str2double(erase(T.Analysis, '#'));
tc.verifyTrue(issorted(a));
% the matrix gives the mono signal one column and each channel of the stereo one its own
M = findobj(fig, 'Tag', 'results_matrix');
tc.verifyEqual(M.ColumnName(:)', {'Analysis', 'Quantity', 'Signal #1, ch1', 'Signal #2, ch1', 'Signal #2, ch2'});
end

function test_gui_all_channels_runs_a_binaural_pair_in_one_call(tc)
% With All on a stereo file, ISO 532-1 loudness gives channels 1 and 2
% interleaved, and ECMA-418-2 loudness gives 1, 2 and Binaural from one call;
% the values equal those of direct calls of the metrics.
fig = SQAT_GUI({tc.TestData.wav_stereo_1s}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Loudness_ECMA418_2'});
il_signal_channel(fig, 1, 'All');
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
iso = T(strcmp(T.Metric, 'Loudness_ISO532_1'), :);
ecma = T(strcmp(T.Metric, 'Loudness_ECMA418_2'), :);
tc.verifyEqual(unique(iso.Channel, 'stable'), {'1'; '2'});
% the channels are interleaved: each quantity of channel 1, then of channel 2
tc.verifyEqual(iso.Channel(1:4), {'1'; '2'; '1'; '2'});
tc.verifyEqual(iso.Quantity(1), iso.Quantity(2));
tc.verifyEqual(iso.Quantity(3), iso.Quantity(4));
tc.verifyEqual(ecma.Channel(1:3), {'1'; '2'; 'Binaural'});
tc.verifyEqual(unique(ecma.Channel, 'stable'), {'1'; '2'; 'Binaural'});
% the values are the ones of the direct calls
ref = audioread(tc.TestData.wav_stereo_1s); fs = tc.TestData.fs;
[~, r2] = evalc('Loudness_ISO532_1(ref(:, 2), fs, 0, 2, 0.5, false)');
row = strcmp(iso.Channel, '2') & strcmp(iso.Quantity, 'N5');
tc.verifyEqual(iso.Value(row), r2.N5);
[~, rb] = evalc('Loudness_ECMA418_2(ref, fs, ''free-frontal'', 0.304, false)');
for c = {'1', '2'}
    row = strcmp(ecma.Channel, c{1}) & strcmp(ecma.Quantity, 'Nmean');
    tc.verifyEqual(ecma.Value(row), rb.Nmean(str2double(c{1})));
end
row = strcmp(ecma.Channel, 'Binaural') & strcmp(ecma.Quantity, 'Nmean');
tc.verifyEqual(ecma.Value(row), rb.Nmean(3));
% the pair was analysed in one call
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyEqual(numel(strfind(log, 'Running Loudness_ECMA418_2 on tone_stereo_1s.wav')), 1);
tc.verifySubstring(log, 'Running Loudness_ECMA418_2 on tone_stereo_1s.wav, both channels');
tc.verifyEqual(numel(strfind(log, 'Running Loudness_ISO532_1 on tone_stereo_1s.wav')), 2);
end

function test_gui_edits_the_parameters_of_a_metric(tc)
% The parameter window of an analysis shows units and limits, Reset restores
% the defaults and the summary in the list, the window stays open while
% another analysis is added, and a stationary method changes what the run
% returns.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
pw = il_open_params(fig, 1);
tc.verifyEqual(findobj(pw, 'Tag', 'unit_time_skip').Text, 's');
tc.verifyEqual(findobj(pw, 'Tag', 'param_time_skip').Limits, [0 Inf]);
c = findobj(pw, 'Tag', 'param_method');
tc.verifyEqual(c.Value, 2);
c.Value = 1; c.ValueChangedFcn(c, []);          % stationary
t = findobj(pw, 'Tag', 'param_time_skip'); t.Value = 1; t.ValueChangedFcn(t, []);
b = findobj(pw, 'Tag', 'params_reset'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(c.Value, 2);
tc.verifyEqual(t.Value, 0.5);
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'free field, time-varying, skip 0.5 s');
c.Value = 1; c.ValueChangedFcn(c, []);          % stationary again, for the run
il_press(fig, 'add_analysis');                  % the window stays open while the list changes
tc.verifyTrue(isvalid(pw));
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyTrue(any(strcmp(T.Quantity, 'Loudness')));
tc.verifyFalse(any(strcmp(T.Quantity, 'N5')));
end

function test_gui_runs_several_metrics_on_several_files(tc)
% Two metrics on two files give rows for every pair, and a value in the table
% equals the one of a direct call of the metric on the same file.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1','Roughness_Daniel1997'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_mono.wav'; 'tone_stereo.wav'});
tc.verifyEqual(unique(T.Metric), {'Loudness_ISO532_1'; 'Roughness_Daniel1997'});
% the value shown is the value of the metric
x_file = audioread(tc.TestData.wav_mono);          % what the interface reads (float32 file)
[~, ref] = evalc('Loudness_ISO532_1(x_file, tc.TestData.fs, 0, 2, 0.5, false)');
row = strcmp(T.File, 'tone_mono.wav') & strcmp(T.Metric, 'Loudness_ISO532_1') & strcmp(T.Quantity, 'N5');
tc.verifyEqual(T.Value(row), ref.N5);
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'Done');
end

function test_gui_uses_the_dbfs_and_channel_of_each_signal(tc)
% The run uses the calibration and the channel of each signal: 100 dBFS gives
% more loudness than 94 dBFS, and channel 2, recorded 10 dB lower, gives less
% than channel 1.
fig = SQAT_GUI({tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_signal_channel(fig, 1, '1');
il_press(fig, 'run');
n94 = il_value(fig, 'tone_stereo.wav', 'Loudness_ISO532_1', 'N5');
il_signal_dbfs(fig, 1, 100);
il_press(fig, 'run');
n100 = il_value(fig, 'tone_stereo.wav', 'Loudness_ISO532_1', 'N5');
tc.verifyGreaterThan(n100, n94);
il_signal_dbfs(fig, 1, 94);
il_signal_channel(fig, 1, '2');
il_press(fig, 'run');
n_ch2 = il_value(fig, 'tone_stereo.wav', 'Loudness_ISO532_1', 'N5');
tc.verifyLessThan(n_ch2, n94);                  % channel 2 is 10 dB lower
end

function test_gui_analyses_the_active_file_first(tc)
% The file shown in the interface is the first one analysed, so its result
% is the first one on screen.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_activate_signal(fig, 2);
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
first = strfind(log, 'Running Loudness_ISO532_1 on tone_1k_60dB.wav');
second = strfind(log, 'Running Loudness_ISO532_1 on tone_mono.wav');
tc.assertNotEmpty(first);
tc.assertNotEmpty(second);
tc.verifyLessThan(first(1), second(1), 'the active file was not analysed first');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'; 'tone_mono.wav'});
end

function test_gui_exports_the_results_to_excel(tc)
% The export writes the results table to a spreadsheet with the values,
% units, calibration, parameters and path of each row, and a Settings sheet
% with the signals, the analyses and the SQAT version.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
xlsx = fullfile(tc.TestData.dir_tmp, 'results.xlsx');
describe = getappdata(fig, 'sqat_run_description');
SQAT_GUI_export(findobj(fig, 'Tag', 'results_table').Data, xlsx, describe());
T = findobj(fig, 'Tag', 'results_table').Data;
R = readtable(xlsx, 'Sheet', 'Results', 'TextType', 'char');
tc.verifyEqual(R.Quantity, T.Quantity);
tc.verifyEqual(R.Value, T.Value, 'AbsTol', 1e-12);
% what is needed to reproduce a value travels with it
tc.verifyEqual(unique(R.Unit), {'asper'});
tc.verifyEqual(unique(R.Cal_dB_SPL), 94);
tc.verifyEqual(unique(R.Parameters), {'Time skip: 0 s'});
tc.verifyEqual(unique(R.Path), {tc.TestData.wav_mono});
S = readtable(xlsx, 'Sheet', 'Settings', 'TextType', 'char');
tc.verifyTrue(any(startsWith(S.Item, 'Signal #1')));
v = S.Value(strcmp(S.Item, 'Signal #1'));
tc.verifySubstring(v{1}, '(default)');
v = S.Value(strcmp(S.Item, 'Analysis #1'));
tc.verifySubstring(v{1}, 'Roughness_Daniel1997');
tc.verifyTrue(any(strcmp(S.Item, 'SQAT version')));
tc.verifyMatches(S.Value{strcmp(S.Item, 'SQAT version')}, '^v?\d');   % a release: v1.3, or 1.3 from citation.cff
end

function test_gui_exported_settings_list_what_ran(tc)
% a run stopped after the first signal, then a removed signal: the settings list only the signals with results
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
stop_run = getappdata(fig, 'sqat_stop');                % the Stop of the progress dialog
tm = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.05, 'TimerFcn', @(t, ~) il_stop_once_running(t, fig, stop_run));
tc.addTeardown(@() delete(tm));
start(tm);
il_press(fig, 'run');
stop(tm);
describe = getappdata(fig, 'sqat_run_description');
S = describe();
T = findobj(fig, 'Tag', 'results_table').Data;
tc.assertEqual(unique(T.File), {'tone_mono.wav'});      % the active signal ran, then the run stopped
tc.verifyEqual(S.Item(startsWith(S.Item, 'Signal')), {'Signal #1'});
il_press(fig, 'run');                                   % both, then the second leaves
il_remove_signal(fig, 2);
S = describe();
tc.verifyEqual(S.Item(startsWith(S.Item, 'Signal')), {'Signal #1'});
end

function test_gui_exports_a_pdf_report(tc)
% Export to .pdf writes the report of the run: the settings and the matrix of
% single values as text, then one page per analysis with its plot.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Do_SLM', 'Loudness_ISO532_1'});
il_press(fig, 'run');
pdf = fullfile(tc.TestData.dir_tmp, 'report.pdf');
write_report = getappdata(fig, 'sqat_write_report');
write_report(pdf);
tc.verifyTrue(isfile(pdf));
tc.verifyGreaterThan(dir(pdf).bytes, 10000);
write_report(pdf);                                             % a second report replaces the first
tc.verifyTrue(isfile(pdf));
end

function test_gui_session_restores_signals_and_analyses(tc)
% A saved session opened in a fresh window gives back the signals with their
% channel, calibration and tick, and the analyses with their parameters.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Do_SLM', 'Loudness_ISO532_1'});
il_signal_channel(fig, 2, '2');
il_signal_dbfs(fig, 2, 100);
il_mark_signal(fig, 1, false);
file = fullfile(tc.TestData.dir_tmp, 'session.mat');
session = getappdata(fig, 'sqat_session');
session('save', file);
fig2 = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig2));
session2 = getappdata(fig2, 'sqat_session');
session2('open', file);
tc.verifyEqual(il_signal_names(fig2), {'tone_mono.wav', 'tone_stereo.wav'});
tc.verifyFalse(findobj(fig2, 'Tag', 'signal_tick_1').Value);
tc.verifyEqual(findobj(fig2, 'Tag', 'signal_channel_2').Value, '2');
tc.verifyEqual(findobj(fig2, 'Tag', 'signal_cal_2').Text, '100 dBFS');
tc.verifyEqual(findobj(fig2, 'Tag', 'analysis_metric_1').Value, 'Do_SLM');
tc.verifyEqual(findobj(fig2, 'Tag', 'analysis_metric_2').Value, 'Loudness_ISO532_1');
end

function test_gui_marks_the_results_when_a_setting_changes(tc)
% The results tab is marked "settings changed" when a parameter or the
% calibration of a signal with results changes, and the mark goes away with
% the next run.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tab = findobj(fig, 'Type', 'uitab', 'Title', 'Results');
tc.assertNumElements(tab, 1);
pw = il_open_params(fig, 1);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
tc.verifySubstring(tab.Title, 'settings changed');
il_press(fig, 'run');
tc.verifyEqual(tab.Title, 'Results');
il_signal_dbfs(fig, 1, 100);                         % the calibration of a signal with results
tc.verifySubstring(tab.Title, 'settings changed');
end

function test_gui_keeps_the_mark_only_for_a_real_change(tc)
% a reset to the values already there is no change; removing a signal keeps the mark of an earlier change
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tab = findobj(fig, 'Type', 'uitab', 'Title', 'Results');
pw = il_open_params(fig, 1);
b = findobj(pw, 'Tag', 'params_reset'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(tab.Title, 'Results');
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
il_remove_signal(fig, 1);
tc.verifySubstring(tab.Title, 'settings changed');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Settings changed');
end

function test_gui_run_keeps_the_results_that_did_not_change(tc)
% an analysis added to the list runs alone; a change of calibration, channel or
% parameters runs the analysis again
fig = SQAT_GUI({tc.TestData.wav_stereo_1s}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
v1 = il_value(fig, 'tone_stereo_1s.wav', 'Loudness_ISO532_1', 'Nmean');
il_select_metrics(fig, {'Loudness_ISO532_1', 'Roughness_Daniel1997'});
il_press(fig, 'run');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifySubstring(log, 'Loudness_ISO532_1 on tone_stereo_1s.wav: kept from the last run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, '1 kept from the last run');
tc.verifyEqual(il_value(fig, 'tone_stereo_1s.wav', 'Loudness_ISO532_1', 'Nmean'), v1);
tc.verifyNotEmpty(il_value(fig, 'tone_stereo_1s.wav', 'Roughness_Daniel1997', 'Rmean'));
il_press(fig, 'run');                                 % nothing changed: nothing runs
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, '2 kept from the last run');
il_signal_dbfs(fig, 1, 104);                          % 10 dB more: every analysis again
il_press(fig, 'run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
tc.verifyEmpty(strfind(findobj(fig, 'Tag', 'status').Text, 'kept'));
tc.verifyGreaterThan(il_value(fig, 'tone_stereo_1s.wav', 'Loudness_ISO532_1', 'Nmean'), v1);
il_signal_channel(fig, 1, '1');                       % one channel instead of both
il_press(fig, 'run');
tc.verifyEmpty(strfind(findobj(fig, 'Tag', 'status').Text, 'kept'));
pw = il_open_params(fig, 2);                          % the parameters of the roughness
c = findobj(pw, 'Tag', 'param_time_skip'); c.Value = 0.5; c.ValueChangedFcn(c, []);
il_press(fig, 'run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, '1 kept from the last run');
end

function test_gui_run_after_a_stop_computes_what_was_left(tc)
% a stopped run keeps what ran; the next run computes the rest, not the kept part again
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
stop_run = getappdata(fig, 'sqat_stop');                % the Stop of the progress dialog
tm = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.05, 'TimerFcn', @(t, ~) il_stop_once_running(t, fig, stop_run));
tc.addTeardown(@() delete(tm));
start(tm);
il_press(fig, 'run');
stop(tm);
tc.assertSubstring(findobj(fig, 'Tag', 'status').Text, 'Stopped');
il_press(fig, 'run');
status = findobj(fig, 'Tag', 'status').Text;
tc.verifySubstring(status, 'Done');
tc.verifySubstring(status, '1 kept from the last run');  % the signal done before the stop
tc.verifyNotEmpty(il_value(fig, 'tone_mono.wav', 'Loudness_ISO532_1', 'Nmean'));
tc.verifyNotEmpty(il_value(fig, 'tone_1k_60dB.wav', 'Loudness_ISO532_1', 'Nmean'));
end

function test_gui_graphs_window_has_no_axes_toolbar(tc)
% Save... keeps the figures; the toolbar of the axes would only move the plots
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
g = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g, 1);
for value = {'sqat', 'all', 'loudness'}
    il_set(g, 'graph_analysis', value{1});
    axs = findall(g, 'Type', 'axes');
    tc.assertNotEmpty(axs, value{1});
    for ax = axs'
        tc.verifyEqual(char(ax.Toolbar.Visible), 'off', value{1});
    end
end
end

function test_gui_calibration_explains_the_full_scale(tc)
% The head of the column explains calibration and its three ways; each
% signal shows its method and, in the tooltip, its full-scale level.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
head = findobj(fig, 'Tag', 'signal_cal_info');
tc.assertNotEmpty(head);
tc.verifySubstring(head.Text, char(9432));
for part = {'between -1 and +1', 'Full-scale level (dBFS)', 'Calibrator recording', 'Relative level'}
    tc.verifySubstring(head.Tooltip, part{1});
end
tc.verifySubstring(findobj(fig, 'Tag', 'signal_cal_1').Tooltip, 'Full scale: 94.00 dB SPL');
tc.verifyFalse(contains(head.Tooltip, 'Greco'));     % no citation the reader cannot look up
end

function test_waveform_and_spectrogram_keep_aligned_plot_areas(tc)
% each plot sits in a box of its own at fixed margins, the same for both, so
% the time axes line up and a zoom moves nothing. A hidden window reports the
% size of its boxes late, so the resizing is left to a look at the screen.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axw = findobj(w, 'Tag', 'waveform_axes');
axs = findobj(w, 'Tag', 'spectrogram');
for lim = {[0 3], [1 1.5]}
    axs.XLim = lim{1}; drawnow
    pw = axw.InnerPosition;
    ps = axs.InnerPosition;
    tc.verifyEqual(pw([1 3]), ps([1 3]));
    tc.verifyEqual(ps(1), 80);
    cb = axs.Colorbar.Position;                       % the colorbar right of the plot
    tc.verifyGreaterThan(cb(1), ps(1) + ps(3));
end
end

function test_gui_calibrates_each_channel_from_a_stereo_calibrator(tc)
% A calibrator of 94 dB recorded in each ear, channel 2 recorded 6 dB lower
% (a less sensitive ear), and a tone that reached both ears with the same
% pressure, so it too is 6 dB lower in channel 2. Each channel gets its own
% full scale (94 dB minus the level of its recording), and the tone reads the
% same loudness in both; the export tells the method.
fs = tc.TestData.fs;
t = (0:2*fs-1)' / fs;
cal = fullfile(tc.TestData.dir_tmp, 'calibrator_stereo.wav');
audiowrite(cal, 0.5 * [sin(2*pi*1000*t), 0.5 * sin(2*pi*1000*t)], fs, 'BitsPerSample', 32);
sig = fullfile(tc.TestData.dir_tmp, 'tone_through_the_ears.wav');
audiowrite(sig, 0.1 * [sin(2*pi*1000*t), 0.5 * sin(2*pi*1000*t)], fs, 'BitsPerSample', 32);
fig = SQAT_GUI({sig}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
set_calibration = getappdata(fig, 'sqat_set_calibration');
tc.verifyTrue(set_calibration(1, 'calibrator', 94, cal));
tc.verifyEqual(findobj(fig, 'Tag', 'signal_cal_1').Text, 'calib. 94 dB');
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
full = @(c) T.Cal_dB_SPL(find(strcmp(T.Channel, c), 1));
tc.verifyEqual(full('1'), 94 - 20*log10(0.5/sqrt(2)), 'AbsTol', 1e-4);
tc.verifyEqual(full('2'), 94 - 20*log10(0.25/sqrt(2)), 'AbsTol', 1e-4);
N = @(c) T.Value(strcmp(T.Channel, c) & strcmp(T.Quantity, 'Nmean'));
tc.verifyEqual(N('2'), N('1'), 'RelTol', 1e-9);
describe = getappdata(fig, 'sqat_run_description');
S = describe();
tc.verifySubstring(S.Value{strcmp(S.Item, 'Signal #1')}, 'calibration calib. 94 dB');
end

function test_gui_calibrates_each_channel_from_its_own_recording(tc)
% The same ears calibrated one at a time, a mono recording per channel: the
% full scale of each channel is the same as from the stereo recording, and
% the export lists both files.
fs = tc.TestData.fs;
t = (0:2*fs-1)' / fs;
cal = fullfile(tc.TestData.dir_tmp, {'calibrator_left.wav', 'calibrator_right.wav'});
audiowrite(cal{1}, 0.5 * sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
audiowrite(cal{2}, 0.25 * sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
sig = fullfile(tc.TestData.dir_tmp, 'tone_through_the_ears.wav');
audiowrite(sig, 0.1 * [sin(2*pi*1000*t), 0.5 * sin(2*pi*1000*t)], fs, 'BitsPerSample', 32);
fig = SQAT_GUI({sig}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
set_calibration = getappdata(fig, 'sqat_set_calibration');
tc.verifyTrue(set_calibration(1, 'calibrator', 94, cal));
tc.verifyEqual(findobj(fig, 'Tag', 'signal_cal_1').Text, 'calib. 94 dB');
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
full = @(c) T.Cal_dB_SPL(find(strcmp(T.Channel, c), 1));
tc.verifyEqual([full('1') full('2')], 94 - 20*log10([0.5 0.25]/sqrt(2)), 'AbsTol', 1e-4);
describe = getappdata(fig, 'sqat_run_description');
S = describe();
tc.verifySubstring(S.Value{strcmp(S.Item, 'Signal #1')}, ['(' cal{1} '; ' cal{2} ')']);
end

function test_gui_reports_a_failing_metric_and_goes_on(tc)
% A metric that raises an error (a time skip longer than the signal) is
% reported as ERROR in the console, and the other metric of the run still
% gives its results.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1','Roughness_Daniel1997'});
pw = il_open_params(fig, 1);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
c = findobj(pw, 'Tag', 'param_time_skip'); c.Value = 10; c.ValueChangedFcn(c, []);  % longer than the signal: the metric raises an error
il_press(fig, 'run');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifySubstring(log, 'ERROR');
tc.verifySubstring(log, 'Loudness_ISO532_1');
st = findobj(fig, 'Tag', 'status');                            % the status bar names the error in red
tc.verifySubstring(st.Text, 'error');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.Metric), {'Roughness_Daniel1997'});
end

function test_gui_run_without_files_or_metrics_only_warns(tc)
% Run with no files, or with no analyses, only writes a warning to the
% console; on the way, the summary of the ECMA-418-2 loudness parameters is
% checked.
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'run');
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No files');
fig2 = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig2));
il_select_metrics(fig2, {});
% the summary of the ECMA-418-2 parameters
il_select_metrics(fig2, {'Loudness_ECMA418_2'});
tc.verifyEqual(findobj(fig2, 'Tag', 'analysis_summary_1').Text, 'free-frontal, skip 0.304 s');
il_select_metrics(fig2, {});
il_press(fig2, 'run');
tc.verifySubstring(strjoin(findobj(fig2, 'Tag', 'console').Value, newline), 'No metrics');
end

function test_gui_player_sits_in_the_waveform_tab(tc)
% The waveform, the sound level and the spectrogram of the signal on screen
% sit in the Waveform tab of the main window; no other window holds them.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
dock = findobj(fig, 'Tag', 'waveform_dock');
tc.verifyNotEmpty(findobj(dock, 'Tag', 'wave_line'));
tc.verifyNotEmpty(findobj(dock, 'Tag', 'spectrogram'));
tc.verifyEmpty(findobj(fig, 'Tag', 'open_waveform'));
tc.verifyEqual(findobj(fig, 'Type', 'uitab', 'Title', 'Waveform').Parent.SelectedTab.Title, 'Waveform');
end

function test_gui_sound_level_follows_the_weightings(tc)
% The sound level below the waveform is Do_SLM of the signal on screen, with
% the frequency weighting of the player and the time weighting above the
% plot; its indicators follow the two exceedance percentages of the spinners.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
[x, fs] = SQAT_GUI_load(tc.TestData.wav_mono, 94, 1);
L = Do_SLM(x, fs, 'Z', 'f', 94);
lbl = findobj(fig, 'Tag', 'level_indicators');
tc.verifySubstring(lbl.Text, sprintf('LZeq %.1f', Get_Leq(L, fs)));
tc.verifySubstring(lbl.Text, sprintf('LZF5 %.1f   LZF90 %.1f dB', get_exceeded_value(L, 5), get_exceeded_value(L, 90)));
tc.verifySubstring(lbl.Text, sprintf('LZE %.1f', Get_Leq(L, fs) + 10*log10(numel(L) / fs)));
il_set(fig, 'wave_weighting', 'A');
il_set(fig, 'level_time_weighting', 's');
L = Do_SLM(x, fs, 'A', 's', 94);
tc.verifySubstring(lbl.Text, sprintf('LASmax %.1f', max(L)));
tc.verifyEqual(findobj(fig, 'Tag', 'level_axes').Title.String, 'Sound pressure level (A-weighted, Slow)');
cbs = findall(findobj(fig, 'Tag', 'spectrogram').Parent, 'Type', 'colorbar');   % one colour bar, relabelled
tc.assertNumElements(cbs, 1);
tc.verifyEqual(cbs.Label.String, 'SPL (dBA)');
il_set(fig, 'level_percentile_1', 10);
il_set(fig, 'level_percentile_2', 50);
tc.verifySubstring(lbl.Text, sprintf('LAS10 %.1f   LAS50 %.1f dB', get_exceeded_value(L, 10), get_exceeded_value(L, 50)));
end

function test_gui_saves_the_three_plots_of_the_waveform_tab(tc)
% Save plots writes the spectrogram, the waveform and the sound level, one
% file each, with the name chosen and a suffix per plot.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_next_file', fullfile(tc.TestData.dir_tmp, 'plots.png'));
il_press(fig, 'save_wave_plots');
for name = {'spectrogram', 'waveform', 'level'}
    tc.verifyTrue(isfile(fullfile(tc.TestData.dir_tmp, ['plots_' name{1} '.png'])), name{1});
end
end

function test_gui_draw_filter_ends_a_zoom_of_the_toolbar(tc)
% A zoom chosen in the toolbar of the axes (the three dots) would take the
% clicks of the filter tool: Draw filter turns it off, and pan with it.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
z = zoom(fig);
z.Enable = 'on';                                   % as the zoom button of the toolbar does
bt = findobj(fig, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
tc.verifyEqual(char(zoom(fig).Enable), 'off');
pan(fig, 'on');
bt.Value = false; bt.ValueChangedFcn(bt, []);
bt.Value = true; bt.ValueChangedFcn(bt, []);
tc.verifyEqual(char(pan(fig).Enable), 'off');
end

function test_gui_overlaid_curves_stay_apart_past_seven_signals(tc)
% Past the seven colours of the axes the curves change line style, so that
% eight or more signals overlaid stay apart: the eighth takes the colour of
% the first, dashed.
wavs = arrayfun(@(k) fullfile(tc.TestData.dir_tmp, sprintf('copy_%d.wav', k)), 1:8, 'UniformOutput', false);
for k = 1:8
    copyfile(tc.TestData.wav_mono_1s, wavs{k});
end
fig = SQAT_GUI(wavs, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Do_SLM'});
il_press(fig, 'run');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'tob_level');
c = flipud(findobj(findobj(g, 'Type', 'axes'), 'Type', 'stair'));   % in the order drawn
tc.assertNumElements(c, 8);
tc.verifyEqual(c(8).Color, c(1).Color);
tc.verifyEqual(char(c(1).LineStyle), '-');
tc.verifyEqual(char(c(8).LineStyle), '--');
tc.verifyEqual(numel(unique(arrayfun(@(h) [mat2str(h.Color) char(h.LineStyle)], c, 'UniformOutput', false))), 8);
end

function test_gui_space_plays_from_the_main_window(tc)
% The space bar in the main window starts the player of the Waveform tab; a
% second press pauses it. Needs an audio output.
il_needs_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
fig.KeyPressFcn(fig, struct('Key', 'space'));
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
b = findobj(fig, 'Tag', 'play');
tc.verifyEqual(b.Text, 'Pause');
pause(0.4);
fig.KeyPressFcn(fig, struct('Key', 'space'));
tc.verifyEqual(b.Text, 'Play');
end

function test_gui_add_metric_adds_an_analysis_with_its_defaults(tc)
% Add metrics opens a list of the metrics: a tick adds one analysis, + and -
% change the count, and OK appends them in the order of the list, each with
% the defaults of the catalogue and the next number; Cancel adds nothing.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'add_metric');
d = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_metric_picker');
tc.assertNumElements(d, 1);
tick = findobj(d, 'Tag', 'picker_tick_Roughness_Daniel1997');
tick.Value = true;
tick.ValueChangedFcn(tick, []);
tc.verifyEqual(findobj(d, 'Tag', 'picker_count_Roughness_Daniel1997').Text, '1x');
il_press(d, 'picker_more_Sharpness_DIN45692');
il_press(d, 'picker_more_Sharpness_DIN45692');
il_press(d, 'picker_more_Sharpness_DIN45692');
il_press(d, 'picker_less_Sharpness_DIN45692');
tc.verifyEqual(findobj(d, 'Tag', 'picker_count_Sharpness_DIN45692').Text, '2x');
tc.verifyTrue(findobj(d, 'Tag', 'picker_tick_Sharpness_DIN45692').Value);
il_press(d, 'picker_ok');
tc.verifyFalse(isvalid(d));
ids = arrayfun(@(k) findobj(fig, 'Tag', sprintf('analysis_metric_%d', k)).Value, 2:4, 'UniformOutput', false);
tc.verifyEqual(sort(ids), sort({'Roughness_Daniel1997', 'Sharpness_DIN45692', 'Sharpness_DIN45692'}));
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_4').Text, '#4');
tc.verifyEqual(findobj(fig, 'Tag', 'run').Text, ['Run 1 signal ' char(215) ' 4 analyses']);
il_press(fig, 'add_metric');
d = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_metric_picker');
il_press(d, 'picker_more_Roughness_Daniel1997');
il_press(d, 'picker_more_Roughness_Daniel1997');
il_set(d, 'picker_all', true);                       % Tick all: every metric once, the 2x stays
counts = findobj(d, '-regexp', 'Tag', '^picker_count_');
tc.verifyTrue(all(~cellfun(@isempty, {counts.Text})));
tc.verifyEqual(findobj(d, 'Tag', 'picker_count_Roughness_Daniel1997').Text, '2x');
il_set(d, 'picker_all', false);                      % and none
tc.verifyTrue(all(cellfun(@isempty, {counts.Text})));
il_press(d, 'picker_cancel');
tc.verifyEqual(findobj(fig, 'Tag', 'run').Text, ['Run 1 signal ' char(215) ' 4 analyses']);
end

function test_gui_compares_one_metric_with_two_sets_of_parameters(tc)
% One metric added twice with different parameters runs as analyses #1 and
% #2, both offered in the graphs window; removing #1 leaves #2 with its own
% number and parameters.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'add_analysis');                       % a copy of the first: Loudness #2
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_metric_2').Value, 'Loudness_ISO532_1');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_1').Text, '#1');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_2').Text, '#2');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'free field, time-varying, skip 0.5 s');
pw = il_open_params(fig, 2);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);   % stationary
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_2').Text, 'free field, stationary, skip 0.5 s');
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.Analysis, 'stable'), {'#1'; '#2'});
tc.verifyEqual(unique(T.Metric), {'Loudness_ISO532_1'});
tc.verifyTrue(any(strcmp(T.Quantity(strcmp(T.Analysis, '#1')), 'N5')));
tc.verifyTrue(any(strcmp(T.Quantity(strcmp(T.Analysis, '#2')), 'Loudness')));
il_press(fig, 'open_graphs');                        % both are on offer in the graphs
g = il_window('SQAT_GUI_graphs');
dm = findobj(g, 'Tag', 'graph_metric');
tc.verifyEqual(dm.Items, {'#1 Loudness (ISO 532-1)', '#2 Loudness (ISO 532-1)'});
il_set(g, 'graph_metric', 'Loudness_ISO532_1#2');
tc.verifyNotEmpty(findobj(g, 'Type', 'axes'));
% removing the first leaves the second with its number: #2 keeps meaning the same analysis
b = findobj(fig, 'Tag', 'analysis_remove_1'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_1').Text, '#2');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'free field, stationary, skip 0.5 s');
tc.verifyEmpty(findobj(fig, 'Tag', 'analysis_metric_2'));
il_press(fig, 'add_analysis');                       % a number is not given twice
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_2').Text, '#3');
end

function test_gui_removing_an_analysis_takes_its_results(tc)
% the results, graphs and exported settings of a removed analysis go with it; the others stay
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'add_analysis');                       % Loudness #1 and #2
pw = il_open_params(fig, 2);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
T2 = T(strcmp(T.Analysis, '#2'), :);
il_press(fig, 'add_analysis');                       % #3 is not in the run: its removal leaves the results
b = findobj(fig, 'Tag', 'analysis_remove_3'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(findobj(fig, 'Tag', 'results_table').Data, T);
b = findobj(fig, 'Tag', 'analysis_remove_1'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(findobj(fig, 'Tag', 'results_table').Data, T2);
tab = findobj(fig, 'Type', 'uitab', 'Title', 'Results');
tc.verifyNumElements(tab, 1);                        % the rest are still the results of the run
describe = getappdata(fig, 'sqat_run_description');
S = describe();
tc.verifyFalse(any(strcmp(S.Item, 'Analysis #1')));
tc.verifyTrue(any(strcmp(S.Item, 'Analysis #2')));
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
tc.verifyEqual(findobj(g, 'Tag', 'graph_metric').Items, {'#2 Loudness (ISO 532-1)'});
end

function test_gui_save_dialog_saves_the_chosen_figures(tc)
% Save in a graphs window opens a dialog: the signals, and per metric the SQAT figure,
% the analyses (the signals overlaid, as in the graphs window) and the statistics (CSV)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
shown = findobj(g, 'Tag', 'graph_analysis').Value;
il_press(g, 'graph_save');
d = il_window('SQAT_GUI_save');
tc.assertNumElements(d, 1);
sig = findobj(d, 'Tag', 'save_signals');
tc.verifyEqual(sort({sig.CheckedNodes.Text}), {'#1 tone_mono.wav', '#2 tone_1k_60dB.wav'});
items = findobj(d, 'Tag', 'save_items');
kids = items.Children(1).Children;
tc.verifyEqual(kids(1).Text, 'SQAT figure');
tc.verifyEqual(kids(end).Text, 'Statistics');
tc.verifyEqual(arrayfun(@(n) n.NodeData.aid, items.CheckedNodes, 'UniformOutput', false), {shown});   % what the window shows
aids = arrayfun(@(n) n.NodeData.aid, kids, 'UniformOutput', false);
items.CheckedNodes = kids(ismember(aids, {'sqat', 'loudness', 'stats'}));
out_dir = fullfile(tc.TestData.dir_tmp, 'saved'); mkdir(out_dir);
h = findobj(d, 'Tag', 'save_folder'); h.Value = out_dir;
il_press(d, 'save_do');
tc.verifyEmpty(il_window('SQAT_GUI_save'));                 % the dialog closes
names = {dir(out_dir).name};
tc.verifyTrue(any(startsWith(names, 'tone_mono_s1_Loudness_ISO532_1_')));
tc.verifyTrue(any(startsWith(names, 'tone_1k_60dB_s2_Loudness_ISO532_1_')));
tc.verifyTrue(ismember('Loudness_ISO532_1_loudness_s1-s2.png', names));
tc.verifyTrue(ismember('Loudness_ISO532_1_statistics_s1-s2.csv', names));
c = readcell(fullfile(out_dir, 'Loudness_ISO532_1_statistics_s1-s2.csv'));
tc.verifyEqual(c(1, :), {'Quantity', 'Signal #1, ch1', 'Signal #2, ch1'});
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), ['saved to ' out_dir]);
end

function test_gui_run_saves_no_figure(tc)
% saving is its own step, after the run: a run draws only the figure of the signal on screen
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
before = dir(pwd);
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifyEqual({dir(pwd).name}, {before.name});
tc.verifyEmpty(findobj(fig, 'Tag', 'figures_folder'));
end

%% Graphs windows ----------------------------------------------------------

function test_gui_plots_the_series_of_the_active_file(tc)
% The graphs window opens at the end of a run and on request, offers the
% metrics that ran, and shows the figure the SQAT function draws itself, with
% the same axes and data, in the inferno colour scale of the toolbox.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1','Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifyNumElements(il_window('SQAT_GUI_graphs'), 1);   % opens at the end of the run
delete(il_window('SQAT_GUI_graphs'));
il_press(fig, 'open_graphs');                            % and again on request
g = il_window('SQAT_GUI_graphs');
tc.verifyNumElements(g, 1);
pm = findobj(g, 'Tag', 'graph_metric');
tc.verifyEqual(pm.ItemsData, {'Loudness_ISO532_1','Roughness_Daniel1997'});
il_set(g, 'graph_metric', 'Roughness_Daniel1997');
% the window holds the figure that Roughness_Daniel1997 draws itself
x_file = audioread(tc.TestData.wav_mono);
before = findall(groot, 'Type', 'figure');
[~, ref] = evalc('Roughness_Daniel1997(x_file, tc.TestData.fs, 0, true)');
after = findall(groot, 'Type', 'figure');
ref_fig = after(~ismember(after, before));
tc.verifyEqual(numel(findobj(g, 'Type', 'axes')), numel(findobj(ref_fig, 'Type', 'axes')));
ln = findobj(g, 'Type', 'line');
tc.verifyTrue(any(arrayfun(@(l) isequal(l.YData(:), ref.InstantaneousRoughness(:)), ln)), ...
    'no line of the window holds the instantaneous roughness');
% with the inferno colour scale of the toolbox
for ax = findobj(g, 'Type', 'axes')'
    tc.verifyEqual(ax.Colormap, load('cmap_inferno.txt'));
end
end

function test_gui_draws_the_figure_of_the_active_file_in_the_analysis(tc)
% The analysis draws the figure of the active file in the same call, so the
% metric runs once and the graphs window only copies it. The figure carries
% the settings of the analysis.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tc.verifyNumElements(il_all_sqat_figures(), 1, 'the analysis kept no figure to plot');
tc.verifyEmpty(il_sqat_figures(), 'the figure of the analysis must stay hidden');
il_signal_dbfs(fig, 1, 100);                            % changed after the run
il_press(fig, 'open_graphs');
il_press(fig, 'open_graphs');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyEqual(numel(strfind(log, 'Drawing the SQAT figure')), 0, ...
    'the graphs window ran the metric again');
x_file = audioread(tc.TestData.wav_mono);                % 94 dBFS, as in the analysis
[~, ref] = evalc('Loudness_ISO532_1(x_file, tc.TestData.fs, 0, 2, 0.5, false)');
ln = findobj(il_window('SQAT_GUI_graphs'), 'Type', 'line');
tc.verifyTrue(any(arrayfun(@(l) isequal(l.YData(:), ref.InstantaneousLoudness(:)), ln)));
end

function test_gui_draws_the_figure_of_another_file_on_request(tc)
% The analysis keeps the figure of the active file alone. The figure of any
% other file is drawn when the graphs window asks for it, once, with the
% settings of the analysis.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tc.verifyNumElements(il_all_sqat_figures(), 1, 'the analysis drew more than the active file');
il_mark_signal(fig, 1, false);                       % the tone is the only signal ticked
il_set(il_window('SQAT_GUI_graphs'), 'graph_analysis', 'sqat');   % the window opened by the run
il_press(fig, 'open_graphs');
il_press(fig, 'open_graphs');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyEqual(numel(strfind(log, 'Drawing the SQAT figure')), 1);
tc.verifyNumElements(il_all_sqat_figures(), 2);
end

function test_gui_graphs_fall_back_when_a_sqat_figure_fails(tc)
% Roughness_ECMA418_2 cannot draw its figure for this signal (CLim error in
% the toolbox); the values stand and the window shows the time series.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_ECMA418_2'});
il_press(fig, 'run');
tc.verifyNotEmpty(findobj(fig, 'Tag', 'results_table').Data);
il_press(fig, 'open_graphs');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeSubstring(log, 'could not be drawn', 'the toolbox now draws this figure');
ln = findobj(il_window('SQAT_GUI_graphs'), 'Type', 'line');
tc.verifyNumElements(ln, 1);
end

function test_gui_close_deletes_the_drawn_figures(tc)
% Closing the main window deletes the graphs window and every figure the
% SQAT functions drew for it.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
fig.CloseRequestFcn(fig, []);
tc.verifyEmpty(il_all_sqat_figures());
tc.verifyEmpty(il_window('SQAT_GUI_graphs'));
end

function test_gui_plot_of_a_nearly_constant_series_is_readable(tc)
% A constant result must not be zoomed into its rounding noise.
fig = SQAT_GUI({tc.TestData.wav_rough}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
il_set(il_window('SQAT_GUI_graphs'), 'graph_analysis', 'roughness');
ax = findobj(il_window('SQAT_GUI_graphs'), 'Type', 'axes');
y = findobj(ax, 'Type', 'line').YData;
tc.assumeLessThan(max(y) - min(y), 5e-3 * mean(y), ...   % the premise: nearly constant (on Linux the
    'the roughness of this signal is not constant within 0.5 % here');   % FFT rounding leaves 0.18 %)
tc.verifyGreaterThanOrEqual(diff(ax.YLim), 0.09 * mean(y));
tc.verifyTrue(ax.YLim(1) <= min(y) && ax.YLim(2) >= max(y));
end

function test_gui_graph_window_offers_the_analyses_of_the_metric(tc)
% With one signal, the graphs window offers the SQAT figure, all the analyses
% and each analysis of the metric; a stationary loudness drops the time
% series from the list.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dd = findobj(g, 'Tag', 'graph_analysis');
% one signal: the SQAT figure and all the analyses come first
tc.verifyEqual(dd.ItemsData, {'sqat','all','loudness','loudness_level','specific_loudness_time','stats'});
tc.verifyEqual(dd.Value, 'sqat');
% the analyses follow the method: a stationary loudness has no time series
c = findobj(il_open_params(fig, 1), 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
il_press(fig, 'run');
dd = findobj(il_window('SQAT_GUI_graphs'), 'Tag', 'graph_analysis');
tc.verifyEqual(dd.ItemsData, {'sqat','all','specific_loudness','stats'});
end

function test_gui_graphs_overlay_the_ticked_signals(tc)
% With two signals, the graphs window leaves out the SQAT figure and overlays
% the loudness of both signals, each curve equal to a direct call.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dd = findobj(g, 'Tag', 'graph_analysis');
% two signals: no SQAT figure (it is one per signal), the analyses that overlay
tc.verifyEqual(dd.ItemsData, {'loudness','loudness_level','specific_loudness_time','stats'});
tc.verifyEqual(dd.Value, 'loudness');
ln = findobj(g, 'Type', 'line');
tc.assertNumElements(ln, 2);
fs = tc.TestData.fs;
wavs = {tc.TestData.wav_mono, tc.TestData.wav_tone};
for k = 1:2
    x_file = audioread(wavs{k});
    [~, ref] = evalc('Loudness_ISO532_1(x_file, fs, 0, 2, 0.5, false)');
    hit = arrayfun(@(l) isequal(l.YData(:), ref.InstantaneousLoudness(:)) && ...
        isequal(l.XData(:), ref.time(:)), ln);
    tc.verifyTrue(any(hit), sprintf('signal %d is not on the plot', k));
end
lg = findobj(g, 'Type', 'legend');
tc.verifyEqual(sort(lg.String(:)), {'Signal #1, ch1'; 'Signal #2, ch1'});   % tags, not file names
% unticking a signal redraws the window that follows the list
il_mark_signal(fig, 2, false);
tc.verifyNumElements(findobj(il_window('SQAT_GUI_graphs'), 'Type', 'line'), 1);
end

function test_gui_graphs_overlay_profiles_over_the_bark_axis(tc)
% The specific roughness of two signals is overlaid against Bark, each curve
% equal to the time-averaged specific roughness of a direct call.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_rough}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'specific_roughness');
ln = findobj(g, 'Type', 'line');
tc.assertNumElements(ln, 2);
fs = tc.TestData.fs;
x_file = audioread(tc.TestData.wav_mono);
[~, ref] = evalc('Roughness_Daniel1997(x_file, fs, 0, false)');
hit = arrayfun(@(l) isequal(l.YData(:), ref.TimeAveragedSpecificRoughness(:)) && ...
    isequal(l.XData(:), ref.barkAxis(:)), ln);
tc.verifyTrue(any(hit));
ax = findobj(g, 'Type', 'axes');
tc.verifySubstring(ax.XLabel.String, 'Bark');
end

function test_gui_graphs_put_maps_side_by_side_on_one_colour_scale(tc)
% The specific roughness maps of two signals are drawn side by side on one
% colour scale, in inferno.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_rough}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'specific_roughness_time');
ax = findobj(g, 'Type', 'axes');
tc.assertNumElements(ax, 2);
tc.verifyNumElements(findobj(g, 'Type', 'surface'), 2);
tc.verifyEqual(ax(1).CLim, ax(2).CLim);
tc.verifyGreaterThan(diff(ax(1).CLim), 0);
for k = 1:2
    tc.verifyEqual(ax(k).Colormap, load('cmap_inferno.txt'));
end
end

function test_gui_graphs_tabulate_the_statistics_of_the_signals(tc)
% The statistics view is a table with one column per signal and channel, and
% its values equal the single values of the results table.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'stats');
st = findobj(g, 'Tag', 'graph_stats');
tc.assertNumElements(st, 1);
tc.verifyEqual(st.ColumnName(:)', {'Quantity', 'Signal #1, ch1', 'Signal #2, ch1'});
q = string(st.Data(:, 1));
names = {'tone_mono.wav', 'tone_1k_60dB.wav'};
for c = 2:3
    row = find(q == "N5");
    tc.verifyEqual(st.Data{row, c}, il_value(fig, names{c - 1}, 'Loudness_ISO532_1', 'N5'));
end
end

function test_gui_graphs_all_analyses_gives_a_tab_each(tc)
% "All analyses" in the graphs window gives one tab per analysis of the
% metric plus one for the statistics: four tabs for Roughness_Daniel1997.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'all');
tabs = findobj(g, 'Type', 'uitab');
tc.verifyEqual(numel(tabs), 4);                 % three analyses and the statistics
tc.verifyNumElements(findobj(g, 'Type', 'surface'), 1);
tc.verifyNumElements(findobj(g, 'Tag', 'graph_stats'), 1);
end

function test_gui_graphs_choose_the_channel_of_a_binaural_analysis(tc)
% For ECMA-418-2 loudness on All, the graphs window offers channels 1, 2,
% Binaural and All in one row, and the curve of channel 2 and of Binaural
% equals the one of a direct call.
fig = SQAT_GUI({tc.TestData.wav_stereo_1s}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ECMA418_2'});
il_signal_channel(fig, 1, 'All');
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dc = findobj(g, 'Tag', 'graph_channel_1');
tc.verifyEqual(dc.Items, {'1','2','Binaural','All'});
tc.verifyTrue(all(arrayfun(@(h) h.Layout.Row, findobj(g, 'Tag', 'graph_channels').Children) == 1));   % one row at the first opening
il_set(g, 'graph_analysis', 'loudness');
ref = audioread(tc.TestData.wav_stereo_1s); fs = tc.TestData.fs;
[~, rb] = evalc('Loudness_ECMA418_2(ref, fs, ''free-frontal'', 0.304, false)');
il_set(g, 'graph_channel_1', '2');
tc.verifyEqual(findobj(g, 'Type', 'line').YData(:), rb.loudnessTDep(:, 2));
il_set(g, 'graph_channel_1', 'Binaural');
tc.verifyEqual(findobj(g, 'Type', 'line').YData(:), rb.loudnessTDepBin(:));
end

function test_gui_pinned_graph_window_keeps_its_signals(tc)
% A pinned graphs window keeps the signals it had when changes are made to the
% list, also when its analysis changes; Open graphs then opens a second
% window, and only that one follows the list.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
il_press(fig, 'open_graphs');                    % a window that follows the list is reused
g1 = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g1, 1);
il_set(g1, 'graph_pin', true);
il_mark_signal(fig, 2, false);
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);       % pinned: as it was
il_set(g1, 'graph_analysis', 'loudness_level');             % drawn again with its own signals
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
il_set(g1, 'graph_analysis', 'loudness');
il_press(fig, 'open_graphs');                    % the pinned window is kept, another opens
g = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g, 2);
g2 = g(g ~= g1);
tc.verifyEqual(findobj(g2, 'Tag', 'graph_analysis').Value, 'sqat');
il_set(g2, 'graph_analysis', 'loudness');
tc.verifyNumElements(findobj(g2, 'Type', 'line'), 1);
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
% a change in the list reaches the window that follows it, not the pinned one
il_mark_signal(fig, 2, true);
tc.verifyNumElements(findobj(g2, 'Type', 'line'), 2);
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
end

function test_gui_pinned_graph_window_keeps_its_results_after_another_run(tc)
% A run replaces the results; a pinned window keeps the ones it was pinned with.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g1 = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g1, 1);
il_set(g1, 'graph_pin', true);
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');                    % the pinned window stays, another opens
g = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g, 2);
g2 = g(g ~= g1);
tc.verifyEqual(findobj(g2, 'Tag', 'graph_metric').ItemsData, {'Roughness_Daniel1997'});
% the pinned window is drawn again from its own results, not from the last run
il_set(g1, 'graph_analysis', 'loudness_level');
il_set(g1, 'graph_analysis', 'loudness');
dm = findobj(g1, 'Tag', 'graph_metric');
tc.verifyEqual(dm.ItemsData, {'Loudness_ISO532_1'});
tc.verifyEqual(dm.Value, 'Loudness_ISO532_1');
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
tc.verifyEqual(findobj(g2, 'Tag', 'graph_metric').Value, 'Roughness_Daniel1997');
end

function test_gui_graphs_compare_a_mono_and_a_stereo_signal(tc)
% each signal has its own channel choice: the mono against channel 1, channel 2,
% the binaural result or every channel of the stereo
fig = SQAT_GUI({tc.TestData.wav_mono_1s, tc.TestData.wav_stereo_1s}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ECMA418_2'});
il_signal_channel(fig, 2, 'All');
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
d1 = findobj(g, 'Tag', 'graph_channel_1');
d2 = findobj(g, 'Tag', 'graph_channel_2');
tc.verifyEqual(d1.Items, {'1'});
tc.verifyEqual(char(d1.Enable), 'off');                    % nothing to choose for the mono
tc.verifyEqual(d2.Items, {'1', '2', 'Binaural', 'All'});
tc.verifyEqual(d2.Value, '1');
il_set(g, 'graph_analysis', 'loudness');
legend_of = @() findobj(g, 'Type', 'legend').String;
tc.verifyEqual(legend_of(), {'Signal #1, ch1', 'Signal #2, ch1'});
il_set(g, 'graph_channel_2', '2');
tc.verifyEqual(legend_of(), {'Signal #1, ch1', 'Signal #2, ch2'});
tc.verifyFalse(contains(findobj(g, 'Type', 'axes').Title.String, 'channel'));   % the legend says it
x2 = audioread(tc.TestData.wav_stereo_1s); fs = tc.TestData.fs;
[~, r2] = evalc('Loudness_ECMA418_2(x2, fs, ''free-frontal'', 0.304, false)');
lines = findobj(g, 'Type', 'line');
tc.verifyEqual(lines(strcmp({lines.DisplayName}, 'Signal #2, ch2')).YData(:), r2.loudnessTDep(:, 2));
il_set(g, 'graph_channel_2', 'Binaural');
tc.verifyEqual(legend_of(), {'Signal #1, ch1', 'Signal #2, binaural'});
il_set(g, 'graph_channel_2', 'All');
tc.verifyEqual(legend_of(), {'Signal #1, ch1', 'Signal #2, ch1', 'Signal #2, ch2', 'Signal #2, binaural'});
il_set(g, 'graph_analysis', 'stats');                      % the choice stays through a redraw
tc.verifyEqual(findobj(g, 'Tag', 'graph_channel_2').Value, 'All');
end

function test_gui_graphs_offer_no_channel_row_for_mono_signals(tc)
% With mono signals only, the channel row of the graphs window is hidden and
% the legend names each signal with ch1.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
row = findobj(g, 'Tag', 'graph_channels');
tc.verifyEqual(row.Parent.RowHeight{2}, 0);
il_set(g, 'graph_analysis', 'loudness');
tc.verifyEqual(findobj(g, 'Type', 'legend').String, {'Signal #1, ch1', 'Signal #2, ch1'});
end

function test_gui_graphs_window_saves_what_it_shows(tc)
% Save in a graphs window opens the dialog set to its signals, metric and analysis
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'loudness');
il_press(g, 'graph_save');
d = il_window('SQAT_GUI_save');
tc.assertNumElements(d, 1);
items = findobj(d, 'Tag', 'save_items');
checked = items.CheckedNodes;
tc.verifyEqual(arrayfun(@(n) n.NodeData.aid, checked, 'UniformOutput', false), {'loudness'});
tc.verifyEqual(checked.NodeData.key, 'Loudness_ISO532_1');
tc.verifyNumElements(findobj(d, 'Tag', 'save_signals').CheckedNodes, 2);
out_dir = fullfile(tc.TestData.dir_tmp, 'saved_graphs'); mkdir(out_dir);
h = findobj(d, 'Tag', 'save_folder'); h.Value = out_dir;
h = findobj(d, 'Tag', 'save_format'); h.Value = 'pdf';
il_press(d, 'save_do');
tc.verifyEqual({dir(fullfile(out_dir, '*.pdf')).name}, {'Loudness_ISO532_1_loudness_s1-s2.pdf'});
end

function test_gui_graphs_window_needs_results(tc)
% Open graphs before any run opens no window and asks in the console to run
% an analysis first.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_graphs');
tc.verifyEmpty(il_window('SQAT_GUI_graphs'));
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'Run an analysis first');
end

%% Waveform window ---------------------------------------------------------

function test_colormap_is_the_inferno_scale_of_the_toolbox(tc)
% The spectrogram of the waveform window uses the 256 colours of the inferno
% scale of the toolbox (cmap_inferno.txt).
c = load('cmap_inferno.txt');
tc.verifySize(c, [256 3]);
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
ax = findobj(fig, 'Tag', 'spectrogram');
tc.verifyEqual(ax.Colormap, c);
end

function test_waveform_window_shows_the_signal_and_follows_playback(tc)
% The waveform window draws the calibrated signal with the playhead at 0; Play
% turns the button into Pause and the playhead moves with the playback. Needs
% an audio output.
il_needs_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
tc.verifyNumElements(w, 1);
ln = findobj(w, 'Tag', 'wave_line');
tc.verifyNumElements(ln, 1);
tc.verifyEqual(ln.YData(:), audioread(tc.TestData.wav_mono));   % calibrated at 94 dBFS
ph = findobj(w, 'Tag', 'playhead');
tc.verifyEqual(ph.Value, 0);
b = findobj(w, 'Tag', 'play');
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyEqual(b.Text, 'Pause');
t_wait = tic;                                      % the audio device may start late
while ph.Value == 0 && toc(t_wait) < 3
    pause(0.05);
end
tc.verifyGreaterThan(ph.Value, 0);
pause(0.2);                                        % a second press within 0.15 s is taken for a repeated key
il_press(w, 'play');                               % pause
tc.verifyEqual(b.Text, 'Play');
il_press(w, 'stop');
tc.verifyEqual(ph.Value, 0);
end

function test_waveform_window_has_the_spectrogram(tc)
% The waveform window has a spectrogram in dB SPL on a log frequency axis, in
% inferno: a 60 dB tone at 1 kHz peaks near 60 dB at 1 kHz, the playhead is
% drawn on it, and the A weighting keeps the colour limits.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
ax = findobj(w, 'Tag', 'spectrogram');
tc.assertNumElements(ax, 1);
tc.verifyEqual(ax.Colormap, load('cmap_inferno.txt'));
tc.verifyEqual(ax.YScale, 'log');
sf = findobj(ax, 'Type', 'surface');
tc.assertNumElements(sf, 1);
% level in dB SPL: the 60 dB tone peaks at 60 dB, less the scalloping of the window
tc.verifyEqual(max(sf.CData(:)), 60, 'AbsTol', 1.5);
[~, i_max] = max(max(sf.CData, [], 2));
tc.verifyEqual(sf.YData(i_max), 1000, 'AbsTol', 48000/1024);
tc.verifyEqual(ax.Colorbar.Label.String, 'SPL (dB SPL)');
tc.verifyNotEmpty(findobj(w, 'Tag', 'playhead_spectrogram'));
cl = ax.CLim;
il_set(w, 'wave_weighting', 'A');
tc.verifyEqual(ax.CLim, cl);                                   % A weighting keeps the colour scale
end

function test_waveform_has_a_tab_per_signal(tc)
% the tabs choose the signal and channel of the waveform window (a stereo signal
% has one tab per channel), and follow the list of signals
fig = SQAT_GUI({tc.TestData.wav_tone, tc.TestData.wav_stereo, tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
tg = findobj(w, 'Tag', 'wave_tabs');
tc.assertNumElements(tg, 1);
titles = @() arrayfun(@(t) t.Title, tg.Children(:)', 'UniformOutput', false);
tc.verifyEqual(titles(), {'#1 tone_1k_60dB.wav', '#2 tone_stereo.wav ch1', '#2 tone_stereo.wav ch2', ...
    '#3 tone_mono.wav'});
tc.verifyEqual(tg.SelectedTab, tg.Children(1));
il_pick_tab(tg, 3);                                        % channel 2 of the stereo
tc.verifySubstring(findobj(w, 'Tag', 'waveform_axes').Title.String, 'tone_stereo.wav, channel 2');
ref = audioread(tc.TestData.wav_stereo);
tc.verifyEqual(findobj(w, 'Tag', 'waveform_axes').Children(end).YData(:), ref(:, 2));
tc.verifyEqual(char(findobj(fig, 'Tag', 'signal_name_2').FontWeight), 'bold');   % the list follows
il_pick_tab(tg, 2);                                        % channel 1 of the same signal
tc.verifySubstring(findobj(w, 'Tag', 'waveform_axes').Title.String, 'tone_stereo.wav, channel 1');
tc.verifyEqual(findobj(w, 'Tag', 'waveform_axes').Children(end).YData(:), ref(:, 1));
il_activate_signal(fig, 3);                                % and the tabs follow the list
tc.verifyEqual(tg.SelectedTab.Title, '#3 tone_mono.wav');
tc.verifySubstring(findobj(w, 'Tag', 'waveform_axes').Title.String, 'tone_mono');
il_activate_signal(fig, 2);                                % back to the stereo: the channel last chosen
tc.verifyEqual(tg.SelectedTab.Title, '#2 tone_stereo.wav ch1');
il_remove_signal(fig, 1);
tc.verifyEqual(titles(), {'#2 tone_stereo.wav ch1', '#2 tone_stereo.wav ch2', '#3 tone_mono.wav'});
tc.verifyEqual(tg.SelectedTab.Title, '#2 tone_stereo.wav ch1');
end

function test_waveform_space_starts_and_pauses_playback(tc)
% The space bar starts, pauses and restarts the playback, other keys do
% nothing, and two space events at once count as one press. Needs an audio
% output.
il_needs_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
b = findobj(w, 'Tag', 'play');
w.KeyPressFcn(w, struct('Key', 'a'));                      % another key does nothing
tc.verifyEqual(b.Text, 'Play');
w.KeyPressFcn(w, struct('Key', 'space'));                  % the first play, with no button pressed
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyEqual(b.Text, 'Pause');
pause(0.4);
w.KeyPressFcn(w, struct('Key', 'space'));
tc.verifyEqual(b.Text, 'Play');
pause(0.4);
w.KeyPressFcn(w, struct('Key', 'space'));                  % and plays again
tc.verifyEqual(b.Text, 'Pause');
% the same key event arriving twice at once is one press (a focused button also answers to space)
pause(0.4);
w.KeyPressFcn(w, struct('Key', 'space'));
w.KeyPressFcn(w, struct('Key', 'space'));
tc.verifyEqual(b.Text, 'Play');
il_press(w, 'stop');
end

function test_waveform_tests_play_in_silence(tc)
% the tests of the waveform window start the player, and must not be heard
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
audio = getappdata(w, 'sqat_audio');
info = getappdata(w, 'sqat_play');
tc.assertGreaterThan(max(abs(audio())), 0.05, 'the file itself is silent');
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyTrue(info().playing);
tc.verifyEqual(double(info().buffer_peak), 0, 'the buffer that goes to the player is not silent');
il_press(w, 'stop');
end

function test_waveform_click_moves_the_playhead(tc)
% A click on the waveform or on the spectrogram moves both playheads to that
% time, kept between the start and the end of the file; the lines and the map
% do not take the click.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
ph = findobj(w, 'Tag', 'playhead'); ph2 = findobj(w, 'Tag', 'playhead_spectrogram');
ax = findobj(w, 'Tag', 'waveform_axes'); axs = findobj(w, 'Tag', 'spectrogram');
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [1.5 0 0]));
tc.verifyEqual(ph.Value, 1.5, 'AbsTol', 1/48000);
tc.verifyEqual(ph2.Value, 1.5, 'AbsTol', 1/48000);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [2.0 1000 0]));   % the spectrogram too
tc.verifyEqual(ph.Value, 2.0, 'AbsTol', 1/48000);
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [99 0 0]));         % beyond the end: the end
tc.verifyEqual(ph.Value, 3.0, 'AbsTol', 2/48000);
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [-4 0 0]));         % before the start: the start
tc.verifyEqual(ph.Value, 0, 'AbsTol', 1/48000);
% what is drawn on the axes does not take the click
tc.verifyTrue(all(strcmp({findobj(ax, 'Type', 'line').PickableParts}, 'none')));
tc.verifyEqual(findobj(axs, 'Type', 'surface').PickableParts, 'none');
end

function test_waveform_play_starts_where_the_click_was_and_stop_goes_back(tc)
% the state of the play is checked, not the time: the audio device starts a second or two late
il_needs_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
ph = findobj(w, 'Tag', 'playhead');
ax = findobj(w, 'Tag', 'waveform_axes');
info = getappdata(w, 'sqat_play');
fs = 48000;
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [1.5 0 0]));
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyTrue(info().playing);
tc.verifyEqual(info().sample, round(1.5*fs) + 1);
il_wait_until(@() ph.Value ~= 1.5);                        % the playhead moves once the device runs
tc.verifyNotEqual(ph.Value, 1.5);
% a click while playing jumps there and goes on
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [0.2 0 0]));
tc.verifyEqual(findobj(w, 'Tag', 'play').Text, 'Pause');
tc.verifyTrue(info().playing);
tc.verifyEqual(info().sample, round(0.2*fs) + 1);
il_press(w, 'stop');
tc.verifyEqual(ph.Value, 0);
tc.verifyFalse(info().playing);
il_press(w, 'play');                                       % after Stop the file starts over
tc.verifyTrue(info().playing);
tc.verifyEqual(info().sample, 1);
il_press(w, 'stop');
end

function test_waveform_loops_by_default(tc)
% The loop is on by default: a short file keeps playing past its end until
% Stop. With the loop off the playback ends at the end of the file and the
% playhead goes back to 0. Needs an audio output that plays in real time.
il_needs_realtime_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_short}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
lp = findobj(w, 'Tag', 'loop');
tc.verifyTrue(lp.Value, 'the loop must be on by default');
b = findobj(w, 'Tag', 'play');
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
pause(1.2);                                                % four times the length of the file
tc.verifyEqual(b.Text, 'Pause', 'the playback ended instead of looping');
il_press(w, 'stop');
tc.verifyEqual(b.Text, 'Play');
pause(0.5);
tc.verifyEqual(b.Text, 'Play', 'the loop restarted after Stop');
% without the loop the playback ends at the end of the file
lp.Value = false;                                          % read when the file ends
il_press(w, 'play');
pause(1.2);
tc.verifyEqual(b.Text, 'Play');
tc.verifyEqual(findobj(w, 'Tag', 'playhead').Value, 0);
end

function test_waveform_box_modes(tc)
% Draw filter with two clicks on the spectrogram draws a box from opposite
% corners, in any order, and disarms the tool; the box modes are isolate (the
% default), remove and loop only (the sound unchanged).
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axs = findobj(w, 'Tag', 'spectrogram');
ph = findobj(w, 'Tag', 'playhead');
audio = getappdata(w, 'sqat_audio');
[ref, fs] = audioread(tc.TestData.wav_two);
tc.verifyEqual(audio(), ref);                              % no box: the file as it is
tc.verifyEqual(findobj(w, 'Tag', 'draw_box').Text, 'Draw filter');
tc.verifyEqual(findobj(w, 'Tag', 'clear_boxes').Text, 'Clear filters');
% a click on the spectrogram moves the playhead until the tool is armed
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [2.5 1500 0]));   % first corner
tc.verifyEqual(ph.Value, 0);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [0.5 2500 0]));   % the opposite one, in any order
box = findobj(w, 'Tag', 'box');
tc.assertNumElements(box, 1);
tc.verifyEqual([min(box.XData) max(box.XData)], [0.5 2.5], 'AbsTol', 1e-9);
tc.verifyEqual([min(box.YData) max(box.YData)], [1500 2500], 'AbsTol', 1e-9);
tc.verifyFalse(bt.Value, 'the box tool stays armed');
% by default only what is inside the box plays
md = findobj(w, 'Tag', 'box_mode');
tc.verifyEqual(md.ItemsData, {'loop', 'isolate', 'remove'});
tc.verifyEqual(md.Items, {'Filter: loop only', 'Filter: isolate', 'Filter: remove'});
tc.verifyEqual(md.Value, 'isolate');
il_set(w, 'box_mode', 'loop');                            % loop only: the sound is not changed
tc.verifyEqual(audio(), ref);
mid = round(1.2*fs):round(1.8*fs); early = 1:round(0.3*fs);
% remove: what is inside goes
il_set(w, 'box_mode', 'remove');
y = audio();
tc.verifyLessThan(il_tone_db(y(mid), fs, 2000) - il_tone_db(ref(mid), fs, 2000), -40);
tc.verifyEqual(il_tone_db(y(mid), fs, 500), il_tone_db(ref(mid), fs, 500), 'AbsTol', 0.1);
tc.verifyEqual(il_tone_db(y(early), fs, 2000), il_tone_db(ref(early), fs, 2000), 'AbsTol', 0.1);
% isolate: only what is inside stays
il_set(w, 'box_mode', 'isolate');
y = audio();
tc.verifyEqual(il_tone_db(y(mid), fs, 2000), il_tone_db(ref(mid), fs, 2000), 'AbsTol', 0.1);
tc.verifyLessThan(il_tone_db(y(mid), fs, 500) - il_tone_db(ref(mid), fs, 500), -40);
tc.verifyLessThan(rms(y(early)), 1e-6);                    % outside the box's time
il_set(w, 'box_mode', 'loop');
tc.verifyEqual(audio(), ref);
% clearing the boxes clears the filter too
il_set(w, 'box_mode', 'remove');
tc.verifyNotEqual(audio(), ref);
il_press(w, 'clear_boxes');
tc.verifyEmpty(findobj(w, 'Tag', 'box'));
tc.verifyEqual(audio(), ref);
% a click on the same axes moves the playhead again
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.0 1000 0]));
tc.verifyEqual(ph.Value, 1.0, 'AbsTol', 1/fs);
end

function test_waveform_box_can_be_dragged(tc)
% A box can also be dragged: a preview follows the mouse, the release draws the
% box and stops following the mouse; a press and release on one spot counts
% as the first of two clicks.
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axs = findobj(w, 'Tag', 'spectrogram');
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [0.5 1500 0]));   % press
tc.assertNotEmpty(w.WindowButtonMotionFcn, 'nothing follows the mouse');
w.WindowButtonMotionFcn(w, struct('IntersectionPoint', [1.5 2000 0]));
prev = findobj(w, 'Tag', 'box_preview');
tc.assertNumElements(prev, 1);                             % the box follows the mouse
tc.verifyEqual([min(prev.XData) max(prev.XData)], [0.5 1.5], 'AbsTol', 1e-9);
w.WindowButtonMotionFcn(w, struct('IntersectionPoint', [2.5 2500 0]));
w.WindowButtonUpFcn(w, struct('IntersectionPoint', [2.5 2500 0]));   % release far away: the box
box = findobj(w, 'Tag', 'box');
tc.assertNumElements(box, 1);
tc.verifyEqual([min(box.XData) max(box.XData)], [0.5 2.5], 'AbsTol', 1e-9);
tc.verifyEqual([min(box.YData) max(box.YData)], [1500 2500], 'AbsTol', 1e-9);
tc.verifyEmpty(findobj(w, 'Tag', 'box_preview'));
tc.verifyEmpty(w.WindowButtonMotionFcn, 'the mouse is still followed');
tc.verifyFalse(bt.Value);
% a press and a release on the same spot is the first click of two
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.0 1000 0]));
w.WindowButtonUpFcn(w, struct('IntersectionPoint', [1.0 1000 0]));
tc.verifyNumElements(findobj(w, 'Tag', 'box'), 1);         % still the first box only
tc.verifyNumElements(findobj(w, 'Tag', 'box_corner'), 1);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [2.0 3000 0]));
tc.verifyNumElements(findobj(w, 'Tag', 'box'), 2);
end

function test_waveform_loops_inside_the_box(tc)
% With a box drawn and the loop on, Play starts at the box and the playhead
% stays inside it for more than three times its length. Needs an audio
% output.
il_needs_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axs = findobj(w, 'Tag', 'spectrogram');
ph = findobj(w, 'Tag', 'playhead');
b = findobj(w, 'Tag', 'play');
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.0 1500 0]));
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.4 2500 0]));
il_press(w, 'play');                                       % starts at the box, not at the start
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
il_wait_until(@() ph.Value > 0);
tc.verifyGreaterThanOrEqual(ph.Value, 1.0 - 0.02);
% it stays in the box for more than three times its length
inside = true; t0 = tic;
while toc(t0) < 1.5
    pause(0.05);
    inside = inside && ph.Value >= 1.0 - 0.02 && ph.Value <= 1.4 + 0.1;
end
tc.verifyTrue(inside, 'the playhead left the box');
tc.verifyEqual(b.Text, 'Pause');
il_press(w, 'stop');
% without the loop the box is played once
lp = findobj(w, 'Tag', 'loop');
lp.Value = false;
il_press(w, 'play');
il_wait_until(@() strcmp(b.Text, 'Play') || ph.Value > 1.0);
il_wait_until(@() strcmp(b.Text, 'Play'));
tc.verifyEqual(b.Text, 'Play', 'the box was played more than once');
end

function test_waveform_box_loop_is_left_and_entered_by_clicks(tc)
% During the playback, a box drawn with the loop on takes the playback into
% it; a click after it leaves the loop and plays on from there, a click inside
% enters it again, a click before it plays on into it and stays, the loop off
% plays the box once, and a box drawn after a pause takes the next play to its
% start. Needs an audio output.
il_needs_audio(tc);
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axs = findobj(w, 'Tag', 'spectrogram'); axw = findobj(w, 'Tag', 'waveform_axes');
ph = findobj(w, 'Tag', 'playhead');
b = findobj(w, 'Tag', 'play');
info = getappdata(w, 'sqat_play');                         % where the play comes from, whatever the audio latency
fs = 48000;
ax_click = @(ax, t) ax.ButtonDownFcn(ax, struct('IntersectionPoint', [t 1000 0]));
draw_box = @(t1, t2) draw_two_corners(w, axs, t1, t2);
% the file plays from the start, outside where the box is going to be
ax_click(axw, 0.05);
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyFalse(info().box_loop);
% a box drawn now, with the loop on, takes the playback into it
draw_box(1.0, 1.4);
tc.verifyTrue(info().playing);
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(1.0*fs) + 1);
il_wait_until(@() ph.Value > 0.5);
tc.verifyTrue(il_stays_between(ph, 0.99, 1.41, 1.5), 'the playhead left the box');
% a click outside the box leaves the loop: the file plays on from there
ax_click(axw, 2.0);
tc.verifyEqual(b.Text, 'Pause');
tc.verifyFalse(info().box_loop);
tc.verifyEqual(info().sample, round(2.0*fs) + 1);
% a click inside the box goes back into the loop
ax_click(axw, 1.2);
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(1.2*fs) + 1);
il_wait_until(@() ph.Value > 1.19 && ph.Value < 1.5);
tc.verifyTrue(il_stays_between(ph, 0.99, 1.41, 1.5), 'the playhead left the box');
% a click before the box plays on into it, and the loop holds it there
ax_click(axw, 0.6);
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(0.6*fs) + 1);
il_wait_until(@() ph.Value > 1.05 && ph.Value < 1.5);
tc.verifyTrue(il_stays_between(ph, 0.99, 1.41, 1.5), 'the play ran past the box');
% the loop tick off: the play ends where the box ends
lp = findobj(w, 'Tag', 'loop');
lp.Value = false; lp.ValueChangedFcn(lp, []);
il_wait_until(@() strcmp(b.Text, 'Play'));
tc.verifyEqual(b.Text, 'Play', 'the box was not played once');
lp.Value = true; lp.ValueChangedFcn(lp, []);
il_press(w, 'stop');
% a box drawn after a pause takes the next play to its start
il_press(w, 'clear_boxes');
ax_click(axw, 2.2);
il_press(w, 'play'); pause(0.3); il_press(w, 'play');     % paused past the box
tc.verifyEqual(b.Text, 'Play');
draw_box(0.5, 0.9);
pause(0.3);
il_press(w, 'play');
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(0.5*fs) + 1);
il_press(w, 'stop');
end

function test_waveform_weighting_filters_the_playback_and_the_spectrogram(tc)
% The Z, A and C weightings filter what plays and the spectrogram: the level
% at 500 Hz drops by the A curve, the colorbar names the weighting, and Z
% gives the file back as it is.
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
wd = findobj(w, 'Tag', 'wave_weighting');
tc.verifyEqual(wd.Items, {'Z', 'A', 'C'});
tc.verifyEqual(wd.Value, 'Z');
audio = getappdata(w, 'sqat_audio');
[ref, fs] = audioread(tc.TestData.wav_two);
sf = findobj(w, 'Tag', 'spectrogram'); sf = findobj(sf, 'Type', 'surface');
z_level = @() max(sf.CData(abs(sf.YData - 500) < 30, :), [], 'all');
level_z = z_level();
il_set(w, 'wave_weighting', 'A');
tc.verifyEqual(audio(), SQAT_GUI_weight(ref, fs, 'A'), 'AbsTol', 1e-12);
sf = findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
tc.verifyEqual(max(sf.CData(abs(sf.YData - 500) < 30, :), [], 'all') - level_z, ...
    SQAT_GUI_weight_curve(500, fs, 'A'), 'AbsTol', 0.6);      % the bin nearest 500 Hz is within 1/2 bin
tc.verifyEqual(findobj(w, 'Tag', 'spectrogram').Colorbar.Label.String, 'SPL (dBA)');
il_set(w, 'wave_weighting', 'C');
tc.verifyEqual(audio(), SQAT_GUI_weight(ref, fs, 'C'), 'AbsTol', 1e-12);
tc.verifyEqual(findobj(w, 'Tag', 'spectrogram').Colorbar.Label.String, 'SPL (dBC)');
il_set(w, 'wave_weighting', 'Z');
tc.verifyEqual(audio(), ref);
end

function test_waveform_spectrogram_options(tc)
% The window type, the FFT degree and the overlap of the plain spectrogram
% start at Hann, 10 and 50 %, change the map as expected, and the degree is a
% spinner limited to 6 to 16.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
wd = findobj(w, 'Tag', 'spec_window');
tc.verifyEqual(wd.ItemsData, {'hann', 'hamming', 'rect', 'blackmanharris'});
tc.verifyEqual(wd.Value, 'hann');
tc.verifyEqual(findobj(w, 'Tag', 'spec_degree').Value, 10);
tc.verifyEqual(findobj(w, 'Tag', 'spec_overlap').Value, 50);
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
c_hann = surf().CData;
[x, fs] = audioread(tc.TestData.wav_tone);
[~, ~, Lh] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 50);
tc.verifyEqual(max(c_hann(:)), max(Lh(:)), 'AbsTol', 1e-9);   % calibrated at 94 dBFS
il_set(w, 'spec_window', 'rect');
tc.verifyNotEqual(surf().CData, c_hann);
il_set(w, 'spec_window', 'hann');
il_set(w, 'spec_degree', 12);
f4096 = (0:2048)' * fs / 4096;
tc.verifyEqual(numel(surf().YData), nnz(f4096 >= 20));
il_set(w, 'spec_overlap', 75);
tc.verifyEqual(numel(surf().XData), floor((numel(x) - 4096) / 1024) + 1);
% values outside the range go back into it
d = findobj(w, 'Tag', 'spec_degree');
tc.verifyClass(d, 'matlab.ui.control.Spinner');           % arrows, no typing and Enter
tc.verifyEqual(d.Limits, [6 16]);
tc.verifyEqual(d.Step, 1);
d.Value = 16; d.ValueChangedFcn(d, []);                   % what the up arrow gives at the top
tc.verifyEqual(d.Value, 16);
tc.verifyEqual(numel(findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface').YData), nnz((0:32768)' * 48000 / 65536 >= 20));
d.Value = 10; d.ValueChangedFcn(d, []);
o = findobj(w, 'Tag', 'spec_overlap');
tc.verifyClass(o, 'matlab.ui.control.Spinner');           % arrows too
tc.verifyEqual(o.Limits, [0 95]);
tc.verifyEqual(o.Step, 5);
o.Value = 95; o.ValueChangedFcn(o, []);                   % what the up arrow gives at the top
tc.verifyEqual(o.Value, 95);
o.Value = 0; o.ValueChangedFcn(o, []);                    % and the down arrow at the bottom: no overlap
[t0, ~, ~] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 0);
tc.verifyEqual(numel(findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface').XData), numel(t0));
end

function test_waveform_spectrogram_takes_a_window_from_a_file(tc)
% Import window reads a window from a text or .mat file: a Hann window from
% a file gives the Hann spectrogram, a window of another length is resampled,
% and a second import replaces the first in the menu.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
c_hann = surf().CData;
wd = findobj(w, 'Tag', 'spec_window');
% a Hann window written to a file gives the Hann spectrogram
file = fullfile(tc.TestData.dir_tmp, 'my_hann.csv');
writematrix(SQAT_GUI_window('hann', 1024), file);
setappdata(w, 'sqat_next_file', file);                     % stands for the file dialog
il_press(w, 'import_window');
tc.verifyEqual(wd.Value, 'custom');
tc.verifySubstring(wd.Items{end}, 'my_hann');
tc.verifyEqual(surf().CData, c_hann, 'AbsTol', 1e-2);     % the text file keeps 5 digits
w_hann = SQAT_GUI_window('hann', 1024);
file = fullfile(tc.TestData.dir_tmp, 'my_hann.mat');
save(file, 'w_hann');
setappdata(w, 'sqat_next_file', file);
il_press(w, 'import_window');
tc.verifyEqual(surf().CData, c_hann, 'AbsTol', 1e-9);
% a window of another length is resampled
file = fullfile(tc.TestData.dir_tmp, 'my_ramp.txt');
writematrix([0 1 3 1 0], file);
setappdata(w, 'sqat_next_file', file);
il_press(w, 'import_window');
tc.verifyEqual(numel(wd.Items), 5, 'a second import must replace the first one');
tc.verifySubstring(wd.Items{end}, 'my_ramp');
tc.verifyNotEqual(surf().CData, c_hann);
% a file that is no window is refused and the window stays
file = fullfile(tc.TestData.dir_tmp, 'not_a_window.csv');
writematrix([1 2; 3 4], file);
before = wd.Value;
setappdata(w, 'sqat_next_file', file);
il_press(w, 'import_window');
tc.verifyEqual(wd.Value, before);
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'not_a_window');
end

%% Enhanced STFT in the waveform window ------------------------------------

function test_waveform_enhanced_stft_switch(tc)
% The enhanced STFT is off by default; switched on, it fades the window, FFT
% degree, overlap and import controls and changes the map.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
sw = findobj(w, 'Tag', 'spec_enhanced');
tc.verifyEqual(sw.Value, 'Off');
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
c_plain = surf().CData;
il_set(w, 'spec_enhanced', 'On');
for tag = {'spec_degree', 'spec_overlap', 'spec_window', 'import_window'}
    tc.verifyEqual(char(findobj(w, 'Tag', tag{1}).Enable), 'off', ['faded: ' tag{1}]);
end
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Title.String, 'Enhanced STFT');
md = findobj(w, 'Tag', 'spec_enhanced_mode');
tc.verifyEqual(char(md.Enable), 'on');
tc.verifyEqual(md.Value, 'readable');
c = surf().CData;
[~, y] = il_map_centres(surf());
tc.verifyEqual(surf().FaceColor, 'texturemap');
il_set(w, 'spec_enhanced_mode', 'sharp');
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Title.String, 'sharp');
tc.verifyNotEqual(surf().CData, c);
il_set(w, 'spec_enhanced_mode', 'readable');
tc.verifyEqual(surf().CData, c);
[~, i] = max(max(c, [], 2));
tc.verifyEqual(y(i), 1000, 'AbsTol', 25);                   % the tone sits at 1 kHz on the log grid
tc.verifyEqual(10*log10(median(sum(10.^(c/10), 1))), 60, 'AbsTol', 0.1);   % the cells of a column add up to the 60 dB SPL of the tone
il_set(w, 'spec_enhanced', 'Off');
for tag = {'spec_degree', 'spec_overlap', 'spec_window', 'import_window'}
    tc.verifyEqual(char(findobj(w, 'Tag', tag{1}).Enable), 'on', ['back: ' tag{1}]);
end
tc.verifyEqual(char(md.Enable), 'off');
tc.verifyEqual(surf().CData, c_plain);
end

function test_enhanced_stft_zoom_recomputes(tc)
% a zoom on the enhanced map recomputes the excerpt with finer columns; zooming out restores the full map
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
t_full = surf().XData;
L_full = surf().CData;
t = il_map_centres(surf());
tc.verifyGreaterThan(t(2) - t(1), 0.0011);                    % 3 s in 2000 columns: 1.5 ms per column
ax.XLim = [1 1.5];
zoom_now = getappdata(w, 'sqat_spec_zoom');
zoom_now();
[t_z, y] = il_map_centres(surf());
tc.verifyEqual(t_z(2) - t_z(1), 0.001, 'AbsTol', 1e-9);        % 1 ms columns in the excerpt
tc.verifyLessThan(t_z(1), 1);
tc.verifyGreaterThan(t_z(end), 1.5);
tc.verifyLessThan(t_z(end) - t_z(1), 1.2);                     % only the excerpt and its margins
[~, i] = max(max(surf().CData, [], 2));
tc.verifyEqual(y(i), 1000, 'AbsTol', 10);
ax.XLim = [0 3];
zoom_now();
tc.verifyEqual(surf().XData, t_full);
tc.verifyEqual(surf().CData, L_full);
end

function test_enhanced_stft_options_keep_the_zoom(tc)
% the maps are kept without weighting: weighting and mode changes keep the zoom, and the
% weighting only adds its curve to each row
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
ax.XLim = [1 1.5];
ax.YLim = [200 5000];
zoom_now = getappdata(w, 'sqat_spec_zoom');
zoom_now();
L_z = surf().CData;
[t_z, f_z] = il_map_centres(surf());
cl = ax.CLim;
il_set(w, 'wave_weighting', 'A');
tc.verifyEqual(ax.CLim, cl);                                   % the colours stay with the weighting
tc.verifyEqual(ax.XLim, [1 1.5]);
tc.verifyEqual(ax.YLim, [200 5000]);
tc.verifyEqual(il_map_centres(surf()), t_z);                   % the same excerpt
d = surf().CData - L_z;
tc.verifyEqual(d, repmat(d(:, 1), 1, numel(t_z)), 'AbsTol', 1e-9);   % one offset per row
a = SQAT_GUI_weight_curve(f_z, 48000, 'A');
tc.verifyEqual(d(:, 1), a(:), 'AbsTol', 1e-6);
il_set(w, 'spec_enhanced_mode', 'sharp');
tc.verifyEqual(ax.CLim, cl);                                   % and with the mode
tc.verifyEqual(ax.XLim, [1 1.5]);
tc.verifyEqual(ax.YLim, [200 5000]);
tc.verifyEqual(il_map_centres(surf()), t_z);
tc.verifySubstring(ax.Title.String, 'sharp');
il_set(w, 'wave_weighting', 'Z');
il_set(w, 'spec_enhanced_mode', 'readable');
tc.verifyEqual(surf().CData, L_z);                             % back where it started
end

function test_enhanced_stft_whole_file_reuses_the_full_map(tc)
% After a zoom of the enhanced spectrogram, going back to the whole file shows
% the full map computed before, with the same times and levels, without
% computing it again.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
x_full = surf().XData;
L_full = surf().CData;
ax.XLim = [1 1.5];
ax.YLim = [500 2000];
zoom_now = getappdata(w, 'sqat_spec_zoom');
zoom_now();
tc.verifyNotEqual(surf().XData, x_full);
ax.XLim = [0 3];                                           % a double click restores the view
ax.YLim = [20 24000];
zoom_now();
tc.verifyEqual(surf().XData, x_full);
tc.verifyEqual(surf().CData, L_full);
end

function test_spectrogram_title_gives_the_frequency_resolution(tc)
% The title of the plain spectrogram gives the window, the points, the
% overlap and the frequency resolution fs/N, as Gil asked (30.09.2026).
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
ax = findobj(fig, 'Tag', 'spectrogram');
t = ax.Title.String;
n = sscanf(regexp(t, '\d+ points', 'match', 'once'), '%d');
tc.verifySubstring(t, sprintf('%cf %.1f Hz', 916, tc.TestData.fs / n));
end

function test_waveform_and_spectrogram_share_the_time_axis(tc)
% a zoom or a pan on either plot moves the other, inside the file and never under 50 ms;
% a new signal or tab goes back to the whole file; the frequency stays apart
fig = SQAT_GUI({tc.TestData.wav_tone, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axw = findobj(w, 'Tag', 'waveform_axes');
axs = findobj(w, 'Tag', 'spectrogram');
both = @() [axw.XLim; axs.XLim];
f_lim = axs.YLim;
axw.XLim = [1 1.5]; drawnow
tc.verifyEqual(both(), [1 1.5; 1 1.5]);
tc.verifyEqual(axs.YLim, f_lim);
axs.XLim = [0.5 2]; axs.YLim = [100 5000]; drawnow
tc.verifyEqual(both(), [0.5 2; 0.5 2]);
axw.XLim = [-1 2]; drawnow
tc.verifyEqual(both(), [0 2; 0 2]);
axs.XLim = [1 1.01]; drawnow
tc.verifyEqual(both(), [0 2; 0 2]);
axs.XLim = [0 3]; drawnow
tc.verifyEqual(both(), [0 3; 0 3]);
axw.XLim = [1 1.5]; drawnow
il_pick_tab(findobj(w, 'Tag', 'wave_tabs'), 2); drawnow
tc.verifyEqual(both(), [0 3; 0 3]);
axw.XLim = [1 1.5]; drawnow
il_activate_signal(fig, 1); drawnow
tc.verifyEqual(both(), [0 3; 0 3]);
end

function test_waveform_zoom_shows_every_sample_of_a_long_signal(tc)
% A 45 s signal is drawn with one sample in two at full view; a zoom to one
% second draws every sample of the view, and the amplitude axis stays where it
% was.
fs = 48000;
wav = fullfile(tc.TestData.dir_tmp, 'long_45s.wav');
audiowrite(wav, 0.1*sin(2*pi*1000*(0:1/fs:45-1/fs)'), fs, 'BitsPerSample', 32);
fig = SQAT_GUI({wav}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
axw = findobj(w, 'Tag', 'waveform_axes');
ln = findobj(axw, 'Tag', 'wave_line');
tc.verifyEqual(ln.XData(2) - ln.XData(1), 2/fs, 'AbsTol', 1e-12);   % 2.16e6 samples: one in two
y_lim = axw.YLim;
tc.verifyEqual(y_lim, 1.05 * max(abs(ln.YData)) * [-1 1], 'RelTol', 1e-3);   % the whole wave, both halves
axw.XLim = [10 11]; drawnow
tc.verifyEqual(ln.XData(2) - ln.XData(1), 1/fs, 'AbsTol', 1e-12);
tc.verifyLessThanOrEqual(ln.XData(1), 10);
tc.verifyGreaterThanOrEqual(ln.XData(end), 11);
tc.verifyEqual(axw.YLim, y_lim);
end

function test_enhanced_stft_follows_a_zoom_of_the_waveform(tc)
% A zoom of the waveform to half a second recomputes the enhanced map of that
% excerpt in 1 ms columns, and the whole file brings the full map back. Needs
% a display.
il_needs_display(tc);
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
axs = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(axs, 'Type', 'surface');
x_full = surf().XData;
axw = findobj(w, 'Tag', 'waveform_axes');
axw.XLim = [1 1.5]; drawnow
pause(1);                                                  % the debounce of 0.3 s
t_z = il_map_centres(surf());
tc.verifyEqual(t_z(2) - t_z(1), 0.001, 'AbsTol', 1e-9);      % the excerpt, in 1 ms columns
axw.XLim = [0 3]; drawnow
pause(1);
tc.verifyEqual(surf().XData, x_full);
end

function test_enhanced_stft_colour_floor_moves_by_5_dB(tc)
% The up and down arrows move the bottom of the colour scale of the enhanced
% map by 5 dB, the floor survives a redraw, and the range stays between 5 dB
% and 60 dB below the default. Needs a display.
il_needs_display(tc);
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
top = ax.CLim(2);
tc.verifyEqual(ax.CLim, top + [-45 0], 'AbsTol', 1e-9);
w.KeyPressFcn(w, struct('Key', 'uparrow'));
tc.verifyEqual(ax.CLim, top + [-40 0], 'AbsTol', 1e-9);
w.KeyPressFcn(w, struct('Key', 'uparrow'));
tc.verifyEqual(ax.CLim, top + [-35 0], 'AbsTol', 1e-9);
il_set(w, 'spec_enhanced_mode', 'sharp');                      % the floor stays through a redraw
tc.verifyEqual(ax.CLim(2) - ax.CLim(1), 35, 'AbsTol', 1e-9);
il_set(w, 'spec_enhanced_mode', 'readable');
w.KeyPressFcn(w, struct('Key', 'downarrow'));
w.KeyPressFcn(w, struct('Key', 'downarrow'));
tc.verifyEqual(ax.CLim, top + [-45 0], 'AbsTol', 1e-9);
for k = 1:20, w.KeyPressFcn(w, struct('Key', 'uparrow')); end
tc.verifyEqual(ax.CLim, top + [-5 0], 'AbsTol', 1e-9);         % at least 5 dB of colour
for k = 1:30, w.KeyPressFcn(w, struct('Key', 'downarrow')); end
tc.verifyEqual(ax.CLim, top + [-105 0], 'AbsTol', 1e-9);       % at most 60 dB below the default
end

function test_enhanced_stft_comes_from_the_background_pool(tc)
% With the parallel pool, the window shows a waiting note until the enhanced
% map arrives, the map equals the direct computation, and the sharp map comes
% with the same job. Needs a display.
il_needs_display(tc);
% with the pool the window waits with a note, and the map that arrives equals the direct computation
il_use_pool(tc);
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
tc.verifyEmpty(surf());
tc.verifyNotEmpty(findobj(ax, 'Tag', 'spec_wait'));
il_wait_until(@() ~isempty(surf()), 60);
tc.assertNotEmpty(surf());
tc.verifyEmpty(findobj(ax, 'Tag', 'spec_wait'));
[x, fs] = SQAT_GUI_load(tc.TestData.wav_tone, 94, 1);
[~, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
keep = f >= 20;
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
tc.verifyEmpty(strfind(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'ERROR'));
il_set(w, 'spec_enhanced_mode', 'sharp');                  % the same job brought the sharp map
tc.verifyEmpty(findobj(ax, 'Tag', 'spec_wait'));
[~, ~, L] = SQAT_GUI_enhanced_stft(x, fs, 'sharp');
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
end

function test_enhanced_stft_of_a_long_signal_shows_a_preview_first(tc)
% For a 61 s signal the enhanced map first appears as a preview, marked in the
% title, and the full map replaces it later. Needs a display.
il_needs_display(tc);
il_use_pool(tc);
fs = 48000;
t = (0:1/fs:61-1/fs)';
wav = fullfile(tc.TestData.dir_tmp, 'long_tone.wav');
audiowrite(wav, 0.1*sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
fig = SQAT_GUI({wav}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
seen_preview = false;
t0 = tic;
while toc(t0) < 120
    pause(0.05);
    seen_preview = seen_preview || (~isempty(surf()) && contains(ax.Title.String, 'preview'));
    if ~isempty(surf()) && ~contains(ax.Title.String, 'preview')
        break
    end
end
tc.verifyTrue(seen_preview);
tc.assertNotEmpty(surf());
tc.verifyFalse(contains(ax.Title.String, 'preview'));
[x, fs] = SQAT_GUI_load(wav, 94, 1);
[~, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
keep = f >= 20;
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
end

function test_enhanced_stft_follows_a_new_signal_while_computing(tc)
% When another signal becomes active while the map of the last one is being
% computed, the window shows the map of the new signal, and a late map of the
% old one does not replace it. Needs a display.
il_needs_display(tc);
% a map still on its way for the last signal never replaces the map of the new one
il_use_pool(tc);
fig = SQAT_GUI({tc.TestData.wav_tone, tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_activate_signal(fig, 1);
w = fig;
il_set(w, 'spec_enhanced', 'On');
il_activate_signal(fig, 2);
tc.verifySubstring(findobj(w, 'Tag', 'waveform_axes').Title.String, 'two_tones');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
il_wait_until(@() ~isempty(surf()), 60);
tc.assertNotEmpty(surf());
pause(3);                                                      % time for a late map of the tone to arrive
tc.verifySubstring(findobj(w, 'Tag', 'waveform_axes').Title.String, 'two_tones');
[x, fs] = SQAT_GUI_load(tc.TestData.wav_two, 94, 1);
[~, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
keep = f >= 20;
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
end

function test_waveform_close_while_a_map_is_computed(tc)
% Closing the main window while enhanced maps are computed leaves no timer
% running. Needs a display.
il_needs_display(tc);
il_use_pool(tc);
timers_before = numel(timerfindall);                           % before the window: it starts its own
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig(isvalid(fig))));
il_set(fig, 'spec_enhanced', 'On');
fig.CloseRequestFcn(fig, []);
tc.verifyFalse(isvalid(fig));
pause(5);                                                      % any map on its way arrives meanwhile
tc.verifyEqual(numel(timerfindall), timers_before);
end

%% Helpers -----------------------------------------------------------------

function draw_two_corners(w, axs, t1, t2)
% arms the filter tool and clicks the two corners of a box from t1 to t2 (s), 1.5 to 2.5 kHz
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [t1 1500 0]));
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [t2 2500 0]));
end

function ok = il_stays_between(ph, lo, hi, seconds)
% true when the playhead stays between lo and hi (s) for the given time
ok = true; t0 = tic;
while toc(t0) < seconds
    pause(0.05);
    ok = ok && ph.Value >= lo && ph.Value <= hi;
end
end

function il_needs_audio(tc)
% the player needs an audio output; a bare CI runner has none ("No audio
% outputs were found"), so the playback never starts there. A virtual one
% (a PulseAudio null sink) is enough
try
    n = numel(audiodevinfo().output);
catch
    n = 0;
end
tc.assumeGreaterThan(n, 0, 'needs an audio output, and this machine has none');
end

function il_needs_realtime_audio(tc)
% the end of a file needs an output that plays in real time: the virtual one
% of the CI runner (a PulseAudio null sink) fills its buffer and stalls
% (02.10.2026), so a 0.2 s sound must end within 1 s here
il_needs_audio(tc);
persistent ok
if isempty(ok)
    fs = 48000;
    p = audioplayer(zeros(round(0.2 * fs), 1), fs);
    play(p);
    t0 = tic;
    while isplaying(p) && toc(t0) < 1
        pause(0.05);
    end
    ok = ~isplaying(p);
    stop(p);
end
tc.assumeTrue(ok, 'the audio output does not play in real time here');
end

function il_needs_display(tc)
% the enhanced spectrogram recomputes on timers and on the background pool;
% on a Linux runner with no display these tests waited forever (30.09.2026),
% so there they need one (a virtual one, Xvfb, is enough)
tc.assumeFalse(isunix && ~ismac && isempty(getenv('DISPLAY')), ...
    'needs a display: the timers of the enhanced spectrogram do not fire without one');
end

function il_wait_until(cond, seconds)
% waits up to 6 s (or the given time) for a condition (the audio device may take seconds to
% start), letting the timers and the background pool callbacks run
if nargin < 2, seconds = 6; end
t0 = tic;
while ~cond() && toc(t0) < seconds
    pause(0.05);
end
end

function [t, f] = il_map_centres(srf)
% time (row) and frequency (column) of the columns and rows of an enhanced map drawn as a
% texture: uniform in time, log-spaced in frequency, between the edges of its one face
[M, F] = deal(size(srf.CData, 2), size(srf.CData, 1));
xe = [min(srf.XData(:)), max(srf.XData(:))];
ye = log([min(srf.YData(:)), max(srf.YData(:))]);
t = xe(1) + ((1:M) - 0.5) * diff(xe) / M;
f = exp(ye(1) + ((1:F)' - 0.5) * diff(ye) / F);
end

function db = il_tone_db(y, fs, f0)
% level of the component of y at f0 (dB re an arbitrary reference), from a windowed FFT
n = numel(y);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/n);
Y = fft(y(:) .* w);
k = round(f0 * n / fs) + 1;
db = 20*log10(max(abs(Y(k-2:k+2))) * 2 / sum(w) + eps);
end

function il_activate_signal(fig, k)
b = findobj(fig, 'Tag', sprintf('signal_name_%d', k));
b.ButtonPushedFcn(b, []);
end

function il_mark_signal(fig, k, tf)
c = findobj(fig, 'Tag', sprintf('signal_tick_%d', k));
c.Value = tf;
c.ValueChangedFcn(c, []);
end

function il_remove_signal(fig, k)
b = findobj(fig, 'Tag', sprintf('signal_remove_%d', k));   % the window is hidden: no question asked
b.ButtonPushedFcn(b, []);
end

function il_signal_channel(fig, k, c)
d = findobj(fig, 'Tag', sprintf('signal_channel_%d', k));
d.Value = c;
d.ValueChangedFcn(d, []);
end

function il_signal_dbfs(fig, k, v)
% the full-scale level of signal k, as the calibration dialog sets it
set_calibration = getappdata(fig, 'sqat_set_calibration');
set_calibration(k, 'dbfs', v);
end

function names = il_signal_names(fig)
names = {};
k = 1;
while ~isempty(findobj(fig, 'Tag', sprintf('signal_name_%d', k)))
    names{end+1} = findobj(fig, 'Tag', sprintf('signal_name_%d', k)).Text; %#ok<AGROW>
    k = k + 1;
end
end

function nums = il_signal_numbers(fig)
nums = {};
for k = 1:numel(il_signal_names(fig))
    nums{end+1} = findobj(fig, 'Tag', sprintf('signal_number_%d', k)).Text; %#ok<AGROW>
end
end

function il_select_metrics(fig, ids)
set_analyses = getappdata(fig, 'sqat_set_analyses');
set_analyses(ids);
end

function w = il_open_params(fig, k)
b = findobj(fig, 'Tag', sprintf('analysis_params_%d', k));
b.ButtonPushedFcn(b, []);
w = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_params');
end

function tf = il_has_theme()
% the theme function, and with it the dark theme of the GUI, came in R2025a
tf = exist('theme', 'file') > 0;
end

function il_pick_tab(tg, k)
% a click on tab k of a tab group
old = tg.SelectedTab;
tg.SelectedTab = tg.Children(k);
tg.SelectionChangedFcn(tg, struct('NewValue', tg.Children(k), 'OldValue', old));
end

function il_use_pool(tc)
% the background pool for this test only: the other tests leave it off, so that
% no enhanced map is computed behind them
rmappdata(groot, 'sqat_no_background');
tc.addTeardown(@() setappdata(groot, 'sqat_no_background', true));
end

function il_set(win, tag, value)
c = findobj(win, 'Tag', tag);
c.Value = value;
c.ValueChangedFcn(c, []);
end

function il_press(fig, tag)
b = findobj(fig, 'Tag', tag);
b.ButtonPushedFcn(b, []);
end

function il_stop_once_running(t, fig, stop_run)
% presses Stop once the first metric has started, so the run ends after it
if any(contains(findobj(fig, 'Tag', 'console').Value, 'Running'))
    stop_run();
    stop(t);
end
end

function w = il_window(tag)
w = findall(groot, 'Type', 'figure', 'Tag', tag);
end

function figs = il_sqat_figures()
% visible figures opened by the metrics (the interface windows carry a SQAT_GUI tag)
figs = il_all_sqat_figures();
figs = figs(strcmp({figs.Visible}, 'on'));
end

function figs = il_all_sqat_figures()
% every figure drawn by a metric, shown or kept hidden by the interface
figs = findall(groot, 'Type', 'figure');
if isempty(figs)
    return
end
figs = figs(~startsWith({figs.Tag}, 'SQAT_GUI'));
end

function v = il_value(fig, file, metric, quantity)
T = findobj(fig, 'Tag', 'results_table').Data;
v = T.Value(strcmp(T.File, file) & strcmp(T.Metric, metric) & strcmp(T.Quantity, quantity));
end
