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

function test_gui_opens_with_one_analysis_and_every_metric(tc)
% The main window opens with one analysis, ISO 532-1 loudness, whose menu
% offers every metric of SQAT_GUI_metrics.
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyClass(fig, 'matlab.ui.Figure');
dd = findobj(fig, 'Tag', 'analysis_metric_1');
m = SQAT_GUI_metrics;
tc.verifyEqual(dd.ItemsData, {m.id});
tc.verifyEqual(dd.Value, 'Loudness_ISO532_1');                 % one analysis at start
tc.verifyEmpty(findobj(fig, 'Tag', 'analysis_metric_2'));
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
tc.verifyNumElements(il_gui_state(fig, 'loaded'), 2);
% each signal has its channel (All by default for a stereo file) and its calibration
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Value, '1');
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_2').Value, 'All');
c = findobj(fig, 'Tag', 'signal_cal_1');
tc.verifySubstring(c.Text, '94');
il_signal_dbfs(fig, 1, 94);
tc.verifyEqual(findobj(fig, 'Tag', 'signal_cal_1').FontAngle, 'normal');   % set, even to the same value
end

function test_gui_signal_list_ticks_and_removes_signals(tc)
% Only the ticked signals are analysed; with none ticked the run warns and
% does nothing; removing a signal takes its results with it and updates the
% file count down to "No files loaded".
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_mark_signal(fig, 2, false);                       % only the first is analysed
tc.verifyEqual(sum([il_gui_state(fig, 'loaded').marked]), 1);
tc.verifyNumElements(il_gui_state(fig, 'analyses'), 1);
il_press(fig, 'run');
T = il_results(fig);
tc.verifyEqual(unique(T.File), {'tone_mono.wav'});
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, '1 file(s)');
il_mark_signal(fig, 2, true);
il_mark_signal(fig, 1, false);
il_press(fig, 'run');
T = il_results(fig);
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'});
% nothing ticked: a warning and no run
il_mark_signal(fig, 2, false);
il_press(fig, 'run');
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No signals ticked');
% removing a signal takes its results away
il_mark_signal(fig, 1, true); il_mark_signal(fig, 2, true);
il_press(fig, 'run');
il_remove_signal(fig, 1);
T = il_results(fig);
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'});
tc.verifyEqual(il_signal_names(fig), {'tone_1k_60dB.wav'});
tc.verifyNumElements(il_gui_state(fig, 'loaded'), 1);
il_remove_signal(fig, 1);
tc.verifyEmpty(il_gui_state(fig, 'loaded'));
tc.verifyEmpty(il_results(fig));
il_press(fig, 'run');                                % the interface goes on
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No files loaded');
end

%% Running the analyses ----------------------------------------------------

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
tc.verifySubstring(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'time-varying');
tc.verifySubstring(findobj(fig, 'Tag', 'analysis_summary_1').Text, '0.5');
c.Value = 1; c.ValueChangedFcn(c, []);          % stationary again, for the run
il_press(fig, 'add_analysis');                  % the window stays open while the list changes
tc.verifyTrue(isvalid(pw));
il_press(fig, 'run');
T = il_results(fig);
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
T = il_results(fig);
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
SQAT_GUI_export(il_results(fig), xlsx, describe());
T = il_results(fig);
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
tc.verifySubstring(findobj(fig2, 'Tag', 'signal_cal_2').Text, '100');
tc.verifyEqual(findobj(fig2, 'Tag', 'analysis_metric_1').Value, 'Do_SLM');
tc.verifyEqual(findobj(fig2, 'Tag', 'analysis_metric_2').Value, 'Loudness_ISO532_1');
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
tc.verifySubstring(findobj(fig, 'Tag', 'signal_cal_1').Text, '94');
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
T = il_results(fig);
full = @(c) T.Cal_dB_SPL(find(strcmp(T.Channel, c), 1));
tc.verifyEqual(full('1'), 94 - 20*log10(0.5/sqrt(2)), 'AbsTol', 1e-4);
tc.verifyEqual(full('2'), 94 - 20*log10(0.25/sqrt(2)), 'AbsTol', 1e-4);
N = @(c) T.Value(strcmp(T.Channel, c) & strcmp(T.Quantity, 'Nmean'));
tc.verifyEqual(N('2'), N('1'), 'RelTol', 1e-9);
describe = getappdata(fig, 'sqat_run_description');
S = describe();
tc.verifySubstring(S.Value{strcmp(S.Item, 'Signal #1')}, 'calibration calib. 94 dB');
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
T = il_results(fig);
tc.verifyEqual(unique(T.Metric), {'Roughness_Daniel1997'});
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
% as PDF each plot is its PNG on a page of its size, written by the GUI
setappdata(fig, 'sqat_next_file', fullfile(tc.TestData.dir_tmp, 'plots.pdf'));
il_press(fig, 'save_wave_plots');
for name = {'spectrogram', 'waveform', 'level'}
    pdf = fullfile(tc.TestData.dir_tmp, ['plots_' name{1} '.pdf']);
    tc.verifyTrue(isfile(pdf), name{1});
    tc.verifySubstring(fileread(pdf), '/Im0 Do', name{1});
end
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
tc.verifyEqual(kids(1).NodeData.aid, 'sqat');
tc.verifyEqual(kids(end).NodeData.aid, 'stats');
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

function test_gui_save_dialog_writes_every_ticked_format(tc)
% The formats are tick boxes: PNG and PDF ticked write both of the SQAT figure
% and of the analysis, the statistics CSV goes once, and the PDF of the SQAT
% figure is one image on a page of its size, written by the GUI. No format
% ticked keeps the dialog open and says so in the console.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
il_press(il_window('SQAT_GUI_graphs'), 'graph_save');
d = il_window('SQAT_GUI_save');
items = findobj(d, 'Tag', 'save_items');
kids = items.Children(1).Children;
aids = arrayfun(@(n) n.NodeData.aid, kids, 'UniformOutput', false);
items.CheckedNodes = kids(ismember(aids, {'sqat', 'loudness', 'stats'}));
out_dir = fullfile(tc.TestData.dir_tmp, 'saved_formats'); mkdir(out_dir);
h = findobj(d, 'Tag', 'save_folder'); h.Value = out_dir;
set(findobj(d, 'Tag', 'save_png'), 'Value', false);
il_press(d, 'save_do');
tc.verifyNotEmpty(il_window('SQAT_GUI_save'));                 % no format: the dialog stays
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'at least one format');
tc.verifyEmpty(setdiff({dir(out_dir).name}, {'.', '..'}));   % nothing written
set(findobj(d, 'Tag', 'save_png'), 'Value', true);
set(findobj(d, 'Tag', 'save_pdf'), 'Value', true);
il_press(d, 'save_do');
tc.verifyEmpty(il_window('SQAT_GUI_save'));
names = {dir(out_dir).name};
for ext = {'png', 'pdf'}
    tc.verifyTrue(any(startsWith(names, 'tone_mono_s1_Loudness_ISO532_1_') & endsWith(names, ['.' ext{1}])), ext{1});
    tc.verifyTrue(ismember(['Loudness_ISO532_1_loudness_s1.' ext{1}], names), ext{1});
end
tc.verifyNumElements(dir(fullfile(out_dir, '*.csv')), 1);       % the statistics once
sqat_pdf = names(startsWith(names, 'tone_mono_s1_Loudness_ISO532_1_') & endsWith(names, '.pdf'));
txt = fileread(fullfile(out_dir, sqat_pdf{1}));
tc.verifySubstring(txt, '/Subtype /Image');
tc.verifySubstring(txt, '/Im0 Do');
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), '5 file(s) saved');
end

%% Graphs windows ----------------------------------------------------------

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

%% Waveform window ---------------------------------------------------------

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
tc.verifyTrue(il_playing(b));
t_wait = tic;                                      % the audio device may start late
while ph.Value == 0 && toc(t_wait) < 3
    pause(0.05);
end
tc.verifyGreaterThan(ph.Value, 0);
pause(0.2);                                        % a second press within 0.15 s is taken for a repeated key
il_press(w, 'play');                               % pause
tc.verifyFalse(il_playing(b));
il_press(w, 'stop');
tc.verifyEqual(ph.Value, 0);
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
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Colorbar.Label.String, 'dBA');
il_set(w, 'wave_weighting', 'C');
tc.verifyEqual(audio(), SQAT_GUI_weight(ref, fs, 'C'), 'AbsTol', 1e-12);
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Colorbar.Label.String, 'dBC');
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
il_set(w, 'spec_scale', 'log');                 % the surface holds the map row by row on a log axis
il_set(w, 'spec_enhanced', 'On');
for tag = {'spec_degree', 'spec_overlap', 'spec_window', 'import_window'}
    tc.verifyEqual(char(findobj(w, 'Tag', tag{1}).Enable), 'off', ['faded: ' tag{1}]);
end
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

function test_waveform_spectrogram_colours(tc)
% The colour scale of the spectrogram is the inferno of the toolbox by default;
% the dropdown beside the overlap swaps it for a MATLAB one of the same length,
% keeps the map as it is, and stays usable with the enhanced spectrogram.
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
w = fig;
dd = findobj(w, 'Tag', 'spec_cmap');
tc.verifyEqual(dd.ItemsData, {'inferno', 'jet', 'parula', 'turbo', 'hot', 'gray'});
tc.verifyEqual(dd.Value, 'inferno');
ax = findobj(w, 'Tag', 'spectrogram');
inferno = load('cmap_inferno.txt');
tc.verifyEqual(ax.Colormap, inferno, 'AbsTol', 1e-6);
surf = @() findobj(ax, 'Type', 'surface');
c = surf().CData;
il_set(w, 'spec_cmap', 'jet');
tc.verifyEqual(ax.Colormap, jet(size(inferno, 1)), 'AbsTol', 1e-6);
tc.verifyEqual(surf().CData, c);                           % the colours change, the map is the same
il_set(w, 'spec_window', 'rect');                          % a new drawing keeps the choice
tc.verifyEqual(ax.Colormap, jet(size(inferno, 1)), 'AbsTol', 1e-6);
il_set(w, 'spec_enhanced', 'On');
tc.verifyEqual(char(dd.Enable), 'on');                     % both maps take the colours
il_set(w, 'spec_cmap', 'gray');
tc.verifyEqual(ax.Colormap, gray(size(inferno, 1)), 'AbsTol', 1e-6);
il_set(w, 'spec_enhanced', 'Off');
il_set(w, 'spec_cmap', 'inferno');
tc.verifyEqual(ax.Colormap, inferno, 'AbsTol', 1e-6);
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

function il_in_calibration_dialog(fig, k, fn)
% opens the calibration dialog of signal k and runs fn(d) on it from a timer,
% since the dialog holds the call until it closes; an error in fn closes the
% dialog and is raised here. fn waits 0.3 s after the dialog appears: a click
% while uiwait is still starting deletes the dialog before uiwait reads it
err = [];
seen = 0;
tm = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.05, 'TasksToExecute', 200, 'TimerFcn', @tick);
start(tm);
try
    il_press(fig, sprintf('signal_cal_%d', k));
catch e
    err = e;
end
stop(tm);
delete(tm);                                         % by hand: tick holds this workspace, so no onCleanup would run
if ~isempty(err)
    rethrow(err);
end

    function tick(~, ~)
        d = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_calibration');
        if isempty(d)
            return
        end
        seen = seen + 1;
        if seen < 6
            return
        end
        stop(tm);
        try
            fn(d);
        catch e
            err = e;
        end
        if isvalid(d)
            delete(d);                                  % fn left it open: Cancel
        end
    end
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

function il_menu(fig, tag)
% the menu item of that tag, as a click on it
m = findobj(fig, 'Tag', tag);
m.MenuSelectedFcn(m, []);
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
T = il_results(fig);
v = T.Value(strcmp(T.File, file) & strcmp(T.Metric, metric) & strcmp(T.Quantity, quantity));
end

function T = il_results(fig)
% the single values of the last run, one row per signal, channel and quantity
T = il_gui_state(fig, 'results');
end

function tf = il_playing(b)
% whether the player of the window of the button b plays (the button shows
% Pause then); the GUI exposes the state of the play for the tests
info = getappdata(ancestor(b, 'figure'), 'sqat_play');
tf = logical(info().playing);
end

function v = il_gui_state(fig, name)
% a variable of the open GUI (results, matrix, stale): the tests read the state
% of a run here, so a change of the layout on screen reaches only this function;
% functions() gives the workspace that SQAT_GUI shares with its callbacks
s = functions(getappdata(fig, 'sqat_run_description'));
v = s.workspace{1}.(name);
end
