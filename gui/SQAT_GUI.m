function varargout = SQAT_GUI(files, varargin)
% function fig = SQAT_GUI(files, varargin)
%
%   Graphical interface to the metrics of SQAT, laid out as the interface of
%   pySQAT. It loads .wav files into a list of signals (a tick marks a signal
%   for use, the bin removes it; each signal has its own channel and calibration),
%   runs the ticked metrics with the chosen parameters on the ticked signals,
%   and lists their single values. The channel to analyse is one channel of
%   the file or All (the default for a stereo file): the ECMA-418-2 metrics take a stereo pair in one call and
%   return the left, the right and (except the tonality) the combined
%   binaural result, so a pair runs once. A run keeps the results of the last
%   run whose signal (channel and calibration) and analysis (metric and
%   parameters) did not change, so an analysis added to the list runs alone.
%
%   A graphs window opens at the end of a run (Open Graphs Window opens it
%   again after it is closed) and plots one metric of the ticked signals. The analysis is
%   chosen in the window (a time series, a profile over the critical bands, a
%   map of band against time, or the statistics): the signals are overlaid,
%   and the maps sit side by side on one colour scale. For one signal the
%   window also offers the figure that the SQAT function draws (with the
%   inferno colour scale of the toolbox, cmap_inferno.txt) and all the analyses at once. Each signal of the list
%   has a number, and the plots tag a curve as Signal #1, ch1. When a signal
%   has several channels, the window offers a channel choice per signal (1, 2,
%   Binaural when the metric gives one, or All), so that a mono signal can be
%   compared with either channel or the binaural result of a stereo signal.
%   Pin keeps a window with its signals and its results, so that a later run
%   leaves it as it is and Open Graphs Window opens another one to compare
%   with. Save in a graphs window opens
%   a dialog to tick the signals and, per metric, the SQAT figure, the analyses
%   (the signals overlaid) and the statistics (CSV), as PNG or PDF. The
%   Waveform tab shows the waveform and the spectrogram of the signal on
%   screen, with its player, and the Results tab the matrix of single
%   values. The results go to a spreadsheet.
%
%   Every value and every plot comes from the output of the SQAT function:
%   the interface builds the call (see SQAT_GUI_metrics) and reads the output
%   (see SQAT_GUI_extract). The analysis draws the figure of the signal on
%   screen in the same call, so each metric runs once; the figure of another
%   signal is drawn when a graphs window asks for it.
%
% USAGE
%   SQAT_GUI                          % opens the interface
%   SQAT_GUI({'a.wav','b.wav'})       % opens it with files loaded
%   fig = SQAT_GUI(files, 'Visible', 'off')
%
% INPUT ARGUMENTS
%   files : cell array of char, paths of .wav files to load (optional)
%   'Visible' : 'on' (default) or 'off'
%
% OUTPUTS
%   fig : the uifigure of the interface
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September and October
% 2026 assisted by Claude Opus 5, Claude Opus 5.5, Claude
% Sonnet 5 and Claude Fable 5.1 (Anthropic). All codes were
% verified by the authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(files)
    files = {};
end
files = cellstr(files);
opts = struct('Visible', 'on');
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

if isempty(which('Loudness_ISO532_1'))
    run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup_SQAT.m'));
end
dir_logos = fullfile(fileparts(mfilename('fullpath')), 'logos');
icon_remove = fullfile(dir_logos, 'trash.png');          % the remove buttons of the lists
green = [0.13 0.55 0.37];

%% State
metrics = SQAT_GUI_metrics;
% the analyses of the list: a number that is never given again, a metric and
% its parameters; the same metric can be there more than once, and a later
% entry of a metric is keyed Metric_id#n (n its number)
analyses = struct('key', 'Loudness_ISO532_1', 'id', 'Loudness_ISO532_1', 'n', 1, ...
    'p', il_default_params(metrics(strcmp({metrics.id}, 'Loudness_ISO532_1'))));
next_analysis = 2;                                    % the number of the next analysis
loaded = struct('path', {}, 'name', {}, 'nch', {}, 'fs', {}, 'marked', {}, 'id', {}, 'channel', {}, 'dBFS', {}, 'cal_set', {}, 'cal', {});
% dBFS: the full-scale level of each channel (dB SPL); cal: how it was set (method, level, file, label)
last_cal = struct('method', 'dbfs', 'level', 94, 'file', '');   % the calibration dialog opens with the last choice
ask_remove = true;                                    % the bin of a signal asks first, until told not to
next_id = 1;                                          % the number of the next signal loaded
cal_help = ['A WAV file stores numbers between -1 and +1, with no unit. To read them as sound ' ...
    'pressure, in Pascal, calibration is necessary. Three possibilities are provided here:' newline newline ...
    '- Full-scale level (dBFS): the level in dB SPL of a sample of value 1, known from the recording ' ...
    'chain or the documentation of the file. SQAT assumes 94 dB (1 = 1 Pa) when nothing is set.' newline newline ...
    '- Calibrator recording: a recording of a sound level calibrator (1 kHz at 94 dB, for example) made ' ...
    'with the same setup; its rms against the level of the calibrator gives the full scale. A stereo ' ...
    'recording calibrates each channel on its own.' newline newline ...
    '- Relative level: with no information, the rms of the file is set to a chosen level. The results ' ...
    'then only compare signals among themselves; they are not absolute levels.'];
active_idx = 0;                                       % the signal on screen in the player
wave_ch_by_id = [];                                    % channel chosen by a tab of the waveform window, by signal number
results = il_empty_results();
store = il_empty_store();                             % analyses of each signal, metric and channel
player = [];
player_fs = [];
play_start = 0;                                       % sample where the next play starts, 0 for the start of the region
playing = false;                                      % the interface's own state: the audio device answers late
play_t0 = tic;                                        % start of the last play, and the time it should last
play_span = 0;
play_y = [];                                          % the sound that plays (filter and weighting applied), until they change
play_map = struct('sample', 1, 'n_first', 0, 'loop_from', 1, 'n_rep', 0);   % where the buffer of the play comes from
play_peak = 0;                                        % peak of the buffer that went to the player
box_loop = false;                                     % the play that runs is inside the filter box
toggle_clock = tic;                                   % time of the last play or pause, to drop a key event that comes twice
last_toggle = -1;
wave_x = [];                                          % the signal of the waveform window, in pascals
wave_y = [];                                          % the same as the file has it, for playback
wave_fs = [];
wave_key = '';                                        % file and channel on screen
boxes = zeros(0, 4);                                  % [t1 t2 f1 f2] removed from what plays
box_corner = [];                                      % first corner of a box being drawn
drag_active = false;                                  % the mouse is being followed to size a box
spec_custom = [];                                     % window read from a file
spec_view = [];                                       % enhanced maps on screen: the full file and the zoomed excerpt
spec_signal = '';                                     % file, channel and calibration on screen: the key of the cache
spec_cache = struct('key', {}, 'map', {});            % full enhanced maps already computed, not weighted
spec_dlg = [];                                        % the wait message of the spectrogram
spec_computing = false;                               % an enhanced map is being computed
spec_jobs = struct('key', {}, 'job', {});             % full maps queued on the background pool
spec_preview = struct('key', {}, 'map', {});          % quick full maps of a long signal, until the exact ones come
preview_min = 60;                                     % a signal longer than this (s) gets a preview first
spec_pool = [];                                       % the background pool, when MATLAB has one
spec_pool_tried = false;
prefetch_timer = [];                                  % starts the background maps of a new signal
spec_busy = false;                                    % the spectrogram limits are being set by the code
spec_timer = [];                                      % waits for the zoom to settle before recomputing
spec_top = [];                                        % top of the colour scale of the spectrogram (dB)
spec_range = 45;                                      % its default span (dB)
spec_black = 0;                                       % dB added to its bottom: more is more black
theme_style = il_if(il_has_theme(), 'dark', 'light');   % no themes before R2025a: the default light look
graph_figs = gobjects(0);                              % the graphs windows
last_metric = '';                                     % the metric the last graphs window showed
win_wave = [];                                        % the figure of the player: the main window
wave_view = [];                                       % the grid of the player
ax_wave = [];
ax_lvl = [];                                          % the sound level vs time, below the waveform
lvl_L = [];                                           % its level at every sample (dB), for the indicators
ax_spec = [];
btn_play = [];
run_settings = struct('signals', loaded, 'analyses', analyses);      % of the last run
stop_requested = false;                                             % the Stop button or the dialog
dlg = [];                                                           % progress dialog of a run
cache = struct('file', {}, 'metric', {}, 'figs', {});               % SQAT figures, hidden
cmap = load('cmap_inferno.txt');                                    % the colour scale of the toolbox (utilities/ECMA418_2)

%% Main window
fig = uifigure('Name', 'SQAT: Sound Quality Analysis Toolbox', ...
    'Position', [40 30 1460 900], 'Visible', opts.Visible, 'Tag', 'SQAT_GUI', ...
    'CloseRequestFcn', @on_close, 'CreateFcn', '', ...   % skips a user default CreateFcn
    'KeyPressFcn', @on_main_key);
m_file = uimenu(fig, 'Text', 'File');
uimenu(m_file, 'Text', 'Open session...', 'MenuSelectedFcn', @(~, ~) on_session('open'));
uimenu(m_file, 'Text', 'Save session...', 'MenuSelectedFcn', @(~, ~) on_session('save'));
main = uigridlayout(fig, [3 2]);
main.RowHeight = {64, '1x', 30};                 % the logo as tall as the Actions panel beside it
main.ColumnWidth = {600, '1x'};                  % the lists get the room, the console the rest

top = uigridlayout(main, [1 2]);
top.Layout.Row = 1; top.Layout.Column = 1;
top.Padding = [20 0 0 0];                       % the logo a little in from the edge
top.ColumnWidth = {84, '1x'};                   % the logo under the row height (1563 x 895 px)
top.ColumnSpacing = 25;                         % from the logo to the title
img_logo = uiimage(top, 'ImageSource', fullfile(dir_logos, 'logo_white.png'), 'Tag', 'logo', ...
    'HorizontalAlignment', 'right');
uilabel(top, 'Text', 'Sound Quality Analysis Toolbox', 'FontSize', 26, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left');             % right after the logo

left = uigridlayout(main, [2 1]);
left.Layout.Row = 2; left.Layout.Column = 1;
left.Padding = [0 0 0 0];
left.RowHeight = {'1x', '1x'};
sig_box = uigridlayout(uipanel(left), [2 1]);          % a box around each list, to set them apart
sig_box.RowHeight = {28, '1x'};
sig_box.Padding = [6 6 6 6];
ana_box = uigridlayout(uipanel(left), [2 1]);
ana_box.RowHeight = {28, '1x'};
ana_box.Padding = [6 6 6 6];
sh = uigridlayout(sig_box, [1 3]);
sh.Padding = [0 0 0 0];
sh.ColumnWidth = {'1x', 100, 130};
uilabel(sh, 'Text', '1   SIGNALS', 'FontWeight', 'bold', 'FontSize', 15);
lbl_files = uilabel(sh, 'Text', 'No files loaded', 'Tag', 'file_count', 'HorizontalAlignment', 'right');
uibutton(sh, 'Text', 'Open WAV files...', 'Tag', 'load_files', 'ButtonPushedFcn', @on_load_files);
signal_list = uigridlayout(sig_box, [1 6], 'Scrollable', 'on', 'Tag', 'signals_list');
signal_list.ColumnWidth = {30, 22, '1x', 62, 96, 26};
signal_list.Padding = [0 0 0 0];
signal_list.RowSpacing = 4;
ah = uigridlayout(ana_box, [1 3]);
ah.Padding = [0 0 0 0];
ah.ColumnWidth = {'1x', 190, 100};
uilabel(ah, 'Text', '2   ANALYSES', 'FontWeight', 'bold', 'FontSize', 15);
uibutton(ah, 'Text', '+ Add metrics...', 'Tag', 'add_metric', 'ButtonPushedFcn', @on_add_metric, ...
    'Tooltip', 'Opens the list of metrics: the ticked ones are added, each with its default parameters');
uibutton(ah, 'Text', 'Copy last', 'Tag', 'add_analysis', 'ButtonPushedFcn', @on_add_analysis, ...
    'Tooltip', 'Adds a copy of the last analysis, to compare the same metric with other parameters');
analysis_list = uigridlayout(ana_box, [1 5], 'Scrollable', 'on', 'Tag', 'analysis_list');
analysis_list.ColumnWidth = {30, 175, '1x', 36, 28};
analysis_list.Padding = [0 0 0 0];

right = uigridlayout(main, [2 1]);
right.Layout.Row = [1 2]; right.Layout.Column = 2;  % Actions up beside the logo, the results below
right.Padding = [0 0 0 0];
right.RowHeight = {64, '1x'};

save_folder = pwd;                                     % the folder of the last save
split_figures = false;                                 % one tab and one file per panel: no control for now

ag = uigridlayout(right, [1 4]);                       % the buttons alone, as tall as the logo row
ag.ColumnWidth = {'1.4x', '1x', '1x', 48};
ag.Padding = [0 4 0 4];
btn_run = uibutton(ag, 'Text', 'Run Analysis', 'Tag', 'run', 'FontWeight', 'bold', 'FontSize', 15, ...
    'BackgroundColor', green, 'FontColor', [1 1 1], 'ButtonPushedFcn', @on_run);
uibutton(ag, 'Text', 'Open Graphs Window', 'Tag', 'open_graphs', 'FontSize', 14, 'ButtonPushedFcn', @on_open_graphs);
uibutton(ag, 'Text', 'Export results...', 'Tag', 'export', 'FontSize', 14, 'ButtonPushedFcn', @on_export);
btn_theme = uibutton(ag, 'Text', char(9788), 'FontSize', 22, 'Tag', 'theme', ...
    'Tooltip', 'Light theme', 'ButtonPushedFcn', @on_theme);   % a sun, or a moon in the light theme

% the waveform and the spectrogram of the signal on screen, with the player,
% take the whole area, as Gil proposed; the single values of the signals side
% by side, the full table and the log are one tab away, and the plots of the
% results live in the graphs window
tabs = uitabgroup(right);
tab_wave = uitab(tabs, 'Title', 'Waveform');
wave_dock = uipanel(uigridlayout(tab_wave, [1 1], 'Padding', [4 4 4 4]), 'BorderType', 'none', ...
    'Tag', 'waveform_dock');                           % the player
tab_results = uitab(tabs, 'Title', 'Results');
matrix = uitable(uigridlayout(tab_results, [1 1]), 'Data', cell(0, 2), ...
    'ColumnName', {'Analysis', 'Quantity'}, 'RowName', {}, 'Tag', 'results_matrix');
tab_table = uitab(tabs, 'Title', 'Table');
tbl = uitable(uigridlayout(tab_table, [1 1]), 'Data', results, 'Tag', 'results_table');
tab_console = uitab(tabs, 'Title', 'Log');
console = uitextarea(uigridlayout(tab_console, [1 1]), 'Value', {''}, 'Editable', 'off', ...
    'Tag', 'console', 'FontName', 'Monospaced');

status_bar = uigridlayout(main, [1 2]);
status_bar.Layout.Row = 3; status_bar.Layout.Column = [1 2];
status_bar.Padding = [0 0 0 0];
status_bar.ColumnWidth = {'1x', 300};
lbl_status = uilabel(status_bar, 'Text', 'Ready', 'Tag', 'status');
gauge = uigauge(status_bar, 'linear', 'Tag', 'progress', 'Limits', [0 100], 'Value', 0, ...
    'MajorTicks', [], 'MinorTicks', []);

%% Start
win_wave = fig;
build_wave_view(wave_dock);
apply_theme();
add_files(files);
if isempty(loaded)
    refresh_signals();                                 % the empty list, with the hint to start
end
refresh_analyses();
setappdata(fig, 'sqat_set_analyses', @set_analyses);   % the list from metric ids, for the tests
setappdata(fig, 'sqat_run_description', @run_description);   % the Settings sheet, for the tests
setappdata(fig, 'sqat_write_report', @write_report);
setappdata(fig, 'sqat_session', @session);           % save or open a session without the file dialog, for the tests   % the PDF report, without its file dialog, for the tests
setappdata(fig, 'sqat_stop', @on_stop_run);           % the Stop of the progress dialog, which a hidden window has not
setappdata(fig, 'sqat_set_calibration', @set_calibration);   % the calibration without its dialog, for the tests
write_log('Ready. Open WAV files, choose metrics and parameters, then press Run Analysis.');
if nargout > 0
    varargout{1} = fig;   % fig itself stays: the callbacks share it
end

%% Callbacks of the main window

    function on_load_files(~, ~)
        [f, p] = uigetfile({'*.wav', 'WAV files (*.wav)'}, 'Open WAV files', 'MultiSelect', 'on');
        focus_gui();
        if isequal(f, 0)
            return
        end
        add_files(fullfile(p, cellstr(f)));
    end

    function on_signal_channel(k, c)
        loaded(k).channel = c;
        if numel(wave_ch_by_id) >= loaded(k).id
            wave_ch_by_id(loaded(k).id) = 0;             % the waveform follows the list again
        end
        signal_changed(k);
    end

    function ok = set_calibration(k, method, level, calfile)
        % calibrates signal k: the full-scale level of each channel from the method
        % (SQAT_GUI_calibration); a failure leaves the signal as it was and returns false
        if nargin < 4
            calfile = '';
        end
        try
            [dBFS, label] = SQAT_GUI_calibration(method, loaded(k).path, level, calfile);
        catch err
            write_log(['Calibration not changed: ' err.message]);
            ok = false;
            return
        end
        loaded(k).dBFS = dBFS;
        loaded(k).cal_set = true;
        loaded(k).cal = struct('method', method, 'level', level, 'file', calfile, 'label', label);
        last_cal = struct('method', method, 'level', level, 'file', calfile);
        write_log(sprintf('%s calibrated: %s (full scale %s dB SPL).', loaded(k).name, label, ...
            strjoin(arrayfun(@(v) sprintf('%.2f', v), dBFS, 'UniformOutput', false), ', ')));
        refresh_signals();
        signal_changed(k);
        ok = true;
    end

    function on_signal_cal(k)
        % the dialog of the calibration of signal k, opened on the last choice
        c = last_cal;
        if loaded(k).cal_set
            c = loaded(k).cal;
        end
        SQAT_GUI_calibration_dialog(fig, loaded(k).name, c, cal_help, ...
            @(method, level, calfile) set_calibration(k, method, level, calfile));
    end

    function signal_changed(k)
        if any(strcmp(results.Path, loaded(k).path))
            mark_stale();
        end
        if k == active_idx && il_is_open(win_wave)
            draw_waveform_window();
        end
    end

    function on_signal_ticked(k, tf)
        loaded(k).marked = tf;
        update_run_label();
        refresh_windows();
    end

    function on_signal_name(k)
        % the name makes the signal the one of the player and the waveform window
        if k ~= active_idx
            active_idx = k;
            show_active();
            refresh_windows();
        end
    end

    function on_signal_remove(k)
        % the bin asks first: the results and the figures of the signal go with it
        if ask_remove && strcmp(fig.Visible, 'on')
            [go, ask_remove] = il_confirm_remove(fig, loaded(k).name);
            if ~go
                return
            end
        end
        remove_signal(k);
    end

    %% The list of analyses

    function refresh_analyses()
        % one row per analysis: its number, the metric, its parameters in short (in full in the tooltip), the gear and the bin
        delete(analysis_list.Children);
        w = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_params');
        if ~isempty(w) && ~ismember(w(1).UserData, [analyses.n])
            close_params();                  % its analysis is gone
        end
        n = numel(analyses);
        analysis_list.RowHeight = repmat({26}, 1, max(n, 1));
        for k = 1:n
            a = analyses(k);
            uilabel(analysis_list, 'Text', sprintf('#%d', a.n), 'Tag', sprintf('analysis_number_%d', k));
            uidropdown(analysis_list, 'Items', {metrics.label}, 'ItemsData', {metrics.id}, ...
                'Value', a.id, 'Tag', sprintf('analysis_metric_%d', k), ...
                'ValueChangedFcn', @(src, ~) on_analysis_metric(k, src.Value));
            uilabel(analysis_list, 'Text', il_param_summary(metrics, a), 'Tag', sprintf('analysis_summary_%d', k), ...
                'Tooltip', il_param_text(metrics(strcmp({metrics.id}, a.id)), a.p));
            uibutton(analysis_list, 'Text', char(9881), 'FontSize', 16, 'Tag', sprintf('analysis_params_%d', k), ...
                'Tooltip', 'Parameters of this analysis', 'ButtonPushedFcn', @(~, ~) on_analysis_params(k));
            uibutton(analysis_list, 'Text', '', 'Icon', icon_remove, 'Tag', sprintf('analysis_remove_%d', k), ...
                'Tooltip', 'Removes this analysis', 'ButtonPushedFcn', @(~, ~) on_remove_analysis(k));
        end
        SQAT_GUI_paint(analysis_list, theme_style);    % the new rows in the colours of the theme
        update_run_label();
    end

    function set_analyses(ids)
        % the list holds these metrics, numbered from 1, with their default parameters
        analyses = analyses([]);
        for k = 1:numel(ids)
            e = metrics(strcmp({metrics.id}, ids{k}));
            analyses(k) = struct('key', '', 'id', e.id, 'n', k, 'p', il_default_params(e));
        end
        next_analysis = numel(ids) + 1;
        assign_keys();
        refresh_analyses();
    end

    function assign_keys()
        % the first analysis of a metric in the list is keyed by its id (it can
        % share computations), a later one by id#n
        for k = 1:numel(analyses)
            analyses(k).key = analyses(k).id;
            if nnz(strcmp({analyses(1:k).id}, analyses(k).id)) > 1
                analyses(k).key = sprintf('%s#%d', analyses(k).id, analyses(k).n);
            end
        end
    end

    function on_add_analysis(~, ~)
        % a copy of the last analysis, whose parameters are then changed to compare
        if isempty(analyses)
            e = metrics(strcmp({metrics.id}, 'Loudness_ISO532_1'));
            analyses = struct('key', '', 'id', e.id, 'n', next_analysis, 'p', il_default_params(e));
        else
            analyses(end+1) = analyses(end);
            analyses(end).n = next_analysis;
        end
        next_analysis = next_analysis + 1;
        assign_keys();
        refresh_analyses();
    end

    function on_add_metric(~, ~)
        % the list of metrics, to tick several at once (see SQAT_GUI_metric_picker)
        SQAT_GUI_metric_picker(fig, metrics, @add_metrics);
    end

    function add_metrics(ids)
        % a new analysis per metric id, with the defaults of the metric
        for id = ids(:)'
            e = metrics(strcmp({metrics.id}, id{1}));
            a = struct('key', '', 'id', e.id, 'n', next_analysis, 'p', il_default_params(e));
            if isempty(analyses)
                analyses = a;
            else
                analyses(end+1) = a; %#ok<AGROW>
            end
            next_analysis = next_analysis + 1;
        end
        assign_keys();
        refresh_analyses();
    end

    function on_analysis_metric(k, id)
        % another metric is another analysis: a new number, and the defaults of the metric
        e = metrics(strcmp({metrics.id}, id));
        write_log(sprintf('#%d became #%d, %s, with its default parameters.', analyses(k).n, next_analysis, e.label));
        analyses(k).id = id;
        analyses(k).n = next_analysis;
        next_analysis = next_analysis + 1;
        analyses(k).p = il_default_params(e);
        mark_stale();
        assign_keys();
        refresh_analyses();
    end

    function on_remove_analysis(k)
        % the analysis leaves the list, and its results and figures of the last run go with it
        n = analyses(k).n;
        k_r = find([run_settings.analyses.n] == n, 1);   % the keys of the run, by number
        analyses(k) = [];
        assign_keys();
        refresh_analyses();
        if isempty(k_r)
            return
        end
        key = run_settings.analyses(k_r).key;
        run_settings.analyses(k_r) = [];
        results = results(~strcmp(results.Analysis, sprintf('#%d', n)), :);
        store = store(~strcmp({store.metric}, key));
        drop = strcmp({cache.metric}, key);
        for k_c = find(drop)
            delete(cache(k_c).figs(isvalid(cache(k_c).figs)));
        end
        cache = cache(~drop);
        show_results();
        refresh_windows();
    end

    function on_analysis_params(k)
        % a small window with the parameters of analysis k: name, value and unit; a change applies at once
        close_params();
        a = analyses(k);
        e = metrics(strcmp({metrics.id}, a.id));
        n = numel(e.params);
        w = uifigure('Name', sprintf('Parameters: #%d %s', a.n, e.label), ...
            'Position', [fig.Position(1) + 440, fig.Position(2) + 300, 440, 60 + 34 * max(n, 1)], ...
            'Visible', fig.Visible, 'Tag', 'SQAT_GUI_params', 'UserData', a.n, 'CreateFcn', '');
        g = uigridlayout(w, [max(n, 1) + 1, 3]);
        g.ColumnWidth = {'1x', 170, 60};
        g.RowHeight = [repmat({26}, 1, max(n, 1)), {28}];
        if n == 0
            uilabel(g, 'Text', 'This metric has no parameters.');
            uilabel(g, 'Text', '');
            uilabel(g, 'Text', '');
        end
        for k_par = 1:n
            q = e.params(k_par);
            [name, unit] = il_label_unit(q.label);
            uilabel(g, 'Text', name);
            current = a.p.(q.name);
            if strcmp(q.type, 'choice')
                c = uidropdown(g, 'Items', q.options(:, 1)', 'ItemsData', q.options(:, 2)', 'Value', current);
            else
                c = uieditfield(g, 'numeric', 'Value', current, 'Limits', [0 Inf], ...
                    'Tooltip', 'Zero or more');
                if strcmp(q.name, 'dt')
                    c.LowerLimitInclusive = 'off';
                    c.Tooltip = 'More than zero';
                end
            end
            c.Tag = ['param_' q.name];
            c.ValueChangedFcn = @(src, ~) set_param(a.n, q.name, src.Value);
            uilabel(g, 'Text', unit, 'Tag', ['unit_' q.name]);
        end
        uibutton(g, 'Text', 'Reset to defaults', 'Tag', 'params_reset', ...
            'ButtonPushedFcn', @(~, ~) reset_params(a.n, w));
        uibutton(g, 'Text', 'Close', 'Tag', 'params_close', 'ButtonPushedFcn', @(~, ~) delete(w));
        uilabel(g, 'Text', '');
        if il_has_theme()
            theme(w, theme_style);
        end
        SQAT_GUI_paint(w, theme_style);
    end

    function reset_params(num, w)
        k = find([analyses.n] == num, 1);
        e = metrics(strcmp({metrics.id}, analyses(k).id));
        for q = e.params
            c = findobj(w, 'Tag', ['param_' q.name]);
            c.Value = q.value;
            set_param(num, q.name, q.value);
        end
    end

    function set_param(num, name, value)
        k = find([analyses.n] == num, 1);
        if ~isempty(k) && ~isequal(analyses(k).p.(name), value)   % a reset to the same value changes nothing
            analyses(k).p.(name) = value;
            mark_stale();
            set(findobj(analysis_list, 'Tag', sprintf('analysis_summary_%d', k)), ...
                'Text', il_param_summary(metrics, analyses(k)), ...
                'Tooltip', il_param_text(metrics(strcmp({metrics.id}, analyses(k).id)), analyses(k).p));
        end
    end

    function close_params()
        delete(findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_params'));
    end

    function e = run_entry(key)
        % the metric of an analysis of the last run, with its key, its label and its parameters
        k = find(strcmp({run_settings.analyses.key}, key), 1);
        a = run_settings.analyses(k);
        e = metrics(strcmp({metrics.id}, a.id));
        e.key = a.key;
        e.number = sprintf('#%d', a.n);           % the number of the analysis, never given again
        e.label = sprintf('#%d %s', a.n, e.label);
        e.p = a.p;
    end

    function [old, same_key] = reusable(f, a)
        % the entries of the last run for signal f and analysis a, when neither
        % changed since: the same channel and calibration, the same metric and
        % parameters, and every channel of the signal done (a run stopped halfway
        % leaves some out). same_key: the analysis kept its key, so its SQAT figures stay
        old = il_empty_store();
        same_key = false;
        k_s = find(strcmp({run_settings.signals.path}, f.path), 1);
        k_a = find([run_settings.analyses.n] == a.n, 1);
        if isempty(k_s) || isempty(k_a)
            return
        end
        s = run_settings.signals(k_s);
        b = run_settings.analyses(k_a);
        if ~strcmp(s.channel, f.channel) || ~isequal(s.dBFS, f.dBFS) || ~strcmp(b.id, a.id) || ~isequal(b.p, a.p)
            return
        end
        prev = store(strcmp({store.file}, f.path) & strcmp({store.metric}, b.key));
        if isempty(prev) || ~all(ismember(arrayfun(@num2str, channel_list(f), 'UniformOutput', false), {prev.channel}))
            return
        end
        [prev.metric] = deal(a.key);          % a removal before it may have changed its key
        old = prev;
        same_key = strcmp(b.key, a.key);
    end

    function mark_stale()
        % a setting changed after the run: the results on screen no longer follow the settings
        if height(results) == 0 || contains(tab_results.Title, 'changed')
            return
        end
        tab_results.Title = 'Results (settings changed: run again)';
        lbl_status.Text = 'Settings changed since the last run: the results are those of the run.';
        write_log('Settings changed since the last run; the results on screen are those of the run.');
    end

    function show_results()
        % the matrix of single values on the Results tab, the full list on Table;
        % the waveform stays on screen, and the graphs window shows the plots
        tbl.Data = results;
        show_matrix();
    end

    function show_matrix()
        % one row per analysis and quantity, one column per signal and channel, in
        % the order of the run; UserData holds the analysis number of each row
        mx_T = results;
        mx_p = str2double(regexp(mx_T.Quantity, '\d+$', 'match', 'once'));   % N10, R95, LAF50: percentiles
        mx_T = mx_T(isnan(mx_p) | ismember(mx_p, [5 90]), :);  % only the 5 and 90 %, as Gil asked; Table keeps all
        mx_col = strcat(mx_T.Signal, ', ch', mx_T.Channel);
        mx_col = strrep(mx_col, ', chBinaural', ', binaural');
        [mx_cols, ~, mx_ic] = unique(mx_col, 'stable');
        mx_row = strcat(mx_T.Analysis, '|', mx_T.Quantity);
        [mx_rows, mx_ir, mx_jr] = unique(mx_row, 'stable');
        mx_data = cell(numel(mx_rows), 2 + numel(mx_cols));
        for mx_k = 1:numel(mx_rows)
            mx_i = mx_ir(mx_k);
            mx_data(mx_k, 1:2) = {il_key_label(store, run_key(mx_T.Analysis{mx_i})), sprintf('%s (%s)', mx_T.Quantity{mx_i}, mx_T.Unit{mx_i})};
        end
        for mx_k = 1:height(mx_T)
            mx_data{mx_jr(mx_k), 2 + mx_ic(mx_k)} = mx_T.Value(mx_k);
        end
        matrix.Data = mx_data;
        matrix.ColumnName = [{'Analysis', 'Quantity'}, strrep(mx_cols', '#', 'Signal #')];
        matrix.UserData = cellfun(@(mx_r) str2double(erase(extractBefore(mx_r, '|'), '#')), mx_rows);
    end

    function mx_key = run_key(mx_a)
        % the key of an analysis of the run ('#2' gives Loudness_ISO532_1#2)
        mx_key = '';
        mx_k = find([run_settings.analyses.n] == str2double(erase(mx_a, '#')), 1);
        if ~isempty(mx_k)
            mx_key = run_settings.analyses(mx_k).key;
        end
    end

    function on_main_key(~, event)
        % the space bar plays and pauses, the arrows move the colour scale of the spectrogram
        if strcmp(event.Key, 'space') && ~isempty(loaded)
            toggle_play();
        else
            on_wave_key([], event);
        end
    end

    function on_stop_run(~, ~)
        % a metric cannot be interrupted, so the run ends at the next step
        if ~stop_requested
            stop_requested = true;
            write_log('Stop requested: the run ends when the metric being computed returns.');
            set_status('Stopping...');
        end
    end

    function on_theme(~, ~)
        if strcmp(theme_style, 'dark')
            theme_style = 'light';
        else
            theme_style = 'dark';
        end
        apply_theme();
    end

    function on_run(~, ~)
        if isempty(loaded)
            write_log('No files loaded. Open WAV files first.');
            return
        end
        use = find([loaded.marked]);
        if isempty(use)
            write_log('No signals ticked. Tick the signals to analyse.');
            return
        end
        if isempty(analyses)
            write_log('No metrics in the list of analyses. Add one with Add metric.');
            return
        end
        sel = {analyses.key};
        show = false;                            % the other figures are drawn on request
        k_active = find(use == active_idx, 1);   % the signal on screen is analysed first,
        if isempty(k_active)                     % and its figure is drawn in the analysis
            k_active = 1;                        % call, so that the metric runs once
        end
        files_order = [use(k_active), use(use ~= use(k_active))];
        active_path = loaded(files_order(1)).path;

        % the results of the last run whose signal and analysis did not change are kept
        kept = cell(1, numel(loaded));                     % entries kept, by signal
        keep_figs = {};                                    % file|key of the SQAT figures that stay
        for i = files_order
            kept{i} = il_empty_store();
            for a = analyses
                [old, same_key] = reusable(loaded(i), a);
                kept{i} = [kept{i}, old];
                if same_key
                    keep_figs{end+1} = [loaded(i).path '|' a.key]; %#ok<AGROW>
                end
            end
        end
        clear_cache(keep_figs);
        run_settings = struct('signals', loaded(use), 'analyses', analyses);
        close_params();

        new_results = il_empty_results();
        new_store = il_empty_store();
        stereo_sel = ismember({analyses.id}, {metrics([metrics.stereo]).id});
        n_total = 0;
        n_kept = 0;
        for i = files_order
            n_ch = numel(channel_list(loaded(i)));
            joint = n_ch == 2;
            todo = ~ismember(sel, {kept{i}.metric});
            n_total = n_total + nnz(~(stereo_sel & joint) & todo) * n_ch + nnz(stereo_sel & joint & todo);
            n_kept = n_kept + nnz(~todo);
        end
        n_done = 0;
        n_errors = 0;
        set_progress(0);
        stop_requested = false;
        dlg = [];
        if strcmp(fig.Visible, 'on')             % the dialog needs a visible window
            try
                dlg = uiprogressdlg(fig, 'Title', 'SQAT analysis', 'Message', 'Starting ...', ...
                    'Cancelable', 'on', 'CancelText', 'Stop', 'Value', 0);
                setappdata(fig, 'sqat_progress', dlg);   % the handle the tests read
            catch err
                write_log(['The progress dialog could not open: ' err.message]);
            end
        end
        dlg_off = onCleanup(@() close_progress()); %#ok<NASGU>
        t_start = tic;
        for i = files_order
            f = loaded(i);
            cl = channel_list(f);
            dBFS = f.dBFS;
            joint = numel(cl) == 2;              % a binaural pair goes in one call
            entries = kept{i};                   % of this file: the kept ones, then in the order of the calls
            todo = ~ismember(sel, {entries.metric});
            for key = sel(~todo)
                write_log(sprintf('%s on %s: kept from the last run, nothing changed.', key{1}, f.name));
            end
            ids_joint = sel(stereo_sel & joint & todo);
            ids_single = sel(~(stereo_sel & joint) & todo);
            for c = cl
                if isempty(ids_single) || stop_requested
                    break
                end
                try
                    [x, fs] = SQAT_GUI_load(f.path, dBFS, c);
                catch err
                    write_log(sprintf('ERROR reading %s: %s', f.name, err.message));
                    n_errors = n_errors + numel(ids_single);
                    n_done = n_done + numel(ids_single);
                    continue
                end
                first = c == cl(1);              % the figure of a signal is drawn once
                plan = share_plan(ids_single, numel(x), fs);
                done = struct();                 % outputs of this channel, by metric id
                for j = 1:numel(plan)
                    drawnow                      % the Stop button gets its turn here
                    poll_cancel();
                    if stop_requested
                        break
                    end
                    e = run_entry(plan(j).id);
                    set_status(sprintf('Running %s on %s (%d of %d)', e.label, f.name, n_done + 1, n_total));
                    try
                        [OUT, new_figs] = run_step(e, plan(j), done, x, fs, f, ...
                            show || (first && strcmp(f.path, active_path)));
                        if strcmp(e.key, e.id)       % a source for the metrics that share its computation
                            done.(e.id) = OUT;
                        end
                        entries(end+1) = il_entry_of(f, e, OUT, num2str(c), 1, 1); %#ok<AGROW>
                        if ~isempty(new_figs)
                            keep_figures(new_figs, f, e.key, first);
                        end
                    catch err
                        write_log(sprintf('ERROR in %s (%s): %s', e.key, f.name, err.message));
                        n_errors = n_errors + 1;
                    end
                    n_done = n_done + 1;
                    set_progress(100 * n_done / n_total);
                end
            end
            if ~isempty(ids_joint) && ~stop_requested
                X = [];
                try
                    [X, fs] = SQAT_GUI_load(f.path, dBFS, cl);
                catch err
                    write_log(sprintf('ERROR reading %s: %s', f.name, err.message));
                    n_errors = n_errors + numel(ids_joint);
                    n_done = n_done + numel(ids_joint);
                end
                for j = 1:numel(ids_joint)
                    drawnow
                    poll_cancel();
                    if isempty(X) || stop_requested
                        break
                    end
                    e = run_entry(ids_joint{j});
                    set_status(sprintf('Running %s on %s (%d of %d)', e.label, f.name, n_done + 1, n_total));
                    write_log(sprintf('Running %s on %s, both channels ...', e.key, f.name));
                    try
                        [OUT, new_figs] = run_metric(e, X, fs, e.p, ...
                            show || strcmp(f.path, active_path));
                        for label = {'1', '2', 'Binaural'}
                            c_out = label{1};
                            if ~strcmp(c_out, 'Binaural')
                                c_out = str2double(c_out);
                            end
                            entry = il_entry_of(f, e, OUT, label{1}, c_out, 2);
                            if ~isempty(entry.analyses) || ~isempty(entry.values)
                                entries(end+1) = entry; %#ok<AGROW>
                            end
                        end
                        if ~isempty(new_figs)
                            keep_figures(new_figs, f, e.key, true);
                        end
                    catch err
                        write_log(sprintf('ERROR in %s (%s): %s', e.key, f.name, err.message));
                        n_errors = n_errors + 1;
                    end
                    n_done = n_done + 1;
                    set_progress(100 * n_done / n_total);
                end
            end
            for j = 1:numel(sel)                 % the table keeps the order of the list
                k_e = find(strcmp({entries.metric}, sel{j}));
                new_store = [new_store, entries(k_e)]; %#ok<AGROW>
                rows = il_empty_results();
                a = run_settings.analyses(strcmp({run_settings.analyses.key}, sel{j}));
                e_j = metrics(strcmp({metrics.id}, a.id));
                for en = entries(k_e)
                    n = height(en.values);
                    unit = cellfun(@(q) il_unit(a.id, q), en.values.Quantity, 'UniformOutput', false);
                    rows = [rows; table(repmat({sprintf('#%d', f.id)}, n, 1), repmat({f.name}, n, 1), ...
                        repmat({en.number}, n, 1), repmat({a.id}, n, 1), repmat({en.channel}, n, 1), ...
                        en.values.Quantity, en.values.Value, unit(:), repmat(il_cal_of(f, en.channel), n, 1), ...
                        repmat({il_param_text(e_j, a.p)}, n, 1), repmat({f.path}, n, 1), ...
                        'VariableNames', il_empty_results().Properties.VariableNames)]; %#ok<AGROW>
                end
                % the channels interleaved: each quantity of channel 1, then of channel 2 (and binaural)
                [~, ~, q] = unique(rows.Quantity, 'stable');
                [~, order] = sortrows([q, (1:height(rows))']);
                new_results = [new_results; rows(order, :)]; %#ok<AGROW>
            end
            if stop_requested
                break                            % what ran so far is kept
            end
        end

        % the list goes analysis by analysis and quantity by quantity, with the
        % signals and their channels side by side (the tab of a signal keeps its own rows)
        a_num = str2double(erase(new_results.Analysis, '#'));
        [~, ~, q] = unique(strcat(new_results.Analysis, '|', new_results.Quantity), 'stable');
        [~, order] = sortrows([a_num, q, (1:height(new_results))']);
        results = new_results(order, :);
        store = new_store;
        if stop_requested                        % the settings list only what ran
            run_settings.signals = run_settings.signals(ismember({run_settings.signals.path}, {store.file}));
            run_settings.analyses = run_settings.analyses(ismember({run_settings.analyses.key}, {store.metric}));
        end
        tab_results.Title = 'Results';           % fresh results: a removal later keeps the mark of a change
        show_results();
        if stop_requested
            msg = sprintf('Stopped: %d value(s) from %d of the %d analysis step(s) in %.1f s', ...
                height(results), n_done, n_total, toc(t_start));
        else
            set_progress(100);
            msg = sprintf('Done: %d value(s) from %d file(s) and %d metric(s) in %.1f s', ...
                height(results), numel(use), numel(sel), toc(t_start));
        end
        if n_kept > 0
            msg = sprintf('%s, %d kept from the last run', msg, n_kept);
        end
        if n_errors > 0
            msg = sprintf('%s, %d error(s)', msg, n_errors);
        end
        lbl_status.Text = msg;
        write_log([msg '.']);
        refresh_graph_windows();
        if ~isempty(store)                       % the results show at once; the button opens them again
            if isempty(live_window())
                on_open_graphs();
            elseif strcmp(fig.Visible, 'on')
                figure(live_window());
            end
        end
        try
            focus(matrix);                       % off the Run button: a space would press it again
        catch
        end
    end

    function on_export(~, ~)
        if height(results) == 0
            write_log('No results to export. Run an analysis first.');
            return
        end
        [f, p] = uiputfile({'*.xlsx', 'Excel workbook (*.xlsx)'; '*.csv', 'CSV file (*.csv)'; ...
            '*.pdf', 'PDF report: settings, single values and plots (*.pdf)'}, ...
            'Export results', 'SQAT_results.xlsx');
        focus_gui();
        if isequal(f, 0)
            return
        end
        try
            if endsWith(f, '.pdf', 'IgnoreCase', true)
                write_report(fullfile(p, f));
                write_log(['Report written to ' fullfile(p, f)]);
                return
            end
            SQAT_GUI_export(results, fullfile(p, f), run_description());
            write_log(['Results exported to ' fullfile(p, f)]);
        catch err
            write_log(['ERROR exporting the results: ' err.message]);
        end
    end

    function on_session(what)
        if strcmp(what, 'save')
            [f, p] = uiputfile('*.mat', 'Save session', 'SQAT_session.mat');
        else
            [f, p] = uigetfile('*.mat', 'Open session');
        end
        focus_gui();
        if ~isequal(f, 0)
            session(what, fullfile(p, f));
        end
    end

    function session(what, file)
        % a session: the signals (path, channel, calibration, tick) and the
        % analyses with their parameters; results are not kept, a run makes them
        if strcmp(what, 'save')
            se = struct('signals', loaded, 'analyses', analyses, 'next_analysis', next_analysis); %#ok<NASGU>
            save(file, '-struct', 'se');
            write_log(['Session saved to ' file]);
            return
        end
        se = load(file);
        while ~isempty(loaded)
            remove_signal(1);                    % the signals of the session replace the list, results too
        end
        add_files({se.signals.path});
        for se_f = se.signals
            se_k = find(strcmp({loaded.path}, se_f.path), 1);
            if isempty(se_k)
                write_log(['ERROR: not found, left out of the session: ' se_f.path]);
                continue
            end
            loaded(se_k).channel = se_f.channel;
            loaded(se_k).marked = se_f.marked;
            loaded(se_k).dBFS = se_f.dBFS;
            loaded(se_k).cal_set = se_f.cal_set;
            loaded(se_k).cal = se_f.cal;
        end
        analyses = se.analyses;
        for se_k = 1:numel(analyses)                 % a parameter added since the session was saved takes its default
            se_d = il_default_params(metrics(strcmp({metrics.id}, analyses(se_k).id)));
            for se_n = setdiff(fieldnames(se_d), fieldnames(analyses(se_k).p))'
                analyses(se_k).p.(se_n{1}) = se_d.(se_n{1});
            end
        end
        next_analysis = se.next_analysis;
        assign_keys();
        refresh_signals();
        refresh_analyses();
        write_log(['Session opened from ' file]);
    end

    function write_report(file)
        % the report: the settings, the matrix of single values and, per analysis
        % of the run, its first plot (a time series, a profile or a map)
        rp_groups = struct('title', {}, 'analyses', {}, 'names', {});
        for rp_n = unique(matrix.UserData, 'stable')'
            rp_e = store(strcmp({store.number}, sprintf('#%d', rp_n)) & ismember({store.file}, {loaded.path}));
            if isempty(rp_e)
                continue
            end
            [~, rp_ids] = il_analysis_items(rp_e);
            rp_ids = rp_ids(~ismember(rp_ids, {'sqat', 'all', 'stats'}));
            if isempty(rp_ids)
                continue
            end
            rp_a = arrayfun(@(e) e.analyses(strcmp({e.analyses.id}, rp_ids{1})), rp_e);
            rp_names = arrayfun(@il_tag, rp_e, 'UniformOutput', false);
            rp_groups(end+1) = struct('title', sprintf('%s: %s', il_key_label(store, rp_e(1).metric), rp_a(1).label), ...
                'analyses', rp_a, 'names', {rp_names}); %#ok<AGROW>
        end
        SQAT_GUI_report(file, run_description(), matrix.ColumnName(:)', matrix.Data, rp_groups);
    end

    function S = run_description()
        % what the run used, for the Settings sheet of the export
        items = {'Exported', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
                 'SQAT version', il_sqat_version();
                 'MATLAB', version};
        for f = run_settings.signals
            items(end+1, :) = {sprintf('Signal #%d', f.id), sprintf('%s | fs %g Hz | %d channel(s) | channel %s | calibration %s%s | full scale %s dB SPL%s', ...
                f.path, f.fs, f.nch, f.channel, f.cal.label, il_if(isempty(f.cal.file), '', [' (' f.cal.file ')']), ...
                strjoin(arrayfun(@(v) sprintf('%.2f', v), f.dBFS, 'UniformOutput', false), ', '), ...
                il_if(f.cal_set, '', ' (default)'))}; %#ok<AGROW>
        end
        for a = run_settings.analyses
            e = metrics(strcmp({metrics.id}, a.id));
            items(end+1, :) = {sprintf('Analysis #%d', a.n), sprintf('%s | %s', a.id, il_param_text(e, a.p))}; %#ok<AGROW>
        end
        S = cell2table(items, 'VariableNames', {'Item', 'Value'});
    end

    function on_open_graphs(~, ~)
        if isempty(store)
            write_log('No results to plot. Run an analysis first.');
            return
        end
        % a window that follows the ticked signals is reused; a pinned one stays as it is
        w = live_window();
        if isempty(w)
            w = new_graph_window();
        end
        draw_window(w);
        if strcmp(fig.Visible, 'on')
            figure(w);
        end
    end

    function build_wave_view(parent)
        % the tabs of the signals, the player, the options of the spectrogram and
        % the plots, in one grid in the Waveform tab
        wave_view = uigridlayout(parent, [4 1]);
        wave_view.RowHeight = {26, 26, 26, '1x'};
        wave_view.Padding = [4 4 4 4];
        uitabgroup(wave_view, 'Tag', 'wave_tabs', 'SelectionChangedFcn', @on_wave_tab, ...
            'Visible', 'off');                          % a tab per signal, shown once there is one
        hw = uigridlayout(wave_view, [1 10]);
        hw.Padding = [0 0 0 0];
        hw.ColumnWidth = {70, 70, 55, 90, 95, 130, 70, 55, '1x', 100};   % narrow enough for the Waveform tab
        hw.ColumnSpacing = 6;
        btn_play = uibutton(hw, 'Text', 'Play', 'Tag', 'play', 'ButtonPushedFcn', @on_play, ...
            'Tooltip', 'Space plays and pauses; a click on the waveform or the spectrogram moves the playhead');
        uibutton(hw, 'Text', 'Stop', 'Tag', 'stop', 'ButtonPushedFcn', @on_stop);
        uicheckbox(hw, 'Text', 'Loop', 'Value', true, 'Tag', 'loop', 'ValueChangedFcn', @on_loop_changed, ...
            'Tooltip', 'Starts again at the end of the file, or of the filter box when there is one');
        uibutton(hw, 'state', 'Text', 'Draw filter', 'Tag', 'draw_box', 'ValueChangedFcn', @on_draw_box, ...
            'Tooltip', 'Drag a box on the spectrogram, or click two opposite corners');
        uibutton(hw, 'Text', 'Clear filters', 'Tag', 'clear_boxes', 'ButtonPushedFcn', @on_clear_boxes);
        uidropdown(hw, 'Items', {'Filter: loop only', 'Filter: isolate', 'Filter: remove'}, ...
            'ItemsData', {'loop', 'isolate', 'remove'}, 'Value', 'isolate', 'Tag', 'box_mode', ...
            'ValueChangedFcn', @on_processing_changed, ...
            'Tooltip', ['Loop only: the sound is not changed, and the loop runs inside the box. ' ...
                        'Isolate: only what is inside the box plays. Remove: what is inside the box is taken out.']);
        uilabel(hw, 'Text', 'Weighting:', 'HorizontalAlignment', 'right');
        uidropdown(hw, 'Items', {'Z', 'A', 'C'}, 'Value', 'Z', 'Tag', 'wave_weighting', ...
            'ValueChangedFcn', @on_weighting_changed, ...
            'Tooltip', 'Frequency weighting of the sound and of the spectrogram (IEC 61672-1)');
        uilabel(hw, 'Text', '');
        uibutton(hw, 'Text', 'Save plots...', 'Tag', 'save_wave_plots', 'ButtonPushedFcn', @on_save_wave_plots, ...
            'Tooltip', 'Saves the spectrogram, the waveform and the sound level, one file each, as PNG or PDF');
        sw = uigridlayout(wave_view, [1 11]);
        sw.Padding = [0 0 0 0];
        sw.ColumnWidth = {50, 115, 70, 72, 52, 78, 52, 62, 70, 85, '1x'};   % the switch needs room for Off and On
        sw.ColumnSpacing = 6;
        uilabel(sw, 'Text', 'Window:', 'HorizontalAlignment', 'right');
        uidropdown(sw, 'Items', {'Hann', 'Hamming', 'Rectangular', 'Blackman-Harris'}, ...
            'ItemsData', {'hann', 'hamming', 'rect', 'blackmanharris'}, 'Value', 'hann', ...
            'Tag', 'spec_window', 'ValueChangedFcn', @on_spec_option);
        uibutton(sw, 'Text', 'Import...', 'Tag', 'import_window', 'ButtonPushedFcn', @on_import_window, ...
            'Tooltip', 'A .txt, .csv, .dat or .mat file with the samples of the window');
        uilabel(sw, 'Text', 'FFT degree:', 'HorizontalAlignment', 'right');
        uispinner(sw, 'Value', 10, 'Limits', [6 16], 'Step', 1, 'Tag', 'spec_degree', ...
            'RoundFractionalValues', 'on', 'ValueChangedFcn', @on_spec_option, ...
            'Tooltip', 'The FFT has 2^degree points (6 to 16); use the arrows');
        uilabel(sw, 'Text', 'Overlap (%):', 'HorizontalAlignment', 'right');
        uispinner(sw, 'Value', 50, 'Limits', [0 95], 'Step', 5, 'Tag', 'spec_overlap', ...
            'ValueChangedFcn', @on_spec_option, 'Tooltip', 'Overlap of the frames (0 to 95); use the arrows');
        uilabel(sw, 'Text', 'Enhanced:', 'HorizontalAlignment', 'right', 'Tooltip', 'Enhanced STFT');
        uiswitch(sw, 'slider', 'Items', {'Off', 'On'}, 'Value', 'Off', 'Tag', 'spec_enhanced', ...
            'ValueChangedFcn', @on_spec_enhanced, ...
            'Tooltip', ['Enhanced STFT (consensus): nine windows of 8 to 512 ms, reassigned and combined, ' ...
                        'so that no window has to be chosen. It replaces the window, the FFT degree and the overlap.']);
        uidropdown(sw, 'Items', {'Readable', 'Sharp'}, 'ItemsData', {'readable', 'sharp'}, 'Value', 'readable', ...
            'Tag', 'spec_enhanced_mode', 'Enable', 'off', 'ValueChangedFcn', @on_spec_enhanced, ...
            'Tooltip', 'Readable: smoothing of 4 ms and 1.45 % of the frequency, continuous lines. Sharp: 1 ms and 1 Hz, the thinnest lines.');
        uilabel(sw, 'Text', '');
        pg = uigridlayout(wave_view, [4 1]);                  % the plots, each with its own tools above it
        pg.RowHeight = {'1.2x', '1x', 24, '1x'};       % the spectrogram on top, then the waveform and the sound level;
                                                        % the up and down arrows move the colour floor of the spectrogram
        pg.Padding = [0 0 0 0];
        pg.RowSpacing = 2;
        % each plot in a box of its own: the axes keep fixed margins for the
        % labels inside it, the same for all, so the time axes stay aligned
        box_wave = uipanel(pg, 'BorderType', 'none', 'AutoResizeChildren', 'off');
        ax_wave = uiaxes(box_wave, 'Tag', 'waveform_axes', 'ButtonDownFcn', @on_wave_click);
        box_wave.Layout.Row = 2;
        lg = uigridlayout(pg, [1 6]);                   % the sound level: its time weighting and indicators
        lg.Layout.Row = 3;
        lg.ColumnWidth = {100, 80, 100, 55, 55, '1x'};
        lg.ColumnSpacing = 6;
        lg.Padding = [0 0 0 0];
        uilabel(lg, 'Text', 'Time weighting:', 'HorizontalAlignment', 'right');
        uidropdown(lg, 'Items', {'Fast', 'Slow', 'Impulse'}, 'ItemsData', {'f', 's', 'i'}, 'Value', 'f', ...
            'Tag', 'level_time_weighting', 'ValueChangedFcn', @(~, ~) draw_level(), ...
            'Tooltip', 'Time weighting of the sound level (IEC 61672-1); the frequency weighting is the one of the player');
        uilabel(lg, 'Text', 'Exceeded (%):', 'HorizontalAlignment', 'right');
        tip = ['The level reached or exceeded during this percentage of the time (1 to 99), ' ...
               'as N5 in ISO 532-1; use the arrows'];
        uispinner(lg, 'Value', 5, 'Limits', [1 99], 'Step', 1, 'RoundFractionalValues', 'on', ...
            'Tag', 'level_percentile_1', 'ValueChangedFcn', @(~, ~) show_level_values(), 'Tooltip', tip);
        uispinner(lg, 'Value', 90, 'Limits', [1 99], 'Step', 1, 'RoundFractionalValues', 'on', ...
            'Tag', 'level_percentile_2', 'ValueChangedFcn', @(~, ~) show_level_values(), 'Tooltip', tip);
        uilabel(lg, 'Text', '', 'Tag', 'level_indicators', 'FontSize', 11, ...   % small enough for five values
            'Tooltip', ['Leq: equivalent level. LE: sound exposure level (SEL), Leq plus 10 lg of the duration in s. ' ...
                        'Lmax: maximum. LN: level reached or exceeded during N % of the time. All in dB.']);
        box_lvl = uipanel(pg, 'BorderType', 'none', 'AutoResizeChildren', 'off');
        ax_lvl = uiaxes(box_lvl, 'Tag', 'level_axes', 'ButtonDownFcn', @on_wave_click);
        box_lvl.Layout.Row = 4;
        box_spec = uipanel(pg, 'BorderType', 'none', 'AutoResizeChildren', 'off');
        box_spec.Layout.Row = 1;
        ax_spec = uiaxes(box_spec, 'Tag', 'spectrogram', 'ButtonDownFcn', @on_wave_click);
        box_wave.SizeChangedFcn = @(src, ~) il_fit_axes(src, ax_wave, 26);
        box_spec.SizeChangedFcn = @(src, ~) il_fit_axes(src, ax_spec, 26);
        box_lvl.SizeChangedFcn = @(src, ~) il_fit_axes(src, ax_lvl);
        il_fit_axes(box_wave, ax_wave, 26);
        il_fit_axes(box_lvl, ax_lvl);
        il_fit_axes(box_spec, ax_spec, 26);
        ax_spec.XAxis.LimitsChangedFcn = @on_spec_limits;
        ax_wave.XAxis.LimitsChangedFcn = @on_wave_limits;
        ax_lvl.XAxis.LimitsChangedFcn = @on_level_limits;
        wave_appdata();
    end

    function on_save_wave_plots(~, ~)
        % the three plots of the Waveform tab as they are on screen (zoom included),
        % one file each: <name>_spectrogram, <name>_waveform and <name>_level
        if isempty(wave_x)
            write_log('No file loaded.');
            return
        end
        f = active_file();
        [~, base] = fileparts(f.name);
        base = sprintf('%s_ch%d', base, channel_of_active());
        if isappdata(fig, 'sqat_next_file')            % a test stands in for the file dialog
            path = getappdata(fig, 'sqat_next_file');
            rmappdata(fig, 'sqat_next_file');
        else
            [name, folder, k] = uiputfile({'*.png', 'PNG image (*.png)'; '*.pdf', 'PDF (*.pdf)'}, ...
                'Save the plots', fullfile(save_folder, [base '.png']));
            focus_gui();
            if isequal(name, 0)
                return
            end
            [~, ~, ext] = fileparts(name);
            if isempty(ext)
                name = [name il_if(k == 2, '.pdf', '.png')];
            end
            path = fullfile(folder, name);
        end
        [folder, base, ext] = fileparts(path);
        save_folder = folder;
        n = 0;
        names = {'spectrogram', 'waveform', 'level'};
        axs = [ax_spec, ax_wave, ax_lvl];
        bottom = [26 26 48];                           % the margin under each plot on screen
        for k = 1:3
            had = axs(k).XLabel.String;
            xlabel(axs(k), 'Time (s)');                % a file stands alone: each carries its time label
            il_fit_axes(axs(k).Parent, axs(k));        % with the room for it
            n = n + il_export(axs(k), folder, sprintf('%s_%s%s', base, names{k}, ext));
            xlabel(axs(k), had);
            il_fit_axes(axs(k).Parent, axs(k), bottom(k));
        end
        write_log(sprintf('%d plot(s) saved to %s', n, folder));
    end

    function wave_appdata()
        setappdata(win_wave, 'sqat_spec_zoom', @apply_spec_zoom);   % the recomputation, for the tests
        setappdata(win_wave, 'sqat_audio', @process_audio);   % what plays, for the tests
        setappdata(win_wave, 'sqat_play', @play_info);
    end

    function on_wave_tab(~, event)
        % a tab makes its signal and channel the ones of the window, as the name
        % of the signal in the list does
        u = event.NewValue.UserData;                    % [signal number, channel]
        k = find([loaded.id] == u(1), 1);
        if isempty(k)
            return
        end
        wave_ch_by_id(u(1)) = u(2);
        if k == active_idx
            show_active();                              % the same signal, another channel
        else
            on_signal_name(k);
        end
    end

    function sync_wave_tabs()
        % one tab per signal of the list, and per channel of a signal with several;
        % the tab on screen selected
        tg = findobj(win_wave, 'Tag', 'wave_tabs');
        titles = {};
        data = zeros(0, 2);
        for k = 1:numel(loaded)
            f = loaded(k);
            for c = 1:f.nch
                titles{end+1} = sprintf('#%d %s', f.id, f.name); %#ok<AGROW>
                if f.nch > 1
                    titles{end} = sprintf('%s ch%d', titles{end}, c);
                end
                data(end+1, :) = [f.id c]; %#ok<AGROW>
            end
        end
        if ~isequal(arrayfun(@(t) t.Title, tg.Children(:)', 'UniformOutput', false), titles)
            delete(tg.Children);
            for k = 1:numel(titles)
                uitab(tg, 'Title', titles{k}, 'UserData', data(k, :));
            end
            SQAT_GUI_paint(tg, theme_style);
        end
        tg.Visible = ~isempty(titles);                   % no empty grey strip before a file is loaded
        if active_idx > 0
            tg.SelectedTab = tg.Children(ismember(data, [active_file().id channel_of_active()], 'rows'));
        end
    end

    function on_wave_key(~, event)
        switch event.Key
            case 'space'
                toggle_play();
            case 'uparrow'
                shift_black(5);
            case 'downarrow'
                shift_black(-5);
        end
    end

    function on_play(~, ~)
        toggle_play();
        if il_is_open(win_wave) && strcmp(win_wave.Visible, 'on')
            figure(win_wave);                % the key events of the window need its focus
        end
    end

    function toggle_play()
        if isempty(loaded) || ~il_is_open(win_wave)
            return
        end
        if toc(toggle_clock) - last_toggle < 0.15   % the same key event reaching the window and a focused button
            return
        end
        if is_playing()
            play_start = file_position(player.CurrentSample);
            stop_player();
            btn_play.Text = 'Play';
            write_log('Paused.');
        else
            start_playback(play_start, false);
        end
        last_toggle = toc(toggle_clock);          % counted from the end: stopping the sound takes a while
    end

    function tf = is_playing()
        tf = playing && ~isempty(player) && isvalid(player);
    end

    function stop_player()
        % stops the sound; the state goes first, so that the stop is not taken for the end of the buffer
        was = playing;
        playing = false;
        if was && ~isempty(player) && isvalid(player)
            stop(player);
        end
    end

    function y = play_audio()
        % the sound that plays: the file with the filter and the weighting, kept until they change
        if isempty(play_y)
            y = process_audio();
            peak = max(abs(y), [], 'all');
            if peak > 1                      % a weighting can lift the level above full scale
                y = y / peak;
                write_log(sprintf('The sound was lowered by %.1f dB so that it does not clip.', 20*log10(peak)));
            end
            play_y = single(y);
        end
        y = play_y;
    end

    function start_playback(sample, quiet)
        % plays from a sample of the file (0 for the start of the box, or of the file). With the
        % loop on, the buffer holds the rest of the play and then the loop repeated for about two
        % minutes, so that the loop has no restart, and so no gap, between two turns
        f = active_file();
        try
            y = play_audio();
            n = numel(y);
            [r1, r2] = play_region();
            if sample == 0
                sample = r1;
            end
            sample = min(max(round(sample), 1), n);
            looping = findobj(win_wave, 'Tag', 'loop').Value;
            % inside the box the box is the loop; with the loop on, a play from
            % before the box runs into it and stays there
            box_loop = ~isempty(boxes) && sample <= r2 && (sample >= r1 || looping);
            after_box = ~isempty(boxes) && sample > r2 && looping;
            if box_loop
                last = r2;
                rep = y(r1:r2);
                loop_from = r1;
            else
                last = n;
                rep = y;
                loop_from = 1;
            end
            first = y(sample:last);
            if looping && ~after_box
                buf = [first; repmat(rep, max(1, ceil(120 * wave_fs / numel(rep))), 1)];
            else                      % past the box: to the end, then from the start (on_player_stopped) into the box
                buf = first;
                rep = zeros(0, 1, 'single');
            end
            play_map = struct('sample', sample, 'n_first', numel(first), 'loop_from', loop_from, ...
                'n_rep', numel(rep));
            if isappdata(groot, 'sqat_gui_mute')       % the tests play in silence
                buf = zeros(size(buf), 'like', buf);
            end
            play_peak = max(abs(buf), [], 'all');
            player = audioplayer(buf, wave_fs);
            player_fs = wave_fs;
            player.TimerPeriod = 0.05;
            player.TimerFcn = @(~, ~) on_player_tick();
            player.StopFcn = @(~, ~) on_player_stopped();
            play(player);
            playing = true;
            play_t0 = tic;
            play_span = numel(buf) / wave_fs;
            btn_play.Text = 'Pause';
            move_playhead(sample);            % now, not at the first tick: the device may start seconds late
            if ~quiet
                write_log(sprintf('Playing %s, channel %d.', f.name, channel_of_active()));
            end
        catch err
            player = [];
            playing = false;
            write_log(['Audio output unavailable: ' err.message]);
        end
    end

    function pos = file_position(i)
        % the sample of the file that sample i of the buffer being played holds
        i = max(i, 1);
        m = play_map;
        if i <= m.n_first || m.n_rep == 0
            pos = m.sample + min(i, m.n_first) - 1;
        else
            pos = m.loop_from + mod(i - m.n_first - 1, m.n_rep);
        end
    end

    function info = play_info()
        % the state of the play, for the tests
        info = struct('playing', playing, 'box_loop', box_loop, 'sample', play_map.sample, ...
            'buffer_peak', play_peak);
    end

    function [r1, r2] = play_region()
        % the samples of the filter boxes' time, or of the whole file
        n = numel(wave_y);
        r1 = 1;
        r2 = n;
        if ~isempty(boxes)
            r1 = min(max(floor(min(boxes(:, 1)) * wave_fs) + 1, 1), n);
            r2 = min(max(ceil(max(boxes(:, 2)) * wave_fs), r1), n);
        end
    end

    function on_stop(~, ~)
        if ~isempty(player)
            stop_player();
            write_log('Stopped.');
        end
        play_start = 0;
        move_playhead(1);
        if il_is_open(win_wave)
            btn_play.Text = 'Play';
        end
    end

    function on_loop_changed(~, ~)
        if is_playing()
            rebuild_audio(false);            % the buffer holds the loop, or not
        end
    end

    function on_player_tick()
        if ~playing || isempty(player) || ~isvalid(player)
            return
        end
        move_playhead(file_position(player.CurrentSample));
    end

    function on_player_stopped()
        % called on every stop. The end of the buffer is a stop that comes while the interface
        % thinks it plays and after the time the play was to last; a stop the interface asked for
        % comes with playing false, or too soon (an old one, after a jump)
        if ~playing || toc(play_t0) < 0.9 * play_span
            return
        end
        if il_is_open(win_wave) && findobj(win_wave, 'Tag', 'loop').Value
            start_playback(play_map.loop_from, true);     % after two minutes of loop: once more
            return
        end
        playing = false;
        play_start = 0;
        if il_is_open(win_wave)
            btn_play.Text = 'Play';
            move_playhead(1);
        end
    end

    function on_wave_click(src, event)
        pt = event.IntersectionPoint;
        if src == ax_spec && findobj(win_wave, 'Tag', 'draw_box').Value
            on_box_click(pt);
        else
            seek(pt(1));
        end
    end

    function seek(t)
        % moves the playhead, and the sound with it when it plays
        if isempty(wave_y)
            return
        end
        sample = min(max(round(t * wave_fs) + 1, 1), numel(wave_y));
        was_playing = is_playing();
        play_start = sample;
        if was_playing
            stop_player();
            start_playback(sample, true);
        end
        move_playhead(sample);
    end

    function rebuild_audio(jump)
        % what plays changed: the sound is made again, from the same position, or from the
        % start of the box when jump is true
        if nargin < 1
            jump = false;
        end
        was_playing = is_playing();
        if was_playing
            play_start = file_position(player.CurrentSample);
            stop_player();
        end
        play_y = [];
        if jump
            play_start = 0;
        end
        if was_playing
            start_playback(play_start, true);
        end
    end

    function y = process_audio()
        % the file with the boxes isolated or removed, and the frequency weighting applied
        y = wave_y;
        mode = findobj(win_wave, 'Tag', 'box_mode').Value;
        if ~isempty(boxes) && ~strcmp(mode, 'loop')
            if strcmp(mode, 'isolate')
                y = SQAT_GUI_spectral_filter(y, wave_fs, boxes, 'keep');
            else
                y = SQAT_GUI_spectral_filter(y, wave_fs, boxes, 'remove');
            end
        end
        y = SQAT_GUI_weight(y, wave_fs, findobj(win_wave, 'Tag', 'wave_weighting').Value);
    end

    function on_processing_changed(~, ~)
        rebuild_audio();
    end

    function on_weighting_changed(~, ~)
        rebuild_audio();
        draw_level();
        draw_spectrogram(true);
    end

    function on_draw_box(src, ~)
        if src.Value                                  % a zoom, pan or data tip mode of the axes
            zoom(win_wave, 'off');                    % toolbar would take the clicks
            pan(win_wave, 'off');
            datacursormode(win_wave, 'off');
        end
        end_drag();
        box_corner = [];
        delete(findobj(ax_spec, 'Tag', 'box_corner'));
        if src.Value
            write_log('Draw filter: drag a box on the spectrogram, or click two opposite corners.');
        end
    end

    function shift_black(d)
        % + moves the bottom of the colour scale up (more black), - down (less black)
        spec_black = min(max(spec_black + d, -60), spec_range - 5);
        apply_black();
    end

    function apply_black()
        if ~isempty(spec_top) && il_is_open(win_wave)
            clim(ax_spec, spec_top + [spec_black - spec_range, 0]);
        end
    end

    function on_box_click(pt)
        if isempty(box_corner)
            box_corner = pt(1:2);
            line(ax_spec, pt(1), pt(2), 1, 'Marker', '+', 'MarkerSize', 12, 'Color', [1 1 1], ...
                'Tag', 'box_corner', 'PickableParts', 'none');
            drag_active = true;              % the mouse may now be dragged to the other corner
            win_wave.WindowButtonMotionFcn = @on_box_motion;
            win_wave.WindowButtonUpFcn = @on_box_release;
            return
        end
        add_box(box_corner, pt(1:2));
    end

    function on_box_motion(~, event)
        if ~drag_active || isempty(box_corner)
            return
        end
        pt = pointer_point(event);
        if isempty(pt)
            return
        end
        delete(findobj(ax_spec, 'Tag', 'box_preview'));
        patch(ax_spec, 'XData', [box_corner(1) pt(1) pt(1) box_corner(1)], ...
            'YData', [box_corner(2) box_corner(2) pt(2) pt(2)], 'ZData', 0.5 * ones(1, 4), ...
            'FaceColor', [1 1 1], 'FaceAlpha', 0.1, 'EdgeColor', [1 1 1], 'LineStyle', ':', ...
            'Tag', 'box_preview', 'PickableParts', 'none');
        drawnow limitrate
    end

    function on_box_release(~, event)
        if ~drag_active
            return
        end
        pt = pointer_point(event);
        end_drag();
        if isempty(pt) || isempty(box_corner)
            return
        end
        d_time = abs(pt(1) - box_corner(1)) / (numel(wave_y) / wave_fs);
        d_freq = abs(log10(max(pt(2), 1) / max(box_corner(2), 1)));
        if d_time > 0.01 || d_freq > 0.05
            add_box(box_corner, pt);         % dragged: the box is done
        end                                  % pressed and released on one spot: the first click
    end

    function end_drag()
        drag_active = false;
        if il_is_open(win_wave)
            win_wave.WindowButtonMotionFcn = '';
            win_wave.WindowButtonUpFcn = '';
            delete(findobj(ax_spec, 'Tag', 'box_preview'));
        end
    end

    function pt = pointer_point(event)
        % where the mouse is on the spectrogram (time, frequency), or [] when it is not known
        pt = [];
        has_point = (isstruct(event) && isfield(event, 'IntersectionPoint')) || ...
            (isobject(event) && isprop(event, 'IntersectionPoint'));
        on_axes = true;                       % the point of an event holds for the object under the mouse
        if isobject(event) && isprop(event, 'HitObject') && ~isempty(event.HitObject)
            on_axes = isequal(event.HitObject, ax_spec);
        end
        if has_point && on_axes && numel(event.IntersectionPoint) >= 2 && all(isfinite(event.IntersectionPoint(1:2)))
            pt = event.IntersectionPoint(1:2);
        else
            cp = ax_spec.CurrentPoint;
            if all(isfinite(cp(1, 1:2)))
                pt = cp(1, 1:2);
            end
        end
    end

    function add_box(c1, c2)
        % the box of two corners, kept inside the spectrogram
        end_drag();
        box_corner = [];
        delete(findobj(ax_spec, 'Tag', 'box_corner'));
        t_lim = [0 numel(wave_y) / wave_fs];
        f_lim = [20 wave_fs / 2];
        c = [c1(1:2); c2(1:2)];
        t12 = sort(min(max(c(:, 1), t_lim(1)), t_lim(2)))';
        f12 = sort(min(max(c(:, 2), f_lim(1)), f_lim(2)))';
        boxes(end+1, :) = [t12 f12];
        set(findobj(win_wave, 'Tag', 'draw_box'), 'Value', false);
        draw_boxes();
        rebuild_audio(findobj(win_wave, 'Tag', 'loop').Value);   % with the loop on, the box takes the play
        write_log(sprintf('Filter box added: %.2f to %.2f s, %.0f to %.0f Hz.', t12, f12));
    end

    function on_clear_boxes(~, ~)
        boxes = zeros(0, 4);
        box_corner = [];
        delete(findobj(ax_spec, 'Tag', 'box_corner'));
        draw_boxes();
        rebuild_audio();
        write_log('Filters cleared.');
    end

    function draw_boxes()
        delete(findobj(ax_spec, 'Tag', 'box'));
        for b = 1:size(boxes, 1)
            patch(ax_spec, 'XData', [boxes(b, 1) boxes(b, 2) boxes(b, 2) boxes(b, 1)], ...
                'YData', [boxes(b, 3) boxes(b, 3) boxes(b, 4) boxes(b, 4)], 'ZData', 0.5 * ones(1, 4), ...
                'FaceColor', [1 1 1], 'FaceAlpha', 0.15, 'EdgeColor', [1 1 1], 'LineStyle', '--', ...
                'Tag', 'box', 'PickableParts', 'none');
        end
    end

    function on_spec_option(~, ~)
        d = findobj(win_wave, 'Tag', 'spec_degree');
        d.Value = round(min(max(d.Value, 6), 16));
        o = findobj(win_wave, 'Tag', 'spec_overlap');
        o.Value = min(max(o.Value, 0), 95);
        draw_spectrogram(true);
    end

    function on_spec_enhanced(~, ~)
        % the enhanced STFT sets its own windows: the controls of the plain spectrogram are faded
        on = strcmp(findobj(win_wave, 'Tag', 'spec_enhanced').Value, 'On');
        state = matlab.lang.OnOffSwitchState(~on);
        set(findobj(win_wave, 'Tag', 'spec_window'), 'Enable', state);
        set(findobj(win_wave, 'Tag', 'import_window'), 'Enable', state);
        set(findobj(win_wave, 'Tag', 'spec_degree'), 'Enable', state);
        set(findobj(win_wave, 'Tag', 'spec_overlap'), 'Enable', state);
        set(findobj(win_wave, 'Tag', 'spec_enhanced_mode'), 'Enable', matlab.lang.OnOffSwitchState(on));
        draw_spectrogram(true);
    end

    function on_import_window(~, ~)
        if isappdata(win_wave, 'sqat_next_file')      % a test stands in for the file dialog
            path = getappdata(win_wave, 'sqat_next_file');
            rmappdata(win_wave, 'sqat_next_file');
        else
            [name, folder] = uigetfile({'*.txt;*.csv;*.dat;*.mat', 'Window files'}, 'Import a window');
            focus_gui();
            if isequal(name, 0)
                return
            end
            path = fullfile(folder, name);
        end
        try
            spec_custom = SQAT_GUI_read_window(path);
        catch err
            write_log(['ERROR importing a window: ' err.message]);
            return
        end
        [~, base, ext] = fileparts(path);
        dd = findobj(win_wave, 'Tag', 'spec_window');
        keep = ~strcmp(dd.ItemsData, 'custom');
        set(dd, 'Items', [dd.Items(keep), {['Custom: ' base ext]}], 'ItemsData', [dd.ItemsData(keep), {'custom'}]);
        dd.Value = 'custom';
        write_log(sprintf('Window %s%s loaded: %d samples, resampled to the FFT size.', base, ext, numel(spec_custom)));
        draw_spectrogram(true);
    end

    function on_close(~, ~)
        if ~isempty(player)
            stop(player);
        end
        clear_cache();
        stop_spec_timer();
        il_delete_timer(prefetch_timer);
        cancel_jobs();
        delete(open_graph_windows());
        delete(fig);
    end

%% Helpers sharing the state

    function focus_gui()
        % a file dialog hands the focus to the MATLAB desktop; take it back
        if strcmp(fig.Visible, 'on')
            figure(fig);
        end
    end

    function add_files(paths)
        % every new file starts at 94 dBFS; its calibration button changes it
        n_before = numel(loaded);
        for k_file = 1:numel(paths)
            path = char(paths{k_file});
            if any(strcmp({loaded.path}, path))
                continue
            end
            try
                info = audioinfo(path);
            catch err
                write_log(sprintf('ERROR opening %s: %s', path, err.message));
                continue
            end
            [~, base, ext] = fileparts(path);
            loaded(end+1) = struct('path', path, 'name', [base ext], ...
                'nch', info.NumChannels, 'fs', info.SampleRate, 'marked', true, 'id', next_id, ...
                'channel', il_default_channel(info.NumChannels), 'dBFS', 94 * ones(1, info.NumChannels), ...
                'cal_set', false, 'cal', struct('method', 'dbfs', 'level', 94, 'file', '', 'label', '94 dBFS')); %#ok<AGROW>
            next_id = next_id + 1;
        end
        if isempty(loaded)
            return
        end
        if active_idx == 0 || active_idx > numel(loaded)
            active_idx = 1;
        end
        refresh_signals();
        write_log(sprintf('%d file(s) loaded.', numel(loaded)));
        if numel(loaded) > n_before
            refresh_windows();
        end
    end

    function remove_signal(k)
        % the signal leaves the list, the results and the drawn figures
        f = loaded(k);
        loaded(k) = [];
        keep = ~strcmp({store.file}, f.path);
        store = store(keep);
        results = results(~strcmp(results.Path, f.path), :);
        run_settings.signals(strcmp({run_settings.signals.path}, f.path)) = [];   % out of the exported settings
        show_results();
        drop = strcmp({cache.file}, f.path);
        for k_c = find(drop)
            delete(cache(k_c).figs(isvalid(cache(k_c).figs)));
        end
        cache = cache(~drop);
        if k < active_idx
            active_idx = active_idx - 1;
        end
        active_idx = min(active_idx, numel(loaded));   % 0 when the list is empty
        write_log(sprintf('%s removed.', f.name));
        refresh_signals();
        refresh_windows();
    end

    function refresh_signals()
        % one row per signal: number, tick, name (click: player), channel, calibration and bin
        delete(signal_list.Children);
        n = numel(loaded);
        signal_list.RowHeight = repmat({24}, 1, n + 1);
        heads = {'', '', 'Signal', 'Channel', ['Calibration ' char(9432)], ''};
        for c = 1:6
            h = uilabel(signal_list, 'Text', heads{c}, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');  % centred over the boxes below
            if c == 5
                set(h, 'Tag', 'signal_cal_info', 'Tooltip', cal_help);
            end
        end
        for k = 1:n
            f = loaded(k);
            uilabel(signal_list, 'Text', sprintf('#%d', f.id), 'Tag', sprintf('signal_number_%d', k));
            uicheckbox(signal_list, 'Text', '', 'Value', f.marked, 'Tag', sprintf('signal_tick_%d', k), ...
                'Tooltip', 'Use this signal in the run and the plots', ...
                'ValueChangedFcn', @(src, ~) on_signal_ticked(k, src.Value));
            uibutton(signal_list, 'Text', f.name, 'HorizontalAlignment', 'left', ...
                'Tag', sprintf('signal_name_%d', k), 'Tooltip', [f.path newline 'Click: show it in the Waveform tab'], ...
                'ButtonPushedFcn', @(~, ~) on_signal_name(k));
            uidropdown(signal_list, 'Items', il_channel_items(f.nch), 'Value', f.channel, ...
                'Tag', sprintf('signal_channel_%d', k), ...
                'Tooltip', ['Channel to analyse, or All (the ECMA-418-2 metrics take a stereo ' ...
                            'file as a binaural pair, in one call)'], ...
                'ValueChangedFcn', @(src, ~) on_signal_channel(k, src.Value));
            angle = 'italic';
            if f.cal_set
                angle = 'normal';
            end
            uibutton(signal_list, 'Text', f.cal.label, 'Tag', sprintf('signal_cal_%d', k), ...
                'FontAngle', angle, ...
                'Tooltip', sprintf(['Full scale: %s dB SPL (%s). Click to change the calibration.' newline ...
                    'In italics while it is the default of SQAT.'], ...
                    strjoin(arrayfun(@(v) sprintf('%.2f', v), f.dBFS, 'UniformOutput', false), ', '), ...
                    il_if(f.nch > 1, 'one value per channel', 'one channel')), ...
                'ButtonPushedFcn', @(~, ~) on_signal_cal(k));
            uibutton(signal_list, 'Text', '', 'Icon', icon_remove, 'Tag', sprintf('signal_remove_%d', k), ...
                'Tooltip', 'Removes this signal and its results', 'ButtonPushedFcn', @(~, ~) on_signal_remove(k));
        end
        if n == 0
            lbl_files.Text = 'No files loaded';
            h = uilabel(signal_list, 'Text', 'Open WAV files to start: the button above, top right.', ...
                'FontAngle', 'italic', 'Tag', 'signals_hint');
            h.Layout.Row = 2;
            h.Layout.Column = [3 6];
        elseif n == 1
            lbl_files.Text = '1 file loaded';
        else
            lbl_files.Text = sprintf('%d files loaded', n);
        end
        mark_active();
        SQAT_GUI_paint(signal_list, theme_style);      % the new rows in the colours of the theme
        update_run_label();
    end

    function update_run_label()
        % the button tells what a run computes: the ticked signals times the analyses
        n_s = nnz([loaded.marked]);
        n_a = numel(analyses);
        if n_s == 0 || n_a == 0
            btn_run.Text = 'Run Analysis';
        else
            btn_run.Text = sprintf('Run %d %s %c %d %s', n_s, il_if(n_s == 1, 'signal', 'signals'), 215, ...
                n_a, il_if(n_a == 1, 'analysis', 'analyses'));
        end
    end

    function mark_active()
        % the name of the signal on screen is shown in bold
        for k = 1:numel(loaded)
            b = findobj(signal_list, 'Tag', sprintf('signal_name_%d', k));
            b.FontWeight = 'normal';
            if k == active_idx
                b.FontWeight = 'bold';
            end
        end
    end

    function show_active()
        mark_active();
        if il_is_open(win_wave)
            draw_waveform_window();
        end
    end

    function f = active_file()
        f = loaded(active_idx);
    end

    function ch = channel_of_active()
        % the channel that the waveform and the player take: the one of its tab, or
        % else the channel of the signal in the list
        id = active_file().id;
        if numel(wave_ch_by_id) >= id && wave_ch_by_id(id) > 0
            ch = wave_ch_by_id(id);
        else
            ch = str2double(active_file().channel);
        end
        if isnan(ch)
            ch = 1;                          % All
        end
        ch = min(ch, active_file().nch);
    end

    function refresh_windows()
        refresh_graph_windows();
        if il_is_open(win_wave)
            if isempty(loaded)
                sync_wave_tabs();
                cla(ax_wave);
                cla(ax_spec);
                wave_x = [];
                draw_level();
            else
                draw_waveform_window();
            end
        end
    end

    %% Graphs windows: signal, metric, analysis, plot

    function ws = open_graph_windows()
        graph_figs = graph_figs(isvalid(graph_figs));
        ws = graph_figs;
    end

    function w = live_window()
        % the first window that follows the ticked signals
        w = [];
        for w_k = open_graph_windows()
            if ~findobj(w_k, 'Tag', 'graph_pin').Value
                w = w_k;
                return
            end
        end
    end

    function refresh_graph_windows()
        for w = open_graph_windows()
            if ~findobj(w, 'Tag', 'graph_pin').Value
                draw_window(w);
            end
        end
    end

    function w = new_graph_window()
        n = numel(open_graph_windows());
        w = uifigure('Name', 'SQAT graphs', 'Position', [160 160 1100 700] + [30 -30 0 0] * mod(n, 6), ...
            'Visible', fig.Visible, 'Tag', 'SQAT_GUI_graphs', 'CreateFcn', '');
        graph_figs(end+1) = w;
        g = uigridlayout(w, [3 1]);
        g.RowHeight = {30, 0, '1x'};                  % the channel row opens when a signal has several
        g.Padding = [8 8 8 8];
        bar = uigridlayout(g, [1 7]);
        bar.Padding = [0 0 0 0];
        bar.ColumnWidth = {50, 230, 60, 260, 70, 80, '1x'};
        uilabel(bar, 'Text', 'Metric:', 'HorizontalAlignment', 'right');
        uidropdown(bar, 'Items', {}, 'Tag', 'graph_metric', 'ValueChangedFcn', @on_graph_control);
        uilabel(bar, 'Text', 'Analysis:', 'HorizontalAlignment', 'right');
        uidropdown(bar, 'Items', {}, 'Tag', 'graph_analysis', 'ValueChangedFcn', @on_graph_control, ...
            'Tooltip', ['SQAT figure and All analyses need one signal; with several signals ' ...
                        'the analyses that can be compared are on offer']);
        uibutton(bar, 'state', 'Text', 'Pin', 'Tag', 'graph_pin', 'ValueChangedFcn', @on_graph_pin, ...
            'Tooltip', ['Keeps this window with its signals and results, whatever runs next; ' ...
                        'Open Graphs Window then opens another one to compare with']);
        uibutton(bar, 'Text', 'Save...', 'Tag', 'graph_save', 'ButtonPushedFcn', @(src, ~) open_save_dialog(ancestor(src, 'figure')), ...
            'Tooltip', 'Save what the window shows, or choose other signals and figures');
        uilabel(bar, 'Text', '');
        row = uigridlayout(g, [1 1], 'Tag', 'graph_channels');   % one channel choice per signal
        row.Padding = [0 0 0 0];
        uipanel(g, 'BorderType', 'none', 'Tag', 'graph_body');
        apply_theme();
    end

    function on_graph_control(src, ~)
        draw_window(ancestor(src, 'figure'));
    end

    function on_graph_pin(src, ~)
        w = ancestor(src, 'figure');
        if src.Value
            paths = {loaded([loaded.marked]).path};
            % the signals and the results of the window stay these, whatever runs next
            src.UserData = struct('paths', {paths}, 'store', store(ismember({store.file}, paths)));
            src.Text = 'Pinned';
        else
            src.UserData = [];
            src.Text = 'Pin';
        end
        draw_window(w);
    end

    function paths = window_paths(w)
        % the signals of a window: the ticked ones, or the ones it was pinned with
        pin = findobj(w, 'Tag', 'graph_pin');
        if pin.Value
            paths = pin.UserData.paths;
        else
            paths = {loaded([loaded.marked]).path};
        end
    end

    function s = window_store(w)
        % the results of a window: the last run, or the ones it was pinned with
        pin = findobj(w, 'Tag', 'graph_pin');
        if pin.Value
            s = pin.UserData.store;
        else
            s = store;
        end
    end

    function draw_window(w)
        % refreshes the selectors of a window and draws what they ask for
        paths = window_paths(w);
        ws = window_store(w);
        pinned = findobj(w, 'Tag', 'graph_pin').Value;
        dm = findobj(w, 'Tag', 'graph_metric');
        da = findobj(w, 'Tag', 'graph_analysis');
        body = findobj(w, 'Tag', 'graph_body');
        delete(body.Children);
        in_window = ismember({ws.file}, paths);
        ids = unique({ws(in_window).metric}, 'stable');   % the keys of the analyses, in the order of the run
        if isempty(ids)
            set([dm, da], 'Items', {});
            channel_row(w, {}, {});
            uilabel(uigridlayout(body, [1 1]), 'HorizontalAlignment', 'center', ...
                'Text', 'No results for the signals of this window. Tick signals and run an analysis.');
            w.Name = 'SQAT graphs';
            return
        end
        wanted = dm.Value;
        if ~il_is_member(wanted, ids)
            wanted = last_metric;
        end
        set(dm, 'Items', cellfun(@(k) il_key_label(ws, k), ids, 'UniformOutput', false), 'ItemsData', ids);
        if il_is_member(wanted, ids)
            dm.Value = wanted;
        end
        id = dm.Value;
        last_metric = id;
        label = il_key_label(ws, id);

        in_metric = in_window & strcmp({ws.metric}, id);
        chans = channel_row(w, paths, ws(in_metric));
        entries = il_empty_store();
        for k_p = 1:numel(paths)
            k_e = find(in_metric & strcmp({ws.file}, paths{k_p}));
            if ~isempty(k_e) && ~strcmp(chans{k_p}, 'All')
                k_e = k_e(strcmp({ws(k_e).channel}, chans{k_p}));
            end
            entries = [entries, ws(k_e)]; %#ok<AGROW>
        end
        chan = unique({entries.channel});
        if isscalar(chan)
            chan = chan{1};                   % the channel of every line, for the titles
        else
            chan = '';                        % mixed channels: the legend names them
        end

        [items, data] = il_analysis_items(entries);
        wanted = da.Value;
        set(da, 'Items', items, 'ItemsData', data);
        if ~il_is_member(wanted, data)
            wanted = data{1};                 % the SQAT figure of one signal, else the first analysis
        end
        da.Value = wanted;
        w.Name = sprintf('SQAT graphs: %s, %s', label, items{strcmp(data, da.Value)});

        switch da.Value
            case 'sqat'
                if ~show_sqat_figure(body, id, entries(1).file, ~pinned)
                    da.Value = data{3};      % the first analysis of the metric, or the statistics
                    if strcmp(da.Value, 'stats')
                        draw_stats(body, entries);
                    else
                        draw_analysis(body, entries, da.Value, chan);
                    end
                end
            case 'all'
                tg = uitabgroup(uigridlayout(body, [1 1], 'Padding', 0));   % the grid sizes it on screen
                for a = entries(1).analyses
                    draw_analysis(uitab(tg, 'Title', a.label), entries, a.id, chan);
                end
                draw_stats(uitab(tg, 'Title', 'Statistics'), entries);
            case 'stats'
                draw_stats(body, entries);
            otherwise
                draw_analysis(body, entries, da.Value, chan);
        end
        % no axes toolbar: Save... keeps the figures, and its tools would move the plots
        for ax = findall(body, 'Type', 'axes')'
            ax.Toolbar.Visible = 'off';
        end
        SQAT_GUI_paint(w, theme_style);             % the new plots in the colours of the theme
    end

    function chans = channel_row(w, paths, results)
        % one channel choice per signal of the window, each with the channels its
        % results hold (and All when there are several); the row shows only when
        % a signal has more than one. Returns the choice of each signal
        row = findobj(w, 'Tag', 'graph_channels');
        before = findobj(row, 'Type', 'uidropdown');
        old_tags = arrayfun(@(d) d.Tag, before, 'UniformOutput', false);
        old_values = arrayfun(@(d) d.Value, before, 'UniformOutput', false);
        delete(row.Children);
        row.ColumnWidth = [repmat({70, 90}, 1, numel(paths)), {'1x'}];   % before the controls: a full grid
        row.RowHeight = {'1x'};                                          % adds a row for each one
        chans = repmat({''}, 1, numel(paths));
        n_c = 0;
        for k_p = 1:numel(paths)
            r = results(strcmp({results.file}, paths{k_p}));
            if isempty(r)
                continue
            end
            c = unique({r.channel}, 'stable');
            c = [sort(c(~strcmp(c, 'Binaural'))), c(strcmp(c, 'Binaural'))];
            if numel(c) > 1
                c{end+1} = 'All'; %#ok<AGROW>
            end
            tag = sprintf('graph_channel_%d', r(1).id);
            old = old_values(strcmp(old_tags, tag));
            chans{k_p} = c{1};
            if ~isempty(old) && il_is_member(old{1}, c)
                chans{k_p} = old{1};          % the choice stays through a redraw
            end
            uilabel(row, 'Text', sprintf('Signal #%d:', r(1).id), 'HorizontalAlignment', 'right', ...
                'Tooltip', r(1).name);
            uidropdown(row, 'Items', c, 'Value', chans{k_p}, 'Tag', tag, 'Enable', numel(c) > 1, ...
                'ValueChangedFcn', @on_graph_control, 'Tooltip', ['Channel of ' r(1).name]);
            n_c = max(n_c, numel(c));
        end
        n = numel(row.Children) / 2;
        row.ColumnWidth = [repmat({70, 90}, 1, n), {'1x'}];
        uilabel(row, 'Text', '');
        row.Parent.RowHeight{2} = il_if(n_c > 1, 30, 0);
    end

    function k_c = sqat_figures(id, path, may_run)
        % the index in the cache of the figures that the SQAT function draws for a
        % signal, drawn now when they are not kept and may_run is true (0 when none)
        f = loaded(strcmp({loaded.path}, path));
        label = il_key_label(store, id);
        k_c = find(strcmp({cache.file}, f.path) & strcmp({cache.metric}, id), 1);
        if ~isempty(k_c) && all(isvalid(cache(k_c).figs))
            return
        end
        k_c = 0;
        if ~may_run
            write_log(sprintf(['The SQAT figure of %s for %s is gone after a new run; ' ...
                'the window shows the first analysis.'], id, f.name));
            return
        end
        write_log(sprintf('Drawing the SQAT figure of %s for %s ...', id, f.name));
        lbl_status.Text = sprintf('Drawing the SQAT figure of %s for %s', label, f.name);
        drawnow limitrate
        try
            k_r = find(strcmp({run_settings.signals.path}, f.path), 1);
            if ~isempty(k_r)
                f = run_settings.signals(k_r);   % the channel and dBFS of the run
            end
            cl = channel_list(f);
            e = run_entry(id);
            if ~e.stereo
                cl = cl(1);                   % one figure per signal: the first channel
            end
            [x, fs] = SQAT_GUI_load(f.path, f.dBFS, cl);
            [~, new_figs] = run_metric(e, x, fs, e.p, true, false);
        catch err
            write_log(sprintf('The SQAT figure of %s could not be drawn: %s', id, err.message));
            lbl_status.Text = 'Ready';
            return
        end
        keep_figures(new_figs, f, id);
        k_c = numel(cache);
        lbl_status.Text = 'Ready';
    end

    function ok = show_sqat_figure(parent, id, path, may_run)
        % the figure that the SQAT function draws, copied into the window; a pinned
        % window does not run the metric again, since the settings may have changed
        label = il_key_label(store, id);
        k_c = sqat_figures(id, path, may_run);
        ok = k_c > 0;
        if ~ok
            return
        end
        delete(parent.Children);
        tg = uitabgroup(uigridlayout(parent, [1 1], 'Padding', 0));   % the grid sizes it on screen
        for k_f = 1:numel(cache(k_c).figs)
            src = cache(k_c).figs(k_f);
            if ~split_figures
                tab = uitab(tg, 'Title', il_figure_title(src, label, k_f));
                copyobj(src.Children, tab);
                set_colormap(tab);
                continue
            end
            axs = flipud(findobj(src, 'Type', 'axes'));   % in the order they were drawn
            for k_ax = 1:numel(axs)
                ax = axs(k_ax);
                tab = uitab(tg, 'Title', il_axes_title(ax, label, k_ax));
                objs = ax;                                 % a legend or colour bar goes with its axes
                if ~isempty(ax.Legend), objs(end+1) = ax.Legend; end %#ok<AGROW>
                if ~isempty(ax.Colorbar), objs(end+1) = ax.Colorbar; end %#ok<AGROW>
                copies = copyobj(objs, tab);
                copies(1).Units = 'normalized';
                copies(1).OuterPosition = [0 0 1 1];
                set_colormap(tab);
            end
        end
        ok = true;
    end

    function draw_analysis(parent, entries, aid, chan)
        % one analysis of the signals: their lines on one axes, or their maps side by side
        A = arrayfun(@(e) e.analyses(strcmp({e.analyses.id}, aid)), entries);
        if numel(entries) > 1
            names = arrayfun(@il_tag, entries, 'UniformOutput', false);   % Signal #1, ch1
        else
            names = {entries.name};
        end
        if strcmp(A(1).kind, 'map')
            n = numel(A);
            n_cols = min(n, 2);
            gl = uigridlayout(parent, [ceil(n / n_cols), n_cols]);
            lo = min(arrayfun(@(a) min(a.z(:), [], 'omitnan'), A));
            hi = max(arrayfun(@(a) max(a.z(:), [], 'omitnan'), A));
            for k = 1:n
                ax = uiaxes(gl);
                surface(ax, A(k).x, A(k).y, zeros(numel(A(k).y), numel(A(k).x)), A(k).z.', ...
                    'EdgeColor', 'none');
                view(ax, 2);
                axis(ax, 'tight');
                ax.Layer = 'top';
                colormap(ax, cmap);
                if isfinite(lo) && hi > lo
                    clim(ax, [lo hi]);            % the same colour scale for every signal
                end
                if strcmp(A(k).bandscale, 'log')
                    ax.YScale = 'log';
                end
                cb = colorbar(ax);
                cb.Label.String = A(k).zlabel;
                xlabel(ax, A(k).xlabel);
                ylabel(ax, A(k).ylabel);
                if numel(entries) > 1
                    title(ax, sprintf('%s: %s', A(k).label, names{k}), 'Interpreter', 'none');
                else
                    title(ax, sprintf('%s: %s, channel %s', A(k).label, names{k}, chan), 'Interpreter', 'none');
                end
            end
            return
        end
        ax = uiaxes(uigridlayout(parent, [1 1]));
        hold(ax, 'on');
        colors = ax.ColorOrder;
        styles = {'-', '--', ':', '-.'};
        for k = 1:numel(A)
            % the colours of the axes, then again dashed, dotted and dash-dotted:
            % from the eighth signal on no two curves look the same
            look = {'Color', colors(mod(k - 1, size(colors, 1)) + 1, :), ...
                    'LineStyle', styles{mod(floor((k - 1) / size(colors, 1)), numel(styles)) + 1}};
            if strcmp(A(k).id, 'tob_level')         % one level per band: a step across its width
                e = A(k).x(:) * 2^(-1/6);
                stairs(ax, [e; A(k).x(end) * 2^(1/6)], [A(k).y(:); A(k).y(end)], look{:});
            else
                plot(ax, A(k).x, A(k).y, look{:});
            end
        end
        hold(ax, 'off');
        y_all = vertcat(A.y);
        y_range = [min(y_all) max(y_all)];
        y_ref = max(abs(y_range));
        if all(isfinite(y_range)) && y_ref > 0 && diff(y_range) <= 5e-3 * y_ref
            % a constant result: show it at +/-5 %, away from its rounding noise
            % (0.5 %: the FFT rounding of Linux leaves 0.18 % on a constant roughness)
            ylim(ax, mean(y_range) + [-0.05 0.05] * y_ref);
        end
        if strcmp(A(1).kind, 'profile') && strcmp(A(1).bandscale, 'log')
            ax.XScale = 'log';
        end
        xlabel(ax, A(1).xlabel);
        ylabel(ax, A(1).ylabel);
        if numel(A) == 1
            what = names{1};
        elseif numel(unique({entries.file})) == 1
            what = entries(1).name;               % the channels of one signal
        else
            what = sprintf('%d signals', numel(unique({entries.file})));
        end
        if isempty(chan)
            title(ax, sprintf('%s: %s', A(1).label, what), 'Interpreter', 'none');
        else
            title(ax, sprintf('%s: %s, channel %s', A(1).label, what, chan), 'Interpreter', 'none');
        end
        if numel(A) > 1
            legend(ax, names, 'Interpreter', 'none', 'Location', 'best');
        end
    end

    function draw_stats(parent, entries)
        % the single values of the signals, one column each
        [names, data] = il_stats_cells(entries);
        uitable(uigridlayout(parent, [1 1]), 'Data', data, 'Tag', 'graph_stats', 'RowName', {}, ...
            'ColumnName', [{'Quantity'}, names], 'ColumnEditable', false);
    end

    function open_save_dialog(w)
        % the signals and the figures of the results to save; from a graphs window it
        % starts with the signals, the metric and the analysis that the window shows
        if nargin == 0
            w = [];
        end
        if isempty(w)
            ws = store;
            paths = {loaded([loaded.marked]).path};
        else
            ws = window_store(w);
            paths = window_paths(w);
        end
        if isempty(ws)
            write_log('No results to save: run an analysis first.');
            return
        end
        delete(findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_save'));
        d = uifigure('Name', 'Save figures', 'Position', [220 160 520 580], 'Visible', fig.Visible, ...
            'Tag', 'SQAT_GUI_save', 'CreateFcn', '');
        g = uigridlayout(d, [6 1]);
        g.RowHeight = {22, '1x', 22, '2x', 30, 30};
        uilabel(g, 'Text', 'Signals', 'FontWeight', 'bold');
        ts = uitree(g, 'checkbox', 'Tag', 'save_signals');
        files = unique({ws.file}, 'stable');
        on = [];
        for k = 1:numel(files)
            en = ws(find(strcmp({ws.file}, files{k}), 1));
            n = uitreenode(ts, 'Text', sprintf('#%d %s', en.id, en.name), 'NodeData', files{k});
            if ismember(files{k}, paths)
                on = [on; n]; %#ok<AGROW>
            end
        end
        if ~isempty(on)
            ts.CheckedNodes = on;
        end
        uilabel(g, 'Text', 'Figures', 'FontWeight', 'bold');
        ti = uitree(g, 'checkbox', 'Tag', 'save_items');
        pre_key = '';
        pre_aid = '';
        if ~isempty(w)
            pre_key = findobj(w, 'Tag', 'graph_metric').Value;
            pre_aid = findobj(w, 'Tag', 'graph_analysis').Value;
        end
        on = [];
        for key = unique({ws.metric}, 'stable')
            m = uitreenode(ti, 'Text', il_key_label(ws, key{1}), 'NodeData', struct('key', key{1}, 'aid', ''));
            A = [ws(strcmp({ws.metric}, key{1})).analyses];
            [ids, i_first] = unique({A.id}, 'stable');
            aids = [{'sqat'}, ids, {'stats'}];
            texts = [{'SQAT figure'}, {A(i_first).label}, {'Statistics'}];
            for j = 1:numel(aids)
                n = uitreenode(m, 'Text', texts{j}, 'NodeData', struct('key', key{1}, 'aid', aids{j}));
                if isempty(w)
                    pick = strcmp(aids{j}, 'sqat');        % from the main window: the SQAT figures
                else
                    pick = strcmp(key{1}, pre_key) && (strcmp(aids{j}, pre_aid) || ...
                        (strcmp(pre_aid, 'all') && ~strcmp(aids{j}, 'sqat')));
                end
                if pick
                    on = [on; n]; %#ok<AGROW>
                end
            end
        end
        if ~isempty(on)
            ti.CheckedNodes = on;
        end
        expand(ti);
        row = uigridlayout(g, [1 3], 'Padding', 0);
        row.ColumnWidth = {60, 90, '1x'};
        uilabel(row, 'Text', 'Format:', 'HorizontalAlignment', 'right');
        uidropdown(row, 'Items', {'PNG', 'PDF'}, 'ItemsData', {'png', 'pdf'}, 'Tag', 'save_format');
        uilabel(row, 'Text', '');
        row = uigridlayout(g, [1 5], 'Padding', 0);
        row.ColumnWidth = {60, '1x', 80, 80, 80};
        uilabel(row, 'Text', 'Folder:', 'HorizontalAlignment', 'right');
        ed = uieditfield(row, 'text', 'Value', save_folder, 'Tag', 'save_folder');
        uibutton(row, 'Text', 'Browse...', 'Tag', 'save_browse', 'ButtonPushedFcn', @(~, ~) browse_save_folder(d, ed));
        uibutton(row, 'Text', 'Save', 'Tag', 'save_do', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) save_chosen(d, ws));
        uibutton(row, 'Text', 'Cancel', 'Tag', 'save_cancel', 'ButtonPushedFcn', @(~, ~) delete(d));
        if il_has_theme()
            theme(d, theme_style);
        end
        SQAT_GUI_paint(d, theme_style);
    end

    function browse_save_folder(d, ed)
        p = uigetdir(ed.Value, 'Folder for the figures');
        figure(d);
        if ~isequal(p, 0)
            ed.Value = p;
        end
    end

    function save_chosen(d, ws)
        % saves the ticked figures of the ticked signals: the SQAT figure of each
        % signal, each analysis with the signals overlaid, the statistics as CSV
        folder = strtrim(findobj(d, 'Tag', 'save_folder').Value);
        if ~isfolder(folder)
            write_log(['ERROR: the folder for the figures does not exist: ' folder]);
            if strcmp(d.Visible, 'on')
                uialert(d, ['The folder does not exist: ' folder], 'Save figures');
            end
            return
        end
        fmt = findobj(d, 'Tag', 'save_format').Value;
        split = split_figures;                          % one file per panel
        n_sig = findobj(d, 'Tag', 'save_signals').CheckedNodes;
        paths = {};
        if ~isempty(n_sig)
            paths = {n_sig.NodeData};
        end
        n_items = findobj(d, 'Tag', 'save_items').CheckedNodes;
        items = struct('key', {}, 'aid', {});
        for k = 1:numel(n_items)
            if ~isempty(n_items(k).NodeData.aid)       % a metric node only groups its figures
                items(end+1) = n_items(k).NodeData; %#ok<AGROW>
            end
        end
        if isempty(paths) || isempty(items)
            write_log('Nothing to save: tick at least one signal and one figure.');
            return
        end
        save_folder = folder;
        delete(d);
        n_saved = 0;
        for it = items
            en = ws(strcmp({ws.metric}, it.key) & ismember({ws.file}, paths));
            [~, o] = sort([en.id]);                     % by signal number, the channels in order
            en = en(o);
            if isempty(en)
                continue
            end
            key_file = strrep(it.key, '#', '_');
            tag = strjoin(arrayfun(@(i) sprintf('s%d', i), unique([en.id]), 'UniformOutput', false), '-');
            switch it.aid
                case 'sqat'
                    for p = unique({en.file}, 'stable')
                        if ~any(strcmp({loaded.path}, p{1}))
                            continue                    % removed from the list since
                        end
                        k_c = sqat_figures(it.key, p{1}, true);
                        if k_c == 0
                            continue
                        end
                        f = loaded(strcmp({loaded.path}, p{1}));
                        [~, base] = fileparts(f.name);
                        base = sprintf('%s_s%d_%s', base, f.id, key_file);
                        figs = cache(k_c).figs;
                        for k_fig = 1:numel(figs)
                            if split
                                axs = flipud(findobj(figs(k_fig), 'Type', 'axes'));
                                for k_ax = 1:numel(axs)
                                    n_saved = n_saved + il_export(axs(k_ax), folder, ...
                                        sprintf('%s_%d_%d.%s', base, k_fig, k_ax, fmt));
                                end
                            else
                                n_saved = n_saved + il_export(figs(k_fig), folder, ...
                                    sprintf('%s_%d.%s', base, k_fig, fmt));
                            end
                        end
                    end
                case 'stats'
                    [names, data] = il_stats_cells(en);
                    writecell([[{'Quantity'}, names]; data], ...
                        il_free_path(fullfile(folder, sprintf('%s_statistics_%s.csv', key_file, tag))));
                    n_saved = n_saved + 1;
                otherwise
                    en = en(arrayfun(@(e) any(strcmp({e.analyses.id}, it.aid)), en));
                    if isempty(en)
                        continue
                    end
                    h = uifigure('Visible', 'off', 'Position', [0 0 1100 650], 'Tag', 'SQAT_GUI_export', ...
                        'CreateFcn', '');
                    chan = unique({en.channel});
                    if ~isscalar(chan)
                        chan = {''};                    % mixed channels: the legend names them
                    end
                    draw_analysis(h, en, it.aid, chan{1});
                    for ax = findobj(h, 'Type', 'axes')'
                        % a hidden window computes the automatic limits only when they are read;
                        % fixed here, they survive the drawing of the export
                        ax.XLim = ax.XLim;
                        ax.YLim = ax.YLim;
                    end
                    n_saved = n_saved + il_export(h, folder, sprintf('%s_%s_%s.%s', key_file, it.aid, tag, fmt));
                    delete(h);
            end
        end
        write_log(sprintf('%d file(s) saved to %s', n_saved, folder));
    end

    function draw_waveform_window()
        sync_wave_tabs();
        cla(ax_wave);
        f = active_file();
        ch = channel_of_active();
        try
            [x, fs] = SQAT_GUI_load(f.path, f.dBFS, ch);
            y = audioread(f.path);
            y = y(:, min(ch, size(y, 2)));
        catch err
            write_log(sprintf('ERROR reading %s: %s', f.name, err.message));
            return
        end
        key = sprintf('%s|%d', f.path, ch);
        if ~strcmp(key, wave_key)               % another signal: what belonged to the last one goes
            stop_player();
            player = [];
            play_y = [];
            play_start = 0;
            boxes = zeros(0, 4);
            box_corner = [];
            set(findobj(win_wave, 'Tag', 'draw_box'), 'Value', false);
            btn_play.Text = 'Play';
            wave_key = key;
        end
        wave_x = x;
        wave_y = y;
        wave_fs = fs;
        spec_signal = sprintf('%s|%s', key, mat2str(f.dBFS));
        cancel_jobs(spec_signal);
        t_end = numel(x) / fs;
        step = max(1, ceil(numel(x) / 2e6));   % display only: at most 2e6 points
        t = (0:numel(x)-1)' / fs;
        plot(ax_wave, t(1:step:end), x(1:step:end), 'PickableParts', 'none', 'Tag', 'wave_line', ...
            'UserData', [step -inf inf]);
        ylim(ax_wave, 1.05 * max(max(abs(x)), eps) * [-1 1]);   % fixed: a zoom redraws the line, the scale stays
        xlim(ax_wave, [0 t_end]);
        ylabel(ax_wave, 'Sound pressure (Pa)');
        title(ax_wave, sprintf('Waveform: %s, channel %d', f.name, ch), 'Interpreter', 'none');   % the signal on screen
        xline(ax_wave, (max(play_start, 1) - 1) / fs, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5, ...
            'Tag', 'playhead', 'PickableParts', 'none');
        draw_level();
        draw_spectrogram();
        il_delete_timer(prefetch_timer);            % the pool starts once the window is drawn, not while
        prefetch_timer = timer('StartDelay', 1, 'ExecutionMode', 'singleShot', ...
            'TimerFcn', @(~, ~) prefetch_maps(), 'ObjectVisibility', 'off');
        start(prefetch_timer);
    end

    function prefetch_maps()
        % both enhanced maps of the signal on screen, computed in one job
        if il_is_open(win_wave)
            request_maps();
        end
    end

    function draw_spectrogram(keep_view)
        % keep_view: an option changed, not the signal, so the zoom of the spectrogram stays
        if isempty(wave_x)
            return
        end
        lims = {};
        if nargin > 0 && keep_view && ~isempty(findobj(ax_spec, 'Type', 'surface'))
            lims = {ax_spec.XLim, ax_spec.YLim};
        end
        zoomed = ~isempty(lims) && diff(lims{1}) < 0.95 * numel(wave_x) / wave_fs;
        dd = findobj(win_wave, 'Tag', 'spec_window');
        degree = findobj(win_wave, 'Tag', 'spec_degree').Value;
        overlap = findobj(win_wave, 'Tag', 'spec_overlap').Value;
        weighting = findobj(win_wave, 'Tag', 'wave_weighting').Value;
        win = dd.Value;
        if strcmp(win, 'custom')
            win = spec_custom;
        end
        enhanced = strcmp(findobj(win_wave, 'Tag', 'spec_enhanced').Value, 'On');
        if enhanced
            mode = findobj(win_wave, 'Tag', 'spec_enhanced_mode').Value;
            full = full_map(mode, false);
            preview = false;
            if isempty(full) && has_pool()
                request_maps();                         % drawn when it arrives; the window stays usable
                full = preview_map(mode);               % meanwhile the preview of a long signal, if it came
                preview = ~isempty(full);
                if ~preview && ~zoomed
                    draw_waiting(mode, weighting);
                    return
                end
            end
            computing = spec_guard();   %#ok<NASGU>
            busy = spec_loading('Drawing the enhanced spectrogram...');   %#ok<NASGU> closes when the map is drawn
            if isempty(full) && ~zoomed
                full = full_map(mode, true);            % no background pool: computed here
            end
            prev = spec_view;
            spec_view = struct('full', full, 'mode', mode, 'weighting', weighting, ...
                'zoomed', zoomed, 'zoom', [], 'preview', preview);
            if zoomed && ~isempty(prev) && ~isempty(prev.zoom) && strcmp(prev.mode, mode)
                spec_view.zoom = prev.zoom;             % the same excerpt, only weighted otherwise
            end
            if zoomed
                m = zoom_map(lims{1});
            else
                m = spec_view.full;
            end
            t_spec = m.t;
            f_spec = m.f;
            L = weighted(m);
            top = enhanced_top(m);
        else
            spec_view = [];
            [t_spec, f_spec, L, info] = SQAT_GUI_spectrogram(wave_x, wave_fs, win, degree, overlap);
            if info.limited
                write_log(sprintf('The spectrogram was limited to %d frames: the overlap is %.0f %%.', ...
                    numel(t_spec), info.overlap));
            end
            keep = f_spec >= 20;
            f_spec = f_spec(keep);
            top = max(L(keep, :), [], 'all');                % without weighting: A or C keep the colours
            L = L(keep, :) + SQAT_GUI_weight_curve(f_spec, wave_fs, weighting);
        end
        spec_busy = true;
        cla(ax_spec);
        if enhanced                  % a log-spaced grid: one textured face, much lighter to draw than a mesh
            [xe, ye] = il_map_edges(t_spec, f_spec);
            surface(ax_spec, xe, ye, zeros(2), 'CData', L, 'FaceColor', 'texturemap', ...
                'EdgeColor', 'none', 'PickableParts', 'none');
        else
            surface(ax_spec, t_spec, f_spec, zeros(numel(f_spec), numel(t_spec)), L, ...
                'EdgeColor', 'none', 'PickableParts', 'none');
        end
        ax_spec.YScale = 'log';
        ax_spec.Layer = 'top';
        if isempty(lims)
            xlim(ax_spec, [0 numel(wave_x) / wave_fs]);
            ylim(ax_spec, [20 wave_fs/2]);
        else
            xlim(ax_spec, lims{1});
            ylim(ax_spec, lims{2});
        end
        colormap(ax_spec, cmap);
        spec_top = top;
        spec_range = il_if(enhanced, 45, 80);
        spec_black = min(spec_black, spec_range - 5);
        apply_black();
        cb = ax_spec.Colorbar;                         % the one colour bar of the plot, kept between drawings
        delete(setdiff(findall(ax_spec.Parent, 'Type', 'colorbar'), cb));   % a stray one would overlap its label
        if isempty(cb)
            cb = colorbar(ax_spec);
        end
        cb.Label.String = sprintf('SPL (%s)', il_level_unit(weighting));   % short: it fits the height of the plot
        il_fit_axes(ax_spec.Parent, ax_spec, 26);      % the time label only on the sound level, the plot at the bottom
        ylabel(ax_spec, 'Frequency (Hz)');
        if enhanced
            spec_title();
        else
            title(ax_spec, sprintf('Spectrogram (%s window, %d points, %.0f %% overlap, %cf %.1f Hz)', ...
                dd.Items{strcmp(dd.ItemsData, dd.Value)}, info.n_fft, info.overlap, 916, wave_fs / info.n_fft), ...
                'Interpreter', 'none');                          % the frequency resolution, fs/N
        end
        xline(ax_spec, (max(play_start, 1) - 1) / wave_fs, 'Color', [1 1 1], 'LineWidth', 1.5, ...
            'Tag', 'playhead_spectrogram', 'PickableParts', 'none');
        draw_boxes();
        spec_busy = false;
        drawnow                                        % the limit events of this drawing arrive now, and not at the next zoom
    end

    function on_wave_limits(~, event)
        % a zoom or a pan of the waveform: the spectrogram follows, and its own
        % LimitsChangedFcn takes the enhanced map to the excerpt
        follow_limits(ax_wave, ax_spec, event);
        draw_wave_line();
        if ~isequal(ax_lvl.XLim, ax_wave.XLim)          % the sound level follows the waveform
            setappdata(ax_lvl, 'sqat_echo', ax_wave.XLim);
            ax_lvl.XLim = ax_wave.XLim;
        end
    end

    function on_level_limits(~, event)
        % a zoom or a pan of the sound level: the waveform and the spectrogram follow
        if isequal(event.NewLimits, getappdata(ax_lvl, 'sqat_echo'))
            setappdata(ax_lvl, 'sqat_echo', []);
            return
        end
        follow_limits(ax_lvl, ax_wave, event);
        follow_limits(ax_lvl, ax_spec);
    end

    function draw_level()
        % the sound level meter of SQAT (Do_SLM) on the signal on screen, with
        % the frequency weighting of the player and the time weighting chosen
        % above the plot; the indicators as the Sound level metric gives them
        cla(ax_lvl);
        lvl_L = [];
        if isempty(wave_x)
            show_level_values();
            return
        end
        fw = findobj(win_wave, 'Tag', 'wave_weighting').Value;
        dd_tw = findobj(win_wave, 'Tag', 'level_time_weighting');
        lvl_L = Do_SLM(wave_x, wave_fs, fw, dd_tw.Value, 94);   % the signal is in Pa: 94 dBFS keeps it
        lvl_L = lvl_L(:);
        step = max(1, round(wave_fs / 1000));                  % the level every millisecond, as the metric
        plot(ax_lvl, (0:step:numel(lvl_L) - 1)' / wave_fs, lvl_L(1:step:end), ...
            'PickableParts', 'none', 'Tag', 'level_line');
        xlim(ax_lvl, ax_wave.XLim);
        ylabel(ax_lvl, sprintf('SPL (%s)', il_level_unit(fw)));   % short: the plot is low
        xlabel(ax_lvl, 'Time (s)');
        title(ax_lvl, sprintf('Sound pressure level (%s-weighted, %s)', fw, dd_tw.Items{strcmp(dd_tw.ItemsData, dd_tw.Value)}));
        xline(ax_lvl, (max(play_start, 1) - 1) / wave_fs, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5, ...
            'Tag', 'playhead_level', 'PickableParts', 'none');
        il_fit_axes(ax_lvl.Parent, ax_lvl);
        show_level_values();
    end

    function show_level_values()
        % Leq, LE (SEL), Lmax and the levels exceeded during the two percentages
        % of the time chosen, in one line above the plot
        lbl = findobj(win_wave, 'Tag', 'level_indicators');
        if isempty(lvl_L)
            lbl.Text = '';
            return
        end
        F = findobj(win_wave, 'Tag', 'wave_weighting').Value;
        T = upper(findobj(win_wave, 'Tag', 'level_time_weighting').Value);
        Leq = Get_Leq(lvl_L, wave_fs);
        txt = sprintf('L%seq %.1f   L%sE %.1f   L%s%smax %.1f', F, Leq, F, Leq + 10*log10(numel(lvl_L) / wave_fs), ...
            F, T, max(lvl_L));
        pct = unique([findobj(win_wave, 'Tag', 'level_percentile_1').Value, ...
            findobj(win_wave, 'Tag', 'level_percentile_2').Value]);
        for p = pct
            txt = sprintf('%s   L%s%s%g %.1f', txt, F, T, p, get_exceeded_value(lvl_L, p));
        end
        lbl.Text = [txt ' ' il_level_unit(F)];
    end

    function follow_limits(src, dst, event)
        % the time axes of the waveform and the spectrogram move together, inside
        % the file; a span under 50 ms is refused. The event arrives at the next
        % drawnow, so the limits one axis gave the other are recognised as its echo
        if isempty(wave_x)
            return
        end
        if nargin > 2 && isequal(event.NewLimits, getappdata(src, 'sqat_echo'))
            setappdata(src, 'sqat_echo', []);
            return
        end
        lim = min(max(src.XLim, 0), numel(wave_x) / wave_fs);
        if diff(lim) < 0.05
            lim = dst.XLim;
        end
        if ~isequal(src.XLim, lim)
            src.XLim = lim;
        end
        if ~isequal(dst.XLim, lim)
            setappdata(dst, 'sqat_echo', lim);
            dst.XLim = lim;
        end
    end

    function draw_wave_line()
        % the waveform of the view and one view to each side, at most 2e6 points:
        % a zoom shows every sample again, a pan inside the drawn part draws nothing
        h = findobj(ax_wave, 'Tag', 'wave_line');
        if isempty(h) || isempty(wave_x)
            return
        end
        lim = ax_wave.XLim;
        n = numel(wave_x);
        i1 = max(1, floor((lim(1) - diff(lim)) * wave_fs) + 1);
        i2 = min(n, ceil((lim(2) + diff(lim)) * wave_fs) + 1);
        step = max(1, ceil((i2 - i1 + 1) / 2e6));
        u = h.UserData;
        if u(1) == step && u(2) <= lim(1) && u(3) >= lim(2)
            return
        end
        i = (i1:step:i2)';
        set(h, 'XData', (i - 1) / wave_fs, 'YData', wave_x(i), ...
            'UserData', [step il_if(i1 == 1, -inf, (i1 - 1) / wave_fs) il_if(i2 == n, inf, (i2 - 1) / wave_fs)]);
    end

    function on_spec_limits(~, event)
        % a zoom or a pan of the spectrogram: the waveform follows, and the enhanced
        % map of the excerpt is recomputed once the limits settle, since the full
        % map has at most 2000 columns and a zoom only stretches them
        if nargin > 1
            follow_limits(ax_spec, ax_wave, event);
        else
            follow_limits(ax_spec, ax_wave);
        end
        if spec_busy || isempty(spec_view)
            return
        end
        stop_spec_timer();
        spec_timer = timer('StartDelay', 0.3, 'ExecutionMode', 'singleShot', ...
            'TimerFcn', @(~, ~) apply_spec_zoom(), 'ObjectVisibility', 'off');
        start(spec_timer);
    end

    function closer = spec_loading(msg)
        % a wait message over the waveform window; it closes when the closer of
        % the step that opened it is cleared, and a later step only changes its text
        if il_is_open(spec_dlg)
            spec_dlg.Message = msg;
            closer = [];
            return
        end
        closer = [];
        if ~strcmp(win_wave.Visible, 'on')                     % the dialog needs a visible window
            return
        end
        spec_dlg = uiprogressdlg(win_wave, 'Title', 'Enhanced STFT', 'Message', msg, 'Indeterminate', 'on');
        drawnow
        d = spec_dlg;
        closer = onCleanup(@() delete(d));
    end

    function stop_spec_timer()
        if ~isempty(spec_timer) && isvalid(spec_timer)
            stop(spec_timer);
            delete(spec_timer);
        end
        spec_timer = [];
    end

    function apply_spec_zoom()
        if isempty(spec_view) || ~il_is_open(win_wave) || spec_computing
            return
        end
        srf = findobj(ax_spec, 'Type', 'surface');
        if isempty(srf)
            return
        end
        computing = spec_guard();   %#ok<NASGU>
        lim = ax_spec.XLim;
        cleanup = onCleanup(@() recheck_zoom(lim));   %#ok<NASGU> after computing: the view may have moved
        span = diff(lim);
        hop = max(round(0.001 * wave_fs), ceil(numel(wave_x) / 2000)) / wave_fs;   % time step of the full map
        if span >= 0.95 * numel(wave_x) / wave_fs || span / 2000 >= 0.999 * hop
            if spec_view.zoomed                                 % back to the full map
                if isempty(spec_view.full)                      % the mode changed during the zoom
                    spec_view.full = full_map(spec_view.mode, ~has_pool());
                    if isempty(spec_view.full)                  % on its way from the background pool
                        request_maps();
                        spec_view.full = preview_map(spec_view.mode);
                        spec_view.preview = ~isempty(spec_view.full);
                        if ~spec_view.preview
                            note_waiting('Computing the whole file in the background...');
                            return
                        end
                    end
                    spec_top = enhanced_top(spec_view.full);
                    apply_black();
                end
                busy = spec_loading('Drawing the enhanced spectrogram...');   %#ok<NASGU>
                show_map(srf(1), spec_view.full);
                spec_view.zoomed = false;
                spec_title();
            end
            return
        end
        busy = spec_loading('Drawing the enhanced spectrogram...');   %#ok<NASGU>
        show_map(srf(1), zoom_map(lim));
        spec_view.zoomed = true;
        spec_title();
    end

    function spec_title()
        % the title of the enhanced map, which says when the map on screen is the preview
        s = sprintf('Enhanced STFT (consensus of 8 to 512 ms windows, %s)', spec_view.mode);
        if spec_view.preview && ~spec_view.zoomed
            s = [s ': preview, the exact map follows'];
        end
        title(ax_spec, s, 'Interpreter', 'none');
    end

    function guard = spec_guard()
        % marks a computation of an enhanced map: a zoom event that arrives
        % meanwhile (during the drawnow of the wait message) does not start another
        spec_computing = true;
        guard = onCleanup(@() end_spec_computing());
    end

    function end_spec_computing()
        spec_computing = false;
    end

    function recheck_zoom(lim)
        if il_is_open(win_wave) && ~isempty(spec_view) && ~isequal(ax_spec.XLim, lim)
            on_spec_limits();
        end
    end

    function m = full_map(mode, compute)
        % the enhanced map of the whole signal, not weighted: from the cache, or
        % computed and kept when compute is true ([] otherwise)
        key = [spec_signal '|' mode];
        k = find(strcmp({spec_cache.key}, key), 1);
        if ~isempty(k)
            m = spec_cache(k).map;
            return
        end
        m = [];
        if compute
            busy = spec_loading('Computing the enhanced spectrogram...');   %#ok<NASGU> closes on return
            ms = enhanced_map(wave_x, {'readable', 'sharp'}, []);   % both modes for the time of one
            store_map([spec_signal '|readable'], ms(1));
            store_map([spec_signal '|sharp'], ms(2));
            m = ms(1 + strcmp(mode, 'sharp'));
        end
    end

    function store_map(key, m)
        spec_cache(strcmp({spec_cache.key}, key)) = [];
        spec_cache(end+1) = struct('key', key, 'map', m);
        if numel(spec_cache) > 6                                % about 15 MB each
            spec_cache(1) = [];
        end
    end

    function tf = has_pool()
        % the background pool of MATLAB (R2021b or newer, no toolbox needed): the
        % full maps are computed there while the window stays usable
        if ~spec_pool_tried
            spec_pool_tried = true;
            if ~isappdata(fig, 'sqat_no_background') && ~isappdata(groot, 'sqat_no_background')   % set by the tests
                try
                    spec_pool = backgroundPool;
                catch
                    spec_pool = [];
                end
            end
        end
        tf = ~isempty(spec_pool);
    end

    function request_maps()
        % queues both full maps of the signal on screen in one job, unless they are
        % kept or already queued; a long signal gets a preview first (see
        % SQAT_GUI_enhanced_stft)
        key = [spec_signal '|both'];
        if isempty(wave_x) || ~has_pool() || all(ismember(strcat(spec_signal, {'|readable', '|sharp'}), {spec_cache.key}))
            return
        end
        if numel(wave_x) > preview_min * wave_fs && ~any(startsWith({spec_preview.key}, [spec_signal '|']))
            queue_job([key '|preview'], {wave_x, wave_fs, {'readable', 'sharp'}, [], [], [], true});
        end
        queue_job(key, {wave_x, wave_fs, {'readable', 'sharp'}});
    end

    function queue_job(key, args)
        if any(strcmp({spec_jobs.key}, key))
            return
        end
        job = parfeval(spec_pool, @SQAT_GUI_enhanced_stft, 4, args{:});
        afterEach(job, @(done) on_map_done(key, done), 0, 'PassFuture', true);
        spec_jobs(end+1) = struct('key', key, 'job', job);
    end

    function cancel_jobs(keep)
        % drops the queued maps of other signals (all of them without keep); the pool may have a
        % single worker (R2024a without the Parallel Computing Toolbox), and then the maps queue
        for k = numel(spec_jobs):-1:1
            if nargin == 0 || ~startsWith(spec_jobs(k).key, [keep '|'])
                cancel(spec_jobs(k).job);
                spec_jobs(k) = [];
            end
        end
        if nargin == 0
            spec_preview(:) = [];
        else
            spec_preview(~startsWith({spec_preview.key}, [keep '|'])) = [];
        end
    end

    function on_map_done(key, job)
        k = find(strcmp({spec_jobs.key}, key), 1);
        if isempty(k) || spec_jobs(k).job.ID ~= job.ID          % cancelled, or asked for again since
            return
        end
        spec_jobs(k) = [];
        if ~isempty(job.Error)
            write_log(['ERROR computing the enhanced spectrogram: ' job.Error.message]);
            return
        end
        [t, f, L, info] = fetchOutputs(job);
        keep = f >= 20;
        preview = endsWith(key, '|preview');
        base = extractBefore(key, '|both');
        for mode = {'readable', 'sharp'}
            k_m = [base '|' mode{1}];
            i_m = 1 + strcmp(mode{1}, 'sharp');
            m = struct('t', t, 'f', f(keep), 'L', single(L{i_m}(keep, :)), 'ref', info.ref(i_m));
            if preview
                if ~any(strcmp({spec_cache.key}, k_m))          % unless the exact map came first
                    spec_preview(end+1) = struct('key', k_m, 'map', m); %#ok<AGROW>
                end
            else
                store_map(k_m, m);
                spec_preview(strcmp({spec_preview.key}, k_m)) = [];
            end
            show_ready_map(k_m);
        end
    end

    function m = preview_map(mode)
        k = find(strcmp({spec_preview.key}, [spec_signal '|' mode]), 1);
        m = [];
        if ~isempty(k)
            m = spec_preview(k).map;
        end
    end

    function show_ready_map(key)
        % a full map arrived: it goes on screen when the spectrogram waits for it,
        % or shows the preview in its place
        if ~il_is_open(win_wave) || isempty(spec_view) || ~strcmp(key, [spec_signal '|' spec_view.mode])
            return
        end
        if spec_computing || spec_busy                          % a step in the foreground: once it ends
            t = timer('StartDelay', 0.2, 'TimerFcn', @(~, ~) show_ready_map(key), ...
                'StopFcn', @(tm, ~) delete(tm), 'ObjectVisibility', 'off');
            start(t);
            return
        end
        m = full_map(spec_view.mode, false);
        preview = isempty(m);
        if preview
            m = preview_map(spec_view.mode);
        end
        if isempty(m) || (~isempty(spec_view.full) && (preview || ~spec_view.preview))
            return                                              % nothing better than what is on screen
        end
        if ~spec_view.zoomed || isempty(findobj(ax_spec, 'Type', 'surface'))
            draw_spectrogram(true);
            return
        end
        spec_view.full = m;                                     % zoomed: for Home, or on screen if at the whole file
        spec_view.preview = preview;
        spec_view.zoom = [];                                    % the excerpt again, on the scale of the new map
        spec_top = enhanced_top(m);
        apply_black();
        stop_spec_timer();
        apply_spec_zoom();
    end

    function draw_waiting(mode, weighting)
        % the enhanced map is on its way from the background pool: a note in its place
        spec_view = struct('full', [], 'mode', mode, 'weighting', weighting, 'zoomed', false, 'zoom', [], ...
            'preview', false);
        spec_busy = true;
        cla(ax_spec);
        ax_spec.YScale = 'log';
        xlim(ax_spec, [0 numel(wave_x) / wave_fs]);
        ylim(ax_spec, [20 wave_fs/2]);
        title(ax_spec, sprintf('Enhanced STFT (%s)', mode), 'Interpreter', 'none');
        note_waiting('Computing the enhanced spectrogram in the background...');
        spec_busy = false;
    end

    function note_waiting(msg)
        delete(findobj(ax_spec, 'Tag', 'spec_wait'));
        text(ax_spec, 0.5, 0.5, msg, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 14, 'Tag', 'spec_wait', 'PickableParts', 'none');
    end

    function m = zoom_map(lim)
        % the enhanced map of the excerpt in view, kept while the view and the mode stay
        z = spec_view.zoom;
        if ~isempty(z) && isequal(z.lim, lim)
            m = z;
            return
        end
        span = diff(lim);
        pad = 0.3;                                             % the longest window reaches 256 ms around each instant
        i1 = max(1, floor((lim(1) - pad) * wave_fs) + 1);
        i2 = min(numel(wave_x), ceil((lim(2) + pad) * wave_fs));
        n_frames = ceil((i2 - i1 + 1) / max(1, round(max(0.001, span / 2000) * wave_fs)));
        busy = spec_loading('Computing the enhanced spectrogram of the zoomed excerpt...');   %#ok<NASGU>
        ref = [];                                              % the scale of the whole map, when it is there
        if ~isempty(spec_view.full)
            ref = spec_view.full.ref;
        end
        m = enhanced_map(wave_x(i1:i2), spec_view.mode, n_frames, ref);
        m.t = m.t + (i1 - 1) / wave_fs;
        m.lim = lim;
        spec_view.zoom = m;
    end

    function m = enhanced_map(x, mode, n_frames, ref)
        % the map of one mode, or a map per mode when mode is a cell array; ref (optional):
        % the references of the whole map, for an excerpt
        if nargin < 4
            ref = [];
        end
        [t, f, L, info] = SQAT_GUI_enhanced_stft(x, wave_fs, mode, n_frames, [], [], [], ref);
        keep = f >= 20;
        if ~iscell(L)
            L = {L};
        end
        for k = numel(L):-1:1
            m(k) = struct('t', t, 'f', f(keep), 'L', single(L{k}(keep, :)), 'ref', info.ref(k));   % single: dB, half the memory
        end
    end

    function top = enhanced_top(m)
        % the top of the colour scale: the loudest cell of the whole file without weighting,
        % in either mode, so that neither the mode, the weighting nor a zoom moves the colours
        top = max(m.L, [], 'all');
        for md = {'readable', 'sharp'}
            o = full_map(md{1}, false);
            if isempty(o)
                o = preview_map(md{1});
            end
            if ~isempty(o)
                top = max(top, max(o.L, [], 'all'));
            end
        end
        top = double(top);
    end

    function L = weighted(m)
        L = double(m.L) + SQAT_GUI_weight_curve(m.f, wave_fs, spec_view.weighting);
    end

    function show_map(srf, m)
        delete(findobj(ax_spec, 'Tag', 'spec_wait'));
        [xe, ye] = il_map_edges(m.t, m.f);
        set(srf, 'XData', xe, 'YData', ye, 'ZData', zeros(2), 'CData', weighted(m));
    end

    function move_playhead(sample)
        if ~il_is_open(win_wave) || isempty(wave_fs)
            return
        end
        t_now = (sample - 1) / wave_fs;
        set(findobj(win_wave, 'Tag', 'playhead'), 'Value', t_now);
        set(findobj(win_wave, 'Tag', 'playhead_spectrogram'), 'Value', t_now);
        set(findobj(win_wave, 'Tag', 'playhead_level'), 'Value', t_now);
    end

    function [OUT, new_figs] = run_step(e, step, done, x, fs, f, show)
        % the result of one metric of the run: taken from another metric that
        % computed it on the way to its own result (see SQAT_GUI_share), or
        % computed here. A result taken this way carries no figure, so the
        % graphs window draws it when it is asked for.
        new_figs = [];
        if ~isempty(step.from) && isfield(done, step.from)
            OUT = SQAT_GUI_take(done.(step.from), step);
            if ~isempty(OUT)
                write_log(sprintf('%s on %s: taken from %s, the same computation.', ...
                    e.key, f.name, step.from));
                return
            end
        end
        write_log(sprintf('Running %s on %s ...', e.key, f.name));
        [OUT, new_figs] = run_metric(e, x, fs, e.p, show);
    end

    function plan = share_plan(keys, n_samples, fs)
        % the first analysis of each metric shares computations (SQAT_GUI_share);
        % a second one of the same metric runs on its own
        first = keys(~contains(keys, '#'));
        P = struct();
        for a = run_settings.analyses
            if ismember(a.key, first)
                P.(a.key) = a.p;
            end
        end
        plan = SQAT_GUI_share(first, P, n_samples, fs);
        for key = keys(contains(keys, '#'))
            plan(end+1) = struct('id', key{1}, 'from', '', 'field', '', 'attach', {{}}, 'restat', []); %#ok<AGROW>
        end
    end

    function [OUT, new_figs] = run_metric(e, x, fs, p, show, fallback)
        % runs one metric; figures are drawn hidden and returned. When the
        % figure fails and fallback is true, the values are computed again
        % without it, so a plotting defect of the toolbox costs no result.
        if nargin < 6
            fallback = true;
        end
        figs_before = findall(groot, 'Type', 'figure');
        default_visible = get(groot, 'DefaultFigureVisible');
        set(groot, 'DefaultFigureVisible', 'off');
        restore = onCleanup(@() set(groot, 'DefaultFigureVisible', default_visible));
        try
            [txt, OUT] = evalc('e.run(x, fs, p, show)');
        catch err
            delete(il_new_figures(figs_before));
            if ~show || ~fallback
                rethrow(err);
            end
            write_log(sprintf('The SQAT figure of %s could not be drawn (%s); values computed without it.', ...
                e.key, err.message));
            [txt, OUT] = evalc('e.run(x, fs, p, false)');
        end
        write_block(txt);
        new_figs = il_new_figures(figs_before);
    end

    function cl = channel_list(f)
        % the channels of file f that its Channel in the list asks for
        if strcmp(f.channel, 'All')
            cl = 1:f.nch;
        else
            cl = min(str2double(f.channel), f.nch);
        end
    end

    function keep_figures(new_figs, f, id, cache_it)
        % inferno colour scale, then kept hidden for the graphs windows and the saving
        if nargin < 4
            cache_it = true;
        end
        for k_fig = 1:numel(new_figs)
            set_colormap(new_figs(k_fig));
        end
        if cache_it
            cache(end+1) = struct('file', f.path, 'metric', id, 'figs', new_figs);
        else
            delete(new_figs(isvalid(new_figs)));
        end
    end

    function clear_cache(keep)
        % the SQAT figures go, except those of the results kept (file|key)
        if nargin < 1
            keep = {};
        end
        stay = false(1, numel(cache));
        for k_c = 1:numel(cache)
            stay(k_c) = ismember([cache(k_c).file '|' cache(k_c).metric], keep);
            if ~stay(k_c)
                delete(cache(k_c).figs(isvalid(cache(k_c).figs)));
            end
        end
        cache = cache(stay);
    end

    function set_colormap(parent)
        for ax = findobj(parent, 'Type', 'axes')'
            colormap(ax, cmap);
        end
    end

    function apply_theme()
        if strcmp(theme_style, 'dark')
            img_logo.ImageSource = fullfile(dir_logos, 'logo_white.png');
            btn_theme.Text = char(9788);         % a sun: the light theme is one click away
            btn_theme.Tooltip = 'Light theme';
        else
            img_logo.ImageSource = fullfile(dir_logos, 'logo.png');
            btn_theme.Text = char(9790);         % a moon
            btn_theme.Tooltip = 'Dark theme';
        end
        if ~il_has_theme()
            btn_theme.Enable = 'off';
            btn_theme.Tooltip = 'Dark theme requires MATLAB R2025a or newer';
        else
            for w = [fig, open_graph_windows()]
                if il_is_open(w) && strcmp(theme_style, 'dark')
                    SQAT_GUI_paint(w, 'dark');   % back to auto before the switch, which then sets them
                end
                if il_is_open(w)
                    theme(w, theme_style);
                end
            end
        end
        if strcmp(theme_style, 'light')
            for w = [fig, open_graph_windows()]
                if il_is_open(w)
                    SQAT_GUI_paint(w, 'light');  % the light theme in white
                end
            end
        end
    end

    function status_colour(is_error)
        % red for an error, else the colour of the theme
        if is_error
            lbl_status.FontColor = [0.85 0.2 0.2];
        elseif isprop(lbl_status, 'FontColorMode')
            lbl_status.FontColorMode = 'auto';
        else
            lbl_status.FontColor = [0 0 0];
        end
    end

    function set_status(text)
        lbl_status.Text = text;
        status_colour(false);
        if il_is_open(dlg)
            dlg.Message = text;
        end
    end

    function poll_cancel()
        % the Stop of the progress dialog ends the run at the next step
        if il_is_open(dlg) && dlg.CancelRequested
            on_stop_run();
        end
    end

    function close_progress()
        if il_is_open(dlg)
            close(dlg);
        end
        dlg = [];
        if isvalid(fig) && isappdata(fig, 'sqat_progress')
            rmappdata(fig, 'sqat_progress');
        end
    end

    function set_progress(value)
        if il_is_open(dlg)
            dlg.Value = min(max(value / 100, 0), 1);
        end
        gauge.Value = value;
        if value > 0
            gauge.ScaleColors = green;
            gauge.ScaleColorLimits = [0 value];
        else
            gauge.ScaleColors = zeros(0, 3);
            gauge.ScaleColorLimits = zeros(0, 2);
        end
        drawnow limitrate
    end

    function write_log(msg)
        stamp = char(datetime('now', 'Format', 'HH:mm:ss'));
        lines = console.Value;
        if isequal(lines, {''})
            lines = {};
        end
        console.Value = [lines(:); {sprintf('[%s] %s', stamp, msg)}];
        if startsWith(msg, {'ERROR', 'No '})       % what needs the eye shows on the status bar too
            lbl_status.Text = msg;
            status_colour(startsWith(msg, 'ERROR'));
        end
        try
            scroll(console, 'bottom');
        catch
        end
        drawnow limitrate
    end

    function write_block(txt)
        lines = strtrim(splitlines(string(txt)));
        lines = lines(lines ~= "");
        if isempty(lines)
            return
        end
        console.Value = [console.Value(:); cellstr("    " + lines)];
    end

end

function new = il_new_figures(figs_before)
% figures created since figs_before, apart from the windows of the interface
figs_after = findall(groot, 'Type', 'figure');
new = figs_after(~ismember(figs_after, figs_before));
if ~isempty(new)
    new = new(~startsWith({new.Tag}, 'SQAT_GUI'));
end
end

function t = il_figure_title(src, label, k)
if ~isempty(src.Name)
    t = src.Name;
elseif k == 1
    t = label;
else
    t = sprintf('%s (%d)', label, k);
end
end

function t = il_axes_title(ax, label, k)
t = ax.Title.String;
if iscell(t)
    t = strjoin(t, ' ');
end
t = strtrim(regexprep(char(t), '[\\{}$^_]', ''));
if isempty(t)
    t = sprintf('%s, panel %d', label, k);
end
end

function tf = il_is_open(h)
tf = ~isempty(h) && isvalid(h);
end

function [names, data] = il_stats_cells(entries)
% the single values of the signals, one column each: the column names and the
% cells (quantity, then the value of each signal)
if numel(entries) > 1
    names = arrayfun(@il_tag, entries, 'UniformOutput', false);   % Signal #1, ch1
else
    names = {entries.name};
end
q = {};
for k = 1:numel(entries)
    q = union(q, entries(k).values.Quantity, 'stable');
end
data = cell(numel(q), 1 + numel(entries));
data(:, 1) = q(:);
for k = 1:numel(entries)
    for r = 1:numel(q)
        i_q = find(strcmp(entries(k).values.Quantity, q{r}), 1);
        if ~isempty(i_q)
            data{r, k + 1} = entries(k).values.Value(i_q);
        end
    end
end
end

function n = il_export(obj, folder, name)
% one figure, axes or window to a file (vector content for a PDF); returns 1
opts = {};
if endsWith(name, '.pdf')
    opts = {'ContentType', 'vector'};
end
exportgraphics(obj, il_free_path(fullfile(folder, name)), opts{:});
n = 1;
end

function [xe, ye] = il_map_edges(t, f)
% outer edges of a map with a uniform time step and log-spaced frequencies, for a
% texture: on a log axis each row of the texture then covers its own band
dt = t(end) - t(end-1);
r = sqrt(f(end) / f(end-1));
xe = [t(1) - dt/2, t(end) + dt/2];
ye = [f(1) / r, f(end) * r];
end

function il_fit_axes(box, ax, bottom)
% the plot area of ax at fixed margins from the edges of its box (pixels):
% room for the ticks and the label of y on the left, for the colorbar and its
% label on the right (kept on all plots so their time axes line up), for the
% title on top and the time ticks and label below
m = [80 100 48 30];                                   % left, right, bottom, top
if nargin > 2
    m(3) = bottom;                                    % a plot with no time label under it
end
p = box.InnerPosition;
w = max(p(3) - m(1) - m(2), 20);
h = max(p(4) - m(3) - m(4), 20);
ax.InnerPosition = [m(1), m(3), w, h];
cb = ax.Colorbar;
if ~isempty(cb)                                       % placed by hand: in its own place it takes width from the plot
    cb.Units = 'pixels';
    cb.Position = [m(1) + w + 12, m(3), 14, h];
end
end

function [go, ask_again] = il_confirm_remove(parent, name)
% asks whether to remove the signal name; its tick box stops the question for
% the rest of the session (closing the window cancels and changes nothing)
go = false;
ask_again = true;
d = uifigure('Name', 'Remove signal', 'WindowStyle', 'modal', 'Tag', 'SQAT_GUI_remove', ...
    'Position', [parent.Position(1) + 200, parent.Position(2) + 200, 420, 140]);
if isprop(parent, 'Theme') && ~isempty(parent.Theme)
    d.Theme = parent.Theme;
end
g = uigridlayout(d, [3 3]);
g.RowHeight = {'1x', 22, 30};
g.ColumnWidth = {'1x', 90, 90};
q = uilabel(g, 'Text', sprintf('Remove %s and its results?', name), 'WordWrap', 'on');
q.Layout.Column = [1 3];
tick = uicheckbox(g, 'Text', 'Do not ask again', 'Tag', 'remove_do_not_ask');
tick.Layout.Column = [1 3];
uilabel(g, 'Text', '');
uibutton(g, 'Text', 'Remove', 'Tag', 'remove_ok', 'ButtonPushedFcn', @(~, ~) answer(true));
uibutton(g, 'Text', 'Cancel', 'Tag', 'remove_cancel', 'ButtonPushedFcn', @(~, ~) delete(d));
uiwait(d);

    function answer(tf)
        go = tf;
        ask_again = ~tick.Value;
        delete(d);
    end
end

function il_delete_timer(t)
if ~isempty(t) && isvalid(t)
    stop(t);
    delete(t);
end
end

function tf = il_has_theme()
% the theme function, and with it the dark theme, came in R2025a
tf = exist('theme', 'file') > 0;
end

function tf = il_is_member(value, list)
% true when value is a non-empty text found in the cell array list
tf = (ischar(value) || isstring(value)) && strlength(string(value)) > 0 ...
    && any(strcmp(list, value));
end

function p = il_default_params(e)
% the default parameters of metric e
p = struct();
for q = e.params
    p.(q.name) = q.value;
end
end

function [id, suffix] = il_split_key(key)
% Metric_id#2 -> Metric_id and '#2'; the first analysis of a metric has no suffix
k = strfind(key, '#');
if isempty(k)
    id = key;
    suffix = '';
else
    id = key(1:k-1);
    suffix = key(k:end);
end
end

function label = il_key_label(S, key)
% the name of an analysis (#2 Loudness (ISO 532-1)), as the results of the run carry it
k = find(strcmp({S.metric}, key), 1);
label = key;
if ~isempty(k)
    label = S(k).label;
end
end

function txt = il_param_summary(metrics, a)
% the parameters of analysis a in short, e.g. free field, time-varying, skip 0.5 s
e = metrics(strcmp({metrics.id}, a.id));
parts = {};
for q = e.params
    v = a.p.(q.name);
    if strcmp(q.type, 'choice')
        k = find(cellfun(@(o) isequal(o, v), q.options(:, 2)), 1);
        parts{end+1} = il_short(q.options{k, 1}); %#ok<AGROW>
    else
        unit = regexp(q.label, '\((\w+)\)', 'tokens', 'once');   % Time skip (s) -> s
        if isempty(unit)
            unit = {''};
        end
        prefix = struct('time_skip', 'skip', 'dt', 'dt', 'threshold', 'threshold');
        if isfield(prefix, q.name)
            name = prefix.(q.name);
        else
            name = q.name;
        end
        parts{end+1} = strtrim(sprintf('%s %g %s', name, v, unit{1})); %#ok<AGROW>
    end
end
txt = strjoin(parts, ', ');
end

function s = il_short(option)
% the short name of an option in the summary of the parameters
names = {'Free field', 'free field'; 'Diffuse field', 'diffuse field'; 'Free-frontal', 'free-frontal'; ...
         'Diffuse', 'diffuse'; 'Stationary', 'stationary'; 'Time-varying', 'time-varying'; ...
         'DIN 45692', 'DIN 45692'; 'Aures', 'Aures'; 'von Bismarck', 'von Bismarck'};
k = find(strcmp(names(:, 1), option), 1);
s = option;
if ~isempty(k)
    s = names{k, 2};
end
end

function u = il_level_unit(w)
% the unit of a level with frequency weighting w: dB SPL unweighted, dBA and
% dBC weighted (the convention of Greco's thesis: the weighting says it is SPL)
u = 'dB SPL';
if ~strcmpi(w, 'Z')
    u = ['dB' upper(w)];
end
end

function u = il_unit(id, q)
% the unit of quantity q of metric id, as the header of the metric states it
switch q
    case 'EPNL',  u = 'EPNdB'; return
    case 'PNLM',  u = 'PNdB';  return
    case 'PNLTM', u = 'TPNdB'; return
    case 'time',  u = 's';     return
    case {'N_ratio', 'ScalarPA'}, u = '-'; return
end
if strcmp(id, 'Do_SLM')
    u = 'dB SPL';                                  % LZeq, LZFmax, TOB: unweighted
    if numel(q) > 1 && ismember(q(2), 'AC')
        u = il_level_unit(q(2));                   % LAeq, LCFmax: dBA, dBC
    end
    return
end
if contains(q, 'Level')
    u = 'phon';
    return
end
units = struct('Loudness_ISO532_1', 'sone', 'Loudness_ECMA418_2', 'sone_HMS', ...
    'Sharpness_DIN45692', 'acum', 'Roughness_Daniel1997', 'asper', 'Roughness_ECMA418_2', 'asper', ...
    'FluctuationStrength_Osses2016', 'vacil', 'Tonality_Aures1985', 't.u.', 'Tonality_ECMA418_2', 'tu_HMS');
u = '-';
if isfield(units, id)
    u = units.(id);
end
end

function txt = il_param_text(e, p)
% the parameters in full: Sound field: Free field; Method: Time-varying; Time skip: 0.5 s
parts = {};
for q = e.params
    [name, unit] = il_label_unit(q.label);
    v = p.(q.name);
    if strcmp(q.type, 'choice')
        k = find(cellfun(@(o) isequal(o, v), q.options(:, 2)), 1);
        parts{end+1} = sprintf('%s: %s', name, q.options{k, 1}); %#ok<AGROW>
    else
        parts{end+1} = strtrim(sprintf('%s: %g %s', name, v, unit)); %#ok<AGROW>
    end
end
txt = strjoin(parts, '; ');
end

function path = il_free_path(path)
% the path, or path_2, path_3, ... when a file of that name is already there
[folder, base, ext] = fileparts(path);
k = 1;
while isfile(path)
    k = k + 1;
    path = fullfile(folder, sprintf('%s_%d%s', base, k, ext));
end
end

function v = il_sqat_version()
% the release of SQAT: its tag (v1.3), or the tag and the commits since it in
% a git checkout (v1.3 + 299 commits (db61b5c)); without git, the version of
% citation.cff (a download of a release)
root = fileparts(fileparts(mfilename('fullpath')));
[status, out] = system(sprintf('git -C "%s" describe --tags --long --match "v*"', root));
tok = regexp(strtrim(out), '^(.*)-(\d+)-g([0-9a-f]+)$', 'tokens', 'once');
if status == 0 && ~isempty(tok)
    v = tok{1};
    if ~strcmp(tok{2}, '0')
        v = sprintf('%s + %s commits (%s)', tok{1}, tok{2}, tok{3});
    end
    return
end
v = 'unknown';
cff = fullfile(root, 'citation.cff');
if isfile(cff)
    k = regexp(fileread(cff), '(?m)^version:\s*(\S+)', 'tokens', 'once');
    if ~isempty(k)
        v = k{1};
    end
end
end

function v = il_cal_of(f, channel)
% the full-scale level (dB SPL) of the channel of a result: '1', '2', ... or
% 'Binaural', which gets the value of the channels when they share it (NaN if not)
c = str2double(channel);
if ~isnan(c)
    v = f.dBFS(min(c, numel(f.dBFS)));
elseif all(f.dBFS == f.dBFS(1))
    v = f.dBFS(1);
else
    v = NaN;
end
end

function out = il_if(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function [name, unit] = il_label_unit(label)
% Time skip (s) -> 'Time skip' and 's'
tok = regexp(label, '^(.*?)\s*\(([^)]*)\)$', 'tokens', 'once');
if isempty(tok)
    name = label;
    unit = '';
else
    name = tok{1};
    unit = tok{2};
end
end

function items = il_channel_items(n)
% the channels a file of n channels offers: each one, and All when there are several
items = arrayfun(@num2str, 1:n, 'UniformOutput', false);
if n > 1
    items{end+1} = 'All';
end
end

function c = il_default_channel(n)
% a stereo file is analysed in all its channels, a mono one in its only channel
c = '1';
if n > 1
    c = 'All';
end
end

function T = il_empty_results()
% one row per value, with what is needed to reproduce it: the signal, the
% analysis, the unit, the calibration (dB SPL of full scale), the parameters and the path
T = table('Size', [0 11], 'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'double', ...
    'cell', 'double', 'cell', 'cell'}, ...
    'VariableNames', {'Signal', 'File', 'Analysis', 'Metric', 'Channel', 'Quantity', 'Value', ...
    'Unit', 'Cal_dB_SPL', 'Parameters', 'Path'});
end

function S = il_empty_store()
S = struct('file', {}, 'name', {}, 'id', {}, 'metric', {}, 'number', {}, 'label', {}, 'channel', {}, 'analyses', {}, 'values', {});
end

function en = il_entry_of(f, e, OUT, label, channel, n_channels)
% the analyses and the single values that OUT holds for one channel
en = struct('file', f.path, 'name', f.name, 'id', f.id, 'metric', e.key, 'number', e.number, ...
    'label', e.label, 'channel', label, ...
    'analyses', SQAT_GUI_extract(OUT, e.id, channel), ...
    'values', SQAT_GUI_single_values(OUT, channel, n_channels));
end

function t = il_tag(en)
% the short name of a signal and channel in the plots: Signal #1, ch1
if strcmp(en.channel, 'Binaural')
    c = 'binaural';
else
    c = ['ch' en.channel];
end
t = sprintf('Signal #%d, %s', en.id, c);
end

function [items, data] = il_analysis_items(entries)
% the analyses on offer for the signals of a window: the SQAT figure and all
% the analyses for one signal, then the analyses that every signal holds
% (they can be overlaid or set side by side), then the statistics
items = {};
data = {};
if numel(entries) == 1
    items = {'SQAT figure', 'All analyses'};
    data = {'sqat', 'all'};
end
ids = {entries(1).analyses.id};
labels = {entries(1).analyses.label};
for k = 2:numel(entries)
    keep = ismember(ids, {entries(k).analyses.id});
    ids = ids(keep);
    labels = labels(keep);
end
items = [items, labels, {'Statistics'}];
data = [data, ids, {'stats'}];
end
