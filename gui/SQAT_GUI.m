function varargout = SQAT_GUI(files, varargin)
% function fig = SQAT_GUI(files, varargin)
%
%   Graphical interface to the metrics of SQAT, laid out as the interface of
%   pySQAT. It loads .wav files into a list of signals (a tick marks a signal
%   for use, the x removes it; each signal has its own channel and dBFS),
%   runs the ticked metrics with the chosen parameters on the ticked signals,
%   and lists their single values. The channel to analyse is one channel of
%   the file or All (the default for a stereo file): the ECMA-418-2 metrics take a stereo pair in one call and
%   return the left, the right and (except the tonality) the combined
%   binaural result, so a pair runs once.
%
%   A graphs window plots one metric of the ticked signals. The analysis is
%   chosen in the window (a time series, a profile over the critical bands, a
%   map of band against time, or the statistics): the signals are overlaid,
%   and the maps sit side by side on one colour scale. For one signal the
%   window also offers the figure that the SQAT function draws (with the
%   heat colour scale of SQAT_GUI_colormap_heat) and all the analyses at once. Each signal of the list
%   has a number, and the plots tag a curve as Signal #1, ch1. The channel All
%   plots every channel of every signal, to compare a mono with a stereo signal.
%   Pin keeps a window with its signals and its results, so that a later run
%   leaves it as it is and Open Graphs Window opens another one to compare
%   with. A player window shows the waveform and the spectrogram, and the
%   results go to a spreadsheet.
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
% AI disclosure: code development in September 2026 assisted
% by Claude Opus 5 (Anthropic). All codes were verified by
% the authors.
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
green = [0.13 0.55 0.37];

%% State
metrics = SQAT_GUI_metrics;
% the analyses of the list: a metric and its parameters; the same metric can
% be there more than once, and its second entry is keyed Metric_id#2
analyses = struct('key', 'Loudness_ISO532_1', 'id', 'Loudness_ISO532_1', ...
    'p', il_default_params(metrics(strcmp({metrics.id}, 'Loudness_ISO532_1'))));
loaded = struct('path', {}, 'name', {}, 'nch', {}, 'fs', {}, 'marked', {}, 'id', {}, 'channel', {}, 'dBFS', {});
next_id = 1;                                          % the number of the next signal loaded
active_idx = 0;                                       % the signal on screen in the player
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
spec_view = [];                                       % full-file enhanced map, kept to redraw a zoomed excerpt
spec_busy = false;                                    % the spectrogram limits are being set by the code
spec_timer = [];                                      % waits for the zoom to settle before recomputing
theme_style = 'dark';
graph_figs = gobjects(0);                              % the graphs windows
last_metric = '';                                     % the metric the last graphs window showed
win_wave = [];
ax_wave = [];
ax_spec = [];
btn_play = [];
run_settings = struct('signals', loaded, 'analyses', analyses);      % of the last run
stop_requested = false;                                             % the Stop button or the dialog
dlg = [];                                                           % progress dialog of a run
cache = struct('file', {}, 'metric', {}, 'figs', {});               % SQAT figures, hidden
cmap = SQAT_GUI_colormap_heat(256);

%% Main window
fig = uifigure('Name', 'SQAT: Sound Quality Analysis Toolbox', ...
    'Position', [40 30 1460 900], 'Visible', opts.Visible, 'Tag', 'SQAT_GUI', ...
    'CloseRequestFcn', @on_close, 'CreateFcn', '');   % skips a user default CreateFcn
main = uigridlayout(fig, [3 2]);
main.RowHeight = {44, '1x', 30};
main.ColumnWidth = {430, '1x'};

top = uigridlayout(main, [1 3]);
top.Layout.Row = 1; top.Layout.Column = [1 2];
top.Padding = [0 0 0 0];
top.ColumnWidth = {70, 320, '1x'};
img_logo = uiimage(top, 'ImageSource', fullfile(dir_logos, 'logo_white.png'), 'Tag', 'logo');
uilabel(top, 'Text', 'Sound Quality Analysis Toolbox', 'FontSize', 16, 'FontWeight', 'bold');
uilabel(top, 'Text', '');

left = uigridlayout(main, [4 1]);
left.Layout.Row = 2; left.Layout.Column = 1;
left.Padding = [0 0 0 0];
left.RowHeight = {28, 170, 28, '1x'};
sh = uigridlayout(left, [1 3]);
sh.Padding = [0 0 0 0];
sh.ColumnWidth = {'1x', 100, 130};
uilabel(sh, 'Text', 'SIGNALS (tick to use, x to remove)', 'FontWeight', 'bold');
lbl_files = uilabel(sh, 'Text', 'No files loaded', 'Tag', 'file_count', 'HorizontalAlignment', 'right');
uibutton(sh, 'Text', 'Open WAV files...', 'Tag', 'load_files', 'ButtonPushedFcn', @on_load_files);
tbl_signals = uitable(left, 'Tag', 'signals_table', 'RowName', {}, ...
    'ColumnName', {'', 'Signal', 'Channel', 'dBFS', ''}, ...
    'ColumnWidth', {28, 180, 70, 60, 28}, 'ColumnEditable', [true false true true false], ...
    'Tooltip', ['Channel: the one to analyse, or All (the ECMA-418-2 metrics take a stereo file ' ...
                'as a binaural pair, in one call). dBFS: dB SPL of a full-scale amplitude ' ...
                '(94: full scale 1.0 is 1 Pa)'], ...
    'CellEditCallback', @on_signal_edited, 'CellSelectionCallback', @on_signal_selected);
ah = uigridlayout(left, [1 2]);
ah.Padding = [0 0 0 0];
ah.ColumnWidth = {'1x', 130};
uilabel(ah, 'Text', 'ANALYSES (gear: parameters)', 'FontWeight', 'bold');
uibutton(ah, 'Text', '+ Add', 'Tag', 'add_analysis', 'ButtonPushedFcn', @on_add_analysis, ...
    'Tooltip', 'Adds a copy of the last analysis, to compare the same metric with other parameters');
analysis_list = uigridlayout(left, [1 5], 'Scrollable', 'on', 'Tag', 'analysis_list');
analysis_list.ColumnWidth = {30, 175, '1x', 36, 28};
analysis_list.Padding = [0 0 0 0];

right = uigridlayout(main, [3 1]);
right.Layout.Row = 2; right.Layout.Column = 2;
right.Padding = [0 0 0 0];
right.RowHeight = {56, 64, '1x'};

opt_panel = uipanel(right, 'Title', 'OPTIONS');
og = uigridlayout(opt_panel, [1 5]);
og.ColumnWidth = {150, 105, 105, '1x', 80};
og.Padding = [6 4 6 4];
cb_show = uicheckbox(og, 'Text', 'Show plots after run', 'Tag', 'show_plots');
cb_split = uicheckbox(og, 'Text', 'Split figures', 'Tag', 'split_figures', ...
    'Tooltip', 'One tab (and one saved file) per panel of a figure', ...
    'ValueChangedFcn', @on_graph_option);
cb_save = uicheckbox(og, 'Text', 'Save figures', 'Tag', 'save_figures');
ed_folder = uieditfield(og, 'text', 'Value', pwd, 'Tag', 'figures_folder', ...
    'Tooltip', 'Folder for the saved figures');
uibutton(og, 'Text', 'Browse...', 'ButtonPushedFcn', @on_browse);

act_panel = uipanel(right, 'Title', 'ACTIONS');
ag = uigridlayout(act_panel, [1 6]);
ag.ColumnWidth = {'1.4x', '1x', '1x', '1x', '1x', 36};
ag.Padding = [6 4 6 4];
uibutton(ag, 'Text', 'Run Analysis', 'Tag', 'run', 'FontWeight', 'bold', ...
    'BackgroundColor', green, 'FontColor', [1 1 1], 'ButtonPushedFcn', @on_run);
btn_stop = uibutton(ag, 'Text', 'Stop', 'Tag', 'stop_run', 'Enable', 'off', ...
    'Tooltip', 'Ends the run after the metric being computed', 'ButtonPushedFcn', @on_stop_run);
uibutton(ag, 'Text', 'Open Graphs Window', 'Tag', 'open_graphs', 'ButtonPushedFcn', @on_open_graphs);
uibutton(ag, 'Text', 'Waveform / Play', 'Tag', 'open_waveform', 'ButtonPushedFcn', @on_open_waveform);
uibutton(ag, 'Text', 'Export results...', 'Tag', 'export', 'ButtonPushedFcn', @on_export);
btn_theme = uibutton(ag, 'Text', char(9788), 'FontSize', 18, 'Tag', 'theme', ...
    'Tooltip', 'Light theme', 'ButtonPushedFcn', @on_theme);   % a sun, or a moon in the light theme

tabs = uitabgroup(right);
tab_console = uitab(tabs, 'Title', 'Console output');
console = uitextarea(uigridlayout(tab_console, [1 1]), 'Value', {''}, 'Editable', 'off', ...
    'Tag', 'console', 'FontName', 'Monospaced');
tab_results = uitab(tabs, 'Title', 'Results');
tbl = uitable(uigridlayout(tab_results, [1 1]), 'Data', results, 'Tag', 'results_table');

status_bar = uigridlayout(main, [1 2]);
status_bar.Layout.Row = 3; status_bar.Layout.Column = [1 2];
status_bar.Padding = [0 0 0 0];
status_bar.ColumnWidth = {'1x', 300};
lbl_status = uilabel(status_bar, 'Text', 'Ready', 'Tag', 'status');
gauge = uigauge(status_bar, 'linear', 'Tag', 'progress', 'Limits', [0 100], 'Value', 0, ...
    'MajorTicks', [], 'MinorTicks', []);

%% Start
apply_theme();
add_files(files);
refresh_analyses();
setappdata(fig, 'sqat_set_analyses', @set_analyses);   % the list from metric ids, for the tests
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

    function on_signal_edited(~, event)
        % the tick, the channel or the dBFS of a signal
        k = event.Indices(1);
        switch event.Indices(2)
            case 1
                loaded(k).marked = logical(event.NewData);
                refresh_windows();
                return
            case 3
                c = char(event.NewData);
                if ~il_is_member(c, il_channel_items(loaded(k).nch))
                    write_log(sprintf('%s has %d channel(s); channel %s is not there.', loaded(k).name, loaded(k).nch, c));
                    refresh_signals();
                    return
                end
                loaded(k).channel = c;
            case 4
                if ~isfinite(event.NewData)
                    refresh_signals();
                    return
                end
                loaded(k).dBFS = event.NewData;
        end
        refresh_signals();                   % the table shows what the signal holds
        if k == active_idx && il_is_open(win_wave)
            draw_waveform_window();
        end
    end

    function on_signal_selected(~, event)
        if isempty(event.Indices)
            return
        end
        k = event.Indices(1, 1);
        if event.Indices(1, 2) == 5
            remove_signal(k);                % the x at the right of the row
            return
        end
        if k ~= active_idx
            active_idx = k;
            show_active();
            refresh_windows();
        end
    end

    function on_graph_option(~, ~)
        % the split option changes how the SQAT figure is laid out
        for w = open_graph_windows()
            if strcmp(findobj(w, 'Tag', 'graph_analysis').Value, 'sqat')
                draw_window(w);
            end
        end
    end

    %% The list of analyses

    function refresh_analyses()
        % one row per analysis: its number, the metric, its parameters in short, the gear and the x
        delete(analysis_list.Children);
        close_params();
        n = numel(analyses);
        analysis_list.RowHeight = repmat({26}, 1, max(n, 1));
        for k = 1:n
            a = analyses(k);
            uilabel(analysis_list, 'Text', sprintf('#%d', k), 'Tag', sprintf('analysis_number_%d', k));
            uidropdown(analysis_list, 'Items', {metrics.label}, 'ItemsData', {metrics.id}, ...
                'Value', a.id, 'Tag', sprintf('analysis_metric_%d', k), ...
                'ValueChangedFcn', @(src, ~) on_analysis_metric(k, src.Value));
            uilabel(analysis_list, 'Text', il_param_summary(metrics, a), 'Tag', sprintf('analysis_summary_%d', k));
            uibutton(analysis_list, 'Text', char(9881), 'FontSize', 16, 'Tag', sprintf('analysis_params_%d', k), ...
                'Tooltip', 'Parameters of this analysis', 'ButtonPushedFcn', @(~, ~) on_analysis_params(k));
            uibutton(analysis_list, 'Text', 'x', 'Tag', sprintf('analysis_remove_%d', k), ...
                'Tooltip', 'Removes this analysis', 'ButtonPushedFcn', @(~, ~) on_remove_analysis(k));
        end
    end

    function set_analyses(ids)
        % the list holds these metrics, with their default parameters
        analyses = analyses([]);
        for k = 1:numel(ids)
            e = metrics(strcmp({metrics.id}, ids{k}));
            analyses(k) = struct('key', '', 'id', e.id, 'p', il_default_params(e));
        end
        assign_keys();
        refresh_analyses();
    end

    function assign_keys()
        % the first analysis of a metric is keyed by its id, the next ones id#2, id#3, ...
        for k = 1:numel(analyses)
            n = nnz(strcmp({analyses(1:k).id}, analyses(k).id));
            analyses(k).key = analyses(k).id;
            if n > 1
                analyses(k).key = sprintf('%s#%d', analyses(k).id, n);
            end
        end
    end

    function on_add_analysis(~, ~)
        % a copy of the last analysis, whose parameters are then changed to compare
        if isempty(analyses)
            e = metrics(strcmp({metrics.id}, 'Loudness_ISO532_1'));
            analyses = struct('key', '', 'id', e.id, 'p', il_default_params(e));
        else
            analyses(end+1) = analyses(end);
        end
        assign_keys();
        refresh_analyses();
    end

    function on_analysis_metric(k, id)
        e = metrics(strcmp({metrics.id}, id));
        analyses(k).id = id;
        analyses(k).p = il_default_params(e);
        assign_keys();
        refresh_analyses();
    end

    function on_remove_analysis(k)
        analyses(k) = [];
        assign_keys();
        refresh_analyses();
    end

    function on_analysis_params(k)
        % a small window with the parameters of analysis k; a change applies at once
        close_params();
        a = analyses(k);
        e = metrics(strcmp({metrics.id}, a.id));
        n = numel(e.params);
        w = uifigure('Name', sprintf('Parameters: #%d %s', k, e.label), ...
            'Position', [fig.Position(1) + 440, fig.Position(2) + 300, 420, 60 + 34 * max(n, 1)], ...
            'Visible', fig.Visible, 'Tag', 'SQAT_GUI_params', 'CreateFcn', '');
        g = uigridlayout(w, [max(n, 1) + 1, 2]);
        g.ColumnWidth = {'1x', 200};
        g.RowHeight = [repmat({26}, 1, max(n, 1)), {28}];
        if n == 0
            uilabel(g, 'Text', 'This metric has no parameters.');
            uilabel(g, 'Text', '');
        end
        for k_par = 1:n
            q = e.params(k_par);
            uilabel(g, 'Text', q.label);
            current = a.p.(q.name);
            if strcmp(q.type, 'choice')
                c = uidropdown(g, 'Items', q.options(:, 1)', 'ItemsData', q.options(:, 2)', 'Value', current);
            else
                c = uieditfield(g, 'numeric', 'Value', current);
            end
            c.Tag = ['param_' q.name];
            c.ValueChangedFcn = @(src, ~) set_param(a.key, q.name, src.Value);
        end
        uilabel(g, 'Text', '');
        uibutton(g, 'Text', 'Close', 'Tag', 'params_close', 'ButtonPushedFcn', @(~, ~) delete(w));
        if exist('theme', 'file')
            theme(w, theme_style);
        end
    end

    function set_param(key, name, value)
        k = find(strcmp({analyses.key}, key), 1);
        if ~isempty(k)
            analyses(k).p.(name) = value;
            set(findobj(analysis_list, 'Tag', sprintf('analysis_summary_%d', k)), ...
                'Text', il_param_summary(metrics, analyses(k)));
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
        e.number = sprintf('#%d', k);             % the row of the analysis in the list
        e.label = sprintf('#%d %s', k, e.label);
        e.p = a.p;
    end

    function show_results()
        % the results as one list, and one tab per signal (Results #1, #2, ...) with its own rows
        tbl.Data = results;
        delete(findobj(tabs, 'Tag', 'results_signal'));
        for f = loaded
            rows = strcmp(results.File, f.name);
            if any(rows)
                t = uitab(tabs, 'Title', sprintf('Results #%d', f.id), 'Tag', 'results_signal');
                uitable(uigridlayout(t, [1 1]), 'Data', results(rows, 2:end), 'RowName', {});
            end
        end
    end

    function on_browse(~, ~)
        d = uigetdir(ed_folder.Value, 'Folder for the figures');
        focus_gui();
        if ~isequal(d, 0)
            ed_folder.Value = d;
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
            write_log('No metrics in the list of analyses. Add one with + Add.');
            return
        end
        sel = {analyses.key};
        save_figs = cb_save.Value;
        split = cb_split.Value;
        folder = strtrim(ed_folder.Value);
        if save_figs && ~isfolder(folder)
            write_log(['ERROR: the folder for the figures does not exist: ' folder]);
            save_figs = false;
        end
        show = save_figs;                        % the SQAT functions draw the figures to save
        k_active = find(use == active_idx, 1);   % the signal on screen is analysed first,
        if isempty(k_active)                     % and its figure is drawn in the analysis
            k_active = 1;                        % call, so that the metric runs once
        end
        files_order = [use(k_active), use(use ~= use(k_active))];
        active_path = loaded(files_order(1)).path;
        clear_cache();
        run_settings = struct('signals', loaded(use), 'analyses', analyses);
        close_params();

        new_results = il_empty_results();
        new_store = il_empty_store();
        stereo_sel = ismember({analyses.id}, {metrics([metrics.stereo]).id});
        n_total = 0;
        for i = files_order
            n_ch = numel(channel_list(loaded(i)));
            joint = n_ch == 2;
            n_total = n_total + nnz(~(stereo_sel & joint)) * n_ch + nnz(stereo_sel & joint);
        end
        n_done = 0;
        n_errors = 0;
        set_progress(0);
        stop_requested = false;
        btn_stop.Enable = 'on';
        stop_off = onCleanup(@() set(btn_stop, 'Enable', 'off')); %#ok<NASGU>
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
            ids_joint = sel(stereo_sel & joint);
            ids_single = sel(~(stereo_sel & joint));
            entries = il_empty_store();          % of this file, in the order of the calls
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
                suffix = '';
                if numel(cl) > 1
                    suffix = sprintf('_ch%d', c);
                end
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
                            keep_figures(new_figs, f, e.key, save_figs, split, folder, suffix, first);
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
                            keep_figures(new_figs, f, e.key, save_figs, split, folder, '', true);
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
                for k_e = find(strcmp({entries.metric}, sel{j}))
                    en = entries(k_e);
                    new_store(end+1) = en; %#ok<AGROW>
                    n = height(en.values);
                    new_results = [new_results; [table(repmat({f.name}, n, 1), repmat({en.number}, n, 1), ...
                        repmat({il_split_key(en.metric)}, n, 1), repmat({en.channel}, n, 1), ...
                        'VariableNames', {'File', 'Analysis', 'Metric', 'Channel'}), en.values]]; %#ok<AGROW>
                end
            end
            if stop_requested
                break                            % what ran so far is kept
            end
        end

        results = new_results;
        store = new_store;
        show_results();
        ran = sel(ismember(sel, {store.metric}));
        if stop_requested
            msg = sprintf('Stopped: %d value(s) from %d of the %d analysis step(s) in %.1f s', ...
                height(results), n_done, n_total, toc(t_start));
        else
            set_progress(100);
            msg = sprintf('Done: %d value(s) from %d file(s) and %d metric(s) in %.1f s', ...
                height(results), numel(use), numel(sel), toc(t_start));
        end
        if n_errors > 0
            msg = sprintf('%s, %d error(s)', msg, n_errors);
        end
        lbl_status.Text = msg;
        write_log([msg '.']);
        if cb_show.Value && ~isempty(ran)
            on_open_graphs();
        else
            refresh_graph_windows();
        end
    end

    function on_export(~, ~)
        if height(results) == 0
            write_log('No results to export. Run an analysis first.');
            return
        end
        [f, p] = uiputfile({'*.xlsx', 'Excel workbook (*.xlsx)'; '*.csv', 'CSV file (*.csv)'}, ...
            'Export results', 'SQAT_results.xlsx');
        focus_gui();
        if isequal(f, 0)
            return
        end
        try
            SQAT_GUI_export(results, fullfile(p, f));
            write_log(['Results exported to ' fullfile(p, f)]);
        catch err
            write_log(['ERROR exporting the results: ' err.message]);
        end
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

    function on_open_waveform(~, ~)
        if isempty(loaded)
            write_log('No file loaded.');
            return
        end
        if ~il_is_open(win_wave)
            win_wave = uifigure('Name', 'Waveform', 'Position', [140 140 1100 700], ...
                'Visible', fig.Visible, 'Tag', 'SQAT_GUI_waveform', 'CloseRequestFcn', @on_close_waveform, ...
                'CreateFcn', '', 'KeyPressFcn', @on_wave_key);
            gw = uigridlayout(win_wave, [4 1]);
            gw.RowHeight = {30, 30, '1x', '1.2x'};
            hw = uigridlayout(gw, [1 9]);
            hw.Padding = [0 0 0 0];
            hw.ColumnWidth = {80, 80, 60, 100, 110, 140, 80, 60, '1x'};
            btn_play = uibutton(hw, 'Text', 'Play', 'Tag', 'play', 'ButtonPushedFcn', @on_play, ...
                'Tooltip', 'Space plays and pauses; a click on the waveform or the spectrogram moves the playhead');
            uibutton(hw, 'Text', 'Stop', 'Tag', 'stop', 'ButtonPushedFcn', @on_stop);
            uicheckbox(hw, 'Text', 'Loop', 'Value', true, 'Tag', 'loop', 'ValueChangedFcn', @on_loop_changed, ...
                'Tooltip', 'Starts again at the end of the file, or of the filter box when there is one');
            uibutton(hw, 'state', 'Text', 'Draw filter', 'Tag', 'draw_box', 'ValueChangedFcn', @on_draw_box, ...
                'Tooltip', 'Drag a box on the spectrogram, or click two opposite corners');
            uibutton(hw, 'Text', 'Clear filters', 'Tag', 'clear_boxes', 'ButtonPushedFcn', @on_clear_boxes);
            uidropdown(hw, 'Items', {'Filter: loop only', 'Filter: isolate', 'Filter: remove'}, ...
                'ItemsData', {'loop', 'isolate', 'remove'}, 'Value', 'loop', 'Tag', 'box_mode', ...
                'ValueChangedFcn', @on_processing_changed, ...
                'Tooltip', ['Loop only: the sound is not changed, and the loop runs inside the box. ' ...
                            'Isolate: only what is inside the box plays. Remove: what is inside the box is taken out.']);
            uilabel(hw, 'Text', 'Weighting:', 'HorizontalAlignment', 'right');
            uidropdown(hw, 'Items', {'Z', 'A', 'C'}, 'Value', 'Z', 'Tag', 'wave_weighting', ...
                'ValueChangedFcn', @on_weighting_changed, ...
                'Tooltip', 'Frequency weighting of the sound and of the spectrogram (IEC 61672-1)');
            uilabel(hw, 'Text', '');
            sw = uigridlayout(gw, [1 11]);
            sw.Padding = [0 0 0 0];
            sw.ColumnWidth = {60, 150, 140, 80, 70, 80, 70, 110, 60, 100, '1x'};
            uilabel(sw, 'Text', 'Window:', 'HorizontalAlignment', 'right');
            uidropdown(sw, 'Items', {'Hann', 'Hamming', 'Rectangular', 'Blackman-Harris'}, ...
                'ItemsData', {'hann', 'hamming', 'rect', 'blackmanharris'}, 'Value', 'hann', ...
                'Tag', 'spec_window', 'ValueChangedFcn', @on_spec_option);
            uibutton(sw, 'Text', 'Import window...', 'Tag', 'import_window', 'ButtonPushedFcn', @on_import_window, ...
                'Tooltip', 'A .txt, .csv, .dat or .mat file with the samples of the window');
            uilabel(sw, 'Text', 'FFT degree:', 'HorizontalAlignment', 'right');
            uispinner(sw, 'Value', 10, 'Limits', [6 16], 'Step', 1, 'Tag', 'spec_degree', ...
                'RoundFractionalValues', 'on', 'ValueChangedFcn', @on_spec_option, ...
                'Tooltip', 'The FFT has 2^degree points (6 to 16); use the arrows');
            uilabel(sw, 'Text', 'Overlap (%):', 'HorizontalAlignment', 'right');
            uispinner(sw, 'Value', 50, 'Limits', [0 95], 'Step', 5, 'Tag', 'spec_overlap', ...
                'ValueChangedFcn', @on_spec_option, 'Tooltip', 'Overlap of the frames (0 to 95); use the arrows');
            uilabel(sw, 'Text', 'Enhanced STFT:', 'HorizontalAlignment', 'right');
            uiswitch(sw, 'slider', 'Items', {'Off', 'On'}, 'Value', 'Off', 'Tag', 'spec_enhanced', ...
                'ValueChangedFcn', @on_spec_enhanced, ...
                'Tooltip', ['Enhanced STFT (consensus): nine windows of 8 to 512 ms, reassigned and combined, ' ...
                            'so that no window has to be chosen. It replaces the window, the FFT degree and the overlap.']);
            uidropdown(sw, 'Items', {'Readable', 'Sharp'}, 'ItemsData', {'readable', 'sharp'}, 'Value', 'readable', ...
                'Tag', 'spec_enhanced_mode', 'Enable', 'off', 'ValueChangedFcn', @on_spec_enhanced, ...
                'Tooltip', 'Readable: smoothing of 4 ms and 1.45 % of the frequency, continuous lines. Sharp: 1 ms and 1 Hz, the thinnest lines.');
            uilabel(sw, 'Text', '');
            ax_wave = uiaxes(gw, 'Tag', 'waveform_axes', 'ButtonDownFcn', @on_wave_click);
            ax_spec = uiaxes(gw, 'Tag', 'spectrogram', 'ButtonDownFcn', @on_wave_click);
            ax_spec.XAxis.LimitsChangedFcn = @on_spec_limits;
            setappdata(win_wave, 'sqat_spec_zoom', @apply_spec_zoom);   % the recomputation, for the tests
            setappdata(win_wave, 'sqat_audio', @process_audio);   % what plays, for the tests
            setappdata(win_wave, 'sqat_play', @play_info);
            apply_theme();
        end
        draw_waveform_window();
        if strcmp(fig.Visible, 'on')
            figure(win_wave);
        end
    end

    function on_wave_key(~, event)
        if strcmp(event.Key, 'space')
            toggle_play();
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
            box_loop = ~isempty(boxes) && sample >= r1 && sample <= r2;   % inside the box: the box is the loop
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
            if findobj(win_wave, 'Tag', 'loop').Value
                buf = [first; repmat(rep, max(1, ceil(120 * wave_fs / numel(rep))), 1)];
            else
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
        draw_spectrogram();
    end

    function on_draw_box(src, ~)
        end_drag();
        box_corner = [];
        delete(findobj(ax_spec, 'Tag', 'box_corner'));
        if src.Value
            write_log('Draw filter: drag a box on the spectrogram, or click two opposite corners.');
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
        draw_spectrogram();
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
        draw_spectrogram();
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
        draw_spectrogram();
    end

    function on_close_waveform(~, ~)
        on_stop();
        stop_spec_timer();
        delete(win_wave);
    end

    function on_close(~, ~)
        if ~isempty(player)
            stop(player);
        end
        clear_cache();
        stop_spec_timer();
        if il_is_open(win_wave), delete(win_wave); end
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
                'channel', il_default_channel(info.NumChannels), 'dBFS', 94); %#ok<AGROW>
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
        results = results(~strcmp(results.File, f.name), :);
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
        % the table of signals, the count and the channels on offer
        n = numel(loaded);
        if n == 0
            tbl_signals.Data = table(false(0, 1), strings(0, 1), categorical(strings(0, 1)), zeros(0, 1), strings(0, 1));
            tbl_signals.RowName = {};
            lbl_files.Text = 'No files loaded';
        else
            tbl_signals.Data = table(logical([loaded.marked]'), string({loaded.name}'), ...
                categorical({loaded.channel}', il_channel_items(max([loaded.nch]))), ...
                [loaded.dBFS]', repmat("x", n, 1));
            tbl_signals.RowName = arrayfun(@(k) sprintf('#%d', k), [loaded.id], 'UniformOutput', false);
            if n == 1
                lbl_files.Text = '1 file loaded';
            else
                lbl_files.Text = sprintf('%d files loaded', n);
            end
        end
        removeStyle(tbl_signals);
        if active_idx > 0
            addStyle(tbl_signals, uistyle('FontWeight', 'bold'), 'row', active_idx);
        end
    end

    function show_active()
        % the row on screen is shown in bold
        removeStyle(tbl_signals);
        if active_idx > 0
            addStyle(tbl_signals, uistyle('FontWeight', 'bold'), 'row', active_idx);
        end
        if il_is_open(win_wave)
            draw_waveform_window();
        end
    end

    function f = active_file()
        f = loaded(active_idx);
    end

    function ch = channel_of_active()
        % the channel that the waveform and the player take
        ch = str2double(active_file().channel);
        if isnan(ch)
            ch = 1;                          % All
        end
        ch = min(ch, active_file().nch);
    end

    function refresh_windows()
        refresh_graph_windows();
        if il_is_open(win_wave)
            if isempty(loaded)
                cla(ax_wave);
                cla(ax_spec);
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
        g = uigridlayout(w, [2 1]);
        g.RowHeight = {30, '1x'};
        g.Padding = [8 8 8 8];
        bar = uigridlayout(g, [1 8]);
        bar.Padding = [0 0 0 0];
        bar.ColumnWidth = {50, 230, 60, 260, 60, 90, 70, '1x'};
        uilabel(bar, 'Text', 'Metric:', 'HorizontalAlignment', 'right');
        uidropdown(bar, 'Items', {}, 'Tag', 'graph_metric', 'ValueChangedFcn', @on_graph_control);
        uilabel(bar, 'Text', 'Analysis:', 'HorizontalAlignment', 'right');
        uidropdown(bar, 'Items', {}, 'Tag', 'graph_analysis', 'ValueChangedFcn', @on_graph_control, ...
            'Tooltip', ['SQAT figure and All analyses need one signal; with several signals ' ...
                        'the analyses that can be compared are on offer']);
        uilabel(bar, 'Text', 'Channel:', 'HorizontalAlignment', 'right');
        uidropdown(bar, 'Items', {'1'}, 'Tag', 'graph_channel', 'ValueChangedFcn', @on_graph_control);
        uibutton(bar, 'state', 'Text', 'Pin', 'Tag', 'graph_pin', 'ValueChangedFcn', @on_graph_pin, ...
            'Tooltip', ['Keeps this window with its signals and results, whatever runs next; ' ...
                        'Open Graphs Window then opens another one to compare with']);
        uilabel(bar, 'Text', '');
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
        dc = findobj(w, 'Tag', 'graph_channel');
        da = findobj(w, 'Tag', 'graph_analysis');
        body = findobj(w, 'Tag', 'graph_body');
        delete(body.Children);
        in_window = ismember({ws.file}, paths);
        ids = unique({ws(in_window).metric}, 'stable');   % the keys of the analyses, in the order of the run
        if isempty(ids)
            set([dm, dc, da], 'Items', {});
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
        chans = unique({ws(in_metric).channel}, 'stable');
        chans = [sort(chans(~strcmp(chans, 'Binaural'))), chans(strcmp(chans, 'Binaural'))];
        wanted = dc.Value;
        if numel(chans) > 1
            dc.Items = [chans, {'All'}];       % every channel of every signal, to compare them
        else
            dc.Items = chans;
        end
        if il_is_member(wanted, dc.Items)
            dc.Value = wanted;
        end
        chan = dc.Value;

        entries = il_empty_store();
        for k_p = 1:numel(paths)
            if strcmp(chan, 'All')
                k_e = find(in_metric & strcmp({ws.file}, paths{k_p}));
            else
                k_e = find(in_metric & strcmp({ws.file}, paths{k_p}) & strcmp({ws.channel}, chan), 1);
            end
            if ~isempty(k_e)
                entries = [entries, ws(k_e)]; %#ok<AGROW>
            end
        end
        n_missing = nnz(ismember(paths, {ws(in_metric).file})) - numel(entries);
        if ~strcmp(chan, 'All') && n_missing > 0
            write_log(sprintf('%d signal(s) have no channel %s of %s and are left out.', n_missing, chan, label));
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
    end

    function ok = show_sqat_figure(parent, id, path, may_run)
        % the figure that the SQAT function draws, copied into the window; a pinned
        % window does not run the metric again, since the settings may have changed
        f = loaded(strcmp({loaded.path}, path));
        label = il_key_label(store, id);
        k_c = find(strcmp({cache.file}, f.path) & strcmp({cache.metric}, id), 1);
        if isempty(k_c) || ~all(isvalid(cache(k_c).figs))
            if ~may_run
                write_log(sprintf(['The SQAT figure of %s for %s is gone after a new run; ' ...
                    'the window shows the first analysis.'], id, f.name));
                ok = false;
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
                ok = false;
                return
            end
            keep_figures(new_figs, f, id, false, false, '');
            k_c = numel(cache);
            lbl_status.Text = 'Ready';
        end
        delete(parent.Children);
        tg = uitabgroup(uigridlayout(parent, [1 1], 'Padding', 0));   % the grid sizes it on screen
        for k_f = 1:numel(cache(k_c).figs)
            src = cache(k_c).figs(k_f);
            if ~cb_split.Value
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
        for k = 1:numel(A)
            plot(ax, A(k).x, A(k).y);
        end
        hold(ax, 'off');
        y_all = vertcat(A.y);
        y_range = [min(y_all) max(y_all)];
        y_ref = max(abs(y_range));
        if all(isfinite(y_range)) && y_ref > 0 && diff(y_range) <= 1e-3 * y_ref
            % a constant result: show it at +/-5 %, away from its rounding noise
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
        title(ax, sprintf('%s: %s, channel %s', A(1).label, what, chan), 'Interpreter', 'none');
        if numel(A) > 1
            legend(ax, names, 'Interpreter', 'none', 'Location', 'best');
        end
    end

    function draw_stats(parent, entries)
        % the single values of the signals, one column each
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
        uitable(uigridlayout(parent, [1 1]), 'Data', data, 'Tag', 'graph_stats', 'RowName', {}, ...
            'ColumnName', [{'Quantity'}, names], 'ColumnEditable', false);
    end

    function draw_waveform_window()
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
        win_wave.Name = sprintf('Waveform: %s, channel %d', f.name, ch);
        t_end = numel(x) / fs;
        step = max(1, ceil(numel(x) / 2e6));   % display only: at most 2e6 points
        t = (0:numel(x)-1)' / fs;
        plot(ax_wave, t(1:step:end), x(1:step:end), 'PickableParts', 'none');
        xlim(ax_wave, [0 t_end]);
        ylabel(ax_wave, 'Sound pressure (Pa)');
        title(ax_wave, 'Waveform');
        xline(ax_wave, (max(play_start, 1) - 1) / fs, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5, ...
            'Tag', 'playhead', 'PickableParts', 'none');
        draw_spectrogram();
    end

    function draw_spectrogram()
        if isempty(wave_x)
            return
        end
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
            [t_spec, f_spec, L, info] = SQAT_GUI_enhanced_stft(wave_x, wave_fs, mode);
        else
            [t_spec, f_spec, L, info] = SQAT_GUI_spectrogram(wave_x, wave_fs, win, degree, overlap);
            if info.limited
                write_log(sprintf('The spectrogram was limited to %d frames: the overlap is %.0f %%.', ...
                    numel(t_spec), info.overlap));
            end
        end
        keep = f_spec >= 20;
        L = L(keep, :) + SQAT_GUI_weight_curve(f_spec(keep), wave_fs, weighting);
        if enhanced
            spec_view = struct('t', t_spec, 'f', f_spec(keep), 'L', L, 'hop', info.hop, ...
                'mode', mode, 'weighting', weighting, 'zoomed', false);
        else
            spec_view = [];
        end
        spec_busy = true;
        cla(ax_spec);
        surface(ax_spec, t_spec, f_spec(keep), zeros(nnz(keep), numel(t_spec)), L, ...
            'EdgeColor', 'none', 'PickableParts', 'none');
        ax_spec.YScale = 'log';
        ax_spec.Layer = 'top';
        xlim(ax_spec, [0 numel(wave_x) / wave_fs]);
        ylim(ax_spec, [20 wave_fs/2]);
        colormap(ax_spec, SQAT_GUI_colormap_heat(256));
        if enhanced
            clim(ax_spec, max(L, [], 'all') + [-45 0]);
        else
            clim(ax_spec, max(L, [], 'all') + [-80 0]);
        end
        cb = colorbar(ax_spec);
        if strcmp(weighting, 'Z')
            cb.Label.String = 'Level (dB SPL)';
        else
            cb.Label.String = sprintf('Level (dB(%s))', weighting);
        end
        xlabel(ax_spec, 'Time (s)');
        ylabel(ax_spec, 'Frequency (Hz)');
        if enhanced
            title(ax_spec, sprintf('Enhanced STFT (consensus of 8 to 512 ms windows, %s, relative colour scale)', mode), 'Interpreter', 'none');
        else
            title(ax_spec, sprintf('Spectrogram (%s window, %d points, %.0f %% overlap)', ...
                dd.Items{strcmp(dd.ItemsData, dd.Value)}, info.n_fft, info.overlap), 'Interpreter', 'none');
        end
        xline(ax_spec, (max(play_start, 1) - 1) / wave_fs, 'Color', [1 1 1], 'LineWidth', 1.5, ...
            'Tag', 'playhead_spectrogram', 'PickableParts', 'none');
        draw_boxes();
        spec_busy = false;
    end

    function on_spec_limits(~, ~)
        % a zoom or a pan of the spectrogram: the enhanced map of the excerpt is
        % recomputed once the limits settle, since the full map has at most
        % 2000 columns and a zoom only stretches them
        if spec_busy || isempty(spec_view)
            return
        end
        stop_spec_timer();
        spec_timer = timer('StartDelay', 0.3, 'ExecutionMode', 'singleShot', ...
            'TimerFcn', @(~, ~) apply_spec_zoom(), 'ObjectVisibility', 'off');
        start(spec_timer);
    end

    function stop_spec_timer()
        if ~isempty(spec_timer) && isvalid(spec_timer)
            stop(spec_timer);
            delete(spec_timer);
        end
        spec_timer = [];
    end

    function apply_spec_zoom()
        if isempty(spec_view) || ~il_is_open(win_wave)
            return
        end
        srf = findobj(ax_spec, 'Type', 'surface');
        if isempty(srf)
            return
        end
        lim = ax_spec.XLim;
        dur = numel(wave_x) / wave_fs;
        span = diff(lim);
        if span >= 0.95 * dur || span / 2000 >= 0.999 * spec_view.hop
            if spec_view.zoomed                                 % back to the full map
                set(srf(1), 'XData', spec_view.t, 'YData', spec_view.f, ...
                    'ZData', zeros(numel(spec_view.f), numel(spec_view.t)), 'CData', spec_view.L);
                spec_view.zoomed = false;
            end
            return
        end
        pad = 0.3;                                             % the longest window reaches 256 ms around each instant
        i1 = max(1, floor((lim(1) - pad) * wave_fs) + 1);
        i2 = min(numel(wave_x), ceil((lim(2) + pad) * wave_fs));
        n_frames = ceil((i2 - i1 + 1) / max(1, round(max(0.001, span / 2000) * wave_fs)));
        [t_z, f_z, L_z] = SQAT_GUI_enhanced_stft(wave_x(i1:i2), wave_fs, spec_view.mode, n_frames);
        keep = f_z >= 20;
        L_z = L_z(keep, :) + SQAT_GUI_weight_curve(f_z(keep), wave_fs, spec_view.weighting);
        t_z = t_z + (i1 - 1) / wave_fs;
        set(srf(1), 'XData', t_z, 'YData', f_z(keep), 'ZData', zeros(nnz(keep), numel(t_z)), 'CData', L_z);
        spec_view.zoomed = true;
    end

    function move_playhead(sample)
        if ~il_is_open(win_wave) || isempty(wave_fs)
            return
        end
        t_now = (sample - 1) / wave_fs;
        set(findobj(win_wave, 'Tag', 'playhead'), 'Value', t_now);
        set(findobj(win_wave, 'Tag', 'playhead_spectrogram'), 'Value', t_now);
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

    function keep_figures(new_figs, f, id, save_figs, split, folder, suffix, cache_it)
        % heat colour scale, optional saving, then kept hidden for the graphs window
        if nargin < 7
            suffix = '';
        end
        if nargin < 8
            cache_it = true;
        end
        for k_fig = 1:numel(new_figs)
            set_colormap(new_figs(k_fig));
        end
        if save_figs
            [~, base] = fileparts(f.name);
            n_saved = 0;
            for k_fig = 1:numel(new_figs)
                if split
                    axs = findobj(new_figs(k_fig), 'Type', 'axes');
                    for k_ax = 1:numel(axs)
                        exportgraphics(axs(k_ax), fullfile(folder, ...
                            sprintf('%s_%s%s_%d_%d.png', base, id, suffix, k_fig, k_ax)));
                        n_saved = n_saved + 1;
                    end
                else
                    exportgraphics(new_figs(k_fig), fullfile(folder, ...
                        sprintf('%s_%s%s_%d.png', base, id, suffix, k_fig)));
                    n_saved = n_saved + 1;
                end
            end
            write_log(sprintf('%d figure(s) saved to %s', n_saved, folder));
        end
        if cache_it
            cache(end+1) = struct('file', f.path, 'metric', id, 'figs', new_figs);
        else
            delete(new_figs(isvalid(new_figs)));
        end
    end

    function clear_cache()
        for k_c = 1:numel(cache)
            delete(cache(k_c).figs(isvalid(cache(k_c).figs)));
        end
        cache = struct('file', {}, 'metric', {}, 'figs', {});
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
        if exist('theme', 'file')                % R2025a or newer
            for w = [fig, open_graph_windows(), win_wave]
                if il_is_open(w)
                    theme(w, theme_style);
                end
            end
        end
    end

    function set_status(text)
        lbl_status.Text = text;
        if il_is_open(dlg)
            dlg.Message = text;
        end
    end

    function poll_cancel()
        % the Stop of the progress dialog ends the run like the Stop button
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
% the parameters of analysis a in short, e.g. FF tv t 0.5s
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
        prefix = struct('time_skip', 't', 'dt', 'dt', 'threshold', 'thr');
        if isfield(prefix, q.name)
            name = prefix.(q.name);
        else
            name = q.name;
        end
        u = unit{1};
        if ~strcmp(u, 's')
            u = [' ' u];
        end
        parts{end+1} = sprintf('%s %g%s', name, v, u); %#ok<AGROW>
    end
end
txt = strjoin(parts, ' ');
end

function s = il_short(option)
% the short name of an option in the summary of the parameters
names = {'Free field', 'FF'; 'Diffuse field', 'DF'; 'Free-frontal', 'FFront'; 'Diffuse', 'Diff'; ...
         'Stationary', 'stat'; 'Time-varying', 'tv'; 'DIN 45692', 'DIN'; 'Aures', 'Aures'; ...
         'von Bismarck', 'Bism'};
k = find(strcmp(names(:, 1), option), 1);
s = option;
if ~isempty(k)
    s = names{k, 2};
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
T = table('Size', [0 6], 'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'cell', 'double'}, ...
    'VariableNames', {'File', 'Analysis', 'Metric', 'Channel', 'Quantity', 'Value'});
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
