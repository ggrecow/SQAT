function varargout = SQAT_GUI(files, varargin)
% function fig = SQAT_GUI(files, varargin)
%
%   Graphical interface to the metrics of SQAT, laid out as the interface of
%   pySQAT. It loads .wav files into a list of signals (a tick marks a signal
%   for use, the x removes it), calibrates them with a dBFS value, runs the
%   selected metrics with the chosen parameters on the ticked signals, and
%   lists their single values. The channel to analyse is one channel of each
%   file or All: the ECMA-418-2 metrics take a stereo pair in one call and
%   return the left, the right and (except the tonality) the combined
%   binaural result, so a pair runs once.
%
%   A graphs window plots one metric of the ticked signals. The analysis is
%   chosen in the window (a time series, a profile over the critical bands, a
%   map of band against time, or the statistics): the signals are overlaid,
%   and the maps sit side by side on one colour scale. For one signal the
%   window also offers the figure that the SQAT function draws (with the
%   heat colour scale of SQAT_GUI_colormap_heat) and all the analyses at once. Pin keeps a window
%   with its signals, so that Open Graphs Window opens another one to compare
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
params = struct();
for k = 1:numel(metrics)
    p = struct();
    for q = metrics(k).params
        p.(q.name) = q.value;
    end
    params.(metrics(k).id) = p;
end
loaded = struct('path', {}, 'name', {}, 'nch', {}, 'fs', {}, 'marked', {});
active_idx = 0;                                       % the signal on screen in the player
results = il_empty_results();
store = il_empty_store();                             % analyses of each signal, metric and channel
player = [];
player_fs = [];
theme_style = 'dark';
graph_figs = gobjects(0);                              % the graphs windows
last_metric = '';                                     % the metric the last graphs window showed
win_wave = [];
ax_wave = [];
ax_spec = [];
btn_play = [];
run_settings = struct('dBFS', 94, 'channel', '1', 'params', params);   % of the last analysis
stop_requested = false;                                             % the Stop button or the dialog
dlg = [];                                                           % progress dialog of a run
cache = struct('file', {}, 'metric', {}, 'figs', {});               % SQAT figures, hidden
cmap = SQAT_GUI_colormap_heat(256);

%% Main window
fig = uifigure('Name', 'SQAT: Sound Quality Analysis Toolbox', ...
    'Position', [100 100 1300 820], 'Visible', opts.Visible, 'Tag', 'SQAT_GUI', ...
    'CloseRequestFcn', @on_close, 'CreateFcn', '');   % skips a user default CreateFcn
main = uigridlayout(fig, [3 2]);
main.RowHeight = {44, '1x', 24};
main.ColumnWidth = {280, '1x'};

top = uigridlayout(main, [1 9]);
top.Layout.Row = 1; top.Layout.Column = [1 2];
top.Padding = [0 0 0 0];
top.ColumnWidth = {70, 110, 130, '1x', 60, 70, 45, 55, 100};
img_logo = uiimage(top, 'ImageSource', fullfile(dir_logos, 'logo_white.png'), 'Tag', 'logo');
lbl_files = uilabel(top, 'Text', 'No files loaded', 'Tag', 'file_count', 'HorizontalAlignment', 'center');
uibutton(top, 'Text', 'Open WAV files...', 'Tag', 'load_files', 'ButtonPushedFcn', @on_load_files);
uilabel(top, 'Text', '');
uilabel(top, 'Text', 'Channel:', 'HorizontalAlignment', 'right');
dd_channel = uidropdown(top, 'Items', {'1'}, 'Tag', 'channel', 'ValueChangedFcn', @on_signal_changed, ...
    'Tooltip', ['Channel to analyse. All: every channel of each file; the ECMA-418-2 ' ...
                'metrics analyse a stereo file as a binaural pair, in one call']);
uilabel(top, 'Text', 'dBFS:', 'HorizontalAlignment', 'right');
ed_dbfs = uieditfield(top, 'numeric', 'Value', 94, 'Tag', 'dbfs', ...
    'Tooltip', 'dB SPL of a full-scale amplitude (94: full scale 1.0 is 1 Pa)', ...
    'ValueChangedFcn', @on_signal_changed);
btn_theme = uibutton(top, 'Text', 'Light theme', 'Tag', 'theme', 'ButtonPushedFcn', @on_theme);

left = uigridlayout(main, [10 1]);
left.Layout.Row = 2; left.Layout.Column = 1;
left.Padding = [0 0 0 0];
left.RowHeight = {20, 150, 20, '1x', 20, 44, 28, 32, 32, 32};
uilabel(left, 'Text', 'SIGNALS (tick to use, x to remove)', 'FontWeight', 'bold');
tbl_signals = uitable(left, 'Tag', 'signals_table', 'RowName', {}, 'ColumnName', {'', 'Signal', ''}, ...
    'ColumnWidth', {28, 190, 28}, 'ColumnEditable', [true false false], ...
    'CellEditCallback', @on_signal_marked, 'CellSelectionCallback', @on_signal_selected);
uilabel(left, 'Text', 'METRICS TO ANALYZE (Ctrl or Cmd + click)', 'FontWeight', 'bold');
lb_metrics = uilistbox(left, 'Items', {metrics.label}, 'ItemsData', {metrics.id}, ...
    'Multiselect', 'on', 'Value', {'Loudness_ISO532_1'}, 'Tag', 'metrics_list', ...
    'ValueChangedFcn', @on_metrics);
uilabel(left, 'Text', 'ACTIONS', 'FontWeight', 'bold');
uibutton(left, 'Text', 'Run Analysis', 'Tag', 'run', 'FontWeight', 'bold', ...
    'BackgroundColor', green, 'FontColor', [1 1 1], 'ButtonPushedFcn', @on_run);
btn_stop = uibutton(left, 'Text', 'Stop', 'Tag', 'stop_run', 'Enable', 'off', ...
    'Tooltip', 'Ends the run after the metric being computed', 'ButtonPushedFcn', @on_stop_run);
uibutton(left, 'Text', 'Open Graphs Window', 'Tag', 'open_graphs', 'ButtonPushedFcn', @on_open_graphs);
uibutton(left, 'Text', 'Waveform / Play', 'Tag', 'open_waveform', 'ButtonPushedFcn', @on_open_waveform);
uibutton(left, 'Text', 'Export results...', 'Tag', 'export', 'ButtonPushedFcn', @on_export);

right = uigridlayout(main, [3 1]);
right.Layout.Row = 2; right.Layout.Column = 2;
right.Padding = [0 0 0 0];
right.RowHeight = {56, 128, '1x'};

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

par_panel = uipanel(right, 'Title', 'PARAMETERS');
pg = uigridlayout(par_panel, [2 1]);
pg.RowHeight = {22, '1x'};
pg.Padding = [6 4 6 4];
ph = uigridlayout(pg, [1 2]);
ph.Padding = [0 0 0 0];
ph.ColumnWidth = {130, 280};
uilabel(ph, 'Text', 'Edit parameters for:');
dd_param = uidropdown(ph, 'Items', {}, 'Tag', 'param_metric', 'ValueChangedFcn', @on_param_metric);
param_grid = uigridlayout(pg, [2 1]);
param_grid.Padding = [0 0 0 0];

tabs = uitabgroup(right);
tab_console = uitab(tabs, 'Title', 'Console output');
console = uitextarea(uigridlayout(tab_console, [1 1]), 'Value', {''}, 'Editable', 'off', ...
    'Tag', 'console', 'FontName', 'Monospaced');
tab_results = uitab(tabs, 'Title', 'Results');
tbl = uitable(uigridlayout(tab_results, [1 1]), 'Data', results, 'Tag', 'results_table');

status_bar = uigridlayout(main, [1 2]);
status_bar.Layout.Row = 3; status_bar.Layout.Column = [1 2];
status_bar.Padding = [0 0 0 0];
status_bar.ColumnWidth = {280, '1x'};
gauge = uigauge(status_bar, 'linear', 'Tag', 'progress', 'Limits', [0 100], 'Value', 0, ...
    'MajorTicks', [], 'MinorTicks', []);
lbl_status = uilabel(status_bar, 'Text', 'Ready', 'Tag', 'status', 'HorizontalAlignment', 'right');

%% Start
apply_theme();
add_files(files);
on_metrics();
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

    function on_signal_marked(~, event)
        k = event.Indices(1);
        loaded(k).marked = logical(event.NewData);
        refresh_windows();
    end

    function on_signal_selected(~, event)
        if isempty(event.Indices)
            return
        end
        k = event.Indices(1, 1);
        if event.Indices(1, 2) == 3
            remove_signal(k);                % the x at the right of the name
            return
        end
        if k ~= active_idx
            active_idx = k;
            show_active();
            refresh_windows();
        end
    end

    function on_signal_changed(~, ~)
        if il_is_open(win_wave)
            draw_waveform_window();
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

    function on_metrics(~, ~)
        sel = lb_metrics.Value;
        [~, idx] = ismember(sel, {metrics.id});
        previous = dd_param.Value;
        set(dd_param, 'Items', {metrics(idx).label}, 'ItemsData', sel);
        if il_is_member(previous, sel)
            dd_param.Value = previous;
        end
        build_params();
    end

    function on_param_metric(~, ~)
        build_params();
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
        sel = lb_metrics.Value;
        if isempty(sel)
            write_log('No metrics selected.');
            return
        end
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
        dBFS = ed_dbfs.Value;
        option = dd_channel.Value;               % a channel number or All
        clear_cache();
        run_settings = struct('dBFS', dBFS, 'channel', option, 'params', params);

        new_results = il_empty_results();
        new_store = il_empty_store();
        stereo_sel = ismember(sel, {metrics([metrics.stereo]).id});
        n_total = 0;
        for i = files_order
            n_ch = numel(channel_list(loaded(i), option));
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
            cl = channel_list(f, option);
            if ~strcmp(option, 'All') && str2double(option) > f.nch
                write_log(sprintf('%s has %d channel(s); channel 1 is used.', f.name, f.nch));
            end
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
                plan = SQAT_GUI_share(ids_single, params, numel(x), fs);
                done = struct();                 % outputs of this channel, by metric id
                for j = 1:numel(plan)
                    drawnow                      % the Stop button gets its turn here
                    poll_cancel();
                    if stop_requested
                        break
                    end
                    e = metrics(strcmp({metrics.id}, plan(j).id));
                    set_status(sprintf('Running %s on %s (%d of %d)', e.label, f.name, n_done + 1, n_total));
                    try
                        [OUT, new_figs] = run_step(e, plan(j), done, x, fs, f, ...
                            show || (first && strcmp(f.path, active_path)));
                        done.(e.id) = OUT;
                        entries(end+1) = il_entry_of(f, e, OUT, num2str(c), 1, 1); %#ok<AGROW>
                        if ~isempty(new_figs)
                            keep_figures(new_figs, f, e.id, save_figs, split, folder, suffix, first);
                        end
                    catch err
                        write_log(sprintf('ERROR in %s (%s): %s', e.id, f.name, err.message));
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
                    e = metrics(strcmp({metrics.id}, ids_joint{j}));
                    set_status(sprintf('Running %s on %s (%d of %d)', e.label, f.name, n_done + 1, n_total));
                    write_log(sprintf('Running %s on %s, both channels ...', e.id, f.name));
                    try
                        [OUT, new_figs] = run_metric(e, X, fs, params.(e.id), ...
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
                            keep_figures(new_figs, f, e.id, save_figs, split, folder, '', true);
                        end
                    catch err
                        write_log(sprintf('ERROR in %s (%s): %s', e.id, f.name, err.message));
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
                    new_results = [new_results; [table(repmat({f.name}, n, 1), ...
                        repmat({en.metric}, n, 1), repmat({en.channel}, n, 1), ...
                        'VariableNames', {'File', 'Metric', 'Channel'}), en.values]]; %#ok<AGROW>
                end
            end
            if stop_requested
                break                            % what ran so far is kept
            end
        end

        results = new_results;
        store = new_store;
        tbl.Data = results;
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
            win_wave = uifigure('Name', 'Waveform', 'Position', [140 140 1100 640], ...
                'Visible', fig.Visible, 'Tag', 'SQAT_GUI_waveform', 'CloseRequestFcn', @on_close_waveform, ...
                'CreateFcn', '');
            gw = uigridlayout(win_wave, [3 1]);
            gw.RowHeight = {30, '1x', '1.2x'};
            hw = uigridlayout(gw, [1 3]);
            hw.Padding = [0 0 0 0];
            hw.ColumnWidth = {90, 90, '1x'};
            btn_play = uibutton(hw, 'Text', 'Play', 'Tag', 'play', 'ButtonPushedFcn', @on_play);
            uibutton(hw, 'Text', 'Stop', 'Tag', 'stop', 'ButtonPushedFcn', @on_stop);
            ax_wave = uiaxes(gw, 'Tag', 'waveform_axes');
            ax_spec = uiaxes(gw, 'Tag', 'spectrogram');
            apply_theme();
        end
        draw_waveform_window();
        if strcmp(fig.Visible, 'on')
            figure(win_wave);
        end
    end

    function on_play(~, ~)
        if isempty(loaded)
            return
        end
        if ~isempty(player) && isplaying(player)
            pause(player);
            btn_play.Text = 'Play';
            write_log('Paused.');
            return
        end
        if ~isempty(player) && player.CurrentSample > 1
            resume(player);
            btn_play.Text = 'Pause';
            write_log('Playing.');
            return
        end
        f = active_file();
        try
            [y, fs] = audioread(f.path);
            ch = min(channel_of_active(), size(y, 2));
            player = audioplayer(y(:, ch), fs);
            player_fs = fs;
            player.TimerPeriod = 0.05;
            player.TimerFcn = @(~, ~) move_playhead(player.CurrentSample);
            player.StopFcn = @(~, ~) on_player_stopped();
            play(player);
            btn_play.Text = 'Pause';
            write_log(sprintf('Playing %s, channel %d.', f.name, ch));
        catch err
            player = [];
            write_log(['Audio output unavailable: ' err.message]);
        end
    end

    function on_stop(~, ~)
        if ~isempty(player)
            stop(player);
            write_log('Stopped.');
        end
        player = [];
        move_playhead(1);
        if il_is_open(win_wave)
            btn_play.Text = 'Play';
        end
    end

    function on_player_stopped()
        % called on pause, on stop and at the end of the file
        if ~isempty(player) && player.CurrentSample == 1 && il_is_open(win_wave)
            btn_play.Text = 'Play';
            move_playhead(1);
        end
    end

    function on_close_waveform(~, ~)
        on_stop();
        delete(win_wave);
    end

    function on_close(~, ~)
        if ~isempty(player)
            stop(player);
        end
        clear_cache();
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
                'nch', info.NumChannels, 'fs', info.SampleRate, 'marked', true); %#ok<AGROW>
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
        tbl.Data = results;
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
            tbl_signals.Data = table(false(0, 1), strings(0, 1), strings(0, 1));
            lbl_files.Text = 'No files loaded';
        else
            tbl_signals.Data = table(logical([loaded.marked]'), string({loaded.name}'), ...
                repmat("x", n, 1));
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
        update_channels();
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
        ch = str2double(dd_channel.Value);
        if isnan(ch)
            ch = 1;                          % All
        end
        ch = min(ch, active_file().nch);
    end

    function update_channels()
        if isempty(loaded)
            return
        end
        n = max([loaded.nch]);
        items = arrayfun(@num2str, 1:n, 'UniformOutput', false);
        if n > 1
            items{end+1} = 'All';
        end
        previous = dd_channel.Value;
        dd_channel.Items = items;
        if il_is_member(previous, items)
            dd_channel.Value = previous;
        end
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

    function build_params()
        delete(param_grid.Children);
        if isempty(dd_param.Items)
            return
        end
        id = dd_param.Value;
        e = metrics(strcmp({metrics.id}, id));
        n = numel(e.params);
        param_grid.RowHeight = {18, 24};
        param_grid.ColumnWidth = repmat({'1x'}, 1, max(n, 3));
        for k_par = 1:n
            q = e.params(k_par);
            lbl = uilabel(param_grid, 'Text', q.label);
            lbl.Layout.Row = 1; lbl.Layout.Column = k_par;
            current = params.(id).(q.name);
            if strcmp(q.type, 'choice')
                c = uidropdown(param_grid, 'Items', q.options(:, 1)', ...
                    'ItemsData', q.options(:, 2)', 'Value', current);
            else
                c = uieditfield(param_grid, 'numeric', 'Value', current);
            end
            c.Tag = ['param_' q.name];
            c.Layout.Row = 2; c.Layout.Column = k_par;
            c.ValueChangedFcn = @(src, ~) set_param(id, q.name, src.Value);
        end
    end

    function set_param(id, name, value)
        params.(id).(name) = value;
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
            'Tooltip', ['Keeps this window with its signals; Open Graphs Window then opens ' ...
                        'another one to compare with']);
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
            src.UserData = {loaded([loaded.marked]).path};   % the signals of the window stay these
            src.Text = 'Pinned';
        else
            src.UserData = {};
            src.Text = 'Pin';
        end
        draw_window(w);
    end

    function paths = window_paths(w)
        % the signals of a window: the ticked ones, or the ones it was pinned with
        pin = findobj(w, 'Tag', 'graph_pin');
        if pin.Value
            paths = pin.UserData;
        else
            paths = {loaded([loaded.marked]).path};
        end
    end

    function draw_window(w)
        % refreshes the selectors of a window and draws what they ask for
        paths = window_paths(w);
        dm = findobj(w, 'Tag', 'graph_metric');
        dc = findobj(w, 'Tag', 'graph_channel');
        da = findobj(w, 'Tag', 'graph_analysis');
        body = findobj(w, 'Tag', 'graph_body');
        delete(body.Children);
        in_window = ismember({store.file}, paths);
        ids = {metrics(ismember({metrics.id}, {store(in_window).metric})).id};
        if isempty(ids)
            set([dm, dc, da], 'Items', {});
            uilabel(uigridlayout(body, [1 1]), 'HorizontalAlignment', 'center', ...
                'Text', 'No results for the signals of this window. Tick signals and run an analysis.');
            w.Name = 'SQAT graphs';
            return
        end
        [~, idx] = ismember(ids, {metrics.id});
        wanted = dm.Value;
        if ~il_is_member(wanted, ids)
            wanted = last_metric;
        end
        set(dm, 'Items', {metrics(idx).label}, 'ItemsData', ids);
        if il_is_member(wanted, ids)
            dm.Value = wanted;
        end
        id = dm.Value;
        last_metric = id;
        label = metrics(strcmp({metrics.id}, id)).label;

        in_metric = in_window & strcmp({store.metric}, id);
        chans = unique({store(in_metric).channel}, 'stable');
        chans = [sort(chans(~strcmp(chans, 'Binaural'))), chans(strcmp(chans, 'Binaural'))];
        wanted = dc.Value;
        dc.Items = chans;
        if il_is_member(wanted, chans)
            dc.Value = wanted;
        end
        chan = dc.Value;

        entries = il_empty_store();
        for k_p = 1:numel(paths)
            k_e = find(in_metric & strcmp({store.file}, paths{k_p}) & strcmp({store.channel}, chan), 1);
            if ~isempty(k_e)
                entries(end+1) = store(k_e); %#ok<AGROW>
            end
        end
        n_missing = nnz(ismember(paths, {store(in_metric).file})) - numel(entries);
        if n_missing > 0
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
                if ~show_sqat_figure(body, id, entries(1).file)
                    da.Value = data{3};      % the first analysis of the metric, or the statistics
                    if strcmp(da.Value, 'stats')
                        draw_stats(body, entries);
                    else
                        draw_analysis(body, entries, da.Value, chan);
                    end
                end
            case 'all'
                tg = uitabgroup(body, 'Units', 'normalized', 'Position', [0 0 1 1]);
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

    function ok = show_sqat_figure(parent, id, path)
        % the figure that the SQAT function draws, copied into the window
        f = loaded(strcmp({loaded.path}, path));
        label = metrics(strcmp({metrics.id}, id)).label;
        k_c = find(strcmp({cache.file}, f.path) & strcmp({cache.metric}, id), 1);
        if isempty(k_c) || ~all(isvalid(cache(k_c).figs))
            write_log(sprintf('Drawing the SQAT figure of %s for %s ...', id, f.name));
            lbl_status.Text = sprintf('Drawing the SQAT figure of %s for %s', label, f.name);
            drawnow limitrate
            try
                cl = channel_list(f, run_settings.channel);
                e = metrics(strcmp({metrics.id}, id));
                if ~e.stereo
                    cl = cl(1);                   % one figure per signal: the first channel
                end
                [x, fs] = SQAT_GUI_load(f.path, run_settings.dBFS, cl);
                [~, new_figs] = run_metric(e, x, fs, run_settings.params.(id), true, false);
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
        tg = uitabgroup(parent, 'Units', 'normalized', 'Position', [0 0 1 1]);
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
        names = {entries.name};
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
                title(ax, sprintf('%s: %s, channel %s', A(k).label, names{k}, chan), 'Interpreter', 'none');
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
        else
            what = sprintf('%d signals', numel(A));
        end
        title(ax, sprintf('%s: %s, channel %s', A(1).label, what, chan), 'Interpreter', 'none');
        if numel(A) > 1
            legend(ax, names, 'Interpreter', 'none', 'Location', 'best');
        end
    end

    function draw_stats(parent, entries)
        % the single values of the signals, one column each
        names = matlab.lang.makeUniqueStrings({entries.name});
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
        cla(ax_spec);
        f = active_file();
        ch = channel_of_active();
        try
            [x, fs] = SQAT_GUI_load(f.path, ed_dbfs.Value, ch);
        catch err
            write_log(sprintf('ERROR reading %s: %s', f.name, err.message));
            return
        end
        win_wave.Name = sprintf('Waveform: %s, channel %d', f.name, ch);
        t_end = numel(x) / fs;
        step = max(1, ceil(numel(x) / 2e6));   % display only: at most 2e6 points
        t = (0:numel(x)-1)' / fs;
        plot(ax_wave, t(1:step:end), x(1:step:end));
        xlim(ax_wave, [0 t_end]);
        ylabel(ax_wave, 'Sound pressure (Pa)');
        title(ax_wave, 'Waveform');
        xline(ax_wave, 0, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5, 'Tag', 'playhead');

        [t_spec, f_spec, L] = il_spectrogram(x, fs);
        keep = f_spec >= 20;
        surface(ax_spec, t_spec, f_spec(keep), zeros(nnz(keep), numel(t_spec)), L(keep, :), ...
            'EdgeColor', 'none');
        ax_spec.YScale = 'log';
        ax_spec.Layer = 'top';
        xlim(ax_spec, [0 t_end]);
        ylim(ax_spec, [20 fs/2]);
        colormap(ax_spec, SQAT_GUI_colormap_heat(256));
        L_max = max(L(keep, :), [], 'all');
        clim(ax_spec, L_max + [-80 0]);
        cb = colorbar(ax_spec);
        cb.Label.String = 'Level (dB SPL)';
        xlabel(ax_spec, 'Time (s)');
        ylabel(ax_spec, 'Frequency (Hz)');
        title(ax_spec, 'Spectrogram');
        xline(ax_spec, 0, 'Color', [1 1 1], 'LineWidth', 1.5, 'Tag', 'playhead_spectrogram');
    end

    function move_playhead(sample)
        if ~il_is_open(win_wave) || isempty(player_fs)
            return
        end
        t_now = (sample - 1) / player_fs;
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
                    e.id, f.name, step.from));
                return
            end
        end
        write_log(sprintf('Running %s on %s ...', e.id, f.name));
        [OUT, new_figs] = run_metric(e, x, fs, params.(e.id), show);
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
                e.id, err.message));
            [txt, OUT] = evalc('e.run(x, fs, p, false)');
        end
        write_block(txt);
        new_figs = il_new_figures(figs_before);
    end

    function cl = channel_list(f, option)
        % the channels of file f that the option of the run analyses
        if strcmp(option, 'All')
            cl = 1:f.nch;
        else
            cl = str2double(option);
            if cl > f.nch
                cl = 1;
            end
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
            btn_theme.Text = 'Light theme';
        else
            img_logo.ImageSource = fullfile(dir_logos, 'logo.png');
            btn_theme.Text = 'Dark theme';
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

function [t, f, L] = il_spectrogram(x, fs)
% Short-time spectrum in dB SPL (x in Pa): periodic Hann window of 1024
% samples, 50 % overlap (widened for long signals, at most 2000 frames),
% amplitude scaled so that a sine on a bin centre reads its RMS level.
n_fft = 1024;
hop = max(n_fft/2, ceil((numel(x) - n_fft) / 2000));
if numel(x) < n_fft
    x(n_fft) = 0;
end
w = 0.5 - 0.5*cos(2*pi*(0:n_fft-1)'/n_fft);
starts = 1:hop:(numel(x) - n_fft + 1);
frames = x(starts + (0:n_fft-1)') .* w;
X = fft(frames);
A = abs(X(1:n_fft/2+1, :)) * 2 / sum(w);        % peak amplitude per bin
A([1 end], :) = A([1 end], :) / 2;
L = 20*log10(max(A / sqrt(2), eps) / 2e-5);
t = (starts - 1 + n_fft/2) / fs;
f = (0:n_fft/2)' * fs / n_fft;
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

function T = il_empty_results()
T = table('Size', [0 5], 'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'double'}, ...
    'VariableNames', {'File', 'Metric', 'Channel', 'Quantity', 'Value'});
end

function S = il_empty_store()
S = struct('file', {}, 'name', {}, 'metric', {}, 'channel', {}, 'analyses', {}, 'values', {});
end

function en = il_entry_of(f, e, OUT, label, channel, n_channels)
% the analyses and the single values that OUT holds for one channel
en = struct('file', f.path, 'name', f.name, 'metric', e.id, 'channel', label, ...
    'analyses', SQAT_GUI_extract(OUT, e.id, channel), ...
    'values', SQAT_GUI_single_values(OUT, channel, n_channels));
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
