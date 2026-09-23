function varargout = SQAT_GUI(files, varargin)
% function fig = SQAT_GUI(files, varargin)
%
%   Graphical interface to the metrics of SQAT, laid out as the interface of
%   pySQAT. It loads .wav files, calibrates them with a dBFS value, runs the
%   selected metrics with the chosen parameters, lists their single values,
%   shows in a graphs window the figure that each SQAT function draws (with
%   the heat colour scale of SQAT_GUI_colormap_heat), shows the waveform and the spectrogram in a
%   player window, and exports the results to a spreadsheet. Every value
%   and every figure comes from the SQAT function itself: the interface only
%   builds the call (see SQAT_GUI_metrics). The analysis draws the figure of
%   the active file in the same call, so each metric runs once; the figure of
%   another file is drawn when the graphs window asks for it.
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
loaded = struct('path', {}, 'name', {}, 'nch', {}, 'fs', {});
results = il_empty_results();
series = struct('file', {}, 'metric', {}, 't', {}, 'y', {}, 'name', {});
player = [];
player_fs = [];
theme_style = 'dark';
win_graphs = [];
win_wave = [];
ax_graphs = [];
ax_wave = [];
ax_spec = [];
btn_play = [];
run_settings = struct('dBFS', 94, 'channel', 1, 'params', params);   % of the last analysis
stop_requested = false;                                             % the Stop button
cache = struct('file', {}, 'metric', {}, 'figs', {});               % SQAT figures, hidden
cmap = SQAT_GUI_colormap_heat(256);

%% Main window
fig = uifigure('Name', 'SQAT: Sound Quality Analysis Toolbox', ...
    'Position', [100 100 1300 820], 'Visible', opts.Visible, 'Tag', 'SQAT_GUI', ...
    'CloseRequestFcn', @on_close, 'CreateFcn', '');   % skips a user default CreateFcn
main = uigridlayout(fig, [3 2]);
main.RowHeight = {44, '1x', 24};
main.ColumnWidth = {280, '1x'};

top = uigridlayout(main, [1 12]);
top.Layout.Row = 1; top.Layout.Column = [1 2];
top.Padding = [0 0 0 0];
top.ColumnWidth = {70, 110, 130, 75, '1x', 60, 50, 45, 55, 80, 230, 100};
img_logo = uiimage(top, 'ImageSource', fullfile(dir_logos, 'logo_white.png'), 'Tag', 'logo');
lbl_files = uilabel(top, 'Text', 'No files loaded', 'Tag', 'file_count', 'HorizontalAlignment', 'center');
uibutton(top, 'Text', 'Open WAV files...', 'Tag', 'load_files', 'ButtonPushedFcn', @on_load_files);
uilabel(top, 'Text', 'Active file:', 'HorizontalAlignment', 'right');
dd_file = uidropdown(top, 'Items', {}, 'Tag', 'active_file', 'ValueChangedFcn', @on_active_file);
uilabel(top, 'Text', 'Channel:', 'HorizontalAlignment', 'right');
dd_channel = uidropdown(top, 'Items', {'1'}, 'Tag', 'channel', 'ValueChangedFcn', @on_signal_changed);
uilabel(top, 'Text', 'dBFS:', 'HorizontalAlignment', 'right');
ed_dbfs = uieditfield(top, 'numeric', 'Value', 94, 'Tag', 'dbfs', ...
    'Tooltip', 'dB SPL of a full-scale amplitude (94: full scale 1.0 is 1 Pa)', ...
    'ValueChangedFcn', @on_signal_changed);
uilabel(top, 'Text', 'Plot metric:', 'HorizontalAlignment', 'right');
dd_plot = uidropdown(top, 'Items', {}, 'Tag', 'plot_metric', 'ValueChangedFcn', @on_plot_changed);
btn_theme = uibutton(top, 'Text', 'Light theme', 'Tag', 'theme', 'ButtonPushedFcn', @on_theme);

left = uigridlayout(main, [8 1]);
left.Layout.Row = 2; left.Layout.Column = 1;
left.Padding = [0 0 0 0];
left.RowHeight = {20, '1x', 20, 44, 28, 32, 32, 32};
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
og = uigridlayout(opt_panel, [1 6]);
og.ColumnWidth = {150, 105, 105, 105, '1x', 80};
og.Padding = [6 4 6 4];
cb_show = uicheckbox(og, 'Text', 'Show plots after run', 'Tag', 'show_plots');
cb_split = uicheckbox(og, 'Text', 'Split figures', 'Tag', 'split_figures', ...
    'Tooltip', 'One tab (and one saved file) per panel of a figure', ...
    'ValueChangedFcn', @on_graph_option);
cb_together = uicheckbox(og, 'Text', 'Plot together', 'Tag', 'plot_together', ...
    'Tooltip', 'The graphs window overlays the selected metric for all analysed files', ...
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

    function on_active_file(~, ~)
        update_channels();
        refresh_windows();
    end

    function on_signal_changed(~, ~)
        if il_is_open(win_wave)
            draw_waveform_window();
        end
    end

    function on_graph_option(~, ~)
        if il_is_open(win_graphs) && ~isempty(series)
            draw_graphs();
        end
    end

    function on_plot_changed(~, ~)
        if il_is_open(win_graphs)
            draw_graphs();
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
            lbl_status.Text = 'Stopping...';
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
        active_path = '';                        % the graphs window plots the active
        if ~isempty(loaded)                      % file, so its figures are drawn in
            active_path = active_file().path;    % the analysis call and the metric
        end                                      % runs once
        dBFS = ed_dbfs.Value;
        channel = str2double(dd_channel.Value);
        clear_cache();
        run_settings = struct('dBFS', dBFS, 'channel', channel, 'params', params);

        new_results = il_empty_results();
        new_series = struct('file', {}, 'metric', {}, 't', {}, 'y', {}, 'name', {});
        n_total = numel(loaded) * numel(sel);
        n_done = 0;
        n_errors = 0;
        set_progress(0);
        stop_requested = false;
        btn_stop.Enable = 'on';
        stop_off = onCleanup(@() set(btn_stop, 'Enable', 'off')); %#ok<NASGU>
        files_order = 1:numel(loaded);           % the file on screen is analysed first
        k_active = find(strcmp({loaded.path}, active_path), 1);
        if ~isempty(k_active)
            files_order = [k_active, files_order(files_order ~= k_active)];
        end
        t_start = tic;
        for i = files_order
            f = loaded(i);
            ch = channel;
            if ch > f.nch
                write_log(sprintf('%s has %d channel(s); channel 1 is used.', f.name, f.nch));
                ch = 1;
            end
            try
                [x, fs] = SQAT_GUI_load(f.path, dBFS, ch);
            catch err
                write_log(sprintf('ERROR reading %s: %s', f.name, err.message));
                n_errors = n_errors + numel(sel);
                n_done = n_done + numel(sel);
                continue
            end
            plan = SQAT_GUI_share(sel, params, numel(x), fs);
            done = struct();                     % outputs of this file, by metric id
            rows = struct();                     % their rows of the results table
            got = struct();                      % their time series
            for j = 1:numel(plan)
                drawnow                          % the Stop button gets its turn here
                if stop_requested
                    break
                end
                e = metrics(strcmp({metrics.id}, plan(j).id));
                lbl_status.Text = sprintf('Running %s on %s (%d of %d)', e.label, f.name, n_done + 1, n_total);
                try
                    [OUT, new_figs] = run_step(e, plan(j), done, x, fs, f, ...
                        show || strcmp(f.path, active_path));
                    done.(e.id) = OUT;
                    T = SQAT_GUI_single_values(OUT);
                    n = height(T);
                    rows.(e.id) = [table(repmat({f.name}, n, 1), repmat({e.id}, n, 1), ...
                        'VariableNames', {'File', 'Metric'}), T];
                    [ts, ys, name] = SQAT_GUI_series(OUT, e);
                    got.(e.id) = struct('file', f.path, 'metric', e.id, ...
                        't', ts, 'y', ys, 'name', name);
                    if ~isempty(new_figs)
                        keep_figures(new_figs, f, e.id, save_figs, split, folder);
                    end
                catch err
                    write_log(sprintf('ERROR in %s (%s): %s', e.id, f.name, err.message));
                    n_errors = n_errors + 1;
                end
                n_done = n_done + 1;
                set_progress(100 * n_done / n_total);
            end
            for j = 1:numel(sel)                 % the table keeps the order of the list
                if isfield(rows, sel{j})
                    new_results = [new_results; rows.(sel{j})]; %#ok<AGROW>
                    new_series(end+1) = got.(sel{j}); %#ok<AGROW>
                end
            end
            if stop_requested
                break                            % what ran so far is kept
            end
        end

        results = new_results;
        series = new_series;
        tbl.Data = results;
        ran = sel(ismember(sel, {series.metric}));
        [~, idx] = ismember(ran, {metrics.id});
        previous = dd_plot.Value;
        set(dd_plot, 'Items', {metrics(idx).label}, 'ItemsData', ran);
        if il_is_member(previous, ran)
            dd_plot.Value = previous;
        end
        if stop_requested
            msg = sprintf('Stopped: %d value(s) from %d of the %d analysis step(s) in %.1f s', ...
                height(results), n_done, n_total, toc(t_start));
        else
            set_progress(100);
            msg = sprintf('Done: %d value(s) from %d file(s) and %d metric(s) in %.1f s', ...
                height(results), numel(loaded), numel(sel), toc(t_start));
        end
        if n_errors > 0
            msg = sprintf('%s, %d error(s)', msg, n_errors);
        end
        lbl_status.Text = msg;
        write_log([msg '.']);
        if cb_show.Value && ~isempty(ran)
            on_open_graphs();
        elseif il_is_open(win_graphs) && ~isempty(dd_plot.Items) && ~isempty(series)
            draw_graphs();       % the signal is the one the player already shows
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
        if isempty(series) || isempty(dd_plot.Items)
            write_log('No results to plot. Run an analysis first.');
            return
        end
        if ~il_is_open(win_graphs)
            win_graphs = uifigure('Name', 'SQAT graphs', 'Position', [160 160 1100 700], ...
                'Visible', fig.Visible, 'Tag', 'SQAT_GUI_graphs', 'CreateFcn', '');
            apply_theme();
        end
        draw_graphs();
        if strcmp(fig.Visible, 'on')
            figure(win_graphs);
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
            ch = min(str2double(dd_channel.Value), size(y, 2));
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
        if il_is_open(win_graphs), delete(win_graphs); end
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
                'nch', info.NumChannels, 'fs', info.SampleRate); %#ok<AGROW>
        end
        if isempty(loaded)
            return
        end
        previous = dd_file.Value;
        set(dd_file, 'Items', {loaded.name}, 'ItemsData', {loaded.path});
        if il_is_member(previous, {loaded.path})
            dd_file.Value = previous;
        end
        if numel(loaded) == 1
            lbl_files.Text = '1 file loaded';
        else
            lbl_files.Text = sprintf('%d files loaded', numel(loaded));
        end
        write_log(sprintf('%d file(s) loaded.', numel(loaded)));
        update_channels();
        refresh_windows();
    end

    function f = active_file()
        f = loaded(strcmp({loaded.path}, dd_file.Value));
    end

    function update_channels()
        if isempty(loaded)
            return
        end
        f = active_file();
        items = arrayfun(@num2str, 1:f.nch, 'UniformOutput', false);
        previous = dd_channel.Value;
        dd_channel.Items = items;
        if il_is_member(previous, items)
            dd_channel.Value = previous;
        end
    end

    function refresh_windows()
        if il_is_open(win_graphs) && ~isempty(dd_plot.Items) && ~isempty(series)
            draw_graphs();
        end
        if il_is_open(win_wave)
            draw_waveform_window();
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

    function draw_graphs()
        id = dd_plot.Value;
        if ~cb_together.Value && show_sqat_figure(id)
            return
        end
        draw_series(id);
    end

    function ok = show_sqat_figure(id)
        % the figure that the SQAT function draws, copied into the graphs window
        f = active_file();
        label = metrics(strcmp({metrics.id}, id)).label;
        k_c = find(strcmp({cache.file}, f.path) & strcmp({cache.metric}, id), 1);
        if isempty(k_c) || ~all(isvalid(cache(k_c).figs))
            if ~any(strcmp({series.file}, f.path) & strcmp({series.metric}, id))
                ok = false;                       % no result for this file
                return
            end
            write_log(sprintf('Drawing the SQAT figure of %s for %s ...', id, f.name));
            lbl_status.Text = sprintf('Drawing the SQAT figure of %s for %s', label, f.name);
            drawnow limitrate
            try
                ch = min(run_settings.channel, f.nch);
                [x, fs] = SQAT_GUI_load(f.path, run_settings.dBFS, ch);
                e = metrics(strcmp({metrics.id}, id));
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
        delete(win_graphs.Children);
        tg = uitabgroup(win_graphs, 'Units', 'normalized', 'Position', [0 0 1 1]);
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
        win_graphs.Name = sprintf('SQAT graphs: %s, %s', label, f.name);
        ok = true;
    end

    function draw_series(id)
        % the time series drawn by the interface: one file, or all of them
        delete(win_graphs.Children);
        ax_graphs = uiaxes(uigridlayout(win_graphs, [1 1]));
        label = metrics(strcmp({metrics.id}, id)).label;
        if cb_together.Value
            s = series(strcmp({series.metric}, id) & ~cellfun(@isempty, {series.t}));
            names = cell(1, numel(s));
            hold(ax_graphs, 'on');
            for k_s = 1:numel(s)
                plot(ax_graphs, s(k_s).t, s(k_s).y);
                names{k_s} = loaded(strcmp({loaded.path}, s(k_s).file)).name;
            end
            hold(ax_graphs, 'off');
            what = sprintf('all files (%d)', numel(s));
        else
            f = active_file();
            s = series(strcmp({series.file}, f.path) & strcmp({series.metric}, id));
            if isempty(s) || isempty(s.t)
                title(ax_graphs, sprintf('%s, %s: no time series (stationary result or no result)', ...
                    label, f.name), 'Interpreter', 'none');
                return
            end
            plot(ax_graphs, s.t, s.y);
            names = {};
            what = f.name;
        end
        if isempty(s)
            title(ax_graphs, sprintf('%s: no time series to overlay', label), 'Interpreter', 'none');
            return
        end
        y_all = vertcat(s.y);
        y_range = [min(y_all) max(y_all)];
        y_ref = max(abs(y_range));
        if all(isfinite(y_range)) && y_ref > 0 && diff(y_range) <= 1e-3 * y_ref
            % a constant result: show it at +/-5 %, away from its rounding noise
            ylim(ax_graphs, mean(y_range) + [-0.05 0.05] * y_ref);
        end
        xlabel(ax_graphs, 'Time (s)');
        ylabel(ax_graphs, s(1).name);
        title(ax_graphs, sprintf('%s, %s', label, what), 'Interpreter', 'none');
        if ~isempty(names)
            legend(ax_graphs, names, 'Interpreter', 'none', 'Location', 'best');
        end
        win_graphs.Name = sprintf('SQAT graphs: %s', label);
    end

    function draw_waveform_window()
        cla(ax_wave);
        cla(ax_spec);
        f = active_file();
        ch = min(str2double(dd_channel.Value), f.nch);
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

    function keep_figures(new_figs, f, id, save_figs, split, folder)
        % heat colour scale, optional saving, then kept hidden for the graphs window
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
                            sprintf('%s_%s_%d_%d.png', base, id, k_fig, k_ax)));
                        n_saved = n_saved + 1;
                    end
                else
                    exportgraphics(new_figs(k_fig), fullfile(folder, sprintf('%s_%s_%d.png', base, id, k_fig)));
                    n_saved = n_saved + 1;
                end
            end
            write_log(sprintf('%d figure(s) saved to %s', n_saved, folder));
        end
        cache(end+1) = struct('file', f.path, 'metric', id, 'figs', new_figs);
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
            for w = [fig, win_graphs, win_wave]
                if il_is_open(w)
                    theme(w, theme_style);
                end
            end
        end
    end

    function set_progress(value)
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
T = table('Size', [0 4], 'VariableTypes', {'cell', 'cell', 'cell', 'double'}, ...
    'VariableNames', {'File', 'Metric', 'Quantity', 'Value'});
end
