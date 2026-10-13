function SQAT_GUI_calibration_dialog(parent, name, nch, c, help_text, apply)
% function SQAT_GUI_calibration_dialog(parent, name, nch, c, help_text, apply)
%
%   Modal dialog that asks how to calibrate one sound file: full-scale level
%   (dBFS), calibrator recording or relative level (see SQAT_GUI_calibration).
%   A file with more than one channel takes one level (and one calibrator
%   recording) for all its channels or, with Same for all channels unticked,
%   one row per channel. It opens on the choice c and waits until it is closed.
%
% INPUT ARGUMENTS
%   parent : the main window of the interface (the dialog opens over it)
%   name : name of the file, for the title
%   nch : number of channels of the file
%   c : struct with method ('dbfs', 'calibrator', 'relative'), level (dB,
%       one value or one per channel) and file (calibrator recording, or a
%       cell with one per channel), the choice the dialog opens on
%   help_text : what calibration is and the three ways, shown on top
%   apply : function handle ok = apply(method, level, calfile); OK closes
%       the dialog when it returns true, and keeps it open otherwise
%
% Cancel, or closing the window, leaves the calibration of the file as it was.
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September and October
% 2026 assisted by Claude Opus 5.5 (Anthropic). All codes
% were verified by the authors.
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

methods = {'dbfs', 'calibrator', 'relative'};
items = {'Full-scale level (dBFS)', 'Calibrator recording', 'Relative level'};
hints = {'Level in dB SPL of a sample of value 1 (94 dB: 1 = 1 Pa).', ...
         ['Level of the calibrator, e.g. 94 dB (1 Pa), 114 dB (10 Pa) or the value with an adapter. ' ...
          'With one recording for all channels, channel k of the file takes channel k of the recording ' ...
          '(a stereo recording of the two ears, for example), or its last channel when it has fewer.'], ...
         ['All channels: the rms of the whole file (all channels together) is set to the level, so the ' ...
          'channels keep their difference of level. One row per channel: the rms of each channel is set ' ...
          'to its own level.' newline ...
          'The rms will be calculated based on the entire signal length. If you desire to have ' ...
          'the rms calculated otherwise, please trim the signal before loading in SQAT.']};
units = {'dB SPL at full scale', 'dB SPL of the calibrator', 'dB SPL rms'};

% the choice to open on: one row per channel when it has a level or a recording per channel
L = c.level(:).';
F = cellstr(c.file);
per_channel = nch > 1 && ((numel(L) == nch) || (iscell(c.file) && numel(F) == nch));

h = 600 + (nch > 1) * (32 + 36 * (min(nch, 6) - 1));   % the checkbox and the rows of the channels
d = uifigure('Name', ['Calibration: ' name], 'WindowStyle', 'modal', ...
    'Position', [parent.Position(1) + 200, parent.Position(2) + 200, 620, h], 'CreateFcn', '');
if isprop(parent, 'Theme') && ~isempty(parent.Theme)
    d.Theme = parent.Theme;
end
g = uigridlayout(d, [7 3]);
g.RowHeight = {'1x', 26, 84, 22, 0, 22, 30};
g.ColumnWidth = {110, '1x', 150};
t = uilabel(g, 'Text', [newline help_text], ...
    'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontSize', 14);
t.Layout.Column = [1 3];
uilabel(g, 'Text', 'Method:', 'HorizontalAlignment', 'right');
dd = uidropdown(g, 'Items', items, 'ItemsData', methods, 'Value', c.method, 'Tag', 'cal_method', ...
    'ValueChangedFcn', @(~, ~) show());
dd.Layout.Column = [2 3];
hint = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontAngle', 'italic', 'VerticalAlignment', 'top');
hint.Layout.Column = [1 3];
same = uicheckbox(g, 'Text', 'Same for all channels', 'Value', ~per_channel, 'Tag', 'cal_same', ...
    'Visible', nch > 1, 'ValueChangedFcn', @(~, ~) show());
same.Layout.Row = 4; same.Layout.Column = [2 3];
if nch == 1
    g.RowHeight{4} = 0;
end

% the rows of the levels: the first for all channels, then one per channel
cg = uigridlayout(g, [nch + 2, 4]);
cg.Layout.Row = 5; cg.Layout.Column = [1 3];
cg.ColumnWidth = {110, 150, '1x', 90};
cg.ColumnSpacing = g.ColumnSpacing;
cg.Padding = [0 0 0 0];
cg.Scrollable = 'on';
uilabel(cg, 'Text', '');
unit = uilabel(cg, 'Text', '', 'FontWeight', 'bold');
uilabel(cg, 'Text', 'Calibrator file', 'FontWeight', 'bold');
uilabel(cg, 'Text', '');
names = [{'All channels:'}, arrayfun(@(k) sprintf('Channel %d:', k), 1:nch, 'UniformOutput', false)];
tags = [{''}, arrayfun(@(k) sprintf('_%d', k), 1:nch, 'UniformOutput', false)];
if nch == 1
    names{1} = 'Level:';
end
lb = gobjects(1, nch + 1);
lv = gobjects(1, nch + 1);
fl = gobjects(1, nch + 1);
br = gobjects(1, nch + 1);
for r = 1:nch + 1
    ch = max(r - 1, 1);                         % the channel of the row (the first row takes the first values)
    lb(r) = uilabel(cg, 'Text', names{r}, 'HorizontalAlignment', 'right');
    lv(r) = uieditfield(cg, 'numeric', 'Value', L(min(ch, end)), 'Tag', ['cal_level' tags{r}]);
    fl(r) = uilabel(cg, 'Text', '', 'Tag', ['cal_file' tags{r}]);
    set_file(r, F{min(ch, end)});
    br(r) = uibutton(cg, 'Text', 'Browse...', 'Tag', ['cal_browse' tags{r}], 'ButtonPushedFcn', @(~, ~) browse(r));
end

msg = uilabel(g, 'Text', '', 'FontColor', [0.76 0.12 0.17], 'Tag', 'cal_message');
msg.Layout.Row = 6; msg.Layout.Column = [1 3];
uilabel(g, 'Text', '');
uibutton(g, 'Text', 'OK', 'Tag', 'cal_ok', 'ButtonPushedFcn', @(~, ~) ok());
uibutton(g, 'Text', 'Cancel', 'Tag', 'cal_cancel', 'ButtonPushedFcn', @(~, ~) delete(d));
SQAT_GUI_paint(d, il_style(parent));             % white in the light theme, as the main window
show();
d.Tag = 'SQAT_GUI_calibration';                  % once complete: the tests find it by its tag
uiwait(d);

    function show()
        % the hint and the unit follow the method, the rows follow Same for all channels
        k = find(strcmp(methods, dd.Value));
        hint.Text = hints{k};
        unit.Text = units{k};
        set([fl br], 'Enable', strcmp(dd.Value, 'calibrator'));
        on = false(1, nch + 1);
        on(rows_in_use()) = true;
        heights = repmat({0}, 1, nch + 1);
        heights(on) = {26};
        cg.RowHeight = [{22}, heights];
        for r2 = 1:nch + 1
            set([lb(r2) lv(r2) fl(r2) br(r2)], 'Visible', on(r2));
        end
        n = min(sum(on), 6);
        g.RowHeight{5} = 22 + n * 26 + n * cg.RowSpacing;
        msg.Text = '';
    end

    function rows = rows_in_use()
        % the rows in use: the first one, or one per channel
        rows = 1;
        if nch > 1 && ~same.Value
            rows = 2:nch + 1;
        end
    end

    function set_file(r, file)
        % the calibrator recording of row r: its name on screen, the full path in the tooltip
        [~, n, e] = fileparts(file);
        fl(r).UserData = file;
        fl(r).Text = [n e];
        fl(r).Tooltip = file;
    end

    function browse(r)
        [f, p] = uigetfile({'*.wav', 'WAV files (*.wav)'}, ['Recording of the calibrator: ' names{r}(1:end-1)]);
        figure(d);
        if ~isequal(f, 0)
            set_file(r, fullfile(p, f));
        end
    end

    function ok()
        rows = rows_in_use();
        level = [lv(rows).Value];
        if ~all(isfinite(level))
            msg.Text = 'Type a level in dB.';
            return
        end
        file = '';
        if strcmp(dd.Value, 'calibrator')
            file = {fl(rows).UserData};
            if any(cellfun(@isempty, file))
                msg.Text = 'Choose the recording of the calibrator for each row.';
                return
            end
            if isscalar(file)
                file = file{1};
            end
        end
        if apply(dd.Value, level, file)
            delete(d);
        else
            msg.Text = 'The calibration could not be set: see the console output.';
        end
    end
end

function style = il_style(parent)
% the theme of the main window: 'dark', or 'light' (also without themes, before R2025a)
style = 'light';
if isprop(parent, 'Theme') && ~isempty(parent.Theme) && strcmp(char(parent.Theme.BaseColorStyle), 'dark')
    style = 'dark';
end
end
