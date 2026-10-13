function d = SQAT_GUI_metric_picker(parent, metrics, add)
% function d = SQAT_GUI_metric_picker(parent, metrics, add)
%
%   Modal dialog that lists the metrics of the interface, each with a tick
%   and a count: the tick adds one analysis of the metric, + and - change
%   how many (to compare one metric with other parameters), Tick all ticks
%   every metric (a count above one stays) or none, and OK adds them all,
%   in the order of the list.
%
% INPUT ARGUMENTS
%   parent : the main window of the interface (the dialog opens over it)
%   metrics : the catalogue of SQAT_GUI_metrics (fields id and label)
%   add : function handle add(ids), ids a cell array of metric ids, one per
%       analysis to add (a metric repeated as many times as its count)
%
% OUTPUTS
%   d : the uifigure of the dialog
%
% Cancel, or closing the window, adds nothing.
%
% Author: Sergio Aguirre and Gil Felix Greco, October 2026
%
% AI disclosure: code development in October 2026 assisted by
% Claude Opus 5.5 (Anthropic). All codes were verified by the
% authors.
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

n = numel(metrics);
count = zeros(1, n);
d = uifigure('Name', 'Add metrics', 'WindowStyle', 'modal', 'Tag', 'SQAT_GUI_metric_picker', ...
    'Position', [parent.Position(1) + 200, parent.Position(2) + 100, 480, min(30 * n + 110, 700)], ...
    'Visible', parent.Visible, 'CreateFcn', '');
if isprop(parent, 'Theme') && ~isempty(parent.Theme)
    d.Theme = parent.Theme;
end
g = uigridlayout(d, [3 1]);
g.RowHeight = {22, '1x', 30};
hg = uigridlayout(g, [1 2]);
hg.ColumnWidth = {'1x', 80};
hg.Padding = [0 0 0 0];
uilabel(hg, 'Text', 'Tick the metrics to analyse; + adds the same metric again, for other parameters.', ...
    'FontAngle', 'italic');
all_tick = uicheckbox(hg, 'Text', 'Tick all', 'Tag', 'picker_all', 'ValueChangedFcn', @(src, ~) tick_all(src.Value));
list = uigridlayout(g, [n 4], 'Scrollable', 'on', 'Tag', 'picker_list');
list.RowHeight = repmat({24}, 1, n);
list.ColumnWidth = {'1x', 26, 36, 26};
list.Padding = [0 0 0 0];
list.RowSpacing = 6;
ticks = gobjects(1, n);
counts = gobjects(1, n);
for k = 1:n
    ticks(k) = uicheckbox(list, 'Text', metrics(k).label, 'Tag', ['picker_tick_' metrics(k).id], ...
        'ValueChangedFcn', @(src, ~) set_count(k, double(src.Value)));
    uibutton(list, 'Text', char(8722), 'Tag', ['picker_less_' metrics(k).id], ...
        'ButtonPushedFcn', @(~, ~) step(k, -1));
    counts(k) = uilabel(list, 'Text', '', 'HorizontalAlignment', 'center', 'Tag', ['picker_count_' metrics(k).id]);
    uibutton(list, 'Text', '+', 'Tag', ['picker_more_' metrics(k).id], ...
        'ButtonPushedFcn', @(~, ~) step(k, 1));
end
bg = uigridlayout(g, [1 3]);
bg.ColumnWidth = {'1x', 90, 90};
bg.Padding = [0 0 0 0];
uilabel(bg, 'Text', '');
uibutton(bg, 'Text', 'OK', 'Tag', 'picker_ok', 'ButtonPushedFcn', @(~, ~) ok());
uibutton(bg, 'Text', 'Cancel', 'Tag', 'picker_cancel', 'ButtonPushedFcn', @(~, ~) delete(d));
SQAT_GUI_paint(d, il_style(parent));             % white in the light theme, as the main window

    function tick_all(on)
        for k_all = 1:n
            if ~on
                set_count(k_all, 0);
            elseif count(k_all) == 0
                set_count(k_all, 1);
            end
        end
    end

    function step(k, delta)
        % read when pressed: an anonymous function would keep the count of its creation
        set_count(k, count(k) + delta);
    end

    function set_count(k, c)
        % the tick shows whether the metric is added, the label how many times
        count(k) = max(c, 0);
        ticks(k).Value = count(k) > 0;
        counts(k).Text = '';
        if count(k) > 0
            counts(k).Text = sprintf('%dx', count(k));
        end
        all_tick.Value = all(count > 0);               % ticked while every metric is
    end

    function ok()
        ids = repelem({metrics.id}, count);
        delete(d);
        if ~isempty(ids)
            add(ids);
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
