function SQAT_GUI_report(filename, S, cols, data, groups)
% function SQAT_GUI_report(filename, S, cols, data, groups)
%
%   Writes the report of a run of SQAT_GUI to a PDF file (A4, landscape),
%   replacing any file of the same name: the settings of the run and the
%   matrix of single values, as text, then one page per analysis with the
%   plot of the signals (lines overlaid, or maps side by side on one colour
%   scale).
%
% INPUT ARGUMENTS
%   filename : char or string, path of the .pdf file to write
%   S : table (Item, Value) with the settings of the run
%   cols : cell array of char, the column names of the matrix
%   data : cell array, the matrix of single values (text in the first two
%       columns, a number or [] in the others)
%   groups : struct array, one element per plot page, with the fields
%       * title    - char, the title of the page
%       * analyses - struct array of analyses of SQAT_GUI_extract, one per
%                    signal, all of the same kind
%       * names    - cell array of char, the name of each signal
%
% Author: Sergio Aguirre and Gil Felix Greco, October 2026
%
% AI disclosure: code development in October 2026 assisted
% by Claude Opus 5.5 (Anthropic). All codes were verified by
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

filename = char(filename);
if isfile(filename)
    delete(filename);                                  % pages are appended to the file
end
lines = [{'SQAT GUI report', ''}, ...
    cellfun(@(i, v) sprintf('%s: %s', i, v), S.Item(:)', S.Value(:)', 'UniformOutput', false), ...
    {'', 'Single values'}, il_matrix_lines(cols, data)];
lines = il_wrap(lines, 150);                           % what the page width holds in Courier 8
per_page = 44;
first = true;
for k = 1:per_page:numel(lines)
    f = il_page();
    ax = axes(f, 'Position', [0.04 0.03 0.92 0.94], 'Visible', 'off');
    page = lines(k:min(k + per_page - 1, numel(lines)));
    for j = 1:numel(page)
        text(ax, 0, 1 - (j - 1) / per_page, page{j}, 'FontName', 'Courier', 'FontSize', 8, ...
            'Interpreter', 'none', 'VerticalAlignment', 'top', ...
            'FontWeight', il_bold(first && j == 1));
    end
    il_write(f, filename, ~first, 'vector');
    first = false;
end
cmap = load('cmap_inferno.txt');
for g = groups(:)'
    f = il_page();
    A = g.analyses;
    if strcmp(A(1).kind, 'map')
        tl = tiledlayout(f, 'flow', 'TileSpacing', 'compact');
        lo = min(arrayfun(@(a) min(a.z(:), [], 'omitnan'), A));
        hi = max(arrayfun(@(a) max(a.z(:), [], 'omitnan'), A));
        for k = 1:numel(A)
            ax = nexttile(tl);
            surface(ax, A(k).x, A(k).y, zeros(numel(A(k).y), numel(A(k).x)), A(k).z.', 'EdgeColor', 'none');
            view(ax, 2);
            axis(ax, 'tight');
            colormap(ax, cmap);
            if isfinite(lo) && hi > lo
                clim(ax, [lo hi]);                     % one colour scale for every signal
            end
            if strcmp(A(k).bandscale, 'log')
                ax.YScale = 'log';
            end
            cb = colorbar(ax);
            cb.Label.String = A(k).zlabel;
            xlabel(ax, A(k).xlabel);
            ylabel(ax, A(k).ylabel);
            title(ax, g.names{k}, 'Interpreter', 'none');
        end
        title(tl, g.title, 'Interpreter', 'none');
    else
        ax = axes(f);
        hold(ax, 'on');
        for k = 1:numel(A)
            plot(ax, A(k).x, A(k).y);
        end
        hold(ax, 'off');
        grid(ax, 'on');
        if strcmp(A(1).kind, 'profile') && strcmp(A(1).bandscale, 'log')
            ax.XScale = 'log';
        end
        xlabel(ax, A(1).xlabel);
        ylabel(ax, A(1).ylabel);
        title(ax, g.title, 'Interpreter', 'none');
        legend(ax, g.names, 'Interpreter', 'none', 'Location', 'best');
    end
    il_write(f, filename, true, 'image');
end
end

function f = il_page()
% an A4 landscape page, white, off screen
f = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [0 0 29.7 21]);
end

function il_write(f, filename, append, content)
exportgraphics(f, filename, 'Append', append, 'ContentType', content, 'Resolution', 200, ...
    'Padding', 'figure');                              % the whole page, also when half of it is text
close(f);
end

function w = il_bold(tf)
w = 'normal';
if tf
    w = 'bold';
end
end

function out = il_wrap(lines, n)
% long lines (paths, parameters) cut into pieces of n characters, indented
out = {};
for k = 1:numel(lines)
    t = lines{k};
    out{end+1} = t(1:min(n, end)); %#ok<AGROW>
    for j = n + 1:n - 4:numel(t)
        out{end+1} = ['    ' t(j:min(j + n - 5, end))]; %#ok<AGROW>
    end
end
end

function lines = il_matrix_lines(cols, data)
% the matrix as text, one padded column each: numbers with 4 significant digits
txt = cell(size(data));
for k = 1:numel(data)
    v = data{k};
    if isnumeric(v) && isscalar(v)
        txt{k} = sprintf('%.4g', v);
    elseif ischar(v) || isstring(v)
        txt{k} = char(v);
    else
        txt{k} = '';
    end
end
txt = [cols(:)'; txt];
w = max(cellfun(@numel, txt), [], 1);
lines = cell(1, size(txt, 1));
for r = 1:size(txt, 1)
    parts = arrayfun(@(c) sprintf('%-*s', w(c), txt{r, c}), 1:size(txt, 2), 'UniformOutput', false);
    lines{r} = strjoin(parts, '  ');
end
end
