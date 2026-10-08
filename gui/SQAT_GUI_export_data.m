function files = SQAT_GUI_export_data(entries, ids, filename, info, unit_of)
% function files = SQAT_GUI_export_data(entries, ids, filename, info, unit_of)
%
%   Writes the analyses of one signal and one metric to a workbook (.xlsx),
%   replacing any file of the same name. The first sheet, Info, describes
%   the signal, the analysis and every other sheet. A series or a profile
%   takes one sheet, with its axis and one column per channel; a map takes
%   one sheet per channel, with its time axis and one column per band; the
%   single values take one sheet, with one column per channel. The maps go
%   last: once a workbook holds a large map, adding even a small sheet takes
%   about as long as writing the map again. A table longer than a sheet
%   holds (1048575 rows under the header) goes to a CSV file next to the
%   workbook instead.
%
%   The name of a sheet is the name of the analysis without ' vs time',
%   after the channel for a map. Excel takes 31 characters, so the long
%   names are shortened ('Time-averaged' to 'Avg.', the weightings of the
%   psychoacoustic annoyance to w_FR and w_S); the column headers keep the
%   full names and units.
%
% INPUT ARGUMENTS
%   entries : the entries of SQAT_GUI for one signal and one metric, one per
%       channel (fields channel, '1', '2' or 'Binaural', analyses, the
%       output of SQAT_GUI_extract, and values, a table Quantity, Value)
%   ids : cell array, the ids of the analyses to write, and 'values' for the
%       single values
%   filename : char or string, path of the workbook to write
%   info : table (Item, Value) describing the signal and the analysis; the
%       channels and the sheets are added to it
%   unit_of : optional function handle, the unit of a single value from the
%       name of its quantity ('-' when not given)
%
% OUTPUTS
%   files : cell array with the paths written, the workbook first, then any
%       CSV file of a table too long for a sheet
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

if nargin < 5
    unit_of = @(q) '-';
end
filename = char(filename);
[folder, base] = fileparts(filename);
names = arrayfun(@(e) il_channel_name(e.channel), entries, 'UniformOutput', false);
sheets = struct('name', {}, 'what', {}, 'table', {});

% the single values, the series and the profiles, in the order asked
for id = ids(:)'
    if strcmp(id{1}, 'values')
        sheets(end+1) = struct('name', 'Single values', ...
            'what', 'Single values, one column per channel', 'table', il_values_table(entries, names, unit_of)); %#ok<AGROW>
        continue
    end
    [A, k] = il_analysis(entries, id{1});
    if isempty(A) || strcmp(A(1).kind, 'map')
        continue
    end
    T = SQAT_GUI_analysis_table(A(1));
    y = T.Properties.VariableNames{2};
    T.Properties.VariableNames{2} = sprintf('%s, %s', y, names{k(1)});
    for j = 2:numel(A)
        if ~isequal(A(j).x(:), A(1).x(:))
            error('SQAT_GUI:export_data', 'The channels of %s have different axes.', A(1).label);
        end
        T.(sprintf('%s, %s', y, names{k(j)})) = A(j).y(:);
    end
    sheets(end+1) = struct('name', il_sheet_name(A(1).label, '', 31), ...
        'what', sprintf('%s, one column per channel', A(1).label), 'table', T); %#ok<AGROW>
end

% the maps, one sheet per channel; the room left by the longest channel name
% is the same for all, so the channels of a map get the same name
for id = ids(:)'
    [A, k] = il_analysis(entries, id{1});
    if isempty(A) || ~strcmp(A(1).kind, 'map')
        continue
    end
    room = 31 - 1 - max(cellfun(@numel, names(k)));
    for j = 1:numel(A)
        sheets(end+1) = struct('name', il_sheet_name(A(j).label, names{k(j)}, room), ...
            'what', sprintf('%s, %s, one column per band', A(j).label, names{k(j)}), ...
            'table', SQAT_GUI_analysis_table(A(j))); %#ok<AGROW>
    end
end

new = matlab.lang.makeUniqueStrings({sheets.name}, {'Info'}, 31);
[sheets.name] = new{:};
fits = arrayfun(@(s) height(s.table) < 1048576 && width(s.table) <= 16384, sheets);
csv = cell(size(sheets));
for j = find(~fits)
    csv{j} = fullfile(folder, sprintf('%s_%s.csv', base, regexprep(sheets(j).name, '[^\w.-]+', '_')));
end

rows = [info.Item, info.Value; {'Channels', strjoin(names, ', ')}];
for j = 1:numel(sheets)
    if fits(j)
        rows(end+1, :) = {['Sheet ' sheets(j).name], sheets(j).what}; %#ok<AGROW>
    else
        [~, f, e] = fileparts(csv{j});
        rows(end+1, :) = {['File ' f e], sprintf('%s (%d rows, more than a sheet holds)', ...
            sheets(j).what, height(sheets(j).table))}; %#ok<AGROW>
    end
end
writetable(cell2table(rows, 'VariableNames', {'Item', 'Value'}), filename, 'Sheet', 'Info', ...
    'WriteMode', 'replacefile');
files = {filename};
for j = 1:numel(sheets)
    if fits(j)
        writetable(sheets(j).table, filename, 'Sheet', sheets(j).name);
    else
        writetable(sheets(j).table, csv{j}, 'WriteMode', 'overwrite');
        files{end+1} = csv{j}; %#ok<AGROW>
    end
end
end

%% -------------------------------------------------------------------------
function [A, k] = il_analysis(entries, id)
% the analysis id in each channel that holds it, and the indices of those channels
A = [];
k = [];
for j = 1:numel(entries)
    a = entries(j).analyses(strcmp({entries(j).analyses.id}, id));
    if ~isempty(a)
        A = [A, a(1)]; %#ok<AGROW>
        k(end+1) = j; %#ok<AGROW>
    end
end
end

function T = il_values_table(entries, names, unit_of)
% the quantities of all channels, in the order they come, one column per channel
q = {};
for e = entries(:)'
    q = [q; e.values.Quantity(~ismember(e.values.Quantity, q))]; %#ok<AGROW>
end
V = nan(numel(q), numel(entries));
for j = 1:numel(entries)
    [in, at] = ismember(entries(j).values.Quantity, q);
    V(at(in), j) = entries(j).values.Value(in);
end
T = [table(q, cellfun(unit_of, q, 'UniformOutput', false), 'VariableNames', {'Quantity', 'Unit'}), ...
    array2table(V, 'VariableNames', strcat({'Value, '}, names))];
end

function name = il_channel_name(channel)
% ch1, ch2 or binaural, as the plots of the GUI name them
if strcmp(channel, 'Binaural')
    name = 'binaural';
else
    name = ['ch' channel];
end
end

function name = il_sheet_name(label, channel, room)
% the label without ' vs time', shortened step by step until it fits in room,
% after the channel of a map; Excel takes 31 characters and none of : \ / ? * [ ]
short = {'Weighting of fluctuation strength and roughness', 'w_FR'
         'Weighting of sharpness and loudness',             'w_S'
         'Time-averaged',                                   'Avg.'
         'Specific',                                        'Spec.'
         'specific',                                        'spec.'
         'fluctuation strength',                            'fluct. strength'};
name = strrep(label, ' vs time', '');
for j = 1:size(short, 1)
    if numel(name) <= room
        break
    end
    name = strrep(name, short{j, 1}, short{j, 2});
end
name = strtrim([channel ' ' name]);
name = regexprep(name, '[:\\/?*\[\]]', '-');
name = name(1:min(end, 31));
end
