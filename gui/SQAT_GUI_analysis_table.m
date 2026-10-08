function T = SQAT_GUI_analysis_table(A)
% function T = SQAT_GUI_analysis_table(A)
%
%   One analysis of SQAT_GUI_extract as a table, ready to be written to a
%   file: a series or a profile gives two columns, its axis and its values;
%   a map gives its time axis and one column per band. The column names
%   carry the units, so the table describes itself: a column of a map
%   reads, for example, 'Specific loudness (sone/Bark) at 12.5 Bark'. The
%   braces of the TeX labels go: sone_{HMS} reads sone_HMS.
%
% INPUT ARGUMENTS
%   A : one element of the output of SQAT_GUI_extract (fields kind, x, y,
%       z, xlabel, ylabel and zlabel)
%
% OUTPUTS
%   T : table with numel(A.x) rows: two columns for a series or a profile,
%       1 + numel(A.y) columns for a map
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

A.xlabel = erase(A.xlabel, {'{', '}'});
A.ylabel = erase(A.ylabel, {'{', '}'});
A.zlabel = erase(A.zlabel, {'{', '}'});
switch A.kind
    case {'series', 'profile'}
        T = table(A.x(:), A.y(:), 'VariableNames', {A.xlabel, A.ylabel});
    case 'map'
        unit = il_unit(A.ylabel);                     % the unit of the bands: Bark, Hz
        bands = arrayfun(@(b) strtrim(sprintf('%s at %g %s', A.zlabel, b, unit)), A.y(:)', ...
            'UniformOutput', false);
        T = array2table([A.x(:), A.z], 'VariableNames', ...
            matlab.lang.makeUniqueStrings([{A.xlabel}, bands]));   % two bands that print alike stay apart
    otherwise
        error('SQAT_GUI:analysis_table', 'Unknown kind of analysis: %s.', A.kind);
end
end

function u = il_unit(label)
% the unit in the last parentheses of an axis label: 'Frequency (Hz)' gives 'Hz'
u = regexp(label, '\(([^()]*)\)\s*$', 'tokens', 'once');
if isempty(u)
    u = '';
else
    u = u{1};
end
end
