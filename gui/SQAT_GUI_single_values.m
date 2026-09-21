function T = SQAT_GUI_single_values(OUT)
% function T = SQAT_GUI_single_values(OUT)
%
%   Collects the single values of a SQAT output struct: its top-level real
%   numeric scalar fields, in field order. Nested structs (results of
%   sub-metrics), vectors, strings and the Bark step <dz> are left out.
%
% INPUT ARGUMENTS
%   OUT : output struct of a SQAT metric
%
% OUTPUTS
%   T : table with the variables Quantity (cell array of char) and Value
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

names = fieldnames(OUT);
keep = false(size(names));
values = nan(size(names));
for i = 1:numel(names)
    v = OUT.(names{i});
    if isnumeric(v) && isscalar(v) && isreal(v) && ~strcmp(names{i}, 'dz')
        keep(i) = true;
        values(i) = double(v);
    end
end
T = table(names(keep), values(keep), 'VariableNames', {'Quantity', 'Value'});
end
