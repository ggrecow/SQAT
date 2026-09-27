function SQAT_GUI_export(T, filename, S)
% function SQAT_GUI_export(T, filename, S)
%
%   Writes the results table of SQAT_GUI to a spreadsheet (.xlsx) or a text
%   file (.csv), replacing any file of the same name. The settings of the
%   run, when given, go to a second sheet named Settings (.xlsx) or to a
%   second file, name_settings.csv (.csv).
%
% INPUT ARGUMENTS
%   T : results table (Signal, File, Analysis, Metric, Channel, Quantity,
%       Value, Unit, Cal_dB_SPL, Parameters, Path)
%   filename : char or string, path of the file to write
%   S : optional table (Item, Value) with the settings of the run
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

if nargin < 3
    S = [];
end
filename = char(filename);
[folder, base, ext] = fileparts(filename);
if strcmpi(ext, '.xlsx')
    writetable(T, filename, 'Sheet', 'Results', 'WriteMode', 'replacefile');
    if ~isempty(S)
        writetable(S, filename, 'Sheet', 'Settings');
    end
else
    writetable(T, filename, 'WriteMode', 'overwrite');
    if ~isempty(S)
        writetable(S, fullfile(folder, [base '_settings' ext]), 'WriteMode', 'overwrite');
    end
end
end
