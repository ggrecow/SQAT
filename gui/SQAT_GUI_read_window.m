function w = SQAT_GUI_read_window(file)
% function w = SQAT_GUI_read_window(file)
%
%   Reads a window from a file: a .txt, .csv or .dat file with the samples
%   in a row or a column, or a .mat file (the first real vector of at
%   least two samples it holds).
%
% INPUT ARGUMENTS
%   file : path of the file
%
% OUTPUTS
%   w : [nx1] samples of the window
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

[~, ~, ext] = fileparts(file);
switch lower(ext)
    case {'.txt', '.csv', '.dat'}
        try
            v = readmatrix(file);
        catch err
            error('SQAT_GUI:window', 'The file %s could not be read: %s', file, err.message);
        end
    case '.mat'
        S = load(file);
        v = [];
        for name = fieldnames(S)'
            c = S.(name{1});
            if isnumeric(c) && isreal(c) && isvector(c) && numel(c) >= 2
                v = c;
                break
            end
        end
    otherwise
        error('SQAT_GUI:window', 'The file type %s is not supported: use .txt, .csv, .dat or .mat.', ext);
end
if isempty(v) || ~isnumeric(v) || ~isreal(v) || ~isvector(v) || numel(v) < 2 || any(~isfinite(v))
    error('SQAT_GUI:window', 'The file %s does not hold a window: it needs one vector of at least 2 finite samples.', file);
end
w = double(v(:));
end
