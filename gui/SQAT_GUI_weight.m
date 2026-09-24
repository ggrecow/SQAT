function y = SQAT_GUI_weight(x, fs, type)
% function y = SQAT_GUI_weight(x, fs, type)
%
%   Applies a frequency weighting to a signal with the filters of
%   SQAT_GUI_weight_filter (IEC 61672-1): the design of Gen_weighting_filters,
%   the filters of the sound level meter of SQAT, written without a toolbox.
%   Z is flat and returns the signal as it is.
%
% INPUT ARGUMENTS
%   x : signal, one column per channel
%   fs : sampling frequency (Hz)
%   type : 'A', 'C' or 'Z'
%
% OUTPUTS
%   y : weighted signal
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

if ~ismember(type, {'A', 'C', 'Z'})
    error('SQAT_GUI:weighting', 'Unknown weighting %s: use A, C or Z.', type);
end
if strcmp(type, 'Z')
    y = x;
    return
end
[b, a] = SQAT_GUI_weight_filter(fs, type);
y = filter(b, a, x);
end
