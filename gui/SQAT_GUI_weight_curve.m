function H_db = SQAT_GUI_weight_curve(f, fs, type)
% function H_db = SQAT_GUI_weight_curve(f, fs, type)
%
%   Response in dB of the weighting filter of SQAT_GUI_weight at the
%   frequencies f, for the level of a spectrogram.
%
% INPUT ARGUMENTS
%   f : frequencies (Hz)
%   fs : sampling frequency (Hz)
%   type : 'A', 'C' or 'Z'
%
% OUTPUTS
%   H_db : response (dB), with the shape of f
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
    H_db = zeros(size(f));
    return
end
[b, a] = SQAT_GUI_weight_filter(fs, type);
z = exp(-1j * 2*pi * f / fs);
H = polyval(fliplr(b(:)'), z) ./ polyval(fliplr(a(:)'), z);
H_db = 20*log10(abs(H));
end
