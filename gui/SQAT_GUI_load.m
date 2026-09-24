function [insig, fs, nch] = SQAT_GUI_load(wavfilename, dBFS, channel)
% function [insig, fs, nch] = SQAT_GUI_load(wavfilename, dBFS, channel)
%
%   Reads one or more channels of a .wav file and calibrates them to pascals with the
%   convention of the SQAT <_from_wavfile> functions: a full-scale amplitude
%   of 1 corresponds to dBFS dB SPL, gain_factor = 10^((dBFS-94)/20).
%
% INPUT ARGUMENTS
%   wavfilename : char or string, path of the .wav file
%   dBFS : number, dB SPL of a full-scale amplitude (94 means 1 = 1 Pa)
%   channel : integer or vector of integers, channels to read (1 for mono files)
%
% OUTPUTS
%   insig : [Nx1] calibrated signal (Pa), or [Nxk] for k channels
%   fs : sampling frequency (Hz)
%   nch : number of channels in the file
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

[y, fs] = audioread(wavfilename);
nch = size(y, 2);
if isempty(channel) || any(channel ~= round(channel)) || any(channel < 1) || any(channel > nch)
    error('SQAT_GUI:channel', 'Channel %s requested, but %s has %d channel(s).', ...
          mat2str(channel), wavfilename, nch);
end
gain_factor = 10^((dBFS-94)/20);
insig = y(:, channel) * gain_factor;
end
