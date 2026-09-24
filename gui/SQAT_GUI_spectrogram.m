function [t, f, L, info] = SQAT_GUI_spectrogram(x, fs, win, degree, overlap)
% function [t, f, L, info] = SQAT_GUI_spectrogram(x, fs, win, degree, overlap)
%
%   Short-time spectrum in dB SPL of a signal in pascals. The amplitude is
%   scaled by the sum of the window, so that a sine on a bin centre reads
%   its RMS level with any window. The number of frames is limited to keep
%   the result in memory (4000 frames, and 8e6 values in all): the overlap
%   is then reduced and info tells it.
%
% INPUT ARGUMENTS
%   x : signal (Pa)
%   fs : sampling frequency (Hz)
%   win : name of a window of SQAT_GUI_window, or a vector with the window
%         (resampled to the FFT size when its length differs)
%   degree : the FFT has 2^degree points
%   overlap : overlap of the frames in percent, from 0 to 95
%
% OUTPUTS
%   t : [1xN] time of the centre of each frame (s)
%   f : [(n/2+1)x1] frequencies (Hz)
%   L : [(n/2+1)xN] level (dB SPL)
%   info : struct with n_fft, hop (samples), overlap (percent, as used) and
%          limited (true when the number of frames was limited)
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

n_fft = 2^degree;
if ischar(win) || isstring(win)
    w = SQAT_GUI_window(char(win), n_fft);
else
    w = SQAT_GUI_window('custom', n_fft, win);
end
x = x(:);
if numel(x) < n_fft
    x(n_fft) = 0;
end
n_bins = n_fft/2 + 1;
max_frames = max(1, min(4000, floor(8e6 / n_bins)));
hop = max(1, round(n_fft * (1 - overlap/100)));
n_frames = floor((numel(x) - n_fft) / hop) + 1;
limited = n_frames > max_frames;
if limited
    hop = ceil((numel(x) - n_fft) / (max_frames - 1));
    n_frames = floor((numel(x) - n_fft) / hop) + 1;
end
starts = 1 + (0:n_frames-1) * hop;
L = zeros(n_bins, n_frames);
chunk = max(1, floor(2e6 / n_fft));                   % frames per FFT call, to bound the memory
for k = 1:chunk:n_frames
    idx = k:min(k + chunk - 1, n_frames);
    frames = x(starts(idx) + (0:n_fft-1)') .* w;
    X = fft(frames);
    A = abs(X(1:n_bins, :)) * 2 / sum(w);             % peak amplitude per bin
    A([1 end], :) = A([1 end], :) / 2;
    L(:, idx) = 20*log10(max(A / sqrt(2), eps) / 2e-5);
end
t = (starts - 1 + n_fft/2) / fs;
f = (0:n_fft/2)' * fs / n_fft;
info = struct('n_fft', n_fft, 'hop', hop, 'overlap', 100 * (1 - hop/n_fft), 'limited', limited);
end
