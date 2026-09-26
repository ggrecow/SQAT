function [t, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, n_frames, f_min, bins_per_octave)
% function [t, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, n_frames, f_min, bins_per_octave)
%
%   Enhanced spectrogram with no window to choose: five Gaussian windows of
%   16, 32, 64, 128 and 256 ms are reassigned (Auger and Flandrin, 1995),
%   each reassigned map is smoothed with a Gaussian (2 cells in time and in
%   frequency) and normalised, and the maps are combined by geometric mean,
%   so that only the energy that all the windows place at the same point of
%   the time-frequency plane remains (the idea of Cheung and Lim, 1991,
%   applied to reassigned maps). Only base MATLAB is used.
%
% INPUT ARGUMENTS
%   x : signal (Pa)
%   fs : sampling frequency (Hz)
%   n_frames : largest number of time columns of the output (default 6000; the step is 2 ms or more)
%   f_min : lowest frequency of the output (Hz, default 20)
%   bins_per_octave : frequency resolution of the output grid (default 96)
%
% OUTPUTS
%   t : [1xN] time of each column (s)
%   f : [Mx1] frequencies (Hz), on a logarithmic grid
%   L : [MxN] level (dB); the maximum is aligned with the maximum of the
%       level of the 64 ms window (dB SPL), so the colour scale is relative
%   info : struct with the windows (ms) and the time step (s)
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September 2026 assisted
% by Claude Sonnet 5 (Anthropic). All codes were verified by
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

if nargin < 3 || isempty(n_frames), n_frames = 6000; end
if nargin < 4 || isempty(f_min), f_min = 20; end
if nargin < 5 || isempty(bins_per_octave), bins_per_octave = 96; end
x = x(:);
n_x = numel(x);
durations = [16 32 64 128 256] * 1e-3;
Ns = 2 * round(durations * fs / 2);                    % even window lengths (samples)
hop_out = max([1, round(0.002 * fs), ceil(n_x / n_frames)]);   % output time step (samples): 2 ms or more
M = numel(0:hop_out:n_x-1);
n_oct = log2((fs/2) / f_min);
F = ceil(n_oct * bins_per_octave) + 1;
sigma = 2;                                             % smoothing (cells)
k = ceil(3 * sigma);
b = exp(-(-k:k).^2 / (2 * sigma^2));
b = b / sum(b);
log_sum = 0;
for w_i = 1:numel(Ns)
    R = il_reassigned(x, fs, Ns(w_i), hop_out, M, F, f_min, bins_per_octave);
    R = conv2(b(:), b(:)', R, 'same');
    R = R / sum(R(:));
    log_sum = log_sum + log(R + 1e-6 * max(R(:)));
end
C = exp(log_sum / numel(Ns));
C = C / sum(C(:));
% level: 10 log10 of the energy density, its peak aligned with the peak of the 64 ms spectrogram (relative colour scale)
[~, ~, L_ref] = SQAT_GUI_spectrogram(x, fs, 'hann', max(6, min(16, round(log2(0.064 * fs)))), 50);
L = (10*log10(C / max(C(:)) + 1e-12))' + max(L_ref, [], 'all');
t = (0:M-1) * hop_out / fs;
f = f_min * 2 .^ ((0:F-1)' / bins_per_octave);
info = struct('windows_ms', 1e3 * Ns / fs, 'hop', hop_out / fs);
end

function R = il_reassigned(x, fs, N, hop_out, M, F, f_min, bpo)
% reassigned spectrogram of one Gaussian window (sigma = N/8 samples, unit energy) accumulated on the output grid
n_x = numel(x);
hop = max(1, round(N/16));
pad = N/2;
xp = [zeros(pad, 1); x; zeros(pad + N, 1)];
n = (-N/2:N/2-1)';
s = N/8;
w = exp(-n.^2 / (2*s^2));
w = w / norm(w);
wt = n .* w;
wd = -(n / s^2) .* w;
h = N/2;
Fk = N/2 + 1;
fbin = (0:Fk-1)' * fs / N;
starts = 0:hop:n_x-1;
n_fr = numel(starts);
R = zeros(M, F);
chunk = max(1, floor(2e6 / N));
for a = 1:chunk:n_fr
    j = a:min(a + chunk - 1, n_fr);
    S = xp((starts(j) + pad) + (-h+1:h)');             % frames centred on sample starts(j)
    X0 = fft([S(h+1:N, :) .* w(h+1:N); S(1:h, :) .* w(1:h)]);
    X1 = fft([S(h+1:N, :) .* wt(h+1:N); S(1:h, :) .* wt(1:h)]);
    X2 = fft([S(h+1:N, :) .* wd(h+1:N); S(1:h, :) .* wd(1:h)]);
    Xh = X0(1:Fk, :);
    P2 = abs(Xh).^2;
    r = conj(Xh) ./ max(P2, eps);
    t_hat = starts(j) + real(X1(1:Fk, :) .* r);        % samples
    f_hat = fbin - (fs / (2*pi)) * imag(X2(1:Fk, :) .* r);   % Hz
    it = round(t_hat / hop_out) + 1;
    jf = round(bpo * log2(max(f_hat, eps) / f_min)) + 1;
    ok = P2 > 1e-10 * max(P2, [], 'all') & it >= 1 & it <= M & jf >= 1 & jf <= F & f_hat > 0;
    R = R + accumarray([it(ok) jf(ok)], P2(ok) * hop, [M F]);
end
end
