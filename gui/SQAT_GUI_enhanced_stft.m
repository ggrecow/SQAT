function [t, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, smoothing, n_frames, f_min, bins_per_octave, preview, ref)
% function [t, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, smoothing, n_frames, f_min, bins_per_octave, preview, ref)
%
%   Enhanced spectrogram with no window to choose: nine Blackman-Harris
%   windows from 8 to 512 ms (geometric spacing) are reassigned (Auger and
%   Flandrin, 1995), the energy of each frame is held over the frame step
%   of its window, each map is smoothed with the chosen smoothing, and
%   the maps are combined by geometric mean, so that only the energy that
%   all the windows place at the same point of the time-frequency plane
%   remains (the idea of Cheung and Lim, 1991, applied to reassigned maps).
%   Only base MATLAB is used.
%
% INPUT ARGUMENTS
%   x : signal (Pa)
%   fs : sampling frequency (Hz)
%   smoothing : 'readable' (default): Gaussian smoothing of 4 ms in time
%               and 1.45 % of the frequency, for continuous lines;
%               'sharp': 1 ms and 1 Hz, for the thinnest lines; a cell
%               array of both gives both maps for the time of one, since
%               the reassignment does not depend on the smoothing
%   n_frames : largest number of time columns of the output (default 2000,
%              about the width of the screen; the time step is 1 ms or more)
%   f_min : lowest frequency of the output (Hz, default 20)
%   bins_per_octave : frequency resolution of the output grid (default 192)
%   preview : true lets each window step by up to one output column (at most
%             half its length) instead of an eighth of its length: on a long
%             signal, whose columns are long, the short windows compute far
%             fewer frames, and each column gets fewer of them, so the map is
%             an approximation (clicks come out weaker). When the columns are
%             1 ms, as in a zoomed excerpt, nothing changes. Default false.
%   ref : info.ref of the map of the whole signal, for an excerpt of it: the
%         noise threshold and the floor of the geometric mean come from the
%         whole signal, so the excerpt reads as the whole map. Default [] (from
%         x itself)
%
% OUTPUTS
%   t : [1xN] time of each column (s)
%   f : [Mx1] frequencies (Hz), on a logarithmic grid
%   L : [MxN] level of each cell (dB SPL): the mean square pressure that the
%       cell holds, so the cells of a column add up to the sound pressure level
%       of that instant. A steady tone spread by the smoothing over several
%       cells reads lower in each of them (about 10 dB with 'readable').
%       With a cell array of smoothings, a cell array of maps in that order
%   info : struct with the windows (ms), the time step (s), the smoothing and
%          ref, the references of this map (one per smoothing), for its excerpts
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

if nargin < 3 || isempty(smoothing), smoothing = 'readable'; end
if nargin < 4 || isempty(n_frames), n_frames = 2000; end
if nargin < 5 || isempty(f_min), f_min = 20; end
if nargin < 6 || isempty(bins_per_octave), bins_per_octave = 192; end
if nargin < 7 || isempty(preview), preview = false; end
if nargin < 8, ref = []; end
if f_min >= fs/2
    error('SQAT_GUI_enhanced_stft:f_min', 'f_min (%g Hz) must be below fs/2 (%g Hz)', f_min, fs/2);
end
x = x(:);
n_x = numel(x);
durations = 8e-3 * 64 .^ ((0:8) / 8);                 % nine windows, 8 to 512 ms
Ns = 2 * round(durations * fs / 2);                    % even window lengths (samples)
% Each window steps by N/step_div, and the energy of each frame is held over
% that step. 8 is the default. N/16 narrows clicks and onsets (about 5 ms
% instead of 6) and separates close tones better, at about twice the time;
% N/32 narrows clicks further (about 4 ms) at about 3.3 times the time of
% N/8 and leaves more spurious points in broadband noise. Below 8 the frames
% are too sparse for the hold to be short.
step_div = 8;
hop_out = max([1, round(0.001 * fs), ceil(n_x / n_frames)]);   % output time step (samples): 1 ms or more
steps = max(1, round(Ns / step_div));                 % frame step of each window (samples)
if preview
    steps = max(steps, min(round(Ns / 2), hop_out));   % up to one output column, at most N/2
end
M = numel(0:hop_out:n_x-1);
n_oct = log2((fs/2) / f_min);
F = ceil(n_oct * bins_per_octave) + 1;
f = f_min * 2 .^ ((0:F-1)' / bins_per_octave);
cell_f = f * (2^(1/bins_per_octave) - 1);              % width of each frequency cell (Hz)
modes = cellstr(smoothing);
gts = cell(size(modes));
sig_fs = cell(size(modes));
for k = 1:numel(modes)
    switch modes{k}
        case 'readable'
            sig_t = 0.004 * fs / hop_out;              % 4 ms, in cells
            sig_f = 2 * (2^(1/96) - 1) * f ./ cell_f;  % 1.45 % of the frequency, in cells of each row
        case 'sharp'
            sig_t = 0.001 * fs / hop_out;              % 1 ms
            sig_f = 1 ./ cell_f;                       % 1 Hz
        otherwise
            error('SQAT_GUI_enhanced_stft:smoothing', 'smoothing must be ''readable'' or ''sharp''');
    end
    sig_t = max(sig_t, 0.3);
    sig_fs{k} = max(sig_f, 0.3);
    kt = ceil(3 * sig_t);
    gt = exp(-(-kt:kt)' .^ 2 / (2 * sig_t^2));
    gts{k} = gt / sum(gt);
end
% noise threshold of the reassignment: the largest frame energy of the whole signal, which
% bounds |X|^2 since the windows have unit energy (a maximum per block of frames would
% depend on the block, and a quiet block would let noise through)
if isempty(ref)
    e_max = arrayfun(@(N) max(movsum(x.^2, N), [], 'all'), Ns);
else
    e_max = ref(1).e_max;
end
log_sum = num2cell(zeros(1, numel(modes)));      % one geometric mean per smoothing
r_max = zeros(numel(Ns), numel(modes));          % the floor of each map
for w_i = 1:numel(Ns)
    R = il_reassigned(x, fs, Ns(w_i), steps(w_i), hop_out, M, F, f_min, bins_per_octave, e_max(w_i));
    R = R / hop_out;                              % energy density: the same scale at any column width
    % the energy of a frame is held over its own step, so that the long windows leave no gaps
    n_hold = 2 * floor(steps(w_i) / hop_out / 2) + 1;
    if n_hold > 1
        R = conv2(ones(n_hold, 1) / n_hold, 1, R, 'same');
    end
    for k = 1:numel(modes)
        Rk = il_smooth(conv2(gts{k}, 1, R, 'same'), sig_fs{k});
        if isempty(ref)
            r_max(w_i, k) = max(Rk(:));
        else
            r_max(w_i, k) = ref(k).r_max(w_i);
        end
        log_sum{k} = log_sum{k} + log(Rk + 1e-6 * r_max(w_i, k));
    end
end
% level: the maps are mean square pressures (Pa^2), so their geometric mean is one too
L = cell(size(modes));
refs = struct('e_max', {}, 'r_max', {});
for k = 1:numel(modes)
    C = exp(log_sum{k} / numel(Ns));
    L{k} = 10*log10(C / (2e-5)^2 + 1e-12)';            % dB SPL; a silent signal sits at -120 dB
    refs(k) = struct('e_max', e_max, 'r_max', r_max(:, k)');
end
if ~iscell(smoothing)
    L = L{1};
end
t = (0:M-1) * hop_out / fs;
info = struct('windows_ms', 1e3 * Ns / fs, 'hop', hop_out / fs, 'smoothing', {smoothing}, 'preview', preview, ...
    'ref', refs);
end

function R = il_smooth(R, sig_f)
% Gaussian smoothing along the frequency (columns) with a width that changes with the row of the output grid
lev = round(log(sig_f) / log(1.25));                   % geometric levels of the width, one kernel per level
out = R;
n_f = size(R, 2);
for l = unique(lev)'
    cols = find(lev == l);
    s = 1.25^l;
    k = ceil(3 * s);
    g = exp(-(-k:k).^2 / (2 * s^2));
    g = g / sum(g);
    c1 = max(1, cols(1) - k);
    c2 = min(n_f, cols(end) + k);
    sm = conv2(1, g, R(:, c1:c2), 'same');
    out(:, cols) = sm(:, cols - c1 + 1);
end
R = out;
end

function R = il_reassigned(x, fs, N, hop, hop_out, M, F, f_min, bpo, e_max)
% reassigned spectrogram of one Blackman-Harris window (N samples, unit energy,
% frame step hop) accumulated on the output grid
n_x = numel(x);
pad = N/2;
xp = [zeros(pad, 1); x; zeros(pad + N, 1)];
n = (-N/2:N/2-1)';
h = N/2;
kk = n + h;
a = [0.35875 0.48829 0.14128 0.01168];                 % 4-term Blackman-Harris
w = a(1) - a(2)*cos(2*pi*kk/N) + a(3)*cos(4*pi*kk/N) - a(4)*cos(6*pi*kk/N);
wd = a(2)*(2*pi/N)*sin(2*pi*kk/N) - a(3)*(4*pi/N)*sin(4*pi*kk/N) + a(4)*(6*pi/N)*sin(6*pi*kk/N);
c = 1 / norm(w);
w = w * c;
wd = wd * c;
wt = n .* w;
Fk = N/2 + 1;
fbin = (0:Fk-1)' * fs / N;
starts = 0:hop:n_x-1;
n_fr = numel(starts);
R = zeros(M, F);
chunk = max(1, floor(2e6 / N));
lin = {};                                              % cells and energies of several blocks, added in one go
val = {};
n_buf = 0;
for b = 1:chunk:n_fr
    j = b:min(b + chunk - 1, n_fr);
    S = xp((starts(j) + pad) + (-h+1:h)');             % frames centred on sample starts(j)
    % the frames are not rotated to put their centre first: that only multiplies
    % every spectrum by (-1)^k, which cancels in the ratios below
    Xh = fft(S .* w);
    Xh = Xh(1:Fk, :);
    Xt = fft(S .* wt);
    Xd = fft(S .* wd);
    Xd = Xd(1:Fk, :);
    P2 = abs(Xh).^2;
    r = conj(Xh) ./ max(P2, eps);
    t_hat = starts(j) + real(Xt(1:Fk, :) .* r);        % samples
    f_hat = fbin - (fs / (2*pi)) * imag(Xd .* r);      % Hz
    it = round(t_hat / hop_out) + 1;
    jf = round(bpo * log2(max(f_hat, eps) / f_min)) + 1;
    ok = P2 > 1e-10 * e_max & it >= 1 & it <= M & jf >= 1 & jf <= F & f_hat > 0;
    lin{end+1} = it(ok) + (jf(ok) - 1) * M;            %#ok<AGROW>
    val{end+1} = P2(ok) * (2 * hop / N);               %#ok<AGROW> one-sided power of a unit-energy window: mean square
    n_buf = n_buf + nnz(ok);
    if n_buf > 1e7 || j(end) == n_fr                   % about 160 MB at most
        R = R + reshape(accumarray(vertcat(lin{:}), vertcat(val{:}), [M*F 1]), M, F);
        lin = {};
        val = {};
        n_buf = 0;
    end
end
end
