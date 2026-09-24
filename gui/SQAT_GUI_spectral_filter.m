function y = SQAT_GUI_spectral_filter(x, fs, boxes)
% function y = SQAT_GUI_spectral_filter(x, fs, boxes)
%
%   Removes rectangles of the time-frequency plane from a signal: the
%   spectrum of the frames whose centre falls in a box's time span is set
%   to zero between the box's frequencies, and the signal is rebuilt by
%   overlap-add (Hann analysis and synthesis windows, 2048 points, hop of a
%   quarter). Frames outside every box pass unchanged, and a signal with no
%   box comes back as it was.
%
% INPUT ARGUMENTS
%   x : signal
%   fs : sampling frequency (Hz)
%   boxes : [Bx4] matrix, one row [t1 t2 f1 f2] per box (s and Hz)
%
% OUTPUTS
%   y : [Nx1] filtered signal
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

x = x(:);
if isempty(boxes)
    y = x;
    return
end
n = 2048;
hop = n / 4;
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/n);
n_x = numel(x);
n_frames = ceil((n_x + n) / hop) + 1;                 % frames that cover the signal with n samples on each side
xp = [zeros(n, 1); x; zeros(n_frames * hop + n - n_x - n, 1)];
yp = zeros(size(xp));
f_bin = (0:n-1)' * fs / n;
f_bin = min(f_bin, fs - f_bin);                       % frequency of each bin, the mirror half included
for k = 1:n_frames
    seg = (k-1)*hop + (1:n);
    t_c = ((k-1)*hop + n/2 - n) / fs;                 % time of the centre of the frame in the signal
    drop = false(n, 1);
    for b = 1:size(boxes, 1)
        if t_c >= boxes(b, 1) && t_c <= boxes(b, 2)
            drop = drop | (f_bin >= boxes(b, 3) & f_bin <= boxes(b, 4));
        end
    end
    fr = xp(seg) .* w;
    if any(drop)
        X = fft(fr);
        X(drop) = 0;
        fr = real(ifft(X));
    end
    yp(seg) = yp(seg) + fr .* w;
end
y = yp(n + (1:n_x)) / 1.5;                            % the squared Hann window sums to 1.5 at this hop
end
