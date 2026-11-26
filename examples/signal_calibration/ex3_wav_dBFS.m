%   Script ex3_wav_dBFS
%
% How to store and read a sound in a .wav file without losing its level:
% the full-scale level (dBFS) convention used in SQAT
%
% A .wav file stores numbers between -1 and +1. In SQAT, the full-scale
% level, dBFS, is the sound pressure level (SPL) of a signal with an rms
% value of 1. A signal p, in Pa, is stored as numbers x, and read back, as
%
%   x = p / 10^((dBFS-94)/20)   (writing)
%   p = x * 10^((dBFS-94)/20)   (reading)
%
% By default, dBFS = 94 dB SPL, i.e. 1 = 1 Pa. Any value outside +-1 is
% clipped when writing the file. Therefore, a sound with peaks above 1 Pa
% (e.g. a tone above 91 dB SPL) must be stored using a higher dBFS. The
% dBFS used to write the file must be known to read it back correctly.
%
% This script:
% - generates a 1-kHz tone with 100 dB SPL (peaks of 2.8 Pa)
% - writes it to a 16-bit .wav file using dBFS = 94 dB and dBFS = 110 dB
% - reads both files back, and compares their levels and spectra with
%   the original signal
%
%  HOW TO RUN THIS CODE: this is a standalone code. Therefore, no additional
%  steps are necessary to run this code. The .wav files are written to the
%  temporary folder of your system and deleted at the end.
%
% Author: Gil Felix Greco, Braunschweig, 26.11.2026
%
% AI disclosure: code development in November 2026 assisted by Claude
% Opus 5.5 (Anthropic). All codes were verified by the authors.
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
clc; clear; close all;

save_figs = 0; % set to 1 to save the figure in the figs folder

%% Sound pressure signal (Pa)

fs = 48000;    % sampling frequency (Hz)
p0 = 20e-6;    % reference sound pressure (Pa)
L_true = 100;  % level of the tone (dB SPL)

t = (0:1/fs:1-1/fs)';                                 % 1 s
p = p0 * 10^(L_true/20) * sqrt(2) * sin(2*pi*1000*t); % 1-kHz tone (Pa)

%% Write and read back using two full-scale levels

dBFS = [94 110];  % full-scale levels (dB SPL)

x = cell(1,2);
p_read = cell(1,2);
for i = 1:2
    fname = fullfile(tempdir, sprintf('ex3_dBFS_%g.wav', dBFS(i)));

    audiowrite(fname, p / 10^((dBFS(i)-94)/20), fs, 'BitsPerSample', 16);  % Pa -> numbers (clipped to +-1)

    [x{i}, fs] = audioread(fname);          % numbers, as stored in the file
    p_read{i} = x{i} * 10^((dBFS(i)-94)/20); % numbers -> Pa
    delete(fname);

    fprintf('dBFS = %g dB: max|x| = %.2f, level read back = %.1f dB SPL (true: %.1f dB SPL)\n', ...
            dBFS(i), max(abs(x{i})), 20*log10(rms(p_read{i})/p0), L_true);
end

% a common mistake: reading the file written with dBFS = 110 dB using the default 94 dB
fprintf('File written with dBFS = 110 dB, read with dBFS = 94 dB: %.1f dB SPL\n', ...
        20*log10(rms(x{2})/p0));

%% Spectra at 1 kHz and its harmonics (rms level per frequency bin)

fh = (1:10)' * 1000;  % 1 kHz and its harmonics up to 10 kHz (Hz)
spec = @(s) 20*log10( abs(fft(s)) * sqrt(2)/length(s) / p0 );  % 1-s signals: bin k+1 is k Hz
S_true = spec(p);         S_true = S_true(fh+1);
S94    = spec(p_read{1}); S94    = S94(fh+1);
S110   = spec(p_read{2}); S110   = S110(fh+1);

% levels below 0 dB SPL are not shown
S_true(S_true < 0) = NaN; S94(S94 < 0) = NaN; S110(S110 < 0) = NaN;

%% Plot

figure('Name', 'Full-scale level (dBFS) and .wav files');

subplot(2,1,1)
idx = t <= 3e-3; % first 3 ms
plot(t(idx)*1000, x{1}(idx), 'r'); hold on;
plot(t(idx)*1000, x{2}(idx), 'b');
yline(1, 'k--'); yline(-1, 'k--');
ylim([-1.5 1.5]);
xlabel('Time, $t$ (ms)', 'Interpreter', 'Latex');
ylabel('Value stored, $x$', 'Interpreter', 'Latex');
legend('dBFS = 94 dB (clipped)', 'dBFS = 110 dB', 'Location', 'southwest');
title('(a) Numbers stored in the .wav file', 'Interpreter', 'Latex');

subplot(2,1,2)
bar(fh/1000, S_true, 0.4, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none'); hold on;
stem(fh/1000, S94, 'r', 'filled', 'LineWidth', 1.5);
plot(fh/1000, S110, 'bo', 'MarkerSize', 9, 'LineWidth', 1.5);
xlim([0 10.5]); ylim([0 150]); xticks(1:10); yticks(0:20:100);
xlabel('Frequency, $f$ (kHz)', 'Interpreter', 'Latex');
ylabel('$L$ (dB re 20 $\mu$Pa)', 'Interpreter', 'Latex');
legend('original', 'dBFS = 94 dB', 'dBFS = 110 dB', 'Location', 'northeast');
title('(b) Spectra of the signals read back', 'Interpreter', 'Latex');

set(gcf, 'color', 'w');

if save_figs == 1
    figs_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end
    saveas(gcf, [figs_dir 'ex3_wav_dBFS'], 'png');
end
