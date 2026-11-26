%   Script ex4_level_adjustment
%
% What to do when there is no information to calibrate a signal, e.g. a
% sound file downloaded from the internet: level adjustment
%
% If neither a calibrator recording nor the full-scale level (dBFS) of a
% file is known, the true level of the sound cannot be recovered. What can
% be done is to assume a level, L_out, and adjust the signal to it. With
% calibrate.m, the signal itself (or a part of it) is used as reference:
%
%   p = calibrate(x, x_ref, L_out)
%
%   x     : signal to be adjusted (numbers)
%   x_ref : part of x whose rms value should correspond to L_out
%   L_out : assumed level (dB SPL)
%   p     : adjusted signal (Pa)
%
% Every result obtained from p is valid only under this assumption. Level
% adjustment is mainly useful to compare signals at the same level, or to
% set the level of synthetic signals.
%
% This script:
% - generates a "downloaded" file: a noise burst with silence before and
%   after it, normalized to a peak value of 0.99, as done by many audio
%   editors (its original level is unknown)
% - adjusts its level to L_out = 65 dB SPL using as reference (i) the
%   whole file and (ii) only the noise burst, and compares the results
%
% Take-home message: calibrate.m sets the rms value of x_ref to L_out. If
% x_ref includes silence, the sound event ends up louder than L_out (here
% by 10*log10(10 s/2.5 s) = 6 dB).
%
%  HOW TO RUN THIS CODE: this is a standalone code. Therefore, no additional
%  steps are necessary to run this code.
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

fs = 48000;  % sampling frequency (Hz)
p0 = 20e-6;  % reference sound pressure (Pa)

%% A "downloaded" file of unknown level

rng(1); % fixed seed, so the script is reproducible
x = [zeros(2.5*fs,1); randn(2.5*fs,1); zeros(5*fs,1)];  % 2.5-s noise burst within a 10-s file
x = 0.99 * x / max(abs(x));                             % peak-normalized to 0.99
idx_burst = round(2.5*fs)+1 : 5*fs;                     % samples of the noise burst

%% Level adjustment

L_out = 65;  % assumed level (dB SPL)

p_whole = calibrate(x, x, L_out);             % (i) reference: whole file, including silence
p_burst = calibrate(x, x(idx_burst), L_out);  % (ii) reference: noise burst only

fprintf('Level of the noise burst:\n');
fprintf('  %-48s %.1f dB SPL\n', 'read with the default dBFS = 94 dB (arbitrary):', 20*log10(rms(x(idx_burst))/p0));
fprintf('  %-48s %.1f dB SPL\n', 'adjusted, reference = whole file:', 20*log10(rms(p_whole(idx_burst))/p0));
fprintf('  %-48s %.1f dB SPL (assumed: %g dB SPL)\n', 'adjusted, reference = noise burst:', ...
        20*log10(rms(p_burst(idx_burst))/p0), L_out);

%% SPL over time (Z-weighting, fast time-weighting)

L_x     = Do_SLM(x, fs, 'Z', 'f', 94);       % 94 dB: input taken as Pa
L_whole = Do_SLM(p_whole, fs, 'Z', 'f', 94);
L_burst = Do_SLM(p_burst, fs, 'Z', 'f', 94);
t = (0:length(x)-1)'/fs;

%% Plot

figure('Name', 'Level adjustment');

subplot(2,1,1)
plot(t, x, 'Color', [0.6 0.6 0.6]); hold on;
plot(t(idx_burst), x(idx_burst), 'b');
xlabel('Time, $t$ (s)', 'Interpreter', 'Latex');
ylabel('Value stored, $x$', 'Interpreter', 'Latex');
legend('whole file', 'noise burst', 'Location', 'northeast');
title('(a) Downloaded file (no unit, level unknown)', 'Interpreter', 'Latex');

subplot(2,1,2)
plot(t, L_x, 'Color', [0.6 0.6 0.6]); hold on;
plot(t, L_whole, 'r');
plot(t, L_burst, 'b');
yline(L_out, 'k--');
xlabel('Time, $t$ (s)', 'Interpreter', 'Latex');
ylabel('$L_\mathrm{Z,F}$ (dB re 20 $\mu$Pa)', 'Interpreter', 'Latex');
legend('dBFS = 94 dB (arbitrary)', 'reference: whole file', 'reference: noise burst', ...
       sprintf('$L_\\mathrm{out}$ = %g dB SPL', L_out), 'Location', 'northeast', 'Interpreter', 'Latex');
ylim([0 100]);
title('(b) SPL after level adjustment', 'Interpreter', 'Latex');

set(gcf, 'color', 'w');

if save_figs == 1
    figs_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end
    saveas(gcf, [figs_dir 'ex4_level_adjustment'], 'png');
end
