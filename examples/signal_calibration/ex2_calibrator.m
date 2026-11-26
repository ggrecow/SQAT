%   Script ex2_calibrator
%
% How to calibrate a recording using a sound level calibrator and the
% calibrate.m function of SQAT
%
% A measurement chain (microphone, amplifier, analog-to-digital converter)
% converts sound pressure into numbers with an unknown overall gain. If a
% calibrator, which emits a tone with a known level (e.g. 94 dB SPL at
% 1 kHz), is recorded with the same chain and settings, the gain can be
% obtained and removed:
%
%   [p, CalFactor, dBFS] = calibrate(x, x_cal, L_cal)
%
%   x     : recorded signal to be calibrated (numbers)
%   x_cal : recording of the calibrator (numbers)
%   L_cal : level of the calibrator (dB SPL)
%   p     : calibrated signal (Pa)
%   dBFS  : full-scale level of the measurement chain (dB SPL), i.e. the
%           SPL of a recorded signal with an rms value of 1
%
% This script:
% - generates the sound pressure (Pa) of a calibrator (94 dB SPL, 1 kHz)
%   and of a sound to be measured (500-Hz tone burst, 70 dB SPL)
%
% - "records" both with a measurement chain of unknown gain
%
% - calibrates the recording, and compares its SPL with the true SPL
%
% - shows a common mistake: keeping the silence at the start of the
%   calibrator recording, which yields a wrong calibration factor
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

%% True sound pressures (Pa), known here only because we generate them

L_cal = 94;  % level of the calibrator (dB SPL), i.e. 1 Pa rms
p_cal = [zeros(2*fs,1); il_tone(1000, L_cal, 8, fs)];  % 2 s before the calibrator is switched on, then 8 s of tone

L_true = 70; % level of the sound to be measured (dB SPL)
p_true = [zeros(2*fs,1); il_tone(500, L_true, 4, fs); zeros(2*fs,1)];  % 500-Hz tone burst

%% "Record" both signals with the same measurement chain of unknown gain

x_cal = il_measurement_chain(p_cal);
x     = il_measurement_chain(p_true);

%% Calibrate

% correct: use only the part of the calibrator recording where the tone is stable
idx_stable = round(3*fs):length(x_cal);  % from 3 s on
[p, CalFactor, dBFS] = calibrate(x, x_cal(idx_stable), L_cal);

% wrong: use the whole calibrator recording, including the initial silence
p_wrong = calibrate(x, x_cal, L_cal);

idx_burst = round(2.5*fs):round(5.5*fs);  % part of the tone burst used to obtain its level
fprintf('%-46s %.1f dB SPL\n', 'True level of the tone burst:', L_true);
fprintf('%-46s %.1f dB SPL\n', 'Recording, not calibrated (as if it were Pa):', 20*log10(rms(x(idx_burst))/p0));
fprintf('%-46s %.1f dB SPL\n', 'Calibrated, stable part of calibrator:', 20*log10(rms(p(idx_burst))/p0));
fprintf('%-46s %.1f dB SPL\n', 'Calibrated, whole calibrator recording:', 20*log10(rms(p_wrong(idx_burst))/p0));
fprintf('\nFull-scale level of this measurement chain: dBFS = %.1f dB SPL\n', dBFS);
fprintf('Calibration factor: CalFactor = %.4f Pa per unit\n', CalFactor);

%% SPL over time (Z-weighting, fast time-weighting)

L_x = Do_SLM(x, fs, 'Z', 'f', 94); % 94 dB: input taken as Pa (here wrong: x is not calibrated)
L_p = Do_SLM(p, fs, 'Z', 'f', 94); % 94 dB: input in Pa
t = (0:length(x)-1)'/fs;
t_cal = (0:length(x_cal)-1)'/fs;

%% Plot

figure('Name', 'Calibration with a sound level calibrator');

subplot(2,1,1)
plot(t_cal, x_cal, 'Color', [0.6 0.6 0.6]); hold on;
plot(t_cal(idx_stable), x_cal(idx_stable), 'b');
xlabel('Time, $t$ (s)', 'Interpreter', 'Latex');
ylabel('Recorded calibrator, $x_\mathrm{cal}$', 'Interpreter', 'Latex');
legend('whole recording', 'part used for calibration', 'Location', 'southeast');
title('(a) Recording of the calibrator (no unit)', 'Interpreter', 'Latex');

subplot(2,1,2)
plot(t, L_x, 'r'); hold on;
plot(t, L_p, 'b');
yline(L_true, 'k--');
xlabel('Time, $t$ (s)', 'Interpreter', 'Latex');
ylabel('$L_\mathrm{Z,F}$ (dB re 20 $\mu$Pa)', 'Interpreter', 'Latex');
legend('not calibrated', 'calibrated', 'true level', 'Location', 'northeast');
ylim([0 L_true+20]);
title('(b) SPL of the tone burst', 'Interpreter', 'Latex');

set(gcf, 'color', 'w');

if save_figs == 1
    figs_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end
    saveas(gcf, [figs_dir 'ex2_calibrator'], 'png');
end

%% Local functions

function x = il_measurement_chain(p)
% Microphone + amplifier + analog-to-digital converter, seen as a black box:
% it converts sound pressure (Pa) into numbers using a gain that, in
% practice, is unknown. Pretend you do not know the value below.
gain = 0.0371;
x = gain * p;
end

function p = il_tone(f, L, dur, fs)
% Sinusoidal tone with frequency f (Hz), rms level L (dB SPL) and duration dur (s)
t = (0:1/fs:dur-1/fs)';
p = 20e-6 * 10^(L/20) * sqrt(2) * sin(2*pi*f*t);
end
