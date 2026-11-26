%   Script ex1_why_calibrate
%
% Why do we need to calibrate a signal?
%
% A digital signal, e.g. read from a .wav file, is a sequence of numbers
% without units. The same numbers can describe a quiet or a loud sound,
% depending on which sound pressure, in Pa, a value of 1 represents. In
% SQAT, this is given by the full-scale level, dBFS: the sound pressure
% level (SPL) of a signal with an rms value of 1. The signal in Pa is
%
%   p = x * 10^((dBFS-94)/20), so that dBFS = 94 dB means 1 = 1 Pa
%
% This script:
% - generates a 1-kHz tone with an rms value of 0.1 (no unit), as it could
%   be stored in a .wav file
% - converts the same numbers to Pa assuming different dBFS values
% - computes the SPL and the loudness (Loudness_ISO532_1) for each
%   assumption, and plots the results
%
% Take-home message: the SPL follows the assumed dBFS (1 dB per dB), but
% the loudness does not, it roughly doubles every 10 dB. The same file
% yields 1 sone or 65 sone, depending only on the assumption. Without
% calibration, the result of a psychoacoustic metric is meaningless.
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

%% A digital signal: numbers without units

fs = 48000;                        % sampling frequency (Hz)
t = (0:1/fs:1-1/fs)';              % 1 s
x = 0.1*sqrt(2)*sin(2*pi*1000*t);  % 1-kHz tone, rms value of 0.1 (no unit)

%% Convert the same numbers to Pa using different full-scale levels

p0 = 20e-6;       % reference sound pressure (Pa)
dBFS = 60:5:120;  % assumed full-scale levels (dB SPL)

L = zeros(size(dBFS));
N = zeros(size(dBFS));
for i = 1:length(dBFS)

    p = x * 10^((dBFS(i)-94)/20);  % numbers -> Pa

    L(i) = 20*log10(rms(p)/p0);    % SPL (dB re 20 uPa)

    OUT = Loudness_ISO532_1(p, fs, ...  % input signal (Pa) and sampling frequency
                            0, ...      % field: free field
                            1, ...      % method: stationary
                            0.5);       % time_skip (s)
    N(i) = OUT.Loudness;           % loudness (sone)

    fprintf('dBFS = %3.0f dB -> L = %5.1f dB SPL, N = %5.1f sone\n', dBFS(i), L(i), N(i));
end

%% Plot

figure('Name', 'Why calibrate?');

subplot(2,1,1)
plot(dBFS, L, 'ko-'); grid on;
xlabel('Assumed full-scale level, dBFS (dB SPL)', 'Interpreter', 'Latex');
ylabel('SPL, $L$ (dB re 20 $\mu$Pa)', 'Interpreter', 'Latex');
title('(a) Same numbers, different SPL', 'Interpreter', 'Latex');

subplot(2,1,2)
plot(dBFS, N, 'ko-'); grid on;
xlabel('Assumed full-scale level, dBFS (dB SPL)', 'Interpreter', 'Latex');
ylabel('Loudness, $N$ (sone)', 'Interpreter', 'Latex');
title('(b) Same numbers, different loudness', 'Interpreter', 'Latex');

set(gcf, 'color', 'w');

if save_figs == 1
    figs_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end
    saveas(gcf, [figs_dir 'ex1_why_calibrate'], 'png');
end
