function [y, t_total] = create_two_tones( fc_1, fc_2, LevelIn, L_before, L_signal, L_after, fs )
% function [y, t_total] = create_two_tones( fc_1, fc_2, LevelIn, L_before, L_signal, L_after, fs )
%
% Generate signal with two sinusoidal tones and silence before and after
% the tones
%
%   Inputs
%   fc_1 : central frequency of the first tone (Hz)
%   fc_2 : central frequency of the second tone (Hz)
%   LevelIn : rms level (of both tones) (dBSPL)
%   L_before : length of silence before the tones (s)
%   L_signal : length of the tones (s)
%   L_after : length of silence after the tones (s)
%   fs : sampling frequency (Hz)
%
%   Output
%   y : output signal (Pa)
%   t_total : time vector of the output signal
%
% Author: Gil Felix Greco
% Date: 17/11/2024
% Date: 29/11/2025 (Alejandro Osses: Creating a separate function and
%   renaming from il_create_signal to create_two_tones)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify signal characteristics

t_before = 0:1/fs:L_before; % Time vector, before signal
t_signal = 0:1/fs:L_signal; % Time vector, during signal
t_after = 0:1/fs:L_after; % Time vector, after signal

%% Make signal

pref = 20e-6; % reference sound pressure for air

A = pref*10^(LevelIn/20)*sqrt(2); % Amplitude
y_signal = A*sin(2*pi*fc_1*t_signal) + A*sin(2*pi*fc_2*t_signal); % Generate signal

% SPL=20.*log10(rms(y)/2e-5); % Verify dBSPL of the generated signal

% make silence parts
y_before = zeros(1,length(t_before));
y_after = zeros(1,length(t_after));

%% make complete signal
y = [y_before y_signal y_after];
t_total = linspace(0,(L_before+L_signal+L_after),length(y));
