% Script ex_Roughness_ECMA418_2
%
% Example: compute Roughness (ECMA 418-2:2024) of reference signal 
% 
% FUNCTION:
%   OUT = Roughness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Roughness_ECMA418_2> for more info
%
% Reference signal: 60 dB 1 kHz tone 100% modulated at 70 Hz should yield 1 asper.
%
% Author: Gil Felix Greco, Braunschweig 14.01.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load .wav RefSignal (mono .wav file)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

% mono signal [Nx1]
[RefSignal, fs]=audioread([dir_ref_sounds 'RefSignal_Roughness_ECMA418_2.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

% make stereo signal [Nx2]
% RefSignal = [RefSignal,RefSignal];

%% Compute roughness (mono signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 320e-3;% time_skip, in seconds for statistical calculations (default: 0.32 seconds - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT = Roughness_ECMA418_2(RefSignal, fs, fieldtype, time_skip, show);
                              
fprintf('\nRoughness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t calculation of reference signal (60 dB 1 kHz tone 100 %% modulated at 70 Hz)\n');
fprintf('\t yields a overall roughness value of %g (asper).\n',OUT.roughness90Pc);