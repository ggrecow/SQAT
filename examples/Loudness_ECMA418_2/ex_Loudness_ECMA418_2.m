% Script ex_Loudness_ECMA418-2
%
% Example: compute Loudness (ECMA 418-2:2024) of reference signal 
%
% FUNCTION:
%   OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Loudness_ECMA418_2> for more info
%
%  Reference signal: pure tone with center frequency of 1 kHz and RMS value
%  of 40 dBSPL equals 1 tu_HMS
%
% Author: Gil Felix Greco, Braunschweig  22.01.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load .wav RefSignal (mono .wav file)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

% mono signal [Nx1]
[RefSignal, fs]=audioread([dir_ref_sounds 'RefSignal_Loudness_ECMA418_2.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

% make stereo signal [Nx2]
% RefSignal = [RefSignal,RefSignal];

%% Compute Loudness (mono signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3;% time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT = Loudness_ECMA418_2(RefSignal, fs, fieldtype, time_skip, show);
                              
fprintf('\t Loudness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t calculation of reference signal (1-kHz pure tone with 40 dBSPL)\n');
fprintf('\t yields a overall loudness value of %g (sone).\n',OUT.loudnessPowAvg);
                                           