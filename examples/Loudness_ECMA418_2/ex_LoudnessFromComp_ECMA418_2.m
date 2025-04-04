% Script ex_LoudnessFromComp_ECMA418_2
%
% Example: compute Loudness (ECMA 418-2:2024) from exemplary stereo signal using 
% pre-calculated specTonalLoudness, specNoiseLoudness as inputs to the 
% <LoudnessFromComp_ECMA418_2.m> function . 
%
% Description: When computing loudness using the <Loudness_ECMA418_2.m> function,
% the <Tonality_ECMA418_2.m> is called inside the <Loudness_ECMA418_2.m>
% function to calculate specTonalLoudness and specNoiseLoudness. Therefore, if
% you are calculating the tonality anyways, it is not necessary to
% calculate those parameters two times. You can use the specTonalLoudness and 
% specNoiseLoudness outputs from <Tonality_ECMA418_2.m> to directly calculate 
% loudness using the <x_LoudnessFromComp_ECMA418_2.m> function.
%
% Stereo signal: binaural audio recording of a 'train station' environment (30 seconds, 2-channel binaural)
% The signal 'TrainStation.7.wav' was extracted from the EigenScape database 
% (https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between
%  01m00s and 01m30s. The EigenScape database, which is described by 
% Green et al (https://doi.org/10.3390/app7111204), is licenced 
% under Creative Commons Attribution 4.0. 
% - The signal is stored in the following folder: <sound_files\reference_signals\>. 
% - Signal label: <ExStereo_TrainStation7-0100-0130.wav>
%
% FUNCTIONS USED:
%
%   OUT = Tonality_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Tonality_ECMA418_2> for more info
%
%   OUT = LoudnessFromComp_ECMA418_2(specTonalLoudness, specNoiseLoudness, time_skip, show)
%   type <help LoudnessFromComp_ECMA418_2> for more info
%
% Author: Gil Felix Greco, Braunschweig 04.04.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load stereo signal

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
stereo_signal_label = 'ExStereo_TrainStation7-0100-0130.wav';

% load stereo signal [Nx2]
[stereo_signal.signal, stereo_signal.fs]=audioread([dir_ref_sounds stereo_signal_label]); 

%% Compute specTonalLoudness, specNoiseLoudness using the <Tonality_ECMA418_2> function

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3;% time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 0; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_tonality = Tonality_ECMA418_2(stereo_signal.signal, stereo_signal.fs, fieldtype, time_skip, show);

%% compute loudness using the <LoudnessFromComp_ECMA418_2.m> function

show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_loudness =  LoudnessFromComp_ECMA418_2(OUT_tonality.specTonalLoudness, OUT_tonality.specNoiseLoudness, time_skip, show);

fprintf('\nLoudness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t- Stereo signal: %s\n', stereo_signal_label );
fprintf('\t- Overall loudness value (channel 1):  %g (sone).\n', OUT_loudness.loudnessPowAvg(1) );
fprintf('\t- Overall loudness value (channel 2):  %g (sone).\n' ,OUT_loudness.loudnessPowAvg(2) );
fprintf('\t- Overall loudness value (combined binaural):  %g (sone).\n' ,OUT_loudness.loudnessPowAvgBin );
                                           