% Script ex_Loudness_ECMA418-2
%
% Example: compute Loudness (ECMA 418-2:2024) of reference (mono) signal 
% and exemplary stereo signal.
%
% Reference signal: Reference signal: pure tone with center frequency of 1 kHz and RMS value of 40 dBSPL equals 1 sone_HMS
% - The signal is stored in the following folder: <sound_files\reference_signals>. 
% - Signal label: <RefSignal_Loudness_ECMA418_2.wav>
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
% FUNCTION:
%   OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Loudness_ECMA418_2> for more info
%
% Author: Gil Felix Greco, Braunschweig 10.02.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load reference signal (mono .wav file)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
mono_signal_label = 'RefSignal_Loudness_ECMA418_2.wav';

% load mono signal [Nx1]
[ref_signal.signal, ref_signal.fs]=audioread([dir_ref_sounds mono_signal_label]); 

%% Compute Loudness (mono signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3;% time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_mono = Loudness_ECMA418_2(ref_signal.signal, ref_signal.fs, fieldtype, time_skip, show);
                              
fprintf('\nLoudness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t- Reference signal (1-kHz pure tone with 40 dBSPL)\n');
fprintf('\t- Overall loudness value: %g (sone).\n', OUT_mono.loudnessPowAvg);

%% Load stereo signal

stereo_signal_label = 'ExStereo_TrainStation7-0100-0130.wav';

% load stereo signal [Nx2]
[stereo_signal.signal, stereo_signal.fs]=audioread([dir_ref_sounds stereo_signal_label]); 

%% Compute loudness (stereo signal) 

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3;% time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_stereo = Loudness_ECMA418_2(stereo_signal.signal, stereo_signal.fs, fieldtype, time_skip, show);
         
fprintf('\nLoudness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t- Stereo signal: %s\n', stereo_signal_label );
fprintf('\t- Overall loudness value (channel 1):  %g (sone).\n', OUT_stereo.loudnessPowAvg(1) );
fprintf('\t- Overall loudness value (channel 2):  %g (sone).\n' ,OUT_stereo.loudnessPowAvg(2) );
fprintf('\t- Overall loudness value (combined binaural):  %g (sone).\n' ,OUT_stereo.loudnessPowAvgBin );
                                           