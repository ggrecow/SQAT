% Script ex_Roughness_ECMA418_2
%
% Example: compute Roughness (ECMA 418-2:2024) of reference (mono) signal 
% and exemplary stereo signal.
% 
% Reference signal: 60 dB 1 kHz tone 100% modulated at 70 Hz should yield 1 asper. 
% - The signal is stored in the following folder: <sound_files\reference_signals\>. 
% - Signal label: <RefSignal_Roughness_ECMA418_2.wav>
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
%   OUT = Roughness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Roughness_ECMA418_2> for more info
%
% Author: Gil Felix Greco, Braunschweig 10.02.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load reference signal (mono .wav file)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
mono_signal_label = 'RefSignal_Roughness_ECMA418_2.wav';

% load mono signal [Nx1]
[ref_signal.signal, ref_signal.fs]=audioread([dir_ref_sounds mono_signal_label]); 

%% Compute roughness (mono signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 320e-3;% time_skip, in seconds for statistical calculations (default: 0.32 seconds - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_mono = Roughness_ECMA418_2(ref_signal.signal, ref_signal.fs, fieldtype, time_skip, show);
                              
fprintf('\nRoughness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t- Reference signal (60 dB 1 kHz tone 100 %% modulated at 70 Hz)\n');
fprintf('\t- Overall roughness value: %g (asper).\n',OUT_mono.roughness90Pc);

%% Load stereo signal 

stereo_signal_label = 'ExStereo_TrainStation7-0100-0130.wav';

% load stereo signal [Nx2]
[stereo_signal.signal, stereo_signal.fs]=audioread([dir_ref_sounds stereo_signal_label]); 

%% Compute roughness (stereo signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 320e-3;% time_skip, in seconds for statistical calculations (default: 0.32 seconds - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_stereo = Roughness_ECMA418_2(stereo_signal.signal, stereo_signal.fs, fieldtype, time_skip, show);
                              
fprintf('\nRoughness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
fprintf('\t- Stereo signal: %s\n', stereo_signal_label );
fprintf('\t- Overall roughness value (channel 1):  %g (asper).\n', OUT_stereo.roughness90Pc(1) );
fprintf('\t- Overall roughness value (channel 2):  %g (asper).\n' ,OUT_stereo.roughness90Pc(2) );
fprintf('\t- Overall roughness value (combined binaural):  %g (asper).\n' ,OUT_stereo.roughness90PcBin );


