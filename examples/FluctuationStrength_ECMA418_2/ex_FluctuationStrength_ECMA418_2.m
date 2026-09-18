% Script ex_FluctuationStrength_ECMA418_2
%
% Example: compute Fluctuation strength (ECMA 418-2:2025) of reference (mono) signal
% and exemplary stereo signal.
%
% Reference signal: 60 dB 1 kHz tone 100% modulated at 4 Hz should yield
% 1 vacilHMS. The reference signal of the Osses et al. fluctuation
% strength model is used (same physical calibration signal):
% - The signal is stored in the following folder: <sound_files\reference_signals\>.
% - Signal label: <RefSignal_FluctuationStrength_Osses2016.wav>
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
%   OUT = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help FluctuationStrength_ECMA418_2> for more info
%
% Authors: Sergio Aguirre & Gil Felix Greco, 17.09.2026
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
clc; clear all; close all;

%% Load reference signal (mono .wav file)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
mono_signal_label = 'RefSignal_FluctuationStrength_Osses2016.wav';

% load mono signal [Nx1]
[ref_signal.signal, ref_signal.fs]=audioread([dir_ref_sounds mono_signal_label]);

%% Compute fluctuation strength (mono signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 700e-3;% time_skip, in seconds for statistical calculations (default: 0.7 seconds - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_mono = FluctuationStrength_ECMA418_2(ref_signal.signal, ref_signal.fs, fieldtype, time_skip, show);

fprintf('\nFluctuation strength (ECMA-418-2:2025 - Hearing Model of Sottek): \n');
fprintf('\t- Reference signal (60 dB 1 kHz tone 100 %% modulated at 4 Hz)\n');
fprintf('\t- Overall fluctuation strength value: %g (vacilHMS).\n',OUT_mono.fluctStrength90Pc);

%% Load stereo signal

stereo_signal_label = 'ExStereo_TrainStation7-0100-0130.wav';

% load stereo signal [Nx2]
[stereo_signal.signal, stereo_signal.fs]=audioread([dir_ref_sounds stereo_signal_label]);

%% Compute fluctuation strength (stereo signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 700e-3;% time_skip, in seconds for statistical calculations (default: 0.7 seconds - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT_stereo = FluctuationStrength_ECMA418_2(stereo_signal.signal, stereo_signal.fs, fieldtype, time_skip, show);

fprintf('\nFluctuation strength (ECMA-418-2:2025 - Hearing Model of Sottek): \n');
fprintf('\t- Stereo signal: %s\n', stereo_signal_label );
fprintf('\t- Overall fluctuation strength value (channel 1):  %g (vacilHMS).\n', OUT_stereo.fluctStrength90Pc(1) );
fprintf('\t- Overall fluctuation strength value (channel 2):  %g (vacilHMS).\n' ,OUT_stereo.fluctStrength90Pc(2) );
fprintf('\t- Overall fluctuation strength value (combined binaural):  %g (vacilHMS).\n' ,OUT_stereo.fluctStrength90PcBin );
