% Script ex_PsychoacousticAnnoyance_Zwicker1999
%
% Example: computes Zwicker & Fastl's psychoacoustic annoyance of signal 14 from the ISO 532-1:2017 
%
% FUNCTION:
%   OUT = PsychoacousticAnnoyance_Zwicker1999(insig,fs,LoudnessField,time_skip,showPA,show)
%   type <help PsychoacousticAnnoyance_Zwicker1999> for more info
%
% FUNCTION:
%   OUT = PsychoacousticAnnoyance_Zwicker_from_percentile(N,S,R,FS)
%   type <help PsychoacousticAnnoyance_Zwicker1999_from_percentile> for more info
%
% Author: Gil Felix Greco, Braunschweig 14.03.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% load .wav RefSignal 

dir_sounds = [basepath_SQAT 'sound_files' filesep 'validation' filesep 'Loudness_ISO532_1' filesep];

% path='SQAT_open_source\sound_files\validation\Loudness_ISO532_1\';   % path of the sound file for reference
[CalSignal,fs]=audioread([dir_sounds 'calibration signal sine 1kHz 60dB.wav']);
lvl_cal_signal = 60;

% path='SQAT_open_source\sound_files\validation\Loudness_ISO532_1\technical_signals_time varying_loudness';  % path of the sound file for reference
insig_fname = [dir_sounds 'Test signal 14 (propeller-driven airplane).wav'];
[RefSignal,fs]=audioread(insig_fname);

[insig_cal, cal_factor, dBFS_out] = calibrate(RefSignal,CalSignal,lvl_cal_signal); % calibrated signal
lvl_rms = 20*log10(rms(insig_cal))+dBFS_out;
fprintf('%s.m: the RMS level of the calibrated input signal is %.1f dB SPL\n',mfilename,lvl_rms);
fprintf('\t(full scale value = %.0f dB SPL)\n', dBFS_out);
fprintf('\t(file being processed: %s)\n',insig_fname);

%% compute psychoacoustic annoyance (from time-varying input signal)

res = PsychoacousticAnnoyance_Zwicker1999(insig_cal, fs,... % input signal and sampling freq.
                                               0,... % field for loudness calculation; free field = 0; diffuse field = 1;
                                             0.2,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                               1,... % show results of PA, 'false' (disable, default value) or 'true' (enable)                                                                    
                                               1);   % show results of loudness, sharpness, roughness and fluctuation strength, 'false' (disable, default value) or 'true' (enable)
PA_from_insig = res.PA5;
                                                               
%% compute psychoacoustic annoyance (from input percentile values)

PA = PsychoacousticAnnoyance_Zwicker1999_from_percentile(res.L.N5,... % loudness percentile
                                                       res.S.S5,... % sharpness percentile
                                                       res.R.R5,... % roughness percentile
                                                       res.FS.FS5); % fluctuation strength percentile

fprintf('%s.m: Psychoacoustic annoyance (PA) from the input signal=%.1f (arbitrary units)\n',mfilename, PA_from_insig);
fprintf('\tPA estimated directly from the percentiles=%.1f (arbitrary units)\n', PA);                                                                                                                                    