% Script ex_PsychoacousticAnnoyance_More2010
%
% Example: compute More's modified psychoacoustic annoyance of signal 14 from the ISO 532-1:2017 
%
% FUNCTION:
%   OUT = PsychoacousticAnnoyance_More2010(insig,fs,LoudnessField,time_skip,showPA,show)
%   type <help PsychoacousticAnnoyance_More2010> for more info
%
% FUNCTION:
%   OUT = PsychoacousticAnnoyance_More2010_from_percentile(N,S,R,FS,K)
%   type <help PsychoacousticAnnoyance_More2010_from_percentile> for more info
%
% Author: Gil Felix Greco, Braunschweig 05.04.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% load .wav RefSignal 

dir_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
insig_fname = [dir_sounds 'RefSignal_Loudness_ISO532_1.wav'];

[insig,fs]=audioread(insig_fname); %'sound_files\reference_signals\' - path of the sound file for reference
lvl_cal_signal = 40;

[insig_cal, cal_factor, dBFS_out] = calibrate(insig,insig,lvl_cal_signal); % calibrate signal

lvl_rms = 20*log10(rms(insig_cal))+dBFS_out;

fprintf('\n%s.m: the RMS level of the calibrated input signal is %.1f dB SPL\n',mfilename,lvl_rms);
fprintf('\n\t(full scale value = %.0f dB SPL)\n', dBFS_out);
fprintf('\n\t(file being processed: %s)\n',insig_fname);

%% compute psychoacoustic annoyance (from time-varying input signal)

res = PsychoacousticAnnoyance_More2010(insig_cal, fs,... % input signal and sampling freq.
                                              0,... % field for loudness calculation; free field = 0; diffuse field = 1;
                                            0.2,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                              1,... % show results of PA, 'false' (disable, default value) or 'true' (enable)                                                                    
                                              1);   % show results of loudness, sharpness, roughness and fluctuation strength, 'false' (disable, default value) or 'true' (enable)
PA_from_insig = res.PA5;

%% compute psychoacoustic annoyance (from input percentile values)

PA = PsychoacousticAnnoyance_More2010_from_percentile(res.L.N5,...   % loudness percentile
                                                      res.S.S5,...   % sharpness percentile
                                                      res.R.R5,...   % roughness percentile
                                                    res.FS.FS5,...   % fluctuation strength percentile
                                                     res.K.K5);      % tonality percentile 

fprintf('\n%s.m: Psychoacoustic annoyance (PA) from the input signal=%.1f (arbitrary units)\n',mfilename, PA_from_insig);
fprintf('\n\tPA estimated directly from the percentiles=%.1f (arbitrary units)\n', PA);      
