clc; clear all; close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Example: compute Di's modified psychoacoustic annoyance of signal 14 from the ISO 532-1:2017 
%
% FUNCTION:
%   OUT = PsychoacousticAnnoyance_Di2016(insig,fs,LoudnessField,time_skip,showPA,show)
%   type <help PsychoacousticAnnoyance_Di2016> for more info
%
% FUNCTION:
%   OUT = PsychoacousticAnnoyance_Di2016_from_percentile(N,S,R,FS,K)
%   type <help PsychoacousticAnnoyance_Di2016_from_percentile> for more info
%
% Gil Felix Greco, Braunschweig 05.04.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load .wav RefSignal 

% path='SQAT_open_source\sound_files\validation\Loudness_ISO532_1\';   % path of the sound file for reference
[CalSignal,fs]=audioread('calibration signal sine 1kHz 60dB.wav');

% path='SQAT_open_source\sound_files\validation\Loudness_ISO532_1\technical_signals_time varying_loudness';  % path of the sound file for reference
[RefSignal,fs]=audioread('Test signal 14 (propeller-driven airplane).wav');

ycal=calibrate(RefSignal,CalSignal,60); % calibrated signal

%% compute psychoacoustic annoyance (from time-varying input signal)

res = PsychoacousticAnnoyance_Di2016(ycal, fs,... % input signal and sampling freq.
                                            0,... % field for loudness calculation; free field = 0; diffuse field = 1;
                                          0.2,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                            1,... % show results of PA, 'false' (disable, default value) or 'true' (enable)                                                                    
                                            1);   % show results of loudness, sharpness, roughness and fluctuation strength, 'false' (disable, default value) or 'true' (enable)
                                        
%% compute psychoacoustic annoyance (from input percentile values)

PA = PsychoacousticAnnoyance_Di2016_from_percentile(res.L.N5,...   % loudness percentile
                                                    res.S.S5,...   % sharpness percentile
                                                    res.R.R5,...   % roughness percentile
                                                    res.FS.FS5,... % fluctuation strength percentile
                                                    res.K.K5);      % tonality percentile 

                                                                                                                                    