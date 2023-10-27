% Script ex_EPNL_FAR_Part36
%
% Example: compute EPNL of an exemplary aircraft flight over
%
% FUNCTION:
%   OUT = EPNL_FAR_Part36( insig, fs, method, dt, threshold, show )
%   type <help EPNL_FAR_Part36> for more info
%
% test signal: an exemplary pseudorecording of an A320-like aircraft during take-off is
% provided. The pseudorecording was auralized by the author of this script
% using inputs provided by research partners, which employs a simulation procedure
% based on semi-empirical models to predict the overall aircraft noise at the ground.  
% The receiver is located at 4.5 km from the runway threshold. This file is used as, 
% to the best of the Author's knowledge, a reference sound file for the EPNL metric calculation
% does not exists. As reference, the EPNL calculated for this very same aircraft
% flight over by another tool is 98.97 EPNdB
% 
% Author: Gil Felix Greco, Braunschweig 27.10.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% load .wav signal and calibrate 

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[ExSignal,fs]=audioread([dir_ref_sounds 'ExSignal_A320_auralized_departure_104dBFS.wav']); % input .wav signal - for reference, the path where the sound file is located is:  'sound_files\reference_signals\' - 

dBFS = 104; % a priori knowledge
ExSignal = ExSignal.*10^( (dBFS-94) /20); % calibrate sound file, which was auralized considering a dBFS = 104 to avoid clipping. This means that the +1/-1 range of the .wav file corresponds to 104 dBSPL (or ~3.17 Pa) 

%% EPNL calculation - method == 1 ( an calibrated sound file is used as input)

% input parameters
method = 1;
dt = 0.5;
threshold = 10;
show = 1;

EPNL_method_1 = EPNL_FAR_Part36(ExSignal, fs,... % input signal and sampling freq.
                                                        method,... % method = 0, insig is a SPL[nTime,nFreq] matrix; method = 1, insig is a sound file
                                                               dt,... % time-step in which the third octave SPLs are averaged, in seconds.
                                                     threshold,... % threshold value used to calculate the PNLT decay from PNLTM during the calculation of the duration correction
                                                        show);     % show results, 'false' (disable, default value) or 'true' (enable)

%% EPNL calculation - method == 0 ( an SPL[nTime,nFreq] matrix is used as input )

ExSignal_method0 = EPNL_method_1.SPL_TOB_spectra; % input signal

% input parameters
method = 0;

EPNL_method_0 = EPNL_FAR_Part36(ExSignal_method0, [],... % input signal and sampling freq.
                                                                     method,... % method = 0, insig is a SPL[nTime,nFreq] matrix; method = 1, insig is a sound file
                                                                            dt,...  % time-step in which the third octave SPLs are averaged, in seconds.
                                                                  threshold,...  % threshold value used to calculate the PNLT decay from PNLTM during the calculation of the duration correction
                                                                      show);     % show results, 'false' (disable, default value) or 'true' (enable)
                                                
