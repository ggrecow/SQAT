% Script ex_Loudness_ISO532_1
%
% Example: compute loudness (ISO 532-1) of stationary and time-varying inputs
%
% FUNCTION:
%   OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
%   type <help Loudness_ISO532_1> for more info
%
% test signal: pure tone with a center frequency of 1kHz,
%              an overall level of 40 dB should yield a loudness value of 1 sone
%
% Author: Gil Felix Greco, Braunschweig 28.02.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;

%% load .wav RefSignal 

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[RefSignal,fs]=audioread([dir_ref_sounds 'RefSignal_Loudness_ISO532_1.wav']);

%% loudness (stationary) calculation 

L_stationary = Loudness_ISO532_1( RefSignal, fs,...   % input signal and sampling freq.
                                              0,...   % field; free field = 0; diffuse field = 1;
                                              1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                              0.5,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                              1);     % show results, 'false' (disable, default value) or 'true' (enable)

%% loudness (time-varying) calculation 

L_time_varying = Loudness_ISO532_1( RefSignal, fs,...  % input signal and sampling freq.
                                               0,...   % field; free field = 0; diffuse field = 1;
                                               2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                               0.5,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                               1);     % show results, 'false' (disable, default value) or 'true' (enable)
      
                                           