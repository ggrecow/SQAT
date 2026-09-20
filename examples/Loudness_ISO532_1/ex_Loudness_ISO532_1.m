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

%% load .wav RefSignal 

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[RefSignal,fs]=audioread([dir_ref_sounds 'RefSignal_Loudness_ISO532_1.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  

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
      
                                           