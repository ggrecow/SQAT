% Script ex_Sharpness_DIN45692_from_loudness
%
% Example: compute stationary and time-varying sharpness (DIN 45692) from specific loudness inputs
%
% FUNCTION:
%   OUT = Sharpness_DIN45692_from_loudness(SpecificLoudness, weight_type, time, time_skip, show_sharpness)
%   type <help Sharpness_DIN45692_from_loudness> for more info
%
% test signal: narrow band noise with a center frequency of 1kHz,
%              a bandwidth of 160 Hz (920 Hz to 1080 Hz) and an overall level of 60 dB
%              should yield a sharpness value of 1 acum
%
% Author: Gil Felix Greco, Braunschweig 28.02.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% load test signal

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[TestSignal,fs]=audioread([dir_ref_sounds 'RefSignal_Sharpness_DIN45692.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  
lvl_cal_signal = 60; % information from the file name...

dBFS_in = lvl_cal_signal-20*log10(rms(TestSignal)); % difference between target and actual full-scale value 
dBFS_out = 94; % dB full scale convention in SQAT: amplitude of 1 = 1 Pa, or 94 dB SPL
dB_correction = dBFS_in - dBFS_out;
insig_cal = TestSignal * 10^(dB_correction/20);

%% loudness (stationary) calculation using SQAT

L_stationary = Loudness_ISO532_1( insig_cal, fs,...   % input signal and sampling freq.
                                               0,...   % field; free field = 0; diffuse field = 1;
                                               1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                             0.5,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                               1);     % show results, 'false' (disable, default value) or 'true' (enable)                    
                                                           
%% sharpness (stationary) calculation using SQAT

S_stationary = Sharpness_DIN45692_from_loudness(L_stationary.SpecificLoudness,  ... % input (stationary) specific loudness
                                                                   'DIN45692'); ... % Weight_Type, type of weighting function used for sharpness calculation
                                                                           
fprintf('\nCalculation of reference signal (60 dBSPL at 1kHz) based on a stationary loudness\n');
fprintf('\tcalculation yields a sharpness value of %g (acum) using the DIN45692 weighting function\n',S_stationary.Sharpness);

%% loudness (time-varying) calculation using SQAT

L_time_varying = Loudness_ISO532_1( insig_cal, fs,...   % input signal and sampling freq.
                                                 0,...   % field; free field = 0; diffuse field = 1;
                                                 2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                               0.5,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                                 1);     % show results, 'false' (disable, default value) or 'true' (enable)                                   
                                                           
%% sharpness (time-varying) calculation using SQAT

S_time_varying = Sharpness_DIN45692_from_loudness(L_time_varying.InstantaneousSpecificLoudness,...  % input (time-varying) specific loudness
                                                                                    'DIN45692',...  % Weight_Type, type of weighting function used for sharpness calculation
                                                                            L_time_varying.time,... % time vector of the loudness calculation
                                                                                            0.5,... % time_skip (second) for statistics calculation
                                                                                             0);    % show sharpness results; true or false

fprintf('\n Calculation of reference signal (60 dBSPL at 1kHz) based on a time-varying loudness\n');
fprintf('\t calculation yields a mean sharpness value over time of %g (acum) using the DIN45692\n',S_time_varying.Smean);
fprintf('\t weighting function.\n');

%% sharpness (time-varying) calculation using SQAT - compare results obtained using different weighting functions

S_time_varying_bismarck = Sharpness_DIN45692_from_loudness(L_time_varying.InstantaneousSpecificLoudness,... % input (time-varying) specific loudness
                                                  'bismarck',...          % Weight_Type, type of weighting function used for sharpness calculation
                                                  L_time_varying.time,... % time vector of the loudness calculation
                                                  0.5,...                 % time_skip (second) for statistics calculation
                                                  0); 
                                              
S_time_varying_aures = Sharpness_DIN45692_from_loudness(L_time_varying.InstantaneousSpecificLoudness,... % input (time-varying) specific loudness
                                                  'aures',...             % Weight_Type, type of weighting function used for sharpness calculation
                                                  L_time_varying.time,... % time vector of the loudness calculation
                                                  0.5,...                 % time_skip (second) for statistics calculation
                                                  0); 
                                              
%% plot sharpness (time-varying) calculation using SQAT - compare results obtained using different weighting functions

figure('NAME','Sharpness analysis (comparison of weighting functions)');  

plot(S_time_varying.time,S_time_varying.InstantaneousSharpness,'k','Linewidth',2); hold on;
plot(S_time_varying_bismarck.time,S_time_varying_bismarck.InstantaneousSharpness,'c--','Linewidth',1.5); 
plot(S_time_varying_aures.time,S_time_varying_aures.InstantaneousSharpness,'r','Linewidth',1.5);  

xlabel('Time, $t$ (s)','Interpreter','Latex'); 
ylabel('Sharpness, $S$ (acum)','Interpreter','Latex'); 
ylim([0.9 1.1]);

legend('DIN45692:2009','von Bismarck','Aures','Location','best','Interpreter','Latex');
legend boxoff

set(gcf,'color','w')
                                           