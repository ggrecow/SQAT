clc;clear all;close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Example: compute stationary and time-varying Sharpness (DIN 45692) from input signal
%
% FUNCTION:
%   OUT = Sharpness_DIN45692(insig, fs, Weight_Type, field, method, time_skip, show_sharpness, show_loudness)
%
% test signal: narrow band noise with a center frequency of 1kHz,
%              a bandwidth of 160 Hz (920 Hz to 1080 Hz) and an overall level of 60 dB
%              should yield a sharpness value of 1 acum
%
% Gil Felix Greco, Braunschweig 28.02.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load test signal

% path=SQAT_open_source\sound_files\reference_signals\sharpness_DIN45692\'; %  path of the sound file for reference
[TestSignal,fs]=audioread('1KHZ60DB.wav');
                                                   
%% sharpness (stationary) calculation using SQAT

S_stationary = Sharpness_DIN45692(TestSignal, fs,...     % input signal and sampling frequency
                                      'DIN45692',...     % Weight_Type, type of weighting function used for sharpness calculation
                                               0,...     % field used for loudness calculation; free field = 0; diffuse field = 1;
                                               1,...     % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;     
                                               0,...     % time_skip (second) for statistics calculation
                                               0,...     % show sharpness results
                                               0);       % show loudness results
% display result
display(sprintf('\nCalculation of reference signal (60 dBSPL at 1kHz) based on a stationary loudness calculation yields\n a sharpness value of %g (acum) using the DIN45692 weighting function.\n',S_stationary.Sharpness));

%% sharpness (time-varying) calculation using SQAT

S_time_varying = Sharpness_DIN45692(TestSignal, fs,...     % input signal and sampling frequency
                                        'DIN45692',...     % Weight_Type, type of weighting function used for sharpness calculation
                                                 0,...     % field used for loudness calculation; free field = 0; diffuse field = 1;
                                                 2,...     % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;     
                                               0.5,...     % time_skip (second) for statistics calculation
                                                 0,...     % show sharpness results
                                                 0);       % show loudness results
                                                                                                        
% display result
display(sprintf('\nCalculation of reference signal (60 dBSPL at 1kHz) based on a time-varying loudness calculation yields\n a mean sharpness value over time of %g (acum) using the DIN45692 weighting function.\n',S_time_varying.Smean));

%% sharpness (time-varying) calculation using SQAT - compare results obtained using different weighting functions

S_time_varying_bismarck = Sharpness_DIN45692(TestSignal, fs,...     % input signal and sampling frequency
                                                 'bismarck',...     % Weight_Type, type of weighting function used for sharpness calculation
                                                          0,...     % field used for loudness calculation; free field = 0; diffuse field = 1;
                                                          2,...     % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;     
                                                        0.5,...     % time_skip (second) for statistics calculation
                                                          0,...     % show sharpness results
                                                         0);       % show loudness results
                                             

S_time_varying_aures = Sharpness_DIN45692(TestSignal, fs,...     % input signal and sampling frequency
                                                 'aures',...     % Weight_Type, type of weighting function used for sharpness calculation
                                                       0,...     % field used for loudness calculation; free field = 0; diffuse field = 1;
                                                       2,...     % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;     
                                                     0.5,...     % time_skip (second) for statistics calculation
                                                       0,...     % show sharpness results
                                                       0);       % show loudness results
                                             

                                              
%% plot sharpness (time-varying) calculation using SQAT - compare results obtained using different weighting functions

figure('NAME','Sharpness analysis (comparison of weighting functions)');  

plot(S_time_varying.time,S_time_varying.InstantaneousSharpness,'k','Linewidth',2); hold on;
plot(S_time_varying_bismarck.time,S_time_varying_bismarck.InstantaneousSharpness,'c--','Linewidth',1.5); 
plot(S_time_varying_aures.time,S_time_varying_aures.InstantaneousSharpness,'r','Linewidth',1.5); 

xlabel('Time, $t$ (s)','Interpreter','Latex'); 
ylabel('Sharpness, $S$ (acum)','Interpreter','Latex'); 
ylim([0.9 1.1]);

legend('DIN45692','von Bismarck','Aures','Location','best','Interpreter','Latex');
legend boxoff

set(gcf,'color','w')

                                           