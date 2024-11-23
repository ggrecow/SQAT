% Script ex_FluctuationStrength_Osses2016_parametrisation
%
%   Example for FluctuationStrength_Osses2016. It computes the fluctuation 
%   strength of the reference signal, using (1) default values (flat tranmission factor), 
%   (2) activating the a0 transfer function (free-field) as defined in the book by Fastl and Zwicker 2007.
%   The signal to be processed is the same as in ex_FluctuationStrength_Osses2016.m:
%   Reference signal: 60 dB SPL, 1 kHz tone that is 100% modulated at 4 Hz.
%   This signal should yield 1 vacil.
%
% Author: Alejandro Osses, 13.11.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load .wav RefSignal 

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[RefSignal,fs]=audioread([dir_ref_sounds 'RefSignal_FluctuationStrength_Osses2016.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

%% Compute fluctuation strength

fprintf('\nFluctuation strength (Osses et al. model): \n');
fprintf('\tcalculation of reference signal (60 dB 1 kHz tone 100 %% modulated at 4 Hz)\n');

method = 1; % method=0, stationary analysis- window size=length(insig); method=1, time_varying analysis - window size=2s
time_skip = 0; % time_skip, in seconds for statistical calculations
bPlot = 0; % boolean variable, internally called 'show'. 1 or true = to plot
OUT = FluctuationStrength_Osses2016(RefSignal,fs,method,time_skip,bPlot);
fprintf('\tUsing default values, the time-averaged fluctuation strength is %g (vacil).\n',OUT.FSmean);

structIn = []; % empty input struct = load all defaults (equivalent to not specifying the struct)
structIn.a0_type = 'Fastl2007'; % a0 curve as defined by Fastl2007, see calculate_a0.m
OUT_a0 = FluctuationStrength_Osses2016(RefSignal,fs,method,time_skip,bPlot,structIn);
fprintf('\tUsing optional parameters, the time-averaged fluctuation strength value of %g (vacil).\n',OUT_a0.FSmean);
               
%% Plot 1 vacil using signal processed with defaults and with the alternative a0 curve:

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% fluctuation strength
plot(OUT.time,OUT.InstantaneousFluctuationStrength,'b-'); hold on; grid on;
plot(OUT_a0.time,OUT_a0.InstantaneousFluctuationStrength,'r--','LineWidth',2);

ylabel('Fluctuation strength, FS (vacil)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 

ylim([0 1.2]);
set(gcf,'color','w')

legend('fluctuationstrength\_osses2016', 'Fastl2007','Location','SouthEast');

%% Plot  the frequency response of the FIR filters generated to represent the a0 transmission factor

K = 2^12; % FIR filter order 

%  Approximate the transfer function of the outer and middle ear by a
%  low-pass filter, where the ear canal resonance is removed from the Fastl's a0 curve. 
%  This filter has a flat response (unit weigthing or 0 dB) for mid and low frequencies.
%  Although not explicitly stated by Osses et al. 2016 (doi: 10.1121/2.0000410), 
%  the simplified a0 transmission curve leads to very similar results during the validation of 
%  their fluctuation strength algorithm. This is the default of the 
% <FluctuationStrength_Osses2016> model,  which is used if no input is given at 
%  all, as defined by the model's author. 

a0_type = 'fluctuationstrength_osses2016'; % Default a0 filter of the <FluctuationStrength_Osses2016> model
[~, freqs_idle, a0_idle] = calculate_a0(fs,K, a0_type);

% Transmission factor for free-field, according to Fig 8.18 (page 226) in 
% Fastl & Zwicker Book, Psychoacoustics: facts and models 3rd edition
% (doi: 10.1007/978-3-540-68888-4)

% a0_type = 'Fastl2007';  % This is the default value of the function <calculate_a0>.
[~, freqs, a0] = calculate_a0(fs,K);  

figure
plot( hz2bark_local( freqs_idle ), mag2db( a0_idle ), 'b-'); hold on;
plot( hz2bark_local( freqs ), mag2db( a0 ), 'r--'); 

ylim([-30 10]);

legend('fluctuationstrength\_osses2016', 'Fastl2007','Location','SouthEast');

xlabel('Critical band rate, $z$ (Bark)','Interpreter','Latex');
ylabel('Magnitude (dB)','Interpreter','Latex'); 

set(gcf,'color','w')
