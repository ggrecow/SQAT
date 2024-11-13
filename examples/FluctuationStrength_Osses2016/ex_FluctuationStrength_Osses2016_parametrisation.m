% Script ex_FluctuationStrength_Osses2016_parametrisation
%
% Example for FluctuationStrength_Osses2016. It computes the fluctuation 
%   strength of the reference signal, using (1) default values, (2) activating
%   the a0 transfer function as defined in the book by Fastl and Zwicker 2007.
%   The signal to be processed is the same as in ex_FluctuationStrength_Osses2016.m:
%   Reference signal: 60 dBSPL 1 kHz tone 100% modulated at 4 Hz. This 
%   signal should yield 1 vacil.
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
structIn.calculate_a0 = 'calculate_a0';
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

legend('Default','Alternative configuration','Location','SouthEast');
