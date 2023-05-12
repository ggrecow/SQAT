% Script ex_FluctuationStrength_Osses2016
%
% Example for FluctuationStrength_Osses2016. It computes the fluctuation 
%   strength of the reference signal.
%   Reference signal: 60 dBSPL 1 kHz tone 100% modulated at 4 Hz. This 
%   signal should yield 1 vacil.
%
% FUNCTION:
%   OUT = FluctuationStrength_Osses2016(insig,fs,time_skip,show)
%   type <help FluctuationStrength_Osses2016> for more info
%
% See also: FluctuationStrength_Osses2016, run_validation_FS_fmod
% Author: Gil Felix Greco, Braunschweig 13.03.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load .wav RefSignal 

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[RefSignal,fs]=audioread([dir_ref_sounds 'RefSignal_FluctuationStrength_Osses2016.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

%% Compute fluctuation strength

OUT = FluctuationStrength_Osses2016(RefSignal,fs,...  % input signal and sampling freq.
                                                    1,...  % method=0, stationary analysis- window size=length(insig); method=1, time_varying analysis - window size=2s
                                                    0,...  % time_skip, in seconds for statistical calculations
                                                    1);    % show results, 'false' (disable, default value) or 'true' (enable)
                                                                                             
display(sprintf('\nFluctuation strength (Osses et al. model): \ncalculation of reference signal (60 dB 1 kHz tone 100 %% modulated at 4 Hz)\nyields a time-averaged fluctuation strength value of %g (vacil).\n',OUT.FSmean));
               
%% Plot 1 vacil and input signal

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% pressure 
yyaxis left
plot(time_insig,RefSignal); hold on;
ylim([-0.1 0.1]);

ylabel('Acoustic pressure, $p$ (Pa)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 

% fluctuation strength
yyaxis right
plot(OUT.time,OUT.InstantaneousFluctuationStrength);

ylabel('Fluctuation strength, FS (vacil)','Interpreter','Latex');

a=plot(OUT.time,OUT.FSmean*ones(length(OUT.time)),'k--');
legend(a,{sprintf('$\\mathrm{FS}_{\\mathrm{mean}}=%.5g$ (vacil)',OUT.FSmean)},'Location','NorthEast','Interpreter','Latex'); %legend boxoff 

axis([0 10 0 1.2]);

set(gcf,'color','w');
