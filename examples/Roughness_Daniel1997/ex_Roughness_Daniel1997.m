% Script ex_Roughness_Daniel1997
%
% Example: compute Roughness (Daniel & Weber model) of reference signal 
%
% FUNCTION:
%   OUT = Roughness_Daniel1997(insig,fs,time_skip,show) 
%   type <help Roughness_Daniel1997> for more info
%
% Reference signal: 60 dB 1 kHz tone 100% modulated at 70 Hz should yield 1 asper.
%
% Author: Gil Felix Greco, Braunschweig 10.03.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;

%% Load .wav RefSignal 
dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep ...
    'Roughness_Daniel1997' filesep];

[RefSignal,fs]=audioread([dir_ref_sounds 'RefSignal_Roughness_1asper_48kHz_32bit.wav']);

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

%% Compute roughness

OUT=Roughness_Daniel1997(RefSignal,fs,...  % input signal and sampling freq.
                                    0,...  % time_skip, in seconds for statistical calculations
                                    1);    % show results, 'false' (disable, default value) or 'true' (enable)  
                                
fprintf('\nRoughness (Daniel & Weber model): \n');
fprintf('\t calculation of reference signal (60 dB 1 kHz tone 100 %% modulated at 70 Hz)\n');
fprintf('\t yields a time-averaged roughness value of %g (asper).\n',OUT.Rmean);

%% Plot 1 asper and input signal

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% pressure 
yyaxis left
plot(time_insig,RefSignal);
axis([0 .11 -0.1 0.1]);
ylabel('Acoustic pressure, $p$ (Pa)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 

% roughness
yyaxis right
plot(OUT.time,OUT.InstantaneousRoughness);
axis([0 0.1 0 1.1]);
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');

set(gcf,'color','w');
