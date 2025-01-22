% Script ex_Roughness_ECMA418_2
%
% Example: compute Roughness (ECMA 418-2) of reference signal 
% 
% FUNCTION:
%   OUT = Roughness_ECMA418_2(insig, fs, fieldtype, binaural, time_skip, show)
%   type <help Roughness_ECMA418_2> for more info
%
% Reference signal: 60 dB 1 kHz tone 100% modulated at 70 Hz should yield 1 asper.
%
% Author: Gil Felix Greco, Braunschweig 14.01.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Load .wav RefSignal (mono .wav file)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

% mono signal [Nx1]
[RefSignal, fs]=audioread([dir_ref_sounds 'RefSignal_Roughness_ECMA418_2.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference  

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

% make stereo signal [Nx2]
% RefSignal = [RefSignal,RefSignal];

%% Compute roughness (mono signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 0.32;% time_skip, in seconds for statistical calculations (default: 0.32 seconds - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT = Roughness_ECMA418_2(RefSignal, fs, fieldtype, time_skip, show);
                              
fprintf('\nRoughness (ECMA-418-2:2024 - Hearing Model of Sottek): \n');
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
plot(OUT.timeOut,OUT.roughnessTDep);
axis([0 5 0 1.1]);
ylabel('Roughness, $R$ (asper$_{\mathrm{HMS}}$)','Interpreter','Latex');

set(gcf,'color','w');
