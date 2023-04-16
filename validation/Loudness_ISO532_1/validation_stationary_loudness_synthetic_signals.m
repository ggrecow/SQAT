clc;clear all; close all;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code computes stationary loudness from the reference signals provided by 
% ISO 532-1:2017 - Annex B.2. (signal 1) and Annex B.3 (signals 2 to 5) 
% using SQAT and plot the comparison against reference values 
%
% Gil Felix Greco, Braunschweig 27.02.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save figs flag

save_figs=0;

%% validation signal 1

signal_1=[-60 -60 78 79 89 72 80 89 75 87 85 79 86 80 71 70 72 71 72 74 69 65 67 77 68 58 45 30]; % 1/3 octave levels provided by ISO 532-1:2017 - Annex B.2.  

[OUT.L{1},OUT.RefScalar{1}]=compute_and_plot(1,...     % insig_num
                                             signal_1,... % insig name str
                                             save_figs,['validation_stationary_loudness_signal_' sprintf('%g',1)]... % savefig inputs
                                                  );
                                              
%% validation signal 2 to 5

signal_str=[ {'Test signal 2 (250 Hz 80 dB).wav'},...
             {'Test signal 3 (1 kHz 60 dB).wav'},...
             {'Test signal 4 (4 kHz 40 dB).wav'},...
             {'Test signal 5 (pinknoise 60 dB).wav'}]; % name of the input signals
 
for i=2:5

[OUT.L{i},OUT.RefScalar{i}]=compute_and_plot(i,...     % insig_num
                                             char(signal_str(1,i-1)),... % insig name str
                                             save_figs,['validation_stationary_loudness_signal_' sprintf('%g',i)]... % savefig inputs
                                             );
end 

%% function (compute loudness and plot comparison

function [OUT,table]=compute_and_plot(insig_num,insig,save_figs,tag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this function computes the loudness using SQAT and plot the comparison
% against the reference values from the ISO 532-1:2017 - Annex B.3. 
%
% INPUTS:
%   insig_num : scalar
%       number of the reference signal to be tested
%
%   insig : string
%       name of the reference signals
%
%   save_figs : scalar
%       1 to save; <else> dont save figures 
%
%   tag : string
%       tag with the name of the figures to be saved
%
% OUTPUTS:
%   OUT : struct
%       contain all outputs from the computed loudness
%
%   table : matrix containing scalar values of N and LN
%           1st col=reference
%           2nd col=computed by SQAT
%           3rd row=relative percentage difference (SQAT minus ref.)
%
% Gil Felix Greco, Braunschweig 27.02.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signals from ISO 532-2:2017

if insig_num~=1 % calculate stationary loudness using SQAT for signals 2 to 5, method == 1
    
    % calibration signal provided in the Annex C of the ISO 532-1:2017

    % path='SQAT_open_source\sound_files\validation\loudness_ISO532_1\';   % path of the sound file for reference
    [RefSignal,~]=audioread('calibration signal sine 1kHz 60dB.wav');

    % Test signal provided in the Annex B.3 of the ISO 532-1:2017
    
    % path='SQAT_open_source\sound_files\validation\loudness_ISO532_1\synthetic_signals_stationary_loudness\';   % path of the sound file for reference
    [signal,fs]=audioread(insig);

    % calibrated .wav signal
    [ycal]=calibrate(signal,RefSignal,60); 

    % Stationary loudness calculation from input audio signal using SQAT
    
    OUT = Loudness_ISO532_1( ycal, fs,...   % input signal and sampling freq.
                                    0,...   % field; free field = 0; diffuse field = 1;
                                    1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                    0,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                    0);     % show results, 'false' (disable, default value) or 'true' (enable)

else % calculate stationary loudness for signal 1 (input from 1/3 octave band) , method == 0
    
    OUT = Loudness_ISO532_1( insig, 1,...   % input signal and sampling freq.
                                    0,...   % field; free field = 0; diffuse field = 1;
                                    0,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                    0,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                    0);     % show results, 'false' (disable, default value) or 'true' (enable)

end

%% calculate difference from reference values given by ISO 532-1:2017

% reference values provided by ISO 532-1:2017 for signals 1 to 5 
reference_loudness=[83.296 14.655 4.019 1.549 10.498];   
reference_loudness_level=[103.802 78.733 60.069 46.317 73.920]; 

reference_loudness=reference_loudness(insig_num); % take ref values from the current signal number
reference_loudness_level=reference_loudness_level(insig_num);

% compute relative percentage difference (SQAT minus ref.)
percentage_difference_loudness=( (OUT.Loudness-reference_loudness)/reference_loudness )*100;
percentage_difference_loudness_level=( (OUT.LoudnessLevel-reference_loudness_level)/reference_loudness_level )*100;

% write results in a table format (1st col=reference; 2nd col=computed by SQAT; 3rd row=relative percentage difference (SQAT minus ref.))
table=[reference_loudness,OUT.Loudness,percentage_difference_loudness;
       reference_loudness_level,OUT.LoudnessLevel,percentage_difference_loudness_level ];

%% plot results (total loudness over time)

h = figure('Name',['Loudness - signal ' sprintf('%g',insig_num)]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

load([pwd '\reference_values\' 'reference_values_ISO532_1_2017_signal_' sprintf('%g',insig_num) '.mat']); % load reference vectors

% reference values 

% plot( reference(:,1), reference(:,2),'b','Linewidth',0.5); % ref N'
a=plot( reference(:,1), reference(:,3),'r:','Linewidth',1); hold on; % ref N'_min
plot( reference(:,1), reference(:,4),'r:','Linewidth',1); % ref N'_max

% SQAT values

b=plot( OUT.barkAxis, OUT.SpecificLoudness,'k','Linewidth',1); % calculated specific loudness

% legend([a,b],'ISO 532-1:2017, 5 \% tolerance','SQAT','Location','Best');
% legend box off

ylabel('Specific loudness, $N^{\prime}$ (sone)','Interpreter','Latex');
xlabel('Critical band rate, $z$ (Bark)','Interpreter','Latex'); 

grid off

set(gcf,'color','w');

if save_figs==1
figuresdir = 'figs\'; 
saveas(gcf,strcat(figuresdir, tag), 'fig');
saveas(gcf,strcat(figuresdir, tag), 'pdf');
saveas(gcf,strcat(figuresdir, tag), 'png');
else
end

end
