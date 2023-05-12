% Script validation_technical_signals_time_varying
%
% This code computes time-varying loudness from the reference signals 
% provided by ISO 532-1:2017 - Annex B.5. using SQAT and plot the 
% comparison against reference values 
%
% Author: Gil Felix Greco, Braunschweig 27.02.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% save figs flag

save_figs=1;

%% input files

signal_str=[ {'Test signal 14 (propeller-driven airplane).wav'},... % free-field
             {'Test signal 15 (vehicle interior 40 kmh).wav'},...   % diffuse field 
             {'Test signal 16 (hairdryer).wav'},...                 % free-field
             {'Test signal 17 (machine gun).wav'},...               % free-field
             {'Test signal 18 (hammer).wav'},...                    % free-field
             {'Test signal 19 (door creak).wav'},...                % free-field
             {'Test signal 20 (shaking coins).wav'},...             % free-field
             {'Test signal 21 (jackhammer).wav'},...                % free-field
             {'Test signal 22 (ratchet wheel (large)).wav'},...     % free-field
             {'Test signal 23 (typewriter).wav'},...                % free-field
             {'Test signal 24 (woodpecker).wav'},...                % free-field
             {'Test signal 25 (full can rattle).wav'}];             % free-field
         
%% validation

for i=14:25 

[OUT.L{i-13},OUT.RefScalar{i-13}]=compute_and_plot(i,...     % insig_num
                                                  char(signal_str(1,i-13)),... % insig name str
                                                  save_figs,['validation_time_varying_loudness_signal_' sprintf('%g',i)]... % savefig inputs
                                                  );
end 


%% function (compute loudness and plot comparison

function [OUT,table]=compute_and_plot(insig_num,fname_insig,save_figs,tag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this function computes the loudness using SQAT and plot the comparison
% against the reference values from the ISO 532-1:2017 - Annex B.5. 
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
%   table : matrix containing scalar values of Nmax and N5 
%           1st col=reference
%           2nd col=computed by SQAT
%           3rd row=relative percentage difference (SQAT minus ref.)
%
% Gil Felix Greco, Braunschweig 27.02.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signals from ISO 532-2:2017

dir_analysis_name = '3_technical_signals_time varying loudness';
dir_out = [fileparts(mfilename('fullpath')) filesep];
  
dir_sounds = get_dir_validation_sounds('Loudness_ISO532_1');
dir_ref_values = get_dir_reference_values('Loudness_ISO532_1',dir_analysis_name);

% calibration signal provided in the Annex C of the ISO 532-1:2017

% path='sound_files\validation\loudness_ISO532_1\';   % path of the sound file for reference
[RefSignal,~]=audioread([dir_sounds 'calibration signal sine 1kHz 60dB.wav']);

% Test signal provided in the Annex B.5 of the ISO 532-1:2017

% path='sound_files\validation\loudness_ISO532_1\';% path of the sound file for reference
[signal,fs]=audioread([dir_sounds fname_insig]);

%% calibrated .wav signal

[ycal]=calibrate(signal,RefSignal,60); 

%% Loudness calculation using SQAT

if insig_num==15 % only signal to be tested considering difuse field
    
OUT = Loudness_ISO532_1( ycal,   fs,...   % input signal and sampling freq.
                                  1,...   % field; free field = 0; diffuse field = 1;
                                  2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                                  0,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                  0);     % show results, 'false' (disable, default value) or 'true' (enable)  
else % free-field
    
OUT = Loudness_ISO532_1( ycal,   fs,...   % input signal and sampling freq.
                                  0,...   % field; free field = 0; diffuse field = 1;
                                  2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                                  0,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                  0);     % show results, 'false' (disable, default value) or 'true' (enable)  
end

%% calculate difference from reference values given by ISO 532-1:2017

% reference values provided by ISO 532-1:2017 for signals 14 to 25
reference_Nmax=[22.640 9.606 38.536 11.211 12.647 10.882 14.880 9.719 8.906 11.186 9.275 7.259];   
reference_N5=[18.081 8.755 36.888 9.722 10.392 10.009 13.093 9.00 8.158 10.421 8.519 5.761]; 

reference_Nmax=reference_Nmax(insig_num-13);
reference_N5=reference_N5(insig_num-13);

% compute relative percentage difference (SQAT minus ref.)
percentage_difference_Nmax=( (OUT.Nmax-reference_Nmax)/reference_Nmax )*100; % max loudness over time
percentage_difference_N5=( (OUT.N5-reference_N5)/reference_N5 )*100; % 5% percentile loudness

table=[reference_Nmax,OUT.Nmax,percentage_difference_Nmax;
       reference_N5,OUT.N5,percentage_difference_N5 ];

%% plot results (total loudness over time)

h = figure('Name',['Loudness - signal ' sprintf('%g',insig_num)]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

reference = []; % to be loaded in the next line...
fname = sprintf('%sreference_values_ISO532_1_2017_signal_%g.mat', dir_ref_values, insig_num);
load(fname); % load reference vectors

% load([pwd '\reference_values\' 'reference_values_ISO532_1_2017_signal_' sprintf('%g',insig_num) '.mat']); % load reference vectors

% reference values

% plot( reference(:,1), reference(:,2),'b','Linewidth',0.5);hold on; % ref N

a=plot( reference(:,1), reference(:,3),'r:','Color',[1 0 0],'Linewidth',1);hold on; % ref N_min (5% tolerance)
plot( reference(:,1), reference(:,4),'r:','Color',[1 0 0],'Linewidth',1); % ref N_max (5% tolerance)

b=plot( reference(:,1), reference(:,5),'Color',[1 0 0],'Linewidth',0.5); % ref N_min (10% tolerance)
plot( reference(:,1), reference(:,6),'Color',[1 0 0],'Linewidth',0.5); % ref N_max (10% tolerance)

% SQAT values
c=plot( OUT.time, OUT.InstantaneousLoudness,'k','Linewidth',1); % calculated specific loudness

leg=legend([a,b,c],'5\% tolerance','10\% tolerance','SQAT','Location','Best');
set(leg,'Orientation','horizontal','Location','northoutside');

ylabel('Loudness, $N$ (sone)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 
grid off

axis([ 0 max(OUT.time) 0 max(reference(:,6)) ])

set(gcf,'color','w');

if save_figs==1
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = tag;
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end
