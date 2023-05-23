% Script run_validation_roughness_FM_level
%
% Verification of DW roughness code for FM tones dependence on the SPL
%   
% - Inputs, signals: FM tones with carrier frequency of 1.6 kHz, freq deviation of 800 Hz
%   fmod=70 Hz as a function of the SPL, from 40 dBSPL to 80dBSPL in 10 dBSPL steps
%
% - Source: fig 11 from (reference results are the experimental data from this fig
%   Daniel, P. and Weber, R.: Psychoacoustical roughness: implementation of an optimized model. 
%   Acustica – Acta Acustica 81, 1–12 (1995).
%
% Roughness computed using:
%   OUT = Roughness_Daniel1997(insig,fs,time_skip,show) 
%   type <help Roughness_Daniel1997> for more info
% 
% In order to run this code, the user needs to download the dataset of 
%  sound files from zenodo (https://doi.org/10.5281/zenodo.7933206).
%  The obtained folder called `validation_SQAT_v1_0` has to be included in 
%  the `sound_files` folder of the toolbox. 
%
% Author: Gil Felix Greco, Braunschweig 02.03.2020 (updated 13.05.2023) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs=0; %% save figs flag

%% path settings 

dir_analysis_name = '4_FM_level';
dir_out = [fileparts(mfilename('fullpath')) filesep];

%% load  reference data

dir_ref_values = get_dir_reference_values('Roughness_Daniel1997',dir_analysis_name); 

load([dir_ref_values 'REF_FM_level.mat']);

%% load  signals to compute roughness (.mat variables in time [s] x sound pressure [Pa] (all with fs=48 kHz)

SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('Roughness_Daniel1997',SQAT_version);

load([dir_sounds 'FM_fmod_fc_1600hz_reference_tone.mat']); % load ref signal 
load([dir_sounds 'FM_fmod70_fc1600hz_fdev800_SPL40-80.mat']);  % load signals with varying level
 
%% compute FS from signals

fs=48e3;

tic
for i=1:size(s,1)
    res{i} = Roughness_Daniel1997(s(i,:)',fs,0,false);
end

for i=1:size(s_ref,1)
    res_ref{i} = Roughness_Daniel1997(s_ref(i,:)',fs,0,false);
end

%% store mean roughness value in vector results[nfmod,nFc]

for i=1:size(s,1)
    results(i)=res{1,i}.Rmean;
end

for i=1:size(s_ref,1)
    results_ref(i)=res_ref{1,i}.Rmean;
end

%% plot

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves
err= (0.17*REF_FM_level(:,2));
errorbar(REF_FM_level(:,1),REF_FM_level(:,2),err,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(REF_FM_level(:,1),REF_FM_level(:,2),'k*-','MarkerSize',8);hold all;

% plot computed results
plot(REF_FM_level(:,1),results.*100./results_ref,'ko:','MarkerSize',8);hold all;

legend('Reference $\pm\:17\:\%\:(\mathrm{JND})$','SQAT','Location','NW','Interpreter','Latex');
legend boxoff

axis([40 80 0 180]);

ax = gca;
set(ax,'XTick',[40 50 60 70 80]);
set(ax,'YTick',[0 20 40 60 80 100 120 140 160 180]);
ax.XAxis.MinorTick = 'on';
%ax.XAxis.MinorTickValues =  0:10:160;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:0.1:3;
     
ylabel('Relative roughness, $R/R_0$ (\%)','Interpreter','Latex');
xlabel('Sound pressure level, $L_{\mathrm{p}}$ (dB re 20 $\mu$Pa)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'validation_FS_fmod_FM_tones_level';
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function used to generate the FM signals with varying level (only for reference, not used here)

function s=il_make_FM_fmod_fc_1600hz
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates signals for roughness algorithm validation
%
% FM tones with carrier frequency of 1.6 kHz, freq deviation of 800 Hz
% fmod=70 Hz as a function of the SPL
%
%
% Gil Felix Greco, Braunschweig 02.03.2020 (updated 10.03.2023)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify signal temporal characteristics

L=5;                             % signal duration (seconds)
fs=48000;                        % sampling frequency
dt=1/fs;                         % time step
t = 0:dt:L;                      % Time vector

%% parameters 

fc=1600;                   % carrier center frequency (Hz)
Am=1;                      % modulation amplitude
freq_dev=800;              % frequency deviation

%% make modulation signal

fmod=70; % modulation frequency (Hz)

s2=sin(2*pi*fc.*t+(Am*freq_dev/fmod).*sin(2*pi.*fmod.*t));

%% modulated signal

LevelIn=[40 50 60 70 80];

A(1:5)=20e-6*10.^(LevelIn(1:5)./20).*sqrt(2);

for i=1:length(A)
    s(i,:)=A(i)*s2; % output signals
end

for i=1:length(A)
    SPL(i)=20.*log10(rms(s(i,:))/2e-5);
end

%% save signal 

% save s s
% 
% %% plot signal
% 
% h  =figure;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% plot(t,s(1,:)); hold on;              
% 
% xlim([0 1]);ylim([-.1 0.1]);
% ylabel('Amplitude (Pa)','Interpreter','Latex'); xlabel('Time, $t$ (s)','Interpreter','Latex'); 
% 
% % box off
% set(gcf,'color','w');

end
%% function used to generate the reference FM signal (only for reference, not used here)

function s=il_make_FM_fmod_fc_1600hz_reference_tone

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates signals for roughness algorithm validation
%
% FM tones with carrier frequency of 1.6 kHz, freq deviation of 800 Hz and L=60dB 
% as a function of modulation frequency
%
%
% Gil Felix Greco, Braunschweig 02.03.2020 (updated 10.03.2023)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify signal temporal characteristics

L=10;                            % signal duration (seconds)
fs=48000;                        % sampling frequency
dt=1/fs;                         % time step
t = 0:dt:L;                      % Time vector

%% parameters 

fc=1600;                   % carrier center frequency (Hz)
Am=1;                      % modulation amplitude
freq_dev=800;              % frequency deviation

%% make modulation signal

fmod=70; % modulation frequency (Hz)

for i=1:length(fmod)
    s2(i,:) =sin(2*pi*fc.*t+(Am*freq_dev/fmod(i)).*sin(2*pi.*fmod(i).*t));
end

%% modulated signal

LevelIn=60;
A=20e-6*10^(LevelIn/20)*sqrt(2);
A=A.*ones(1,length(fmod));

for i=1:length(fmod)
    s(i,:)=A(i)*s2(i,:); % output signals
end

for i=1:length(fmod)
    SPL(i)=20.*log10(rms(s(i,:))/2e-5);
end

%% save generated signal 

% s_ref=s;
% 
% save s_ref s_ref
% 
% 
% %% plot signal
% 
% h  =figure;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% plot(t,s(1,:)); hold on;                 
% 
% xlim([0 1]);ylim([-.1 0.1]);
% ylabel('Amplitude [Pa]','Interpreter','Latex'); xlabel('Time, $t$ [s]','Interpreter','Latex'); 
% 
% 
% % box off
% set(gcf,'color','w');

end


