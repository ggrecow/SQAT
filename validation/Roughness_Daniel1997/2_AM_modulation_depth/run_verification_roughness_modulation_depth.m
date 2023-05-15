% Script run_verification_roughness_modulation_depth
%  
% Verification of DW roughness code for roughness dependence on the 
%  modulation depth for AM tones
%
% - Inputs, signals: 1 kHz sinosoidal tone modulated at 70 Hz (60 dBSPL) 
%   with varying modulation depth, md (from 0 to 1, in 0.05 increments).
%   Refer to the function in the end of this script to see how the
%   signals were originally generated.
%
% - Ref. values are the power law md^(~1.6), from 
%   Zwicker, E. and Fastl, H. Second ed, Psychoacoustics, Facts and Models, page 258
%
% - An increment of roughness becomes audible for an increment in the degree
%   of modulation of about 10%, which corresponds to an increment of about
%   17% in roughness. Source: Zwicker, E. and Fastl, H. Second ed,
%   Psychoacoustics, Facts and Models, page 260
%
% Roughness computed using:
%   OUT = Roughness_Daniel1997(insig,fs,time_skip,show) 
%   type <help Roughness_Daniel1997> for more info
%
% Author: Gil Felix Greco, Braunschweig 17.02.2020 (updated in 13.05.2023)
% 
% In order to run this code, the user needs to download the dataset of 
%  sound files from zenodo (https://doi.org/10.5281/zenodo.7933206).
%  The obtained folder called `validation_SQAT_v1_0` has to be included in 
%  the `sound_files` folder of the toolbox. 
%
% This code is part of SQAT v1.0, released 14.05.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% save figs flag

save_figs=0;

%% path settings 

SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('Roughness_Daniel1997',SQAT_version);
dir_out = [fileparts(mfilename('fullpath')) filesep];
 
%% load signals (not .wav)

load([dir_sounds 'vary_modulation_depth.mat']);

%% compute roughness using SQAT

res=cell(1,size(s,1));

for i=1:size(s,1)
    res{i} = Roughness_Daniel1997(s(i,:)',fs,0,false);
end

%% plot results

h  =gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% create a vector of results with time-averaged roughness values
for i=1:size(s,1)
    results(i)=res{1,i}.Rmean;
end

md=0:.05:1;  % modulation depth vector 

% ref curve, and jnd of 17% for roughness
a=1;
power_law=a.*md.^(1.6);

err= (0.17*power_law);
errorbar(md,power_law,err,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot results SQAT
plot(md,results,'ko:','MarkerSize',8);hold all % results from SQAT

legend('Ref. - $m_{\mathrm{d}}^{1.6}\pm17\:\%\:(\mathrm{JND})$','SQAT','Location','NW','Interpreter','Latex');
legend boxoff

axis([0 1 0 1.2]);
ax = gca;
set(ax,'XTick',[0 0.2 0.4 0.6 .8 1]);
set(ax,'YTick',[0 0.2 0.4 0.6 .8 1 1.2 1.4 1.6]);
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =  0:0.05:1;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:0.1:25;
     
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
xlabel('Modulation depth, $m_{\mathrm{d}}$','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'verification_roughness_dependence_md';
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function used to generate the signals (only for reference, not used here)

function [s,fs] = make_AM_varying_md

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates amplitude modulated (AM) signals for roughness algorithm validation
%
% Roughness dependence on the modulation depth for AM tonesGenerates amplitude modulated (AM) signals for roughness algorithm validation
%
%   - 1 kHz sinosoidal tone modulated at 70 Hz (60 dBSPL) with varying
%       modulation depth, md (from 0 to 1, in 0.05 increments)%
%
% Gil Felix Greco, Braunschweig 17.02.2020 (updated in 10.03.2023)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_figs=0;  % save fig flag

%% Specify signal temporal characteristics

L=10;                            % Length (seconds)
fs=48000;                        % Sampling frequency
dt=1/fs;                         % Time step
t = 0:dt:L;                      % Time vector

%% Make carrier signal 

fc=1000;                   % carrier center frequency (Hz)

s1=cos(2*pi*fc.*t); % carrier signal 70 dB

%% make modulation signal

m=[0 0.05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95 1];    % modulation depth
fm=70;                                 % modulation frequency (Hz)

for i=1:length(m)
    
    s2(i,:) = (1-m(i)./2+m(i)./2.*cos(2*pi*fm.*t));

end

%% modulated signal

s=s1.*s2;

levelOut = 60;                % desired signal level (dBSPL)

for i=1:length(m)
    
    levelIn = 20*log10(rms(s(i,:))/2e-5);
    s(i,:) = s(i,:) * 10^((levelOut-levelIn)/20);

end

for i=1:length(m)
    
    SPL(i)=20.*log10(rms(s(i,:))/2e-5);

end

end
