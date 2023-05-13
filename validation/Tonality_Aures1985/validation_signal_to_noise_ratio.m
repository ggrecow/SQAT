% Script validation_signal_to_noise_ratio
%
%  This routine plot the validation of Aures' tonality code for pure tone 
%    signals (85 dBSPL @ 1 kHz) and varying emergence level [0 10 20 30 40 
%    50 60 70 80] dBSPL over a bandpass (white noise) in the critical band 
%    centered around the tone 
%   
%  The reference data is taken from:
%    Aaron Hastings, Kyoung Hoon Lee, Patricia Davies, and Aimée M. Surprenant
%    "Measurement of the attributes of complex tonal components commonly 
%    found in product sound" Noise Control Eng. J. 51, 2003.
%
% Author: Gil Felix Greco, Braunschweig, 22.03.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% save settings

save_figs=0;

%% define input/output paths

SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('Tonality_Aures1985',SQAT_version); 

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% reference data 
% Source: Hastings et al. (2003), fig. 1

ref=[
0  0.015823397764229252
10  0.07533719859035015
20  0.20912727004956122
30  0.40794580496985177
40  0.5923422380983235
50  0.7627954783966964
60  0.8820348816074157
70  0.9558121564573688
80  0.9825276908062728
     ];

%% load signals

tag_str={'0','10','20','30','40','50','60','70','80'};

N_signals = length(tag_str);
for i=1:N_signals
    
    % dir_sounds_local = [dir_sounds 'MakeToneEmergence' filesep];
    fname_insig = sprintf('%s1Bark_tone_prominence_%sdB_fc_1khz_44khz_64bit.wav',dir_sounds,tag_str{i}); 
    [insigs(:,i),fs]=audioread(fname_insig);

end

%% compute tonality
for i=1:N_signals
    
    bPlot_figure = 0; % set to 1 to plot the raw outputs
    T{i} = Tonality_Aures1985(insigs(:,i),fs,0,0,bPlot_figure);
    results(i)= [T{i}.Kmean]; % create vector with time-averaged tonality values
    
end

%% plot

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 
yyaxis left

% reference results
a=plot( ref(:,1),ref(:,2),'k*-','MarkerSize',8);hold on;

% SQAT
b=plot(0:10:80,results,'ko:','MarkerSize',8);

% dummy plot to make x-axis larger
plot(1,0.1);hold on;
plot(5,0.1);hold on;

ylim([0 1.1]);
xlim([0 80]);

ylabel('Time-averaged tonality, $K_{\mathrm{mean}}$ (t.u.)','Interpreter','Latex');
xlabel(sprintf('Signal to noise ratio in the critical band\n centered about the tone (dBSPL)'),'Interpreter','Latex'); 
grid off

legend([a,b],{ 'Literature',...
        'SQAT'},...
        'Location','SE','Interpreter','Latex');

legend boxoff

yyaxis right

plot(0:10:80,results'-ref(:,2),'o-','MarkerSize',4,'MarkerFaceColor',[0.85,0.33,0.10],'HandleVisibility','off');hold on;

ylim([-0.2 0.2]);

ylabel('SQAT minus Reference (t.u.)','Interpreter','Latex');

ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =  0:10:80;

set(ax,'XTick',[0 20 40 60 80]);     
set(ax,'YTick',[0 .2 .4 .6 .8 1]);

yticks([-0.2 -0.1 0 0.1 0.2])

ax.YAxis(1).Color = 'k';

set(gcf,'color','w');

%%

if save_figs==1
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'tonality_validation_SNR_tone_85dBSPL_1khz';
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf, figname_out, 'fig');
%     saveas(gcf, figname_out, 'pdf');
    saveas(gcf, figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
