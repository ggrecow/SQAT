% Script run_validation_FS_FM_freq_dev
%
% This routine plot the validation of fluctuation strength (Osses et al. model)
%  for FM tones as a function of the frequency deviation
%
% - Inputs, signals: FM tones with carrier frequency of 1.5 kHz (fmod=4 Hz and L=70dB)
%   as a function of frequency deviation fdev=[16 32 100 300 700];  
%   
% - Reference values taken from: 
%   Zwicker, E. and Fastl, H. Second ed, Psychoacoustics, Facts and Models, 
%   ed. M.R. Schroeder. Springer-Verlag, Berlin, 1999. pg. 251
%
% Fluctuation strength computed using:
%   OUT = FluctuationStrength_Osses2016(insig,fs,method,time_skip,show)
%   type <help FluctuationStrength_Osses2016> for more info
%
% In order to run this code, the user needs to download the dataset of 
%  sound files from zenodo (https://doi.org/10.5281/zenodo.7933206).
%  The obtained folder called `validation_SQAT_v1_0` has to be included in 
%  the `sound_files` folder of the toolbox. 
%
% Author: Gil Felix Greco, Braunschweig, 02/03/2020 (updated in 13.03.2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% reference data from Zwicker, E. and Fastl, Psychoacoustics, Facts and Models, pg. 251 - data taken using webplotdigitizer

ref=[
13.943586989337412 0.05163743230551754
16.02959831142325 0.06320556537909061
18.751442438248745 0.08614646380636914
21.746617245248554 0.111202658854761
25.221863174913047 0.14448930163983897
29.25377718053758 0.18337264888586446
33.93200114606875 0.22884035432123762
39.36021507997048 0.2802339821270263
45.65906942672327 0.33788275021269554
52.9677384030226 0.39981135112144206
61.44984873051068 0.46898274603847234
71.29343423203059 0.5437508454164497
82.71601867720712 0.6218111238891013
95.32547763911934 0.7000643204933077
109.1176246001727 0.774823194643034
124.0675703223092 0.8492482382691673
141.06820999613305 0.9258461200977859
160.39716901309913 1.001478296058633
182.3747286864928 1.077231185252952
208.77300164572932 1.1588707137105478
237.38621574650546 1.2384461886314608
268.10193316248115 1.3135835830907534
302.7997491881304 1.3919399971092834
341.99035344538163 1.4711011660176223
383.64715183365263 1.5448082514159571
430.3893989069232 1.621834950734755
486.09360822117134 1.700996119643094
549.0039587372956 1.779352533661624
629.2563298237576 1.8618268718721231
709.7979531055141 1.9376182531246597
780.6788600495487 1.9986030530240373
844.8392838725672 2.051412675389174];

%% load  signals to compute roughness (.mat variables in time [s] x sound pressure [Pa] (all with fs=48 kHz)

SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('FluctuationStrength_Osses2016',SQAT_version);

load ([dir_sounds 'vary_modulation_freq_fmod_1khz.mat']); % variable s_1khz

%% compute FS from signals

fs=48e3;

res=cell(1,size(s_1khz,1));  % declaring for memory allocation
results_1khz=zeros(1,size(s_1khz,1));  % declaring for memory allocation

tic
for i=1:size(s_1khz,1)
    
    res{i} = FluctuationStrength_Osses2016(s_1khz(i,:),fs,...  % input signal and sampling freq.
                                                             0,...  % method, stationary analysis =0 - window size=length(insig), time_varying analysis - window size=2s
                                                             0,...  % time_skip, in seconds for statistical calculations
                                                             0);    % show results, 'false' (disable, default value) or 'true' (enable)
                                                   
    results_1khz(i)=res{1,i}.FSmean;  % store mean roughness value in vector results[nfmod,nFc] 
end
t=toc/60;  % time to compute FS, in minutes


%% plot

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves
hAx=axes;                     % new axes; save handle
N=1;

err= (0.1*ref(1:N:length(ref),2));
errorbar(ref(1:N:length(ref),1),ref(1:N:length(ref),2),err,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% N=1;
% semilogx(ref(1:N:length(ref),1),ref(1:N:length(ref),2),'k*-','MarkerSize',8);hold all

% plot computed results
freq_dev=[16 32 100 300 700]; 

semilogx(freq_dev,results_1khz,'ko:','MarkerSize',8);hold all;

hAx.XScale='log';              % turn to semilogx form

legend('Reference $\pm\:10\:\%\:(\mathrm{JND})$','SQAT','Location','NW','Interpreter','Latex');
legend boxoff

% axis([16 700 0 2.5]);
ax = gca;
set(ax,'XTick',[16 32 100 300 700]);
set(ax,'YTick',[0 .5 1 1.5 2 2.5 ]);
ax.XAxis.MinorTick = 'off';
ax.XAxis.MinorTickValues =  0:10:160;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:0.1:3;
     
ylabel('Fluctuation strength, FS (vacil)','Interpreter','Latex');
xlabel('Frequency deviation, $f_{\mathrm{dev}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'validation_FS_FM_tones_fdev';
    figname_out = [figures_dir figname_short];
    % figures_dir = 'figs\'; % backslash is the fileseparator on Windows only 
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
