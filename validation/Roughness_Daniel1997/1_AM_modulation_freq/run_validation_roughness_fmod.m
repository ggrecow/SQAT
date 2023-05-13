% Script run_validation_roughness_fmod
%
%  This routine plot the validation of DW roughness code for 100% AM tones
%  at different frequencies in dependence of the modulation frequency
%   
%   - inputs loaded: The reference data is is taken from jury-test from the 
%       paper (FIG3): Daniel, Peter, and Reinhard Weber. "Psychoacoustical 
%      roughness: Implementation of an optimized model." Acta Acustica united 
%       with Acustica 83 (1997): 113-123.
%   
%   - load signals - AM tones with different carrier frequencies (md=1 and L=60dB) 
%     as a function of fmod [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160]
%     (refer to the function in the end of this script to see how the
%     signals were originally generated)
%
% Author: Gil Greco, Braunschweig, 21/02/2020 (updated 13.05.2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs=0; %% save figs flag

%% path settings 

dir_analysis_name = '1_AM_modulation_freq';
dir_out = [fileparts(mfilename('fullpath')) filesep];
 
% Figure where the figures (and the results) will be stored:
figures_dir = [dir_out 'figs' filesep];
if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

fname_res = 'Rmean_results.mat'; % name for the result file
fname_res_full = [figures_dir fname_res];

bCalculation = ~exist(fname_res_full,'file');
if bCalculation == 0
    fprintf('%s.m: Results file found on disk, confirm that you want to load those prestored\n');
    fprintf('\t results or whether you want to re-run the calculations.\n');
    bCalculation = input('Enter your choice (1=re-run; 0=read stored results): ');
end
bLoad = ~bCalculation;

%% load  reference data

dir_ref_values = get_dir_reference_values('Roughness_Daniel1997',dir_analysis_name); 

fmod_125hz = [];
load([dir_ref_values 'fmod_125hz.mat']);

fmod_250hz = [];
load([dir_ref_values 'fmod_250hz.mat']);

fmod_500hz = [];
load([dir_ref_values 'fmod_500hz.mat']);

fmod_1khz = [];
load([dir_ref_values 'fmod_1khz.mat']);

fmod_2khz = [];
load([dir_ref_values 'fmod_2khz.mat']);

fmod_4khz = [];
load([dir_ref_values 'fmod_4khz.mat']);

fmod_8khz = [];
load([dir_ref_values 'fmod_8khz.mat']);

%% load signals to compute roughness (.mat variables in time [s] x sound pressure [Pa] (all with fs=48 kHz)

SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('Roughness_Daniel1997',SQAT_version);

% variable s_125hz
s_125hz = load ([dir_sounds 'vary_modulation_freq_fmod_125hz.mat']); 
s_125hz=transpose(cell2mat(struct2cell(s_125hz)));

% variable s_250hz
s_250hz = load ([dir_sounds 'vary_modulation_freq_fmod_250hz.mat']); 
s_250hz=transpose(cell2mat(struct2cell(s_250hz)));

% variable s_500hz
s_500hz = load ([dir_sounds 'vary_modulation_freq_fmod_500hz.mat']); 
s_500hz=transpose(cell2mat(struct2cell(s_500hz)));

% variable s_1khz
s_1khz = load ([dir_sounds 'vary_modulation_freq_fmod_1khz.mat']); 
s_1khz=transpose(cell2mat(struct2cell(s_1khz)));

% variable s_2khz
s_2khz = load ([dir_sounds 'vary_modulation_freq_fmod_2khz.mat']); 
s_2khz=transpose(cell2mat(struct2cell(s_2khz)));

% variable s_4khz
s_4khz = load ([dir_sounds 'vary_modulation_freq_fmod_4khz.mat']); 
s_4khz=transpose(cell2mat(struct2cell(s_4khz)));

% variable s_8khz
s_8khz = load ([dir_sounds 'vary_modulation_freq_fmod_8khz.mat']); 
s_8khz=transpose(cell2mat(struct2cell(s_8khz)));

%% compute roughness from signals

fs=48e3; % sampling frequency (Hz)

N_signals = size(s_125hz,2); % all test frequencies have the same number of fmod

% signal with fc=125 hz
res_125hz=cell([N_signals 1]);  % declaring for memory allocation
res_250hz=cell([N_signals 1]);  % declaring for memory allocation
res_500hz=cell([N_signals 1]);  % declaring for memory allocation
res_1khz=cell([N_signals 1]);  % declaring for memory allocation
res_4khz=cell([N_signals 1]);  % declaring for memory allocation
res_8khz=cell([N_signals 1]);  % declaring for memory allocation
     
if bCalculation
    fprintf('%s.m: testing f=125 Hz\n',mfilename);
    tic
    for i=1:N_signals
        res_125hz{i} = Roughness_Daniel1997(s_125hz(:,i),fs,0,false);
    end
    t_125=toc;

    fprintf('%s.m: testing f=250 Hz\n',mfilename);
    % signal with fc=250 hz
 
    tic
    for i=1:N_signals
        res_250hz{i} = Roughness_Daniel1997(s_250hz(:,i),fs,0,false);
    end
    t_250=toc;

    fprintf('%s.m: testing f=500 Hz\n',mfilename);
    % signal with fc=500 hz
   
    tic
    for i=1:N_signals
        res_500hz{i} = Roughness_Daniel1997(s_500hz(:,i),fs,0,false);
    end
    t_500=toc;

    fprintf('%s.m: testing f=1 kHz\n',mfilename);
    % signal with fc=1 kHz
     
    tic
    for i=1:N_signals
        res_1khz{i} = Roughness_Daniel1997(s_1khz(:,i),fs,0,false);
    end
    t_1khz=toc;

    fprintf('%s.m: testing f=2 kHz\n',mfilename);
    % signal with fc=2 khz
    res_2khz=cell([N_signals 1]);  % declaring for memory allocation

    tic
    for i=1:N_signals
        res_2khz{i} = Roughness_Daniel1997(s_2khz(:,i),fs,0,false);
    end
    t_2khz=toc;

    fprintf('%s.m: testing f=4 kHz\n',mfilename);
    % signal with fc=4 khz
     
    tic
    for i=1:N_signals
        res_4khz{i} = Roughness_Daniel1997(s_4khz(:,i),fs,0,false);
    end
    t_4khz=toc;

    fprintf('%s.m: testing f=8 kHz\n',mfilename);
    % signal with fc=8 khz
    
    tic
    for i=1:N_signals
        res_8khz{i} = Roughness_Daniel1997(s_8khz(:,i),fs,0,false);
    end
    t_8khz=toc;

    t_calculation=(t_125+t_250+t_500+t_1khz+t_2khz+t_4khz+t_8khz)./60; % total time to compute roughness

    %% saving results so is not need to run roughness calculation again
    
    save(fname_res_full,'res_125hz','res_250hz','res_500hz','res_1khz', ...
                        'res_2khz','res_4khz','res_8khz','t_calculation');
                    
end
if bLoad
    load(fname_res_full);
end
%% store mean roughness value in vector results[nfmod,nFc] 

% format - results_125hz[1,fmod=0:10:100]
for i=1:N_signals
    results_125hz(i)=res_125hz{i}.Rmean;
    results_250hz(i)=res_250hz{i}.Rmean;
    results_500hz(i)=res_500hz{i}.Rmean;
    results_1khz(i)=res_1khz{i}.Rmean;
    results_2khz(i)=res_2khz{i}.Rmean;
    results_4khz(i)=res_4khz{i}.Rmean;
    results_8khz(i)=res_8khz{i}.Rmean;
end

%% plot 125 Hz and 500 Hz
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves
fmod_125=20:10:100;  % modulation freq vector for reference values 
fmod=20:10:160;  % modulation freq vector for reference values 

err_125hz= (0.17*fmod_125hz);
errorbar(fmod_125,fmod_125hz,err_125hz,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod_125,fmod_125hz,'k-*','MarkerSize',6);hold all;

err_500hz= (0.17*fmod_500hz);
errorbar(fmod,fmod_500hz,err_500hz,'b-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod,fmod_500hz,'b-*','MarkerSize',6);hold all;

% plot computed results
fmod_125_computed=0:10:100;  % modulation freq vector for reference values 
fmod_computed=0:10:160;      % modulation freq vector for reference values 

plot(fmod_125_computed,results_125hz(1:11),'ko:','MarkerSize',8);hold all;
plot(fmod_computed,results_500hz,'bo:','MarkerSize',8);

legend('Ref. - $f_{\mathrm{c}}$=125 Hz$\pm17\:\%\:(\mathrm{JND})$',...
       'Ref. - $f_{\mathrm{c}}$=500 Hz$\pm17\:\%\:(\mathrm{JND})$',...  
       'SQAT - $f_{\mathrm{c}}$=125 Hz',...
       'SQAT - $f_{\mathrm{c}}$=500 Hz',... 
       'Location','NorthEast','Interpreter','Latex');
   
legend boxoff

axis([0 160 0 1.2]);
 ax = gca;
     set(ax,'XTick',[0 20 40 60 80 100 120 140 160]);
     set(ax,'YTick',[0 0.2 0.4 0.6 .8 1]);
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues =  0:10:160; 
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = 0:0.1:1.2; 
     
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

if save_figs==1
    figname_short = 'validation_roughness_fmod_125_500';
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% plot 2 kHz and 8 kHz
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves

fmod=20:10:160;  % modulation freq vector for reference values 

err_2khz= (0.17*fmod_2khz);
errorbar(fmod,fmod_2khz,err_2khz,'b-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod,fmod_2khz,'b-*','MarkerSize',6);hold all;

err_8khz= (0.17*fmod_8khz);
errorbar(fmod,fmod_8khz,err_8khz,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod,fmod_8khz,'k-*','MarkerSize',6);hold all;

% plot computed results
fmod_computed=0:10:160;      % modulation freq vector for reference values 

plot(fmod_computed,results_2khz,'bo:','MarkerSize',8);
plot(fmod_computed,results_8khz,'ko:','MarkerSize',8);hold all;

legend('Ref. - $f_{\mathrm{c}}$=2 kHz$\pm17\:\%\:(\mathrm{JND})$',...
       'Ref. - $f_{\mathrm{c}}$=8 kHz$\pm17\:\%\:(\mathrm{JND})$',...  
       'SQAT - $f_{\mathrm{c}}$=2 kHz',...
       'SQAT - $f_{\mathrm{c}}$=8 kHz',... 
       'Location','NW','Interpreter','Latex');
   
legend boxoff

axis([0 160 0 1.4]);
 ax = gca;
     set(ax,'XTick',[0 20 40 60 80 100 120 140 160]);
     set(ax,'YTick',[0 0.2 0.4 0.6 .8 1 1.2]);
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues =  0:10:160; 
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = 0:0.1:1.4; 
     
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    figname_short = 'validation_roughness_fmod_2k_8k';
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot 1 kHz
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves

fmod=20:10:160;  % modulation freq vector for reference values 

err_1khz= (0.17*fmod_1khz);
errorbar(fmod,fmod_1khz,err_1khz,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod,fmod_1khz,'k-*','MarkerSize',6);hold all;

% plot computed results
fmod_computed=0:10:160;      % modulation freq vector for reference values 

plot(fmod_computed,results_1khz,'ko:','MarkerSize',8);hold all;

legend('Ref. - $f_{\mathrm{c}}$=1 kHz$\pm17\:\%\:(\mathrm{JND})$',...
       'SQAT - $f_{\mathrm{c}}$=1 kHz',...  
       'Location','NW','Interpreter','Latex');
   
legend boxoff

axis([0 160 0 1.6]);
 ax = gca;
     set(ax,'XTick',[0 20 40 60 80 100 120 140 160]);
     set(ax,'YTick',[0 0.2 0.4 0.6 .8 1 1.2 1.4]);
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues =  0:10:160; 
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = 0:0.1:1.6; 
     
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

if save_figs==1
    figname_short = 'validation_roughness_fmod_1k';
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% plot 250 Hz and 4 kHz
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves

fmod_250=20:10:130;  % modulation freq vector for reference values 
fmod=20:10:160;  % modulation freq vector for reference values 

err_250hz= (0.17*fmod_250hz);
errorbar(fmod_250,fmod_250hz,err_250hz,'k-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod_250,fmod_250hz,'k-*','MarkerSize',6);hold all;

err_4khz= (0.17*fmod_4khz);
errorbar(fmod,fmod_4khz,err_4khz,'b-*','MarkerSize',6,'Linewidth',.5);hold all % results from SQAT

% plot(fmod,fmod_4khz,'b-*','MarkerSize',6);hold all;

% plot computed results
fmod_computed=0:10:160;      % modulation freq vector for reference values 

plot(fmod_computed,results_250hz,'ko:','MarkerSize',8);hold all;
plot(fmod_computed,results_4khz,'bo:','MarkerSize',8);

legend('Ref. - $f_{\mathrm{c}}$=250 Hz$\pm17\:\%\:(\mathrm{JND})$',...
       'Ref. - $f_{\mathrm{c}}$=4 kHz$\pm17\:\%\:(\mathrm{JND})$',...  
       'SQAT - $f_{\mathrm{c}}$=250 Hz',...
       'SQAT - $f_{\mathrm{c}}$=4 kHz',... 
       'Location','NW','Interpreter','Latex');
   
legend boxoff

axis([0 160 0 1.2]);
 ax = gca;
     set(ax,'XTick',[0 20 40 60 80 100 120 140 160]);
     set(ax,'YTick',[0 0.2 0.4 0.6 .8 1 1.2]);
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues =  0:10:160; 
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = 0:0.1:1.2; 
     
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

if save_figs==1
    figname_short = 'validation_roughness_fmod_250_4k';
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% function used to generate the signals (only for reference, not used here)

function insig = il_make_AM_fmod(fc)
% Generates signals for roughness algorithm validation
%
% AM tones with carrier frequency of 1 kHz (md=1 and L=60dB) as a function 
%   of modulation frequency. The sounds are generated such that amplitudes 
%   are expressed in Pascals (i.e., amplitude 1 is equal to 94 dB SPL).
%
% Author: Gil Felix Greco, Braunschweig 17.02.2020 (updated 10.03.2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify signal temporal characteristics

L=5;                             % signal duration (seconds)
fs=48000;                        % sampling frequency
dt=1/fs;                         % time step
t = 0:dt:L;                      % Time vector
t = t(:); % Added by AO, to convert the time vector (and the later 'insig')
          % to a column array.
%% make carrier signal 

% fc=1000;                   % carrier center frequency (Hz)

s1=cos(2*pi*fc.*t); % carrier signal 

%% make modulation signal

m=1;    % modulation depth
fm=[0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160]; % modulation frequency (Hz)

for i=1:length(fm)
    
    s2(:,i) = (1-m/2+m/2.*cos(2*pi*fm(i).*t));
    
end

%% modulated signal

s=s1.*s2;

levelOut = 60;  

for i=1:length(fm)
    
    levelIn = 20*log10(rms(s(:,i))/2e-5);
    s(:,i) = s(:,i) * 10^((levelOut-levelIn)/20);

    SPL(i)=20.*log10(rms(s(:,i))/2e-5); % Checking the level
end

insig=s; % output signal

end
