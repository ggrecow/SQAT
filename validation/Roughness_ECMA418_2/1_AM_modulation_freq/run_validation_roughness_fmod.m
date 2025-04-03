% Script run_validation_roughness_fmod
%
%  Verification of the roughness implementation according to ECMA-418-2:2024
%  Verification case: 100% AM tones at different center frequencies in dependence 
% of the modulation frequency
%
% - Inputs, ref. data: the reference data is is taken from jury-test from the 
%   paper (FIG3, open symbols): Daniel, Peter, and Reinhard Weber. "Psychoacoustical 
%   roughness: Implementation of an optimized model." Acta Acustica united 
%   with Acustica 83 (1997): 113-123.
%   
% - Inputs, signals: AM tones with different carrier frequencies (md=1 and L=60dB) 
%   as a function of fmod [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160]
%   (refer to the function in the end of this script to see how the
%   signals were originally generated)
%
% Roughness computed using:
%   OUT = Roughness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Roughness_ECMA418_2> for more info
%
%  In order to run this code, the user needs to download the dataset of 
%  sound files from zenodo (https://doi.org/10.5281/zenodo.7933206).
%  The obtained folder called `validation_SQAT_v1_0` has to be included in 
%  the `sound_files` folder of the toolbox. 
%
% Author: Gil Felix Greco, Braunschweig, 22.01.2025
% Updated results: Gil Felix Greco, Cala Rajada, 21.03.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

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

    fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse');

    tic
    for i=1:N_signals
        res_125hz{i} = Roughness_ECMA418_2(s_125hz(:,i), fs, fieldtype);
    end
    t_125=toc;
 
    fprintf('%s.m: testing f=250 Hz\n',mfilename);
    % signal with fc=250 hz
 
    tic
    for i=1:N_signals
        res_250hz{i} = Roughness_ECMA418_2(s_250hz(:,i), fs, fieldtype);
    end
    t_250=toc;

    fprintf('%s.m: testing f=500 Hz\n',mfilename);
    % signal with fc=500 hz
   
    tic
    for i=1:N_signals
        res_500hz{i} = Roughness_ECMA418_2(s_500hz(:,i), fs, fieldtype);
    end
    t_500=toc;

    fprintf('%s.m: testing f=1 kHz\n',mfilename);
    % signal with fc=1 kHz
     
    tic
    for i=1:N_signals
        res_1khz{i} = Roughness_ECMA418_2(s_1khz(:,i), fs, fieldtype);
    end
    t_1khz=toc;

    fprintf('%s.m: testing f=2 kHz\n',mfilename);
    % signal with fc=2 khz
    res_2khz=cell([N_signals 1]);  % declaring for memory allocation

    tic
    for i=1:N_signals
        res_2khz{i} = Roughness_ECMA418_2(s_2khz(:,i), fs, fieldtype);
    end
    t_2khz=toc;

    fprintf('%s.m: testing f=4 kHz\n',mfilename);
    % signal with fc=4 khz
     
    tic
    for i=1:N_signals
        res_4khz{i} = Roughness_ECMA418_2(s_4khz(:,i), fs, fieldtype);
    end
    t_4khz=toc;

    fprintf('%s.m: testing f=8 kHz\n',mfilename);
    % signal with fc=8 khz
    
    tic
    for i=1:N_signals
        res_8khz{i} = Roughness_ECMA418_2(s_8khz(:,i), fs, fieldtype);
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
    results_125hz(i)=res_125hz{i}.R90;
    results_250hz(i)=res_250hz{i}.R90;
    results_500hz(i)=res_500hz{i}.R90;
    results_1khz(i)=res_1khz{i}.R90;
    results_2khz(i)=res_2khz{i}.R90;
    results_4khz(i)=res_4khz{i}.R90;
    results_8khz(i)=res_8khz{i}.R90;
end

%% plot - fc=125 Hz

title_tag = '$f_{\mathrm{c}}=125$~Hz';
save_tag = 'validation_roughness_fmod_125hz';

xRef = 20:10:100;  
yRef = fmod_125hz; 

xSQAT = 0:10:100; 
ySQAT = results_125hz(1:11);

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot - fc=250 Hz

title_tag = '$f_{\mathrm{c}}=250$~Hz';
save_tag = 'validation_roughness_fmod_250hz';

xRef = 20:10:130;  
yRef = fmod_250hz; 

xSQAT = 0:10:160;
ySQAT = results_250hz;

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot - fc=500 Hz

title_tag = '$f_{\mathrm{c}}=500$~Hz';
save_tag = 'validation_roughness_fmod_500hz';

xRef = 20:10:160;  
yRef = fmod_500hz; 

xSQAT = 0:10:160;
ySQAT = results_500hz;

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot - fc=1 kHz

title_tag = '$f_{\mathrm{c}}=1$~kHz';
save_tag = 'validation_roughness_fmod_1khz';

xRef = 20:10:160;  
yRef = fmod_1khz; 

xSQAT = 0:10:160;
ySQAT = results_1khz;

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot - fc=2 kHz

title_tag = '$f_{\mathrm{c}}=2$~kHz';
save_tag = 'validation_roughness_fmod_2khz';

xRef = 20:10:160;  
yRef = fmod_2khz; 

xSQAT = 0:10:160;
ySQAT = results_2khz;

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot - fc=4 kHz

title_tag = '$f_{\mathrm{c}}=4$~kHz';
save_tag = 'validation_roughness_fmod_4khz';

xRef = 20:10:160;  
yRef = fmod_4khz; 

xSQAT = 0:10:160;
ySQAT = results_4khz;

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot - fc=8 kHz

title_tag = '$f_{\mathrm{c}}=8$~kHz';
save_tag = 'validation_roughness_fmod_8khz';

xRef = 20:10:160;  
yRef = fmod_8khz; 

xSQAT = 0:10:160;
ySQAT = results_8khz;

il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir);

%% plot function 

function il_plot(xRef, yRef, xSQAT, ySQAT, title_tag, save_figs, save_tag, figures_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate plots of  the roughness verification results 
%
%   INPUT:
%   xRef - ref. data (x-axis) 
%   xRef - ref. data (x-axis) 
%   xSQAT - results of <Roughness_ECMA418_2.m> (x-axis) 
%   ySQAT - results of <Roughness_ECMA418_2.m> (y-axis) 
%   title_tag : string (figure title)
%   save_figs : boolean - save figure? (0=no; 1=yes)
%   save_figs : string (figure name)
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot settings

xLimits = [10 300];
yLimits = [0.07 1.2];
xTikz = [10 20 50 100 200 350];
yTikz = [0.1 0.2 0.5 1];

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves
a = loglog(xRef, yRef, 'b--', 'Linewidth', .5); hold all; % reference results
b = loglog(xRef, (yRef)+0.1, 'b:', 'Linewidth', 1); % reference results - +0.1 (asper) tolerance
loglog(xRef, (yRef)-0.1, 'b:', 'Linewidth', 1); % reference results - - 0.1 (asper) tolerance

% plot results from SQAT
c = loglog(xSQAT, ySQAT,'ko','MarkerSize',8); 

lgd = legend([a, b, c], 'Fastl \& Zwicker', ...
    'Tolerance: $\pm0.1$ (asper)', ...
    'SQAT', ...
    'Location', 'Best', 'Interpreter', 'Latex' );

% lgd.Title.String = '$f_{\mathrm{c}}=2$~kHz';
title(title_tag, 'Interpreter', 'Latex' );

xlim(xLimits);
ylim(yLimits);
grid on;

ax = gca;
set(ax, 'XTick', xTikz);
set(ax, 'YTick', yTikz);
% ax.XAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues =  0:10:160;
% ax.YAxis.MinorTick = 'on';
% ax.YAxis.MinorTickValues = 0:0.1:1.2;

ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    figname_short = save_tag;
    figname_out = [figures_dir figname_short];

    % saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
