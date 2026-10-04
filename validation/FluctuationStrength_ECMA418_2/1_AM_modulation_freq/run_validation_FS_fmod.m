% Script run_validation_FS_fmod
%
%  Verification of the fluctuation strength implementation according to
%  ECMA-418-2:2025
%  Verification case: 100% AM tones (fc = 1 kHz) in dependence of the
%  modulation frequency
%
% - Inputs, ref. data: fluctuation strength of AM tones from Fastl & Zwicker
%   (same reference data used in validation/FluctuationStrength_Osses2016/1_AM_tones_fmod):
%   Fastl, H., & Zwicker, E. (2007). Psychoacoustics: facts and models,
%   Third edition. Springer-Verlag.
%
% - Inputs, signals: AM tones, fc = 1 kHz, md = 1, SPL = 70 dB,
%   fmod = [1 2 4 8 16 32] Hz (same signals used in
%   validation/FluctuationStrength_Osses2016/1_AM_tones_fmod)
%
% Fluctuation strength computed using:
%   OUT = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help FluctuationStrength_ECMA418_2> for more info
%
%  In order to run this code, the user needs to download the dataset of
%  sound files from zenodo (https://doi.org/10.5281/zenodo.7933206).
%  The obtained folder called `validation_SQAT_v1_0` has to be included in
%  the `sound_files` folder of the toolbox.
%
% Authors: Sergio Aguirre & Gil Felix Greco, 17.09.2026
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];

% Figure where the figures (and the results) will be stored:
figures_dir = [dir_out 'figs' filesep];
if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

fname_res = 'FS_results.mat'; % name for the result file
fname_res_full = [figures_dir fname_res];

bCalculation = ~exist(fname_res_full,'file');
if bCalculation == 0
    fprintf('%s.m: Results file found on disk!\n', fname_res);
    fprintf('Do you want to load those prestored results or re-run the calculations?\n');
    bCalculation = input('Enter your choice (1=re-run; 0=read stored results): ');
end
bLoad = ~bCalculation;

%% reference data (Fastl & Zwicker)

ref=[1    2    4    8    16    32    ; % fmod (Hz)
     0.39 0.84 1.25 1.30  0.36  0.06]; % FS values (vacil)

%% compute fluctuation strength from signals

SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('FluctuationStrength_Osses2016',SQAT_version);

dBFS_in  = 100; % dB full scale convention from the input sounds
dBFS_out =  94; % dB full scale convention in SQAT: amplitude of 1 = 1 Pa, or 94 dB SPL
dB_correction = dBFS_in - dBFS_out;

N_signals = size(ref,2);

res=cell([N_signals 1]);  % declaring for memory allocation

if bCalculation

    fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
    time_skip = 700e-3; % time_skip, in seconds for statistical calculations (default: 0.7 seconds)

    tic
    for i=1:N_signals

        fname = sprintf('%sAM-tone-fc-1000_fmod-%.0f_mdept-100-SPL-70-dB.wav',dir_sounds,ref(1,i));
        [insig,fs]=audioread(fname);

        insig  = insig * 10^(dB_correction/20); % correct rms SPL to desired levelOut

        res{i} = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip);
    end
    t_calculation=toc/60; % time to compute fluctuation strength in minutes

    %% saving results so is not need to run fluctuation strength calculation again

    save(fname_res_full,'res','t_calculation');

end
if bLoad
    load(fname_res_full);
end

%% store fluctuation strength values (90th percentile) in vector results[1,nfmod]

results = zeros(1,N_signals);
for i=1:N_signals
    results(i)=res{i}.fluctStrength90Pc;
end

%% plot results

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves
hAx=axes;                     % new axes; save handle
err= (0.1*ref(2,:));
errorbar(ref(1,:),ref(2,:),err,'k-*','MarkerSize',6,'Linewidth',.5);hold all % reference data +/- JND

% plot computed results
semilogx(ref(1,:),results,'ko:','MarkerSize',8);hold all;

hAx.XScale='log';              % turn to semilogx form

legend('Fastl \& Zwicker $\pm\:10\:\%\:(\mathrm{JND})$','SQAT','Location','NW','Interpreter','Latex');
legend boxoff

axis([0 32 0 1.6]);

ax = gca;
set(ax,'XTick',[0 1 2 4 8 16 32]);
set(ax,'YTick',[0 0.4 .8 1.2 1.6]);
ax.XAxis.MinorTick = 'off';
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:0.1:1.8;

title('$f_{\mathrm{c}}=1$~kHz, $L_{\mathrm{p}}=70$~dB~SPL', 'Interpreter', 'Latex' );
ylabel('Fluctuation strength, $F$ (vacil$_{\mathrm{HMS}}$)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    figname_short = 'validation_FS_fmod_1k';
    figname_out = [figures_dir figname_short];

    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');

    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
