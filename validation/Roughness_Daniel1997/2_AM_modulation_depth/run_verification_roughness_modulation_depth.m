% Script run_verification_roughness_modulation_depth
%  
% Verification of DW roughness code for roughness dependence on the 
%  modulation depth for AM tones
%
% - Inputs, signals: 1 kHz sinusoidal tone modulated at 70 Hz (70 dBSPL) 
%   with varying modulation depth, md (from 0 to 1, in 0.05 increments).
%   The signals are generated on the fly by the local function at the end 
%   of this script. The AM tone follows Eq. (1) of Daniel & Weber (1997):
%       p(t) = p0*[1 + md*cos(2*pi*fm*t)]*cos(2*pi*fc*t)
%
% - Ref. values are the power law R = 1.36*md^(1.6), from 
%   Daniel, P. and Weber, R. Psychoacoustical Roughness: Implementation of
%   an Optimized Model, Acustica 83 (1997) 113-123, Fig. 5 (fc = 1 kHz,
%   fmod = 70 Hz, L = 70 dBSPL). The exponent 1.6 is from Zwicker, E. and
%   Fastl, H. Second ed, Psychoacoustics, Facts and Models, page 258.
%
%   NOTE: the test conditions of this reference curve (70 dBSPL, prefactor
%   1.36) differ from the conditions defining the asper unit (1 kHz, 70 Hz,
%   md = 1, 60 dBSPL -> 1 asper). Changing <SPL> below therefore requires
%   changing <a> accordingly.
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
% Unlike the other roughness validation cases, this script does NOT require 
% the dataset of sound files from zenodo 
% (https://doi.org/10.5281/zenodo.7933206). The signals are generated 
% locally by this script; see the correction note below.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Log:
%
% Author: Gil Felix Greco, Braunschweig 17.02.2020 (updated in 13.05.2023)
%
% Corrected by Gil Felix Greco, Braunschweig 08.09.2026
%   The version of this script released with SQAT v1.0 loaded pre-generated
%   signals from `vary_modulation_depth.mat` and contained two errors in the
%   test conditions:
%
%   (1) Wrong level and reference anchor. The signals were generated at
%       60 dBSPL and compared against the power law md^1.6 with a unity
%       prefactor. The reference condition for the modulation-depth
%       dependence is 70 dBSPL with R = 1.36*md^1.6 (Daniel & Weber 1997,
%       Fig. 5; see also Schrader, J.E., A MATLAB implementation of a model
%       of auditory roughness, TU Eindhoven, 2002, Fig. 13). The 60 dBSPL /
%       1 asper condition defines the asper unit, but is not the condition
%       of this curve.
%
%   (2) Wrong modulation index. The AM tones were generated as
%       [1 - md/2 + (md/2)*cos(2*pi*fm*t)]*cos(2*pi*fc*t), whose envelope
%       spans [1-md, 1] and whose modulation index is therefore md/(2-md),
%       not md. The two definitions coincide only at md = 0 and md = 1, so
%       all intermediate points were plotted against the wrong abscissa.
%
%   Both errors are corrected here. The signals are now generated inline
%   following Eq. (1) of Daniel & Weber (1997) at 70 dBSPL, and the
%   reference curve uses a = 1.36. Results published from the v1.0 version
%   of this script are superseded.
%
%   The zenodo dataset (https://doi.org/10.5281/zenodo.7933206) has 
%   deliberately NOT been updated. Its `vary_modulation_depth.mat` file 
%   still holds the erroneous 60 dBSPL signals and is kept unchanged so 
%   that the SQAT v1.0 record stays reproducible as published. That file is 
%   no longer read by this script and must not be used for this test case.
%
%   Scope: only this verification case is affected. The remaining roughness
%   validation cases use either md = 1, where the two AM definitions are
%   proportional and therefore identical after level calibration, or FM and
%   noise signals, and remain valid as published.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% save figs flag

save_figs = 0;

%% path settings 

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% signal settings

fs  = 48000;      % sampling frequency, Hz
L   = 10;         % signal length, s
fc  = 1000;       % carrier frequency, Hz
fm  = 70;         % modulation frequency, Hz
SPL = 70;         % signal level, dBSPL (Daniel & Weber 1997, Fig. 5)
md  = 0:0.05:1;   % modulation depth vector

%% generate signals

s = il_make_AM_varying_md(md,fc,fm,SPL,fs,L);

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
results = zeros(1,size(s,1));

for i=1:size(s,1)
    results(i)=res{1,i}.Rmean;
end

% ref curve, and jnd of 17% for roughness
a=1.36;
power_law=a.*md.^(1.6);

err= (0.17*power_law);
errorbar(md,power_law,err,'k-*','MarkerSize',6,'Linewidth',.5);hold all % ref. values

% plot results SQAT
plot(md,results,'ko:','MarkerSize',8);hold all % results from SQAT

legend('Ref. - $1.36\,m_{\mathrm{d}}^{1.6}\pm17\:\%\:(\mathrm{JND})$','SQAT','Location','NW','Interpreter','Latex');
legend boxoff

axis([0 1 0 1.8]);
ax = gca;
set(ax,'XTick',[0 0.2 0.4 0.6 .8 1]);
set(ax,'YTick',[0 0.2 0.4 0.6 .8 1 1.2 1.4 1.6 1.8]);
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

%% local function used to generate the signals

function s = il_make_AM_varying_md(md,fc,fm,SPL,fs,L)
% function s = il_make_AM_varying_md(md,fc,fm,SPL,fs,L)
%
% Generates amplitude modulated (AM) tones for roughness verification,
% following Eq. (1) of Daniel & Weber (1997):
%
%   p(t) = p0*[1 + md*cos(2*pi*fm*t)]*cos(2*pi*fc*t)
%
% so that <md> is the modulation index used in the reference literature.
% Each signal is scaled to the desired overall (rms-based) level.
%
% INPUTS:
%   md  : vector, modulation depth (0 to 1)
%   fc  : scalar, carrier frequency, Hz
%   fm  : scalar, modulation frequency, Hz
%   SPL : scalar, signal level, dBSPL
%   fs  : scalar, sampling frequency, Hz
%   L   : scalar, signal length, s
%
% OUTPUT:
%   s   : matrix, length(md) x length(t), one AM tone per row
%
% Gil Felix Greco, Braunschweig 17.02.2020 (updated in 08.09.2026)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1/fs;        % time step
t  = 0:dt:L;      % time vector

pref = 2e-5;      % reference sound pressure, Pa

carrier = cos(2*pi*fc.*t);

s = zeros(length(md),length(t));

for i = 1:length(md)

    modulator = 1 + md(i).*cos(2*pi*fm.*t);

    s(i,:) = carrier.*modulator;

    % calibrate each signal to the desired SPL (rms-based)
    levelIn = 20*log10( sqrt(mean(s(i,:).^2)) / pref );
    s(i,:)  = s(i,:) .* 10^((SPL-levelIn)/20);

end

end