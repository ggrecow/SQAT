% Script validation_Sharpness_DIN45692_narrowband_and_broadband_signals
%
% - Validation of sharpness calculation of narrowband and broadband using 
%   test signals from DIN 45692(2009)
%
% - Sharpness computed using:
%   OUT = Sharpness_DIN45692_from_loudness(SpecificLoudness, weight_type, time, time_skip, show_sharpness)
%   type <help Sharpness_DIN45692_from_loudness> for more info
%
% - Specific Loudness is computed within this code according to 
%   ISO 532:1-2017, stationary, free-field
%
% - Loudness computed using:
%   OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
%   type <help Loudness_ISO532_1> for more info
%
% - Test signals generated within this code, following DIN 45692 [1]: ideal
%   bandpass noise (flat spectrum, random phase), one critical band wide
%   (Table A.2) or from a variable lower cutoff frequency to 10 kHz
%   (Table A.3). The band edges are those of the critical bands of Zwicker
%   [2], whose centre frequencies are the ones of Table A.2. Each signal is
%   set to a loudness of 4 sone (ISO 532-1, stationary), as Tables A.2 and
%   A.3 and Section 6 of [1] require. The tolerance of [1] is 5 % or
%   0.05 acum; the stricter 0.05 acum is used here.
%
% - Before November 2026, this code used the test signals provided in the
%   Zenodo dataset of SQAT v1.0, doi:10.5281/zenodo.7933206,
%   all calibrated against the 1 kHz tone of the standard at 60 dB SPL,
%   which gives them 4.1 to 4.5 sone instead of 4 sone (issue #90). That
%   dataset is kept unchanged so that the SQAT v1.x results remain
%   reproducible (see README.md).
%
% - References:
%   [1] DIN 45692:2009. Messtechnische Simulation der Hoerempfindung
%       Schaerfe.
%   [2] Zwicker, E. (1961). Subdivision of the audible frequency range into
%       critical bands (Frequenzgruppen). J. Acoust. Soc. Am. 33(2), 248.
%
% Author: Gil Felix Greco, Braunschweig 01.03.2023
%
% Modified: Gil Felix Greco, November 2026 (test signals generated within
%   the code and set to 4 sone, instead of the DIN 45692 files calibrated
%   against a 1 kHz tone; see issues #89 and #90)
%
% AI disclosure: modifications performed in November 2026 were assisted
% by Claude Opus 5.5 (Anthropic). All codes were verified by the authors.
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
clear all; close all; clc;

save_figs = 0; % save figure flag

%% test signal settings

fs       = 48000;      % sampling frequency (Hz)
N        = 5*fs;       % 5 s per signal
N_target = 4;          % loudness of the test signals (sone), Tables A.2 and A.3 of [1]
f_upper  = 10000;      % upper cutoff frequency of the broadband signals (Hz), Table A.3 of [1]

% critical band edges (Hz) from Zwicker [2], at z = 0, 1, ..., 24 Bark
z_edges = [0 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 ...
           2700 3150 3700 4400 5300 6400 7700 9500 12000 15500];

f  = (0:N-1).'*fs/N; f(f > fs/2) = f(f > fs/2) - fs;
fa = abs(f);           % absolute frequency of each FFT bin
rng(1);                % the noise is random, so that the run is repeatable

%% Sharpness of narrowband test sounds

% Tabelle A.2 - Sollwerte der Schaerfe fuer frequenzgruppenbreites Schmalbandrauschen mit 
% der Lautheit 4 sone (entsprechend dem Standard-Pruefsignal bei 1 kHz und 60 dB)

% Reference values from the standard DIN45692 for the narrowband test sounds
freq_narrowband = [250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400 4000 4800 5800 7000 8500 10500];
bark_narrowband = 2.5:22.5;
acum_narrowband = [.38 .49 .60 .71 .82 .93 1 1.13 1.26 1.35 1.49 1.64 1.78 2.06 2.4 2.82 3.48 4.43 5.52 6.81 8.55];

% generate narrowband signals, set to 4 sone, compute sharpness 
N_narrow = zeros(1,length(freq_narrowband));
for i=1:length(freq_narrowband)

    % one critical band around the centre frequency, e.g. 920 to 1080 Hz at 1 kHz
    f_lo = z_edges(floor(bark_narrowband(i))+1);
    f_hi = z_edges(floor(bark_narrowband(i))+2);
    insig = il_bandpass_noise(f_lo, f_hi, fa, N);

    % stationary loudness, at 4 sone
    [~, L_narrow] = il_set_loudness(insig, fs, N_target);
    N_narrow(i) = L_narrow.Loudness;
                                                                       
    % sharpness calculation (standard weighting function)                                         
    s_narrow(i) = Sharpness_DIN45692_from_loudness(L_narrow.SpecificLoudness,... % input (stationary) specific loudness
                                                                 'DIN45692');    % type of weighting function used for sharpness calculation
                                                 
    % sharpness calculation (aures weighting function)                                         
    s_narrow_aures(i) = Sharpness_DIN45692_from_loudness(L_narrow.SpecificLoudness,... % input (stationary) specific loudness
                                                                          'aures');    % type of weighting function used for sharpness calculation
                                                  
    % sharpness calculation (von bismarck weighting function)                                         
    s_narrow_bismarck(i) = Sharpness_DIN45692_from_loudness(L_narrow.SpecificLoudness,... % input (stationary) specific loudness
                                                                          'bismarck');    % type of weighting function used for sharpness calculation
                                                  
end

%% Sharpness of broadband test sounds

% Tabelle A.3: Sollwerte der Schaerfe fuer Breitbandrauschen mit fester oberer Eckfrequenz 
% von 10 kHz und variabler unterer Eckfrequenz fu mit der Lautheit 4 sone (entsprechend dem Standard-Pruefsignal bei 1 kHz und 60 dB)

% Reference values from the standard DIN45692 for the broadband test sounds
freq_broadband= [250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400 4000 4800 5800 7000 8500];
bark_broadband=2.5:21.5;
acum_broadband=[2.7 2.74 2.78 2.85 2.91 2.96 3.05 3.12 3.20 3.3 3.42 3.53 3.69 3.89 4.12 4.49 5.04 5.69 6.47 7.46];

% generate broadband signals, set to 4 sone, compute sharpness 
N_broad = zeros(1,length(freq_broadband));
for i=1:length(freq_broadband)

    % from the lower cutoff frequency fu to 10 kHz
    insig = il_bandpass_noise(freq_broadband(i), f_upper, fa, N);

    % stationary loudness, at 4 sone
    [~, L_broad] = il_set_loudness(insig, fs, N_target);
    N_broad(i) = L_broad.Loudness;
                                                    
    % sharpness calculation (standard weighting function)                                         
    s_broad(i) = Sharpness_DIN45692_from_loudness(L_broad.SpecificLoudness,... % input (stationary) specific loudness
                                                               'DIN45692');    % type of weighting function used for sharpness calculation

    % sharpness calculation (aures weighting function)                                         
    s_broad_aures(i) = Sharpness_DIN45692_from_loudness(L_broad.SpecificLoudness,... % input (stationary) specific loudness
                                                                        'aures');    % type of weighting function used for sharpness calculation
                                                  
    % sharpness calculation (von bismarck weighting function)                                         
    s_broad_bismarck(i) = Sharpness_DIN45692_from_loudness(L_broad.SpecificLoudness,... % input (stationary) specific loudness
                                                                        'bismarck');    % type of weighting function used for sharpness calculation
                                                  
end

%% report

err_narrow = vertcat(s_narrow.Sharpness)' - acum_narrowband;
err_broad  = vertcat(s_broad.Sharpness)'  - acum_broadband;
fprintf('Loudness of the test signals: %.3f to %.3f sone (target %g sone)\n', ...
        min([N_narrow N_broad]), max([N_narrow N_broad]), N_target);
fprintf('Narrowband signals (Table A.2): max |error| = %.3f acum, %d of %d within 0.05 acum\n', ...
        max(abs(err_narrow)), sum(abs(err_narrow) <= 0.05), numel(err_narrow));
fprintf('Broadband signals (Table A.3):  max |error| = %.3f acum, %d of %d within 0.05 acum\n', ...
        max(abs(err_broad)), sum(abs(err_broad) <= 0.05), numel(err_broad));

%% plot validation for narrowband and broadband

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% reference curves
plot(bark_narrowband,acum_narrowband,'k','Linewidth',1);hold all
plot(bark_broadband,acum_broadband,'b','Linewidth',1);hold all

% computed values
plot(bark_narrowband,vertcat(s_narrow.Sharpness),'ok','Linewidth',1,'MarkerSize',7);
plot(bark_broadband,vertcat(s_broad.Sharpness),'ob','Linewidth',1,'MarkerSize',7);

% title('DIN 45692 - Sharpness of narrowband noise');

legend('DIN 45692:2009 - ref. for narrowband signals','DIN 45692:2009 - ref. for broadband signals','SQAT - narrowband signals','SQAT - broadband signals','Location','NW','Interpreter','Latex');
legend boxoff
 
grid off
axis([0 24 0 11]);
ylabel('Sharpness, $S$ (acum)','Interpreter','Latex');
xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');

ax = gca;
set(ax,'XTick',[0 4 8 12 16 20 24]);
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =  0:0.5:24; 
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:0.5:20; 
     
% hold off
set(gcf,'color','w');

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = 'sharpness_validation_narrowband_and_broadband';
    figname_out = [figures_dir figname_short];
    
    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%%  plot error for narrowband and broadband
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

a=plot(bark_narrowband,-0.05*ones(1,21),'r--','Linewidth',1);hold on; % error tolerance of -5%
plot(bark_narrowband,0.05*ones(1,21),'r--','Linewidth',1);  % error tolerance of +5%

b=plot(bark_narrowband,vertcat(s_narrow.Sharpness)'-acum_narrowband,'k','Linewidth',1);hold all
c=plot(bark_broadband,vertcat(s_broad.Sharpness)'-acum_broadband,'b','Linewidth',1);hold all

legend([a,b,c],'DIN 45692:2009 tolerance','Narrowband signals','Broadband signals','Location','NW','Interpreter','Latex');
legend boxoff
 
% title('DIN 45692 - Sharpness of narrowband noise');

grid off
axis([0 24 -.12 .12]);
ylabel('SQAT minus DIN 45692:2009 (acum)','Interpreter','Latex');
% ylabel('Relative error [acum]','Interpreter','Latex');

xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');

set(gcf,'color','w');

ax = gca;
set(ax,'XTick',[0 4 8 12 16 20 24]);
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =  0:0.5:24; 
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.2:0.01:0.2;
     
if save_figs==1
    figname_short = 'sharpness_validation_narrowband_broadband_error';
    figname_out = [figures_dir figname_short];
    
    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% narrow band comparison between models 

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

plot(bark_narrowband,acum_narrowband,'k','Linewidth',1);hold all
plot(bark_narrowband,vertcat(s_narrow.Sharpness),'ob','Linewidth',1,'MarkerSize',7);
plot(bark_narrowband,vertcat(s_narrow_aures.Sharpness),'+r','Linewidth',1,'MarkerSize',8);
plot(bark_narrowband,vertcat(s_narrow_bismarck.Sharpness),'xc','Linewidth',1,'MarkerSize',8);

% title('DIN 45692 - Sharpness of narrowband noise');

legend('Ref. values from DIN 45692:2009','SQAT - DIN 45692:2009','SQAT - Aures','SQAT - von Bismarck','Location','NorthWest','Interpreter','Latex');
legend boxoff
 
grid off
axis([0 24 0 10]);
ylabel('Sharpness, $S$ (acum)','Interpreter','Latex');
xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');

 ax = gca;
 set(ax,'XTick',[0 4 8 12 16 20 24]);
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues =  0:0.5:24; 
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = 0:0.5:10; 
     
% hold off
set(gcf,'color','w');

if save_figs==1
    figname_short = 'sharpness_validation_narrowband_compare_models';
    figname_out = [figures_dir figname_short];
    
    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% narrow band comparison
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

plot(bark_broadband,acum_broadband,'k','Linewidth',1);hold all

plot(bark_broadband,vertcat(s_broad.Sharpness),'ob','Linewidth',1,'MarkerSize',7);
plot(bark_broadband,vertcat(s_broad_aures.Sharpness),'+r','Linewidth',1,'MarkerSize',8);
plot(bark_broadband,vertcat(s_broad_bismarck.Sharpness),'xc','Linewidth',1,'MarkerSize',8);

% title('DIN 45692 - Sharpness of broadband noise');

legend('Ref. values from DIN 45692:2009','SQAT - DIN 45692:2009','SQAT - Aures','SQAT - von Bismarck','Location','NorthWest','Interpreter','Latex');
legend boxoff
 
legend boxoff
 
grid off
axis([0 24 0 10]);
ylabel('Sharpness, $S$ (acum)','Interpreter','Latex');
xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');

ax = gca;
 set(ax,'XTick',[0 4 8 12 16 20 24]);
     ax.XAxis.MinorTick = 'on';
     ax.XAxis.MinorTickValues =  0:0.5:24; 
     ax.YAxis.MinorTick = 'on';
     ax.YAxis.MinorTickValues = 0:0.5:10; 
     
% hold off
set(gcf,'color','w');

if save_figs==1
    figname_short = 'sharpness_validation_broadband_model_comparison';
    figname_out = [figures_dir figname_short];
    
    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% local functions

function x = il_bandpass_noise(f_lo, f_hi, fa, N)
% function x = il_bandpass_noise(f_lo, f_hi, fa, N)
%
% Generates ideal bandpass noise: a flat magnitude spectrum between f_lo and
% f_hi (Hz), zero elsewhere, with a random phase in each frequency bin. The
% spectrum is built directly in the frequency domain, so the band edges are
% exact (no filter slopes) and the signal is periodic in N samples (no
% onset or offset transients).
%
% INPUT ARGUMENTS
%   f_lo, f_hi : scalars
%   lower and upper edges of the band (Hz), both included
%
%   fa : [Nx1] array
%   absolute frequency of each bin of an N-point FFT (Hz), i.e. abs(f) with
%   f = (0:N-1)*fs/N mapped to [-fs/2, fs/2). The frequency resolution is
%   fs/N (0.2 Hz for 5 s at 48 kHz)
%
%   N : integer
%   number of samples (must be even)
%
% OUTPUT
%   x : [Nx1] array
%   real-valued noise with unit rms (to be scaled to the desired level)
%
% The phase is drawn with rand, so the noise depends on the state of the
% random number generator (set with rng before the first call).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mag = zeros(N,1);
mag(fa >= f_lo & fa <= f_hi) = 1;                % unit magnitude inside the band
ph  = exp(1j*2*pi*rand(N,1)); ph(1) = 1;         % random phase, uniform in [0, 2*pi); DC bin real

X   = mag.*ph; X(N/2+2:end) = conj(flipud(X(2:N/2))); % Hermitian symmetry: negative frequencies are the complex conjugate of the positive ones
x   = real(ifft(X));                             % real signal (real() removes round-off imaginary parts)
x   = x/rms(x);                                  % unit rms

end % end of function <il_bandpass_noise>

function [x, L] = il_set_loudness(x, fs, N_target)
% function [x, L] = il_set_loudness(x, fs, N_target)
%
% Scales the signal x so that its stationary loudness according to
% ISO 532-1 (free field, time_skip of 0.5 s, the same settings used for the
% validation) equals N_target, as required by DIN 45692 for its test
% signals.
%
% The gain is found by bisection, which relies on the loudness increasing
% monotonically with the level. The search interval is 30 to 110 dB SPL; the
% geometric mean of the interval limits is used as the midpoint, so that
% the interval is halved in dB. After 14 iterations the interval is
% 80/2^14 = 0.005 dB wide, which gives the target loudness to about 0.001
% sone at 4 sone.
%
% INPUT ARGUMENTS
%   x : [Nx1] array
%   signal of arbitrary amplitude (it is normalised to unit rms here)
%
%   fs : integer
%   sampling frequency (Hz)
%
%   N_target : scalar
%   target loudness (sone)
%
% OUTPUTS
%   x : [Nx1] array
%   scaled signal (Pa)
%
%   L : struct
%   output of Loudness_ISO532_1 for the scaled signal (type <help
%   Loudness_ISO532_1> for more info)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x  = x/rms(x);                                   % unit rms, so that the gain g is the rms value in Pa
lo = 20e-6*10^(30/20); hi = 20e-6*10^(110/20);   % gain range: 30 to 110 dB SPL

for it = 1:14
    g = sqrt(lo*hi);                             % midpoint of the interval, in dB
    if Loudness_ISO532_1(x*g, fs, 0, 1, 0.5, 0).Loudness < N_target, lo = g; else, hi = g; end
end

x = x*sqrt(lo*hi);                               % midpoint of the final interval
L = Loudness_ISO532_1(x, fs, 0, 1, 0.5, 0);      % loudness of the scaled signal

end % end of function <il_set_loudness>