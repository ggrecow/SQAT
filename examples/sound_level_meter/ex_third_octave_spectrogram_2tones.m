%   Script ex_third_octave_spectrogram_2tones
%
% - create one signal containing two sinusoidal tones, with rms SPL level of 60 dB SPL each
%   and central frequencies of 100 Hz and 1 kHz. The signals are created with silence before
%   and after the tones in order to verify the temporal response of the fast and slow
%   time-weightings
%
% - Decompose the sound pressure signal into 1/3 octave-bands
%
%   1/3 octave-bands are computed using the following function:
%   [outsig, fc] = Do_OB13_ISO532_1(insig, fs, fmin, fmax)
%   type <help Do_OB13_ISO532_1> for more info
%
% - Compute SPL time-histories of each 1/3 octave-bands using the Z- and A-weighting and the
%   Slow- and Fast time-weightings using the <sound_level_meter> of SQAT
%
%   SPL is computed using the following function:
%   [outsig_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS)
%   type <help Do_SLM> for more info
%
% - Plot results
%
%  HOW TO RUN THIS CODE: this is a standalone code. Therefore, no additional steps are
%  necessary to run this code.
%
%   Gil Felix Greco, 17.11.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Create sound pressure signal containing two sinusoidal tones and
% silence before and after the tones

fc = [100; 1000];      % Center frequency (Hz)
LevelIn = 60;           % Level (dB SPL)
L_before = 5;          % Length of silence before the sinusoidal signal (seconds)
L_signal = 15;         % Length of the sinusoidal signal (seconds)
L_after = 10;           % Length of silence before the sinusoidal signal (seconds)
fs = 48000;             % Sampling frequency

[insig, t_total] = create_two_tones( fc(1), fc(2), LevelIn, L_before, L_signal, L_after, fs );
insig = insig'; % insig must be a row vector

%% Decompose the signal into 1/3 octave-bands

  % filter insig to get 1/3-OB
  fmin = 50; % min freq of 1/3-OB is 50 Hz
  fmax = 10000; % max freq of 1/3-OB is 10 kHz

  % get 1/3-OB spectra from insig - output is p [nTime,nFreq]. OBS: if
  % <fmin, fmax> are not provided, then the default values are fmin=12.5 Hz
  % and fmax=12.5 kHz
  [insig_P_TOB, fnom] = Do_OB13_ISO532_1(insig, fs, fmin, fmax);

%% compute SPL of each third-octave band: Z-weighted

weight_freq = 'Z'; % Z-frequency weighting
weight_time = 's'; % slow leak - time constant = 1s

SPL_Z_grid = zeros( size(insig_P_TOB,1) ,size(insig_P_TOB,2) ); % declare variable for memory allocation

for i = 1:size(insig_P_TOB,2)
    SPL.Z_weighted(:,i) = Get_SPL(insig_P_TOB(:,i), fs, weight_freq, weight_time);
    SPL_Z_grid(:,i) = SPL.Z_weighted(i).lvl_dB; % [nTime x nBands]
end

%% compute SPL of each third-octave band: A-weighted

weight_freq = 'A'; % A-frequency weighting
weight_time = 'f'; % slow leak - time constant = 1s

SPL_A_grid = zeros( size(insig_P_TOB,1) ,size(insig_P_TOB,2) ); % declare variable for memory allocation

for i = 1:size(insig_P_TOB,2) % for each freq. band
    SPL.A_weighted(:,i) = Get_SPL(insig_P_TOB(:,i), fs, weight_freq, weight_time);
    SPL_A_grid(:,i) = SPL.A_weighted(i).lvl_dB; % [nTime x nBands]
end

%% plot spectrogram (1/3 octave bands in dt time steps)

figure('name', '1/3-octave band spectrogram',...
          'units', 'normalized', 'outerposition', [0 0 1 1]); % plot fig in full screen

xmax = t_total(end); % used to define the x-axis on the plots

% plot input signal
subplot( 2, 2, [1,2] )
plot( t_total, insig ); hold on;
try
    a=yline(rms(insig),'k--');
catch
    % In Octave:
    XL = get(gca,'XLim');
    a = plot(XL,rms(insig)*[1 1],'k--');
end
legend(a,sprintf('$p_{\\mathrm{rms}}=$%g Pa',rms(insig)),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('Sound pressure, $p$ (Pa)','Interpreter','Latex'); %grid on;
ax = axis; axis([0 xmax max(insig)*-2 max(insig)*2]);
title('Input signal','Interpreter','Latex');

% spectrogram - Z-weighted, Slow time-weighting
subplot( 2, 2, 3 )
[xx,yy] = meshgrid( t_total, fnom );
pcolor( xx, yy, SPL_Z_grid' );
shading interp; colorbar; axis tight;
colormap jet;

% freq labels
ax=gca;
set(ax,'YScale', 'log'); % set y-axis to log-scale
ylabel('Center frequency, $f$ (Hz)','Interpreter','Latex');
try
    ax.YAxis.MinorTick = 'on';
    ax.YAxis.MinorTickValues = fnom;
end

xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel(colorbar,'SPL, $L_{\mathrm{Z,S}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex');
clim([0 max(max(SPL_Z_grid))]);
set(ax,'layer','top');
box on
title('Spectrogram (1/3 oct. bands, Z-weighting, Slow time-weighting)', 'Interpreter', 'Latex');

% spectrogram - A-weighted, Fast time-weighting
subplot( 2, 2, 4 )
[xx,yy] = meshgrid( t_total, fnom );
pcolor( xx, yy, SPL_A_grid' );
shading interp; colorbar; axis tight;
colormap jet;

% freq labels
ax=gca;
set(ax,'YScale', 'log'); % set y-axis to log-scale
ylabel('Center frequency, $f$ (Hz)','Interpreter','Latex');
try
    ax.YAxis.MinorTick = 'on';
    ax.YAxis.MinorTickValues = fnom;
end

xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel(colorbar,'SPL, $L_{\mathrm{A,F}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex');
clim([0 max(max(SPL_A_grid))]);
set(ax,'layer','top');
box on
title('Spectrogram (1/3 oct. bands, A-weighting, Fast time-weighting)', 'Interpreter', 'Latex');

set(gcf,'color','w');
