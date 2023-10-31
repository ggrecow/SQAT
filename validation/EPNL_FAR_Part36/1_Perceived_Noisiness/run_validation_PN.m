% Script run_validation_PN
%
% Verification of the Perceived Noisiness calculation used
% in the EPNL code
%
%   - The Effective Perceived Noise Level (EPNL) is a metric used mainly in
%   the context of environmental aircraft noise, mainly for aircraft noise
%   certification. To the best of knowledge of the Author of this code, a
%   reference signal to validate/verify the implementation of the EPNL
%   calculation does not exist. However, it is possible to verify
%   particular parts of the EPNL calculation.
%
%   In this code, we use a test signal in order to verify 1) the whole
%   pre-processing stage of the input signal, and 2) the conversion of SPLs
%   to Perceived Noisiness.
%
%   - Reference of the Perceived Noisiness, PN (Noys): the numerical value of 1
%     Noy was assigned to the PN of a third octave of random noise-centered
%     at 1 kHz, and with a rms level of 40 dB SPL.
%
%       * in this code, we verify if the PN is being correctly calculated by generating
%          this reference signal as a white noise, which is band-passed between
%          a lower freq 891.3 Hz and a upper freq. 1122 Hz, which corresponds to the
%          third-octave band centered around fc=1 kHz, and a rms SPL of 40 dB SPL.
%          PLEASE NOTE: computing the EPNL from this signal does not make any sense as this is not
%          an aircraft flyover nor (more generally) a signal with time-varying amplitude.
%          It is only used here for the purpose of verifying the PN implementation.
%
%       Source: Smith, M. (1989). Human reaction to aircraft noise.
%       In Aircraft Noise (Cambridge Aerospace Series, pp. 1-19). Cambridge University Press.
%       doi:10.1017/CBO9780511584527.002 (relevant info can be found in page 7)
%
%       Source: Bennett, R.L. & Persons, K.S. (1981) Handbook of aircraft noise metrics.
%       NASA - Contractor Report no. NASA CR-3406 (permanent link:
%       https://ntrs.nasa.gov/citations/19810013341 - Last viewed 27 Oct.
%       2023) (relevant info can be found in page 53)
%
% EPNL (and/or associated quantities) computed using the following function:
%   OUT = EPNL_FAR_Part36( insig, fs, method, dt, threshold, show )
%   type <help EPNL_FAR_Part36> for more info
%
% Author: Gil Felix Greco, Braunschweig, 30.10.2023
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

%% generate band pass white noise

L = 60;          % signal length (seconds)
fs = 48000;     % sampling freq
levelOut = 40;  % desired rms level (dBSPL)
fc= 1000; % center freq of the 1/3 oct. band
flow = fc.*2^(-1/6); % lower freq of the toct band centered at fc
fhigh = fc.*2^(1/6); % upper freq of the toct band centered at fc

RefSignal_PN = il_make_bandpass_white_noise(L, fs, levelOut, flow,fhigh);

%% compute EPNL of  RefSignal_PN to get the Perceived Noise, in Noys

compute_PN = EPNL_FAR_Part36(RefSignal_PN, fs);

%% plot results of EPNL calculation from RefSignal_PN

figure('name','EPNL calculation - verification of PN implementation',...
    'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen

xmax = compute_PN.time(end); % used to define the x-axis on the plots

% plot input signal
subplot( 2, 6, [1,2] )
plot( compute_PN.time_insig, RefSignal_PN );
a=yline(rms(RefSignal_PN),'k--');
legend(a,sprintf('$p_{\\mathrm{rms}}=$%g Pa',rms(RefSignal_PN)),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('Sound pressure, $p$ (Pa)','Interpreter','Latex'); %grid on;
ax = axis; axis([0 xmax max(RefSignal_PN)*-2 max(RefSignal_PN)*2]);
title('Input signal','Interpreter','Latex');

% plot instantaneous sound pressure level (dBSPL) from original signal and time-averaged over a given dt value
subplot( 2, 6, [3,4] )
plot( compute_PN.time_insig, compute_PN.InstantaneousSPL_insig ); hold on;
plot( compute_PN.time, compute_PN.InstantaneousSPL,'Linewidth',2 );
legend(sprintf('dt=%g sec',1/fs),sprintf('dt=%g sec - Mean value: %.4g dB SPL',0.5,mean(compute_PN.InstantaneousSPL)),'Location','SouthWest','Interpreter','Latex'); %legend boxoff;
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex'); grid on;
ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
title('Instantaneous overall SPL (1/3 oct. bands)','Interpreter','Latex');

% plot spectrogram (1/3 octave bands in dt time steps)

subplot( 2, 6, [5,6] )
fnom = compute_PN.TOB_freq./1000; % convert center freq to kHz to plot ?
[xx,yy] = meshgrid( compute_PN.time, fnom );
pcolor( xx,yy,compute_PN.SPL_TOB_spectra' );
shading interp; colorbar; axis tight;
colormap jet;

% freq labels
ax=gca;
%         set(ax,'YScale', 'log');
%     set(ax,'YTick', fnom);
yticks( [fnom(1) fnom(14) fnom(17:24)] );
ylabel('Center frequency, $f$ (kHz)','Interpreter','Latex');
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = fnom;

xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel(colorbar,'SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex');
caxis([0 max(max(compute_PN.SPL_TOB_spectra))]);
set(ax,'layer','top');
box on
title(sprintf('Spectrogram (1/3 oct. bands, dt=%.g sec)',0.5),'Interpreter','Latex');

% plot perceived noisiness (noys vs. time)
subplot( 2, 6, [7,8,9,10,11,12] )
plot( compute_PN.time, compute_PN.PN ); hold on;
a = yline(mean(compute_PN.PN),'k--');
legend(a,sprintf('Mean PN = %.4g noys',mean(compute_PN.PN)));
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('PN (noys)','Interpreter','Latex'); grid on;
ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
title('Perceived noisiness','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'validation_Perceived_Noisiness';
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function il_make_bandpass_white_noise

function OutputSignal= il_make_bandpass_white_noise(L, fs, levelOut, flow,fhigh)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates a bandpass white noise with a desired bandwidth and rms SPL
%
% INPUT
%  L: signal's length (seconds)
%  fs: sampling frequency (Hz)
%  levelOut: desired rms level (dBSPL)
%  flow: lower freq. used for the bandpass (Hz)
%  fhigh: upper freq. used for the bandpass (Hz)
%
% OUTPUT:
%  OutputSignal: generated signal
%
% This function uses two subfunctions (il_bpnoise.m and il_cut.m) from the AFC toolbox
%
% Alejandro Osses, 15.06.2023
% Gil Felix Greco, Braunschweig 20.10.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save = 0;

%% band pass noise

% L = 30;              % signal length (seconds)
% fs = 48000;         % sampling freq
dt = 1/fs;          % time resolution
t = 0:dt:L-dt;    % Time vector (trick: -dt to have a nice number of samples)
N_samples = length(t);

p0 = 2e-5; % Pa, reference pressure

% flow  = 891.3;  % lower freq of the 1-khz toct octave band
% fhigh = 1122;  % upper freq of the 1-khz toct octave band

bpass = il_bpnoise(N_samples,flow,fhigh,fs); % generates bandpass noise. It has an RMS of '1'

% Calibrate level of the band passed signal
% levelOut=40;  % desired level (dBSPL)

levelIn = 20*log10(rms(bpass)/p0);
cal_factor = 10^((levelOut-levelIn)/20);
bpass = cal_factor*bpass;
bpass = transpose(bpass);

%% final signal

OutputSignal=bpass; % generate final signal


%% save wave file

if save==1
    FileName='RefSignal_PerceivedNoisiness_40dBSPL.wav';
    audiowrite(FileName,OutputSignal,fs,'BitsPerSample',64)
else
end

end

%% bpnoise.m function

function out = il_bpnoise(len,flow,fhigh,fs)
% bpnoise.m - generates spectrally rectangular shaped band-pass noise.
%
% Usage: out = bpnoise(len,flow,fhigh,fs)
%
% len    = output length in samples
% flow   = lower cutoff frequency in Hz
% fhigh  = upper cutoff frequency in Hz
% fs     = sampling rate in Hz
%
% out    = output vector

% out = real(ifft(scut(fft(randn(len,1)),flow,fhigh,fs)));
% out = real(ifft(scut(1*exp(1j*rand(len,1)),flow,fhigh,fs)));

singleSidedSpectrum=1*exp(1j*(2*pi*(rand(len/2+1,1)-0.5)));
doubleSidedSpectrum=[singleSidedSpectrum;conj(flipud(singleSidedSpectrum(2:end-1)))];
cutSpectrum = il_scut(doubleSidedSpectrum,flow,fhigh,fs);
out = real(ifft(cutSpectrum));
out = out/(norm(out,2)/sqrt(len));

end

%% scut.m function

function cut = il_scut(in,flow,fhigh,fs)
% scut.m - cuts rectangular portion of a signal's fourier transform.
%			  Fourier coefficients outside the passband specified by flow
%			  and fhigh are set to zero.
%
% Usage: cut = scut(in,flow,fhigh,fs)
%
% in     = input vector containing a signal's discrete fourier transform
% flow   = lower cutoff frequency in Hz
% fhigh  = upper cutoff frequency in Hz
% fs     = sampling rate in Hz
%
% cut    = output vector

len = length(in);
flow = round(flow*len/fs);
fhigh = round(fhigh*len/fs);
cut = zeros(len,1);
cut(flow+1:fhigh+1) = in(flow+1:fhigh+1);

% HACK: if lowpass ( flow = 0) index would be greater than len (len +1)
if flow == 0
    flow = 1;
end

cut(len-fhigh+1:len-flow+1) = in(len-fhigh+1:len-flow+1);

end
