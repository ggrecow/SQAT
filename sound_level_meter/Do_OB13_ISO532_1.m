function [outsig, fc] = Do_OB13_ISO532_1(insig, fs, fmin, fmax)
% [outsig, fc] = Do_OB13_ISO532_1(insig, fs, fmin, fmax)
%
% This code employs a hard-coded implementation of a one-third octave band
%   filterbank for a sampling frequency fs=48000 Hz. For this reason, if
%   the input fs adopts a value different from 48 kHz, the signal is
%   appropriately resampled. 
%
% This code was extracted from the loudness implementation from the AARAE
%   toolbox refactored by the SQAT toolbox team. Therefore, it complies with
%   ISO 532-1:2017 (and IEC 61260-1:2014). However, the freq range computed
%   is only between 25 Hz - 12.5 kHz
%
% INPUTS:
%
%   insig : double
%   [Nx1] array, corresponding to a monophonic audio signal with N time samples 
%
%   fs : integer
%     sampling frequency (Hz) of insig
%
%   fmin (optional, default=25 Hz): integer
%     minimun centre frequency in the range between 25 Hz and 12.5 kHz in 
%     one-third octave-band steps
%
%   fmax (optional, default=12500 Hz): integer
%     maximum centre frequency in the range between 25 Hz and 12.5 kHz
%
% OUTPUT
%
%   outsig : double  
%   [N x nFreq] array containing nFreq (default, nFreq=28) time signals of 
%     length N. Each nFreq column contains the filtered waveform of the 
%     corresponding one-third octave-band signal centred at fc(nFreq).
%
%   fc: double
%   [nFreq x 1] array containing the corresponding centre frequencies of 
%     each one-third octave band with a maximum of nFreq = 28 elements for
%     the default fmin and fmax values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_Loudness_ISO532_1.wav'];
%   [insig,fs] = audioread(fname);
%   dBFS = 94; % A priori knowledge
%   [OB_filt,fc] = Do_OB13_ISO532_1(insig,fs);
%   % Sub example: obtaining the calibrated RMS level of each band:
%   OB_lvls = 20*log10(rms(OB_filt))+dBFS;
%
%   figure;
%   semilogx(fc,OB_lvls,'bo-'); hold on;
%   set(gca,'XLim',[min(fc)-1 max(fc)+1000]);
%   set(gca,'XTick',fc);
%   set(gca,'XTickLabel',round(fc));
%   xlabel('Frequency (Hz)');
%   ylabel('Band level (dB SPL)');
%
%   % Total level of the input signal:
%   lvl_orig = 20*log10(rms(insig))+dBFS;
%   fprintf('Level as obtained directly from the input signal=%.1f dB SPL\n',lvl_orig);
%
%   % The sum of the power of each band leads to the total level:
%   lvl_from_outsig = 10*log10(sum(10.^(OB_lvls/10))); % '/10' means squared, therefore it is 10*log10()
%   fprintf('Level obtained from the filtered signals=%.1f dB SPL\n',lvl_from_outsig);
%
% Author: Ella Manor - MATLAB implementation for AARAE (2015)
% Author (modifications): Gil Felix Greco (22/02/2023)
% Author (stand-alone function): Alejandro Osses (20/10/2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    bLimit_range = 0; 
else
    % Then nargin is 3 or 4
    if nargin < 4
        % Then only fmin has been specified as input
        fmax = 12500; % Hz
    end
    bLimit_range = 1;
    if fmax > 12500
        warning('The input parameter fmax cannot be greater than 12500 Hz, setting fmax to this value');
        fmax = 12500; % Hz, maximum possible frequency of the filterbank
    end
    if fmin < 25
        warning('The input parameter fmin cannot be lower than 25 Hz, setting fmin to this value');
        fmin = 25; % Hz, minimum possible frequency of the filterbank
    end
end

%%  resample to 48 kHz if necessary

if fs ~= 48000
    insig = resample(insig,48000,fs);
    fs = 48000;
    fprintf('%s.m: This script has only been validated at a sampling frequency fs=48 kHz, resampling to this fs value\n',mfilename);
end

len = size(insig,1);

%% Create filter bank and filter the signal

% reference
br = [1,2,1;1,0,-1;1,-2,1];
ar = [1,-2,1;1,-2,1;1,-2,1];

% filter 'a' coefficient offsets TABLES A.1 A.2
ad = cat(3,[0,-6.70260e-004,6.59453e-004;...
            0,-3.75071e-004,3.61926e-004;...
            0,-3.06523e-004,2.97634e-004],... % 25 Hz
           [0,-8.47258e-004,8.30131e-004;...
            0,-4.76448e-004,4.55616e-004;...
            0,-3.88773e-004,3.74685e-004],... % 31.5 Hz
           [0,-1.07210e-003,1.04496e-003;...
            0,-6.06567e-004,5.73553e-004;...
            0,-4.94004e-004,4.71677e-004],... % 40 Hz
           [0,-1.35836e-003,1.31535e-003;...
            0,-7.74327e-004,7.22007e-004;...
            0,-6.29154e-004,5.93771e-004],... % 50 Hz
           [0,-1.72380e-003,1.65564e-003;...
            0,-9.91780e-004,9.08866e-004;...
            0,-8.03529e-004,7.47455e-004],... % 63 Hz
           [0,-2.19188e-003,2.08388e-003;...
            0,-1.27545e-003,1.14406e-003;...
            0,-1.02976e-003,9.40900e-004],... % 80 Hz
           [0,-2.79386e-003,2.62274e-003;...
            0,-1.64828e-003,1.44006e-003;...
            0,-1.32520e-003,1.18438e-003],... % 100 Hz
           [0,-3.57182e-003,3.30071e-003;...
            0,-2.14252e-003,1.81258e-003;...
            0,-1.71397e-003,1.49082e-003],... % 125 Hz
           [0,-4.58305e-003,4.15355e-003;...
            0,-2.80413e-003,2.28135e-003;...
            0,-2.23006e-003,1.87646e-003],... % 160 Hz
           [0,-5.90655e-003,5.22622e-003;...
            0,-3.69947e-003,2.87118e-003;...
            0,-2.92205e-003,2.36178e-003],... % 200 Hz
           [0,-7.65243e-003,6.57493e-003;...
            0,-4.92540e-003,3.61318e-003;...
            0,-3.86007e-003,2.97240e-003],... % 250 Hz
           [0,-1.00023e-002,8.29610e-003;...
            0,-6.63788e-003,4.55999e-003;...
            0,-5.15982e-003,3.75306e-003],... % 315 Hz
           [0,-1.31230e-002,1.04220e-002;...
            0,-9.02274e-003,5.73132e-003;...
            0,-6.94543e-003,4.71734e-003],... % 400 Hz
           [0,-1.73693e-002,1.30947e-002;...
            0,-1.24176e-002,7.20526e-003;...
            0,-9.46002e-003,5.93145e-003],... % 500 Hz
           [0,-2.31934e-002,1.64308e-002;...
            0,-1.73009e-002,9.04761e-003;...
            0,-1.30358e-002,7.44926e-003],... % 630 Hz
           [0,-3.13292e-002,2.06370e-002;...
            0,-2.44342e-002,1.13731e-002;...
            0,-1.82108e-002,9.36778e-003],... % 800 Hz
           [0,-4.28261e-002,2.59325e-002;...
            0,-3.49619e-002,1.43046e-002;...
            0,-2.57855e-002,1.17912e-002],... % 1000 Hz
           [0,-5.91733e-002,3.25054e-002;...
            0,-5.06072e-002,1.79513e-002;...
            0,-3.69401e-002,1.48094e-002],... % 1250 Hz
           [0,-8.26348e-002,4.05894e-002;...
            0,-7.40348e-002,2.24476e-002;...
            0,-5.34977e-002,1.85371e-002],... % 1600 Hz
           [0,-1.17018e-001,5.08116e-002;...
            0,-1.09516e-001,2.81387e-002;...
            0,-7.85097e-002,2.32872e-002],... % 2000 Hz
           [0,-1.67714e-001,6.37872e-002;...
            0,-1.63378e-001,3.53729e-002;...
            0,-1.16419e-001,2.93723e-002],... % 2500 Hz
           [0,-2.42528e-001,7.98576e-002;...
            0,-2.45161e-001,4.43370e-002;...
            0,-1.73972e-001,3.70015e-002],... % 3150 Hz
           [0,-3.53142e-001,9.96330e-002;...
            0,-3.69163e-001,5.53535e-002;...
            0,-2.61399e-001,4.65428e-002],... % 4000 Hz
           [0,-5.16316e-001,1.24177e-001;...
            0,-5.55473e-001,6.89403e-002;...
            0,-3.93998e-001,5.86715e-002],... % 5000 Hz
           [0,-7.56635e-001,1.55023e-001;...
            0,-8.34281e-001,8.58123e-002;...
            0,-5.94547e-001,7.43960e-002],... % 6300 Hz
           [0,-1.10165e+000,1.91713e-001;...
            0,-1.23939e+000,1.05243e-001;...
            0,-8.91666e-001,9.40354e-002],... % 8000 Hz
           [0,-1.58477e+000,2.39049e-001;...
            0,-1.80505e+000,1.28794e-001;...
            0,-1.32500e+000,1.21333e-001],... % 10000 Hz
           [0,-2.50630e+000,1.42308e-001;...
            0,-2.19464e+000,2.76470e-001;...
            0,-1.90231e+000,1.47304e-001]); % 12500 Hz

% filter gains
filtgain = [4.30764e-011;... % 25 Hz
        8.59340e-011;... % 31.5 Hz
        1.71424e-010;... % 40 Hz
        3.41944e-010;... % 50 Hz
        6.82035e-010;... % 63 Hz
        1.36026e-009;... % 80 Hz
        2.71261e-009;... % 100 Hz
        5.40870e-009;... % 125 Hz
        1.07826e-008;... % 160 Hz
        2.14910e-008;... % 200 Hz
        4.28228e-008;... % 250 Hz
        8.54316e-008;... % 315 Hz
        1.70009e-007;... % 400 Hz
        3.38215e-007;... % 500 Hz
        6.71990e-007;... % 630 Hz
        1.33531e-006;... % 800 Hz
        2.65172e-006;... % 1000 Hz
        5.25477e-006;... % 1250 Hz
        1.03780e-005;... % 1600 Hz
        2.04870e-005;... % 2000 Hz
        4.05198e-005;... % 2500 Hz
        7.97914e-005;... % 3150 Hz
        1.56511e-004;... % 4000 Hz
        3.04954e-004;... % 5000 Hz
        5.99157e-004;... % 6300 Hz
        1.16544e-003;... % 8000 Hz
        2.27488e-003;... % 10000 Hz
        3.91006e-003];   % 12500 Hz

%%  Calculate centre frequencies (25 Hz - 12.5 kHz)

N_bands = 28; % number of bands
CenterFrequency = zeros(N_bands,1);

for i = 1:N_bands
    CenterFrequency(i) = 10^(((i-1)-16)/10.) * 1000; % calculate centre frequencies
end

% limit freq range if required
if bLimit_range
    % find idx of fmin and fmax and fix freq range before filtering
    % fmin and fmax are only used if they are specified as inputs:
    idx_fmin = find(CenterFrequency>=fmin,1);
    idx_fmax = find(CenterFrequency>=fmax,1);
    filtgain = filtgain(idx_fmin:idx_fmax); % adjust filtgain to fmin and fmax 
    ad = ad(:,:,idx_fmin:idx_fmax); % adjust ad to fmin and fmax 	
end

%% Filter signal

N_bands = size(filtgain,1);
outsig = zeros(len,N_bands);

for n = 1:N_bands
    % Three cascaded filters:
    outsig(:,n) = filtgain(n) * filter(br(3,:),ar(3,:)-ad(3,:,n),...
        filter(br(2,:),ar(2,:)-ad(2,:,n),...
        filter(br(1,:),ar(1,:)-ad(1,:,n),insig)));
end

fc = CenterFrequency;

if bLimit_range
    fc = CenterFrequency(idx_fmin:idx_fmax);
end