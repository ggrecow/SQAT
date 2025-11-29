% Script ex_Do_OB13_ISO532_1
%
% This script produces the impulse response of the one-third octave band
% module included in the sound level meter of the SQAT toolbox.
%
% IMPORTANT INFO: the 1/3 octave band filter bank implemented complies with
%   ISO 532-1:2017 (and IEC 61260-1:2014). However, the possible freq range computed
%   is only between 25 Hz - 12.5 kHz.
%
% The calculation is done for a unit impulse (amplitude 1) of 1 second of
% duration at a sampling frequency of 48000 Hz, because the filter bank is
% hard coded to work on this sampling frequency.
% If you want to use as an input any other sound, just obtain your insig
% and fs using audioread, loading any wavfilename available on your disk, i.e.,
% use [insig,fs] = audioread(wavfilename).
%
% Author: Alejandro Osses
% Date: October 2023
% Date: 29/11/2025 (Small visualisation improvements)
% MATLAB / GNU Octave compatible: Yes / Yes
% See also: Do_OB13_ISO532_1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all ;close all;

fs = 48000;

% Impulse, to obtain an impulse response:
insig = [zeros(fs/2-1,1); 1; zeros(fs/2,1)]; dBFS = 100;
YL = [-12 2]; % expected maximum = 0 dB

%%% Checking the average level based on the calibration:
rms_val = 20*log10(rms(insig))+dBFS;

[OB_filt,fc] = Do_OB13_ISO532_1(insig,fs);

hfreq = [];
K = 8192; % K = N/2 - arbitrary number of points for the FFT
for i = 1:length(fc)
    BW(i) = fc(i).*(2^(1/6)-2^(-1/6)); % just extra information, the theoretical bandwidth of each 1/3 octave band

    [hfreq_here,w] = freqz(OB_filt(:,i),1,K);
    if i == 1
        f = (w/pi)*fs/2; % obtaining the frequency in Hz, only for the first fc
    end
    hfreq(:,end+1) = 20*log10(abs(hfreq_here)); % appending the frequency response (dB) in the last column
end

% Setting the first frequency bin to a non-null value:
f(1) = f(2)*.5;

figure;
semilogx(f,hfreq,'k'); grid on;
set(gca,'XTick',fc);
set(gca,'XTickLabel',round(fc));
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
XL = [min(fc) max(fc)];
xlim(XL);
ylim(YL);

hold on; plot(XL,-3*[1 1],'r--'); % theoretical amplitude where the filters should overlap

% Sub example: obtaining the calibrated RMS level of each band:
OB_lvls = 20*log10(rms(OB_filt))+dBFS;

flow = fc.*2^(-1/6);
fhigh = fc.*2^(1/6);
freqs = [flow fc fhigh]; % just if you want to see a list of low, centre, and high frequencies for each filter

% Total level of the input signal:
lvl_orig = 20*log10(rms(insig))+dBFS;
fprintf('Level as obtained directly from the input signal=%.1f dB SPL\n',lvl_orig);

% The sum of the power of each band leads to the total level:
lvl_from_outsig = 10*log10(sum(10.^(OB_lvls/10))); % '/10' means squared, therefore it is 10*log10()
fprintf('Level obtained from the filtered signals=%.1f dB SPL\n',lvl_from_outsig);

figure;
semilogx(fc,OB_lvls,'bo-'); hold on; grid on
set(gca,'XLim',[min(fc)-1 max(fc)+1000]);
set(gca,'XTick',fc);
set(gca,'XTickLabel',round(fc));
xlabel('Frequency (Hz)');
ylabel('Band level (dB SPL)');
