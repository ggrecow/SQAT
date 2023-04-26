% Script ex_sound_level_meter
%
% Example: compute loudness (ISO 532-1) of stationary and time-varying inputs
%
% FUNCTION:
%   OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
%   type <help Loudness_ISO532_1> for more info
%
% test signal: narrow band noise with a center frequency of 1kHz,
%              an overall level of 40 dB should yield a loudness value of 1 sone
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;

dir_sounds = get_dir_reference_sounds('Sharpness_DIN45692');
wavfilename = '1KHZ60DB.WAV';
dBFS = 90; 

fname = [dir_sounds wavfilename];
[insig, fs] = audioread(fname);

%%% Checking the average level based on the calibration:
rms_val = 20*log10(rms(insig))+dBFS;

%%% Running the sound level meter using A-weighting curve:
weight_freq = 'A'; % A-frequency weighting
weight_time = 'f'; % Time weighting

[lvl_dBA, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS);
t = (1:length(lvl_dBA))/fs;

LAeq = Get_Leq(lvl_dBA,fs); % Make sure you enter only mono signals
T = length(lvl_dBA)/fs;

LAFmax = max(lvl_dBA);
SEL  = LAeq + 10*log10(T);

%%% Running the sound level meter again, now only changing the weighting to
%     'Z' (no weighting)
weight_freq = 'Z';
[lvl_dBZ, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS);
LZeq = Get_Leq(lvl_dBZ,fs); % Make sure you enter only mono signals

LZFmax = max(lvl_dBZ);
SEL_Z  = LAeq + 10*log10(T);

figure;
plot(t,lvl_dBA,'b-'); hold on;
plot(t,lvl_dBZ,'r--');
xlabel('Time (s)');
ylabel('Sound pressure level (dB)');

legend({'dB(A)','dB(Z)'});

fprintf('\n%s.m: The level information for sound in %s, analysed using a %s time weighting\n',mfilename,wavfilename,weight_time);
fprintf('\tAverage (RMS) level using dBFS=%.0f: lvl_rms=%.1f dB(Z)\n',dBFS,rms_val);
fprintf('\tDuration of the sound: dur=%.1f s\n\n',T);

fprintf('\tA-weighted maximum sound pressure level Lmax=%.1f dB(A)\n',LAFmax);
fprintf('\tA-weighted equivalent sound pressure level Leq=%.1f dB(A)\n',LAeq);
fprintf('\tA-weighted sound exposure level SEL=%.1f dB(A)\n\n',SEL);

fprintf('\tZ-weighted maximum sound pressure level Lmax=%.1f dB(A)\n',LZFmax);
fprintf('\tZ-weighted equivalent sound pressure level Leq=%.1f dB(A)\n',LZeq);
fprintf('\tZ-weighted sound exposure level SEL=%.1f dB(A)\n',SEL_Z);
