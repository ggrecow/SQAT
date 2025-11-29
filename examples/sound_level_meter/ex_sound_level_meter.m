% Script ex_sound_level_meter
%
% It analyses the reference sound for sharpness:
% Narrowband noise, fc = 1 kHz, BW = 160 Hz and Lp = 60 dBSPL
%
% using the sound_level_meter module.
% This example is meant to show how to use the functions Do_SLM and Get_Leq,
% and the selected sound file can be replaced by any other wavfilename
% available on disk.
%
% Author: Alejandro Osses
% Date: 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all ;close all;

dir_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
wavfilename='RefSignal_Sharpness_DIN45692.wav';
[insig,fs]=audioread([dir_sounds 'RefSignal_Sharpness_DIN45692.wav']); % 'sound_files\reference_signals\' -  path of the sound file for reference

dBFS = 90;

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
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('Sound pressure level, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex');

legend({'dB(A)','dB(Z)'});

set(gcf,'color','w');

fprintf('\n%s.m: The level information for sound in %s, analysed using a %s time weighting\n',mfilename,wavfilename,weight_time);
fprintf('\tAverage (RMS) level using dBFS=%.0f: lvl_rms=%.1f dB(Z)\n',dBFS,rms_val);
fprintf('\tDuration of the sound: dur=%.1f s\n\n',T);

fprintf('\tA-weighted maximum sound pressure level Lmax=%.1f dB(A)\n',LAFmax);
fprintf('\tA-weighted equivalent sound pressure level Leq=%.1f dB(A)\n',LAeq);
fprintf('\tA-weighted sound exposure level SEL=%.1f dB(A)\n\n',SEL);

fprintf('\tZ-weighted maximum sound pressure level Lmax=%.1f dB(Z)\n',LZFmax);
fprintf('\tZ-weighted equivalent sound pressure level Leq=%.1f dB(Z)\n',LZeq);
fprintf('\tZ-weighted sound exposure level SEL=%.1f dB(Z)\n',SEL_Z);
