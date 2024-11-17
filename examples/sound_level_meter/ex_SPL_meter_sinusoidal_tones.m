%   Script ex_SPL_meter_sinusoidal_tones
%
% - create two signals: sinusoidal tones, with rms SPL level of 60 dB SPL and central
%   frequencies of 100 Hz and 1 kHz. The signals are created with silence before and after 
%   the tones in order to verify the temporal response of the fast and slow
%   time-weightings
%
% - Compute SPL of the tones using the Z- and A-weighting and the 
%   Slow- and Fast time-weightings using the sound_level_meter of SQAT
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

%% Create tones

fc = [100; 1000];      % Center frequency (Hz)
Level = 60;              % rms Level (dB SPL)
L_before = 5;          % Length of silence before the sinusoidal signal (seconds)
L_signal = 15;         % Length of the sinusoidal signal (seconds)
L_after = 10;            % Length of silence before the sinusoidal signal (seconds)
fs = 48000;              % Sampling frequency

for i = 1:size(fc,1)
    [insig(i,:), t_total] = il_create_tones( fc(i,1), Level, L_before, L_signal, L_after, fs );
end

%% compute SPL: Z-weighted

weight_freq = 'Z'; % Z-frequency weighting
weight_time = 's'; % slow leak - time constant = 1s  

for i = 1:size(insig,1)
    SPL.Z_weighted(:,i) = il_get_SPL(insig(i,:), fs, weight_freq, weight_time);
end

%% compute SPL:  A-weighted

weight_freq = 'A'; % A-frequency weighting
weight_time = 'f'; % slow leak - time constant = 1s  

for i = 1:size(insig,1)
    SPL.A_weighted(:,i) = il_get_SPL(insig(i,:), fs, weight_freq, weight_time);
end
 
%% plot insig

h  = figure;
set(h, 'name', 'Input sound pressure signals' );
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

xmax = t_total(end); % used to define the x-axis on the plots

% plot input signal - 100-Hz tone
subplot(2,1,1)
plot( t_total, insig(1,:) );  
a=yline(rms(insig(1,:)),'k--');
legend(a,sprintf('$p_{\\mathrm{rms}}=$%g Pa',rms(insig(1,:))),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('$p$ (Pa)','Interpreter','Latex'); %grid on;
ax = axis; axis([0 xmax max(insig(1,:))*-2 max(insig(1,:))*2]);
title('Input signal - 100-Hz tone','Interpreter','Latex');

% plot input signal - 1000-Hz tone
subplot(2,1,2)
plot( t_total, insig(2,:) );  
a=yline(rms(insig(2,:)),'k--');
legend(a,sprintf('$p_{\\mathrm{rms}}=$%g Pa',rms(insig(2,:))),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('$p$ (Pa)','Interpreter','Latex'); %grid on;
ax = axis; axis([0 xmax max(insig(2,:))*-2 max(insig(2,:))*2]);
title('Input signal - 1000-Hz tone','Interpreter','Latex');

set( gcf, 'color', 'w' );

%% plot SPL versus time

h  = figure;
set(h, 'name', 'SPL of sinusoidal tones' );
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

orange = [1.00,0.41,0.16];
L_window = fs*0.01; % lenght of window used to apply a moving average on the input signal for display purposes
pref = 20e-6; % reference sound pressure for air

subplot(2,1,1)
% plot(t_total, 10*log10( movmean(abs(insig(1,:)).^2, L_window)./(pref).^2) ,'b'); hold on;  % input signal - fc = 100 Hz    
plot(SPL.Z_weighted(1).t, SPL.Z_weighted(1).lvl_dB,'k'); hold on;     % fc = 100 Hz
plot(SPL.Z_weighted(2).t, SPL.Z_weighted(2).lvl_dB, '--o', 'Color', orange ,'MarkerIndices',1:50000:length(SPL.Z_weighted(2).t)); % fc = 1 kHz

xlabel( 'Time, $t$ (s)', 'Interpreter', 'Latex' );
ylabel('$L_\mathrm{Z,S}$ (dB re 20 $\mu$Pa)', 'Interpreter','Latex');

legend( sprintf('$f_{\\mathrm{c}}=%g$ Hz', fc(1)), ...
             sprintf('$f_{\\mathrm{c}}=%g$ kHz', fc(2)/1000), ...
                        'Location', 'south', 'Interpreter','Latex') ;

legend Box off

ytikz = [0 20 40 60]; % tiks of y-axis
yticks(ytikz);

% specify limit of the axis
x2 = (L_before + L_signal + L_after) ; y1 = 0; y2=Level+10;
axis([0 x2 y1 y2]);
  
title('Z-weighting and Slow time-weighting', 'Interpreter','Latex');

subplot(2,1,2)
% plot(t_total, 10*log10( movmean(abs(insig(1,:)).^2, L_window)./(pref).^2) ,'b'); hold on;  % input signal - fc = 100 Hz     
plot(SPL.A_weighted(1).t, SPL.A_weighted(1).lvl_dB, 'k'); hold on; % fc = 100 Hz
plot(SPL.A_weighted(2).t, SPL.A_weighted(2).lvl_dB , '--o', 'Color', orange ,'MarkerIndices',1:50000:length(SPL.Z_weighted(1).t)); % fc = 1 kHz

xlabel( 'Time, $t$ (s)', 'Interpreter', 'Latex' );
ylabel('$L_\mathrm{A,F}$ (dB re 20 $\mu$Pa)', 'Interpreter','Latex');

yticks(ytikz);
axis([0 x2 y1 y2]);
 
title('A-weighting and Fast time-weighting', 'Interpreter','Latex');

set( gcf, 'color', 'w' );

%% function : calculate SPL using the sound_level_meter in SQAT

function output = il_get_SPL(insig, fs, weight_freq, weight_time)
% function output = il_get_SPL(insig, fs, weight_freq, weight_time)
%
% Get SPL using the  <Do_SLM> function of SQAT
%
%   Inputs
%   insig : input signal (Pa)
%   fs : sampling frequency (Hz)
%   weight_freq : frquency-weighting
%   weight_time : frquency-weighting
%
% type <Help Do_SLM> for more info about the input format
%
%   Output
%   output : Struct containing the following results
%   t : time vector
%   lvl_dB : sound pressure level over time (vector) 
%   Lmax : maximum sound pressure level (scalar)
%   SEL : sound exposure level (scalar)
%
%   Gil Felix Greco, 15.11.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dBFS = 94;

[output.lvl_dB, ~] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS);
output.t = (1:length(output.lvl_dB))/fs;

output.Leq = Get_Leq(output.lvl_dB,fs); % Make sure you enter only mono signals
T = length(output.lvl_dB)/fs;

output.Lmax = max(output.lvl_dB);
output.SEL  = output.Leq + 10*log10(T);

end

%% function : create tones

function [y, t_total] = il_create_tones( fc, LevelIn, L_before, L_signal, L_after, fs )
% function [y, t_total] = il_create_tones( fc, LevelIn, L_before, L_signal, L_after, fs )
%
% Generate sinusoidal tone
%
%   Inputs
%   fc : central frequency (Hz)
%   LevelIn : rms level (dBSPL)
%   L_before : length of silence before the tones (s)
%   L_signal : length of the tones (s)
%   L_after : length of silence after the tones (s)
%   fs : sampling frequency (Hz)
%
%   Output
%   y : sinousidal signal (Pa)
%
%   Gil Felix Greco, 15.11.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify signal characteristics

t_before = 0:1/fs:L_before; % Time vector, before signal
t_signal = 0:1/fs:L_signal; % Time vector, during signal
t_after = 0:1/fs:L_after; % Time vector, after signal

%% Make signal

pref = 20e-6; % reference sound pressure for air

A = pref*10^(LevelIn/20)*sqrt(2); % Amplitude
y_signal = A*sin(2*pi*fc*t_signal);     % Generate signal

% SPL=20.*log10(rms(y)/2e-5); % Verify dBSPL of the generated signal

% make silence parts
y_before = zeros(1,length(t_before));
y_after = zeros(1,length(t_after));

%% make complete signal
y = [y_before y_signal y_after];
t_total = linspace(0,(L_before+L_signal+L_after),length(y));

end