clc;clear all;close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% Example: compute tonality of a pure tone signal using the Aures Tonality metric
%
% FUNCTION:
%   OUT = Tonality_Aures1985(insig,fs)
%   type <help Tonality_Aures1985> for more info
%
% test signal: a pure tone with frequency 1000 Hz and level 60 dB has a tonality of 1 t.u. 
%
% Gil Felix Greco, Braunschweig 20.03.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_figs=0;
%% load .wav RefSignal 

% path='sound_files\reference_signals\Tonality_Aures1985'; %  % path of the sound file for reference
[RefSignal,fs]=audioread('RefSignal_Tonality_Aures1985_1kHz_60dBSPL_44100hz_64bit.wav');

time_insig=(0 : length(RefSignal)-1) ./ fs;  % time vector of the audio input, in seconds

%% Aures' tonality

OUT = Tonality_Aures1985(RefSignal,fs,...  % input signal and sampling freq.
                                    0,...  % field for loudness calculation; free field = 0; diffuse field = 1;
                                    0,...  % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                    1);... % show results, 'false' (disable, default value) or 'true' (enable)
   
display(sprintf('\nTonality (Aures model): \ncalculation of reference signal (60 dBSPL 1 kHz tone)\nyields a time-averaged tonality value of %g (t.u.).\n',OUT.Kmean));
                       
%% Plot 1 t.u. and input signal

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% pressure 
yyaxis left
plot(time_insig,RefSignal);
axis([0 .1 -0.1 0.1]);
ylabel('Acoustic pressure, $p$ (Pa)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 

% roughness
yyaxis right
plot(OUT.time,OUT.InstantaneousTonality); hold on;
ylabel('Aures tonality, $K$ (t.u.)','Interpreter','Latex');

a=plot(OUT.time,OUT.Kmean*ones(length(OUT.time)),'k--','Linewidth',.5);
legend(a,{sprintf('$\\mathrm{K}_{\\mathrm{mean}}=%.5g$ (t.u.)',OUT.Kmean)},'Location','NorthEast','Interpreter','Latex'); %legend boxoff 

axis([0 .1 0 1.3]);

set(gcf,'color','w');

if save_figs==1
figuresdir = 'figs\'; 
saveas(gcf,strcat(figuresdir, 'aures_tonality_reference_signal'), 'fig');
saveas(gcf,strcat(figuresdir, 'aures_tonality_reference_signal'), 'pdf');
saveas(gcf,strcat(figuresdir, 'aures_tonality_reference_signal'), 'png');
else
end
