% Script run_validation_FS_fmod
%  This routine plot the validation of fluctuation strength (Osses et al. model )
%
%  SIGNALS: Am tones, fc=1 kHz, m=100%, SPL=70 dB, fmod=[1 2 4 8 16 32]
%
%  reference values taken from : 
%  Osses, Alejandro, Rodrigo García, and Armin Kohlrausch (2016). "Modelling 
%      the sensation of fluctuation strength." Proceedings of Meetings on 
%      Acoustics. Vol. 28. 050005. doi: 10.1121/2.0000410
%
%  Author: Gil Felix Greco, Braunschweig, 02/03/2020 (updated in 13.03.2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;

save_figs=0; %% save figs flag

%% reference data from Osses et al 2019

ref=[1    2    4    8    16    32    ; % fmod (Hz)
     0.39 0.84 1.25 1.30  0.36  0.06]; % FS values (vacil)

%% compute FS from signals

levelOut = 70;  % target signal's SPL
j=1; % init counter

res=cell(1,size(ref,2));  % declaring for memory allocation

tic
for i=1:size(ref,2)
    
    %path='SQAT_open_source\sound_files\validation\fluctuation_strength_Ossesetal2016' % path of the sound files for reference
    [insig,fs]=audioread(['AM-tone-fc-1000_fmod-' sprintf('%.0f',ref(1,j)) '_mdept-100-SPL-70-dB.wav']);
    
    levelIn = 20*log10(rms(insig)/2e-5); % SPL of the signal
    insig = insig * 10^((levelOut-levelIn)/20); % correct rms SPL to desired levelOut
    SPL(i)=20.*log10(rms(insig)/2e-5); % verify final SPL of the signal
    
    res{i} = FluctuationStrength_Osses2016(insig,fs,...  % input signal and sampling freq.
                                                       0,...  % method, stationary analysis =0 - window size=length(insig), time_varying analysis - window size=2s
                                                       0,...  % time_skip, in seconds for statistical calculations
                                                       0);    % show results, 'false' (disable, default value) or 'true' (enable)
    
    results(i)=res{1,i}.FSmean; % store mean FS value in vector results[nfmod,nFc] 
    
    j=j+1;
    
end
t=toc/60; % time to run FS calculation in minutes

%% plot results

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% plot reference curves
semilogx(ref(1,:),ref(2,:),'k*-','MarkerSize',8);hold all;

% plot computed results
semilogx(ref(1,:),results,'ko:','MarkerSize',8);hold all;

legend('Reference','SQAT','Location','NorthEast','Interpreter','Latex');
legend boxoff

axis([0 32 0 1.4]);

ax = gca;
set(ax,'XTick',[0 1 2 4 8 16 32]);
set(ax,'YTick',[0 0.2 0.4 0.6 .8 1 1.2 1.4]);
ax.XAxis.MinorTick = 'off';
ax.XAxis.MinorTickValues =  0:10:160;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0:0.1:1.4;
     
ylabel('Fluctuation strength, FS (vacil)','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
figuresdir = 'figs\'; 
saveas(gcf,strcat(figuresdir, 'validation_FS_fmod_1k'), 'fig');
saveas(gcf,strcat(figuresdir, 'validation_FS_fmod_1k'), 'pdf');
saveas(gcf,strcat(figuresdir, 'validation_FS_fmod_1k'), 'png');
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
