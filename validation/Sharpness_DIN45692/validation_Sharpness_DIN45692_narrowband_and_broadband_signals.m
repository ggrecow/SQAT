% Script validation_Sharpness_DIN45692_narrowband_and_broadband_signals
%
% - Validation of sharpness calculation of narrowband and broadband using 
%   test signals from DIN 45692(2009)
%   type <help Sharpness_DIN45692_from_loudness> for more info
% 
% - Specific Loudness is computed within this code (SQAT) according to 
%   ISO 532:1-2017 (type <help Loudness_ISO532_1> for more info)
%   stationary, free-field, 
%
%  signals are calibrated with a ref. tone of 1kHz@60dB
%
% Author: Gil Felix Greco, Braunschweig 01.03.2023 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

save_figs=0; % save figure flag

%% Sharpness of narrowband test sounds

% Tabelle A.2 — Sollwerte der Schärfe für frequenzgruppenbreites Schmalbandrauschen mit 
% der Lautheit 4 sone (entsprechend dem Standard-Prüfsignal bei 1 kHz und 60 dB)

% Reference values from the standard DIN45692 for the narrowband test sounds
freq_narrowband = [250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400 4000 4800 5800 7000 8500 10500];
bark_narrowband = 2.5:22.5;
acum_narrowband = [.38 .49 .60 .71 .82 .93 1 1.13 1.26 1.35 1.49 1.64 1.78 2.06 2.4 2.82 3.48 4.43 5.52 6.81 8.55];

dir_sounds     = get_dir_validation_sounds('Sharpness_DIN45692');
dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep ];

% load calibration signal 

% path='sound_files\reference_signals\sharpness_DIN45692\'; % path of the sound file for reference
[RefSignal,~]=audioread([dir_ref_sounds 'RefSignal_Sharpness_DIN45692.wav']);

%%% Calibration using the concept of dB full scale:
lvl_cal_signal = 60;                               % from file name: 'calibration signal sine 1kHz 60dB.wav'
    
[~,~,dBFS]=calibrate(RefSignal,RefSignal,lvl_cal_signal); % just to obtain the dBFS value

% load narrowband signals, calibrate, compute loudness and sharpness 
for i=1:length(freq_narrowband)

    % load test signal 
    fname_insig = ['narrowband_' sprintf( '%02d',freq_narrowband(i) ) '.wav'];
    [insig,fs]=audioread([dir_sounds fname_insig]); % 'sound_files\validation\Sharpness_DIN45692\'; % path of the sound file for reference

    % calibrate .wav signal
    [insig_cal]=calibrate(insig,RefSignal,lvl_cal_signal); 
    
    % stationary loudness calculation
    L_narrow = Loudness_ISO532_1( insig_cal, fs,...   % input signal and sampling freq.
                                  0,...   % field; free field = 0; diffuse field = 1;
                                  1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                  0.5,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                  0);     % show results, 'false' (disable, default value) or 'true' (enable)
                                                                       
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

% Tabelle A.3: Sollwerte der Schärfe für Breitbandrauschen mit fester oberer Eckfrequenz 
% von 10 kHz und variabler unterer Eckfrequenz fu mit der Lautheit 4 sone (entsprechend dem Standard-Prüfsignal bei 1 kHz und 60 dB)

% Reference values from the standard DIN45692 for the broadband test sounds
freq_broadband= [250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400 4000 4800 5800 7000 8500];
bark_broadband=2.5:21.5;
acum_broadband=[2.7 2.74 2.78 2.85 2.91 2.96 3.05 3.12 3.20 3.3 3.42 3.53 3.69 3.89 4.12 4.49 5.04 5.69 6.47 7.46];

% load broadband signals, calibrate, compute loudness and sharpness 
for i=1:length(freq_broadband)

    % load test signal 
    [insig,fs]=audioread([dir_sounds 'broadband_' sprintf( '%02d',freq_broadband(i) ) '.wav']);  'sound_files\validation\Sharpness_DIN45692\'; % path of the sound file for reference

    % calibrated .wav signal
    insig_cal=calibrate(insig,RefSignal,lvl_cal_signal); 
    
    % stationary loudness calculation
    L_broad = Loudness_ISO532_1( insig_cal, fs,...   % input signal and sampling freq.
                                  0,...   % field; free field = 0; diffuse field = 1;
                                  1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                                  0.5,... % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                  0);     % show results, 'false' (disable, default value) or 'true' (enable)
                                                    
    % sharpness calculation (standard weighting function)                                         
    s_broad(i) = Sharpness_DIN45692_from_loudness(L_broad.SpecificLoudness,... % input (stationary) specific loudness
                                                               'DIN45692');... % type of weighting function used for sharpness calculation
       

    % sharpness calculation (aures weighting function)                                         
    s_broad_aures(i) = Sharpness_DIN45692_from_loudness(L_broad.SpecificLoudness,... % input (stationary) specific loudness
                                                                        'aures');... % type of weighting function used for sharpness calculation
                                                  
    % sharpness calculation (von bismarck weighting function)                                         
    s_broad_bismarck(i) = Sharpness_DIN45692_from_loudness(L_broad.SpecificLoudness,... % input (stationary) specific loudness
                                                                        'bismarck');... % type of weighting function used for sharpness calculation
                                                  
end

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
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
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
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
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

title('DIN 45692 - Sharpness of narrowband noise');

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
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
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

title('DIN 45692 - Sharpness of broadband noise');

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
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

