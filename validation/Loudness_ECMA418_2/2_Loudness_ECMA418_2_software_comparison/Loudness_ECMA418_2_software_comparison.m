% Script Loudness_ECMA418_2_software_comparison
%
% Compute Loudness (ECMA 418-2:2024) and compare with results from
% commercial software 
% 
% Signal: binaural audio recording of a 'train station' environment (30 seconds, 2-channel binaural)
% The signal 'TrainStation.7.wav' was extracted from the EigenScape database 
% (https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between
%  01m00s and 01m30s. The EigenScape database, which is described by 
% Green et al (https://doi.org/10.3390/app7111204), is licensed 
% under Creative Commons Attribution 4.0.
%
% FUNCTION:
%   OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Loudness_ECMA418_2> for more info
%
% Author: Gil Felix Greco, Braunschweig 12.02.2025
% Modified: 19.03.2025 Mike Lotinga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Load  stereo .wav file
 
dir_sound = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

fileTag = 'ExStereo_';
wav_file = 'TrainStation7-0100-0130';

[insig, fs] = audioread([dir_sound fileTag wav_file '.wav']); 

%% Compute loudness (stereo signal)

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3; % time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 1; % show results, 'false' (disable, default value) or 'true' (enable)

OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show);

%% load reference results 

% referece results - main path
ref_path = [pwd '\reference_results\'];

% reference results - file name
AvgSpec_fileName = [ wav_file '.Specific Loudness (Hearing Model).asc']; % channel 1 and 2
AvgSpecCombBinaural_fileName = [ wav_file '.Specific Loudness (Hearing Model)_combined_binaural.asc']; % combined binaural

TDep_fileName = [ wav_file '.Loudness (Hearing Model) vs. Time.asc']; % channel 1 and 2
TDepCombBinaural_fileName = [ wav_file '.Loudness (Hearing Model) vs. Time_combined_binaural.asc']; % combined binaural

Spec_TDep_fileName_channel_1 = [ wav_file '.Specific Loudness (Hearing Model) vs. Time_channel_1.asc']; % channel 1 
Spec_TDep_fileName_channel_2 = [ wav_file '.Specific Loudness (Hearing Model) vs. Time_channel_2.asc']; % channel 2  
Spec_TDep_fileName_combined_binaural = [ wav_file '.Specific Loudness (Hearing Model) vs. Time_combined_binaural.asc']; % combined binaural

 % reference results - load files
ref_results.AvgSpec = readmatrix( [ref_path AvgSpec_fileName] , 'FileType', 'text');
ref_results.AvgSpecCombBinaural = readmatrix( [ref_path AvgSpecCombBinaural_fileName] , 'FileType', 'text');

ref_results.TDep = readmatrix( [ref_path TDep_fileName] , 'FileType', 'text');
ref_results.TDepCombBinaural = readmatrix( [ref_path TDepCombBinaural_fileName] , 'FileType', 'text');

ref_results.Spec_TDep_channel_1 = readmatrix( [ref_path Spec_TDep_fileName_channel_1] , 'FileType', 'text');
ref_results.Spec_TDep_channel_2 = readmatrix( [ref_path Spec_TDep_fileName_channel_2] , 'FileType', 'text');
ref_results.Spec_TDep_combined_binaural = readmatrix( [ref_path Spec_TDep_fileName_combined_binaural] , 'FileType', 'text');

% single values
% Analysis Name    Channel Name    N/soneHMS    AnalysisRange
% Loudness (Hearing Model) vs. Time    Ch1    8,10    [0,304 - 29,9 s]
% Loudness (Hearing Model) vs. Time    Ch2    6,43    [0,304 - 29,9 s]
% Loudness (Hearing Model) vs. Time    CombinedBinaural Ch1&Ch2    7,32    [0,304 - 29,9 s]

ref_results.single_values = [8.10 6.43 7.32]; % [Ch1 Ch2 combBinaural]

%% plot - time-varying loudness 

% plot -  Channel 1
xRef = ref_results.TDep(:,1);
yRef = ref_results.TDep(:,2);
xSQAT = OUT.timeOut;
ySQAT = OUT.loudnessTDep(:,1);
label_title = [wav_file '.wav (Channel 1)'];
label_fig = [wav_file ' (Channel 1)' '_TDep_Loudness'];

il_plt_tDep(xRef, yRef, xSQAT, ySQAT, label_title, save_figs, label_fig);

% plot - Channel 2

yRef = ref_results.TDep(:,3);
ySQAT = OUT.loudnessTDep(:,2);
label_title = [wav_file '.wav (Channel 2)'];
label_fig = [wav_file ' (Channel 2)' '_TDep_Loudness'];

il_plt_tDep(xRef, yRef, xSQAT, ySQAT, label_title, save_figs, label_fig);

% plot - combined binaural
yRef = ref_results.TDepCombBinaural(:,2);
ySQAT = OUT.loudnessTDepBin(:,1);
label_title = [wav_file '.wav (Combined binaural)'];
label_fig = [wav_file ' (Combined binaural)' '_TDep_Loudness'];

il_plt_tDep(xRef, yRef, xSQAT, ySQAT, label_title, save_figs, label_fig);

%% plot - avg specific loudness  

% plot -  Channel 1
yRef = ref_results.AvgSpec(:,2);
ySQAT = OUT.specLoudnessPowAvg(:,1);
label_title = [wav_file '.wav (Channel 1)'];
label_fig = [wav_file ' (Channel 1)' '_avgSpecific_Loudness'];

il_plt_avgSpecific( yRef, ySQAT, label_title, save_figs, label_fig)

% plot -  Channel 2
yRef = ref_results.AvgSpec(:,3);
ySQAT = OUT.specLoudnessPowAvg(:,2);
label_title = [wav_file '.wav (Channel 2)'];
label_fig = [wav_file ' (Channel 2)' '_avgSpecific_Loudness'];

il_plt_avgSpecific( yRef, ySQAT, label_title, save_figs, label_fig)

% plot -  combined binaural
yRef = ref_results.AvgSpecCombBinaural(:,2);
ySQAT = OUT.specLoudnessPowAvgBin(:,1);
label_title = [wav_file '.wav (Combined binaural)'];
label_fig = [wav_file ' (Combined binaural)' '_avgSpecific_Loudness'];

il_plt_avgSpecific( yRef, ySQAT, label_title, save_figs, label_fig)

%% plot - overall loudness

ch1_reference = ref_results.single_values(1);
ch2_reference = ref_results.single_values(2);
combBinaural_reference = ref_results.single_values(3);
ch1_SQAT = OUT.loudnessPowAvg(1);
ch2_SQAT = OUT.loudnessPowAvg(2);
combBinaural_SQAT = OUT.loudnessPowAvgBin(1);

single_values =[ch1_reference ...
                          ch1_SQAT; ...
                          ch2_reference ...
                          ch2_SQAT; ...
                          combBinaural_reference ...
                          combBinaural_SQAT ];

label_title = [wav_file '.wav'];
label_fig = [wav_file '_singleValues_Loudness'];

il_plt_singleValues(single_values, label_title, save_figs, label_fig)

%% plot - time-dependent specific loudness

% plot -  commercial software (Channel 1)
xAxis =  ref_results.Spec_TDep_channel_1(2:end,1) ; % time vector
yAxis =  ref_results.Spec_TDep_channel_1(1, 2:end) ; % freq vector
zAxis =  ref_results.Spec_TDep_channel_1(2:end, 2:end) ; % specific loudness
label_fig = [wav_file ' (Channel 1)' '_tDep_Specific_Loudness_ref'];

il_plt_spectrogram(xAxis, yAxis, zAxis', save_figs, label_fig)

% plot -  commercial software (Channel 2)
xAxis =  ref_results.Spec_TDep_channel_2(2:end,1) ; % time vector
yAxis =  ref_results.Spec_TDep_channel_2(1, 2:end) ; % freq vector
zAxis =  ref_results.Spec_TDep_channel_2(2:end, 2:end) ; % specific loudness
label_fig = [wav_file ' (Channel 2)' '_tDep_Specific_Loudness_ref'];

il_plt_spectrogram(xAxis, yAxis, zAxis', save_figs, label_fig)

% plot -  implementation (Channel 1)
xAxis =  OUT.timeOut; % time vector
yAxis =  OUT.bandCentreFreqs; % freq vector
zAxis =  OUT.specLoudness(:,:,1) ; % specific loudness
label_fig = [wav_file ' (Channel 1)' '_tDep_Specific_Loudness_implementation'];

il_plt_spectrogram(xAxis, yAxis, zAxis', save_figs, label_fig)

% plot -  implementation (Channel 2)
% xAxis =  OUT.timeOut; % time vector
% yAxis =  OUT.bandCentreFreqs; % freq vector
zAxis =  OUT.specLoudness(:,:,2) ; % specific loudness
label_fig = [wav_file ' (Channel 2)' '_tDep_Specific_Loudness_implementation'];

il_plt_spectrogram(xAxis, yAxis, zAxis', save_figs, label_fig)

%% function / plot - time-varying quantity

function il_plt_tDep(xRef, yRef, xSQAT, ySQAT, label_title, save_figs, label_fig)

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

cmap  = load('cmap_viridis.txt');
cmap1 = 166;
cmap2 = 34;

ax = nexttile;
% reference values 
plot(xRef, yRef, 'Color', cmap(cmap1, :), 'Linewidth',1); hold all;

% SQAT values 
plot(xSQAT, ySQAT, ':', 'Color', cmap(cmap2, :), 'Linewidth', 1.5);  

xlabel('Time (s)');
ylabel('Loudness (sone_{HMS})');

legend('Reference', 'Implementation', 'Location', 'NE');
set(gcf,'color','w');

ylim([0 12]);

% title(label_title);

ax.FontName = 'Times';
ax.FontSize = 16;
% ax.Title.FontWeight = 'normal';
% ax.Title.FontSize = 16;

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];
    
    resolution = '-r600'; % Default resolution of 600 DPI

    % saveas(gcf,figname_out, 'fig');
    % print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end

%% function / plot - avg specific quantity

function il_plt_avgSpecific( yRef, ySQAT, label_title, save_figs, label_fig)

cmap  = load('cmap_viridis.txt');
cmap1 = 166;
cmap2 = 34;

h = figure('Position', [200, 200, 1500, 550]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
 
barkAxis = linspace(0.5, 26.5, 53);

left_color = [0 0 0];
right_color = [0 0 0];
set(h,'defaultAxesColorOrder',[left_color; right_color]);

ax = nexttile;

yyaxis left;

a = bar(barkAxis,...
    yRef,...
    'EdgeColor', cmap(cmap1, :),...
    'EdgeAlpha', 0.2, ...
    'FaceColor', cmap(cmap1, :),...
    'FaceAlpha', 0.2, ...
    'LineWidth', 1, 'LineStyle', '-',...
    'BarWidth',0.5, ...
    'DisplayName', "Reference");

hold on;

b = bar(barkAxis,...
    ySQAT,...
    'EdgeColor', cmap(cmap2, :),...
    'FaceColor',  'none',...
    'LineWidth', 0.5, 'LineStyle', '-',...
    'BarWidth',0.5, ...
    'DisplayName', "Implementation");

xlabel( 'Critical band rate (Bark_{HMS})' );
ylabel( 'Specific loudness (sone_{HMS}/Bark_{HMS})' );

xtickangle(90);
xticks(barkAxis);

ylim([0 0.6]);

% title(label_title);

ax.FontName = 'Times';
ax.FontSize = 16;
% ax.Title.FontWeight = 'normal';
% ax.Title.FontSize = 16;

yyaxis right;
c = plot( barkAxis, ySQAT-yRef, 'k*-');

ylabel( 'Implementation-Reference (sone_{HMS}/Bark_{HMS})' );

ylim([-0.002 0.002]);

legend([a b c], {'Reference', 'Implementation', 'Implementation - Reference'});

set(h,'color','w');

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];
    
    resolution = '-r600'; % Default resolution of 600 DPI

    % saveas(gcf,figname_out, 'fig');
    % print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end

%% function / single values 

function il_plt_singleValues(single_values, label_title, save_figs, label_fig)

% Single values
% h  =figure('Position', [200, 200, 550, 550]);
h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
movegui(h, 'center');

str_xlabel = {'Ch1';'Ch2';'Combined binaural'};

br=bar(str_xlabel, single_values);
ylabel("Loudness (sone_{HMS})")

cMap  = load('cmap_viridis.txt');
colororder([cMap(166, :); cMap(34, :)])

for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, string(round(br(bb).YData, 2)),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 12)
end

ylim([0 round(max(single_values(1))+5)]);

legend('Reference', 'Implementation', 'Location', 'NE');

% t=title(label_title);

% t.FontWeight = 'normal';
% t.FontSize = 16;
ax = gca;
ax.FontName = 'Times';
ax.FontSize = 16;

set(gcf,'color','w');

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];

    resolution = '-r600'; % Default resolution of 600 DPI

    % saveas(gcf,figname_out, 'fig');
    % print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end

%% function / plot spectrogram

function il_plt_spectrogram(xAxis, yAxis, zAxis, save_figs, label_fig)

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
movegui(h, 'center');


[xx,yy] = meshgrid(xAxis, yAxis);
xx(1) = 0;  % truncating first idx of time to zero (just a visual thing because ref results may not start exactly on zero)
pcolor(xx, yy, zAxis);
shading interp; colorbar; axis tight;
cMap  = load('cmap_viridis.txt'); colormap(cMap);

ax = gca;
ax.FontName = 'Times';
ax.FontSize = 16;
ax.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3];
ax.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k","8k", "16k"];
ax.YScale = 'log';
ax.YLabel.String = 'Frequency (Hz)';
ax.XLabel.String = 'Time (s)';

 h = colorbar;
 zString = 'Specific loudness (sone_{HMS}/Bark_{HMS})';
 set(get(h,'label'),'string', zString);

 box on
 set(gca,'layer','top');

 set(gcf,'color','w');

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];

    resolution = '-r600'; % Default resolution of 600 DPI

    % saveas(gcf,figname_out, 'fig');
    % print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end