% Script FluctuationStrength_ECMA418_2_software_comparison
%
% Compute fluctuation strength (ECMA-418-2:2025) and compare with results
% from commercial software
%
% Signal 1: binaural audio recording of a 'train station' environment
% (30 seconds, 2-channel binaural). The signal 'TrainStation.7.wav' was
% extracted from the EigenScape database
% (https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between
% 01m00s and 01m30s. The EigenScape database, which is described by
% Green et al (https://doi.org/10.3390/app7111204), is licensed
% under Creative Commons Attribution 4.0. This signal ships with SQAT.
%
% Signal 2: ambisonic recording of a 'park' environment with unmanned
% aircraft system (UAS / drone) flight overhead (25 seconds, 2-channel
% binaural), as described in the folder of
% pub_Lotinga2025_Forum_Acusticum_ECMA418_2. This signal is distributed
% through zenodo (DOI: 10.5281/zenodo.15132460) and has to be downloaded,
% see the README of this folder. The script reports the instructions and
% moves on to the next signal when the file is absent.
%
% FUNCTION:
%   OUT = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help FluctuationStrength_ECMA418_2> for more info
%
% Authors: Sergio Aguirre & Gil Felix Greco, 18.09.2026
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];
ref_path = [dir_out 'reference_results' filesep];

fileTag = 'ExStereo_';

% signal 1 ships with SQAT, signal 2 is downloaded from zenodo
signals(1).wav_file = 'TrainStation7-0100-0130';
signals(1).dir_sound = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
signals(1).label = 'train station';

signals(2).wav_file = 'Park3-0002-0027_UAS';
signals(2).dir_sound = [basepath_SQAT 'publications' filesep ...
    'pub_Lotinga2025_Forum_Acusticum_ECMA418_2' filesep 'data' filesep 'Audio' filesep];
signals(2).label = 'park with UAS flight overhead';

%% analysis parameters

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 700e-3; % time_skip, in seconds for statistical calculations (default: 0.7 seconds)
show = 0; % show results, 'false' (disable, default value) or 'true' (enable)

%% loop over the signals

for i = 1:length(signals)

    wav_path = [signals(i).dir_sound fileTag signals(i).wav_file '.wav'];

    if ~exist(wav_path, 'file')
        fprintf('\n%s.m: the signal %s%s.wav was not found in\n\t%s\n', ...
            mfilename, fileTag, signals(i).wav_file, signals(i).dir_sound);
        fprintf('\tThis signal is distributed through zenodo (DOI: 10.5281/zenodo.15132460).\n');
        fprintf('\tDownload the pub_Lotinga2025_Forum_Acusticum_ECMA418_2 file, unzip it and run\n');
        fprintf('\tthe provided copy_data_for_pub_Lotinga2025_to_SQAT.m script, which places the\n');
        fprintf('\tfiles in the correct folders. Skipping this signal.\n');
        continue;
    end

    fprintf('\n%s.m: %s (%s)\n', mfilename, signals(i).wav_file, signals(i).label);

    %% compute fluctuation strength (stereo signal)

    [insig, fs] = audioread(wav_path);

    OUT = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip, show);

    %% load reference results

    ref_results = il_load_reference_results(ref_path, signals(i).wav_file);

    %% plot - time-dependent fluctuation strength
    % the lower tile of each figure carries the absolute difference between
    % the implementation and the reference, following the figures of
    % Lotinga, Torjussen and Felix Greco (2025), Forum Acusticum

    il_plt_tDep(ref_results.TDep(:,1), ref_results.TDep(:,2), ...
        OUT.timeOut, OUT.fluctStrengthTDep(:,1), ...
        [signals(i).wav_file ' (Channel 1)_TDep_FluctuationStrength'], save_figs);

    il_plt_tDep(ref_results.TDep(:,1), ref_results.TDep(:,3), ...
        OUT.timeOut, OUT.fluctStrengthTDep(:,2), ...
        [signals(i).wav_file ' (Channel 2)_TDep_FluctuationStrength'], save_figs);

    il_plt_tDep(ref_results.TDepCombBinaural(:,1), ref_results.TDepCombBinaural(:,2), ...
        OUT.timeOut, OUT.fluctStrengthTDepBin, ...
        [signals(i).wav_file ' (Combined binaural)_TDep_FluctuationStrength'], save_figs);

    %% plot - time-averaged specific fluctuation strength

    il_plt_avgSpecific(ref_results.AvgSpec(:,2), OUT.specFluctStrengthAvg(:,1), ...
        [signals(i).wav_file ' (Channel 1)_avgSpecific_FluctuationStrength'], save_figs);

    il_plt_avgSpecific(ref_results.AvgSpec(:,3), OUT.specFluctStrengthAvg(:,2), ...
        [signals(i).wav_file ' (Channel 2)_avgSpecific_FluctuationStrength'], save_figs);

    il_plt_avgSpecific(ref_results.AvgSpecCombBinaural(:,2), OUT.specFluctStrengthAvgBin, ...
        [signals(i).wav_file ' (Combined binaural)_avgSpecific_FluctuationStrength'], save_figs);

    %% plot - single values, 90th percentile of the time-dependent quantity
    % the reference files carry the time series alone, so the same statistic
    % is applied to both, on the time samples of each

    ref_90Pc = [il_get_90Pc(ref_results.TDep(:,1), ref_results.TDep(:,2), time_skip), ...
                il_get_90Pc(ref_results.TDep(:,1), ref_results.TDep(:,3), time_skip), ...
                il_get_90Pc(ref_results.TDepCombBinaural(:,1), ref_results.TDepCombBinaural(:,2), time_skip)];

    sqat_90Pc = [OUT.fluctStrength90Pc(1), OUT.fluctStrength90Pc(2), OUT.fluctStrength90PcBin];

    il_plt_singleValues([ref_90Pc; sqat_90Pc], ...
        [signals(i).wav_file '_singleValues_FluctuationStrength'], save_figs);

    %% plot - time-dependent specific fluctuation strength

    il_plt_spectrogram(ref_results.Spec_TDep_channel_1(2:end,1), ...
        ref_results.Spec_TDep_channel_1(1,2:end), ref_results.Spec_TDep_channel_1(2:end,2:end), ...
        [signals(i).wav_file ' (Channel 1)_tDep_Specific_FluctuationStrength_ref'], save_figs);

    il_plt_spectrogram(OUT.timeOut, OUT.bandCentreFreqs, OUT.specFluctStrength(:,:,1), ...
        [signals(i).wav_file ' (Channel 1)_tDep_Specific_FluctuationStrength_implementation'], save_figs);

    il_plt_spectrogram(ref_results.Spec_TDep_channel_2(2:end,1), ...
        ref_results.Spec_TDep_channel_2(1,2:end), ref_results.Spec_TDep_channel_2(2:end,2:end), ...
        [signals(i).wav_file ' (Channel 2)_tDep_Specific_FluctuationStrength_ref'], save_figs);

    il_plt_spectrogram(OUT.timeOut, OUT.bandCentreFreqs, OUT.specFluctStrength(:,:,2), ...
        [signals(i).wav_file ' (Channel 2)_tDep_Specific_FluctuationStrength_implementation'], save_figs);

    %% printed comparison

    fprintf('\t90th percentile (vacil_HMS)   reference : %.4f %.4f %.4f\n', ref_90Pc);
    fprintf('\t                         implementation : %.4f %.4f %.4f\n', sqat_90Pc);

    labels = {'Channel 1', 'Channel 2', 'Combined binaural'};
    refTDep = {ref_results.TDep(:,[1 2]), ref_results.TDep(:,[1 3]), ref_results.TDepCombBinaural(:,[1 2])};
    sqatTDep = {OUT.fluctStrengthTDep(:,1), OUT.fluctStrengthTDep(:,2), OUT.fluctStrengthTDepBin};

    for k = 1:3
        [rmsDiff, maxDiff] = il_difference(refTDep{k}(:,1), refTDep{k}(:,2), OUT.timeOut, sqatTDep{k});
        fprintf('\ttime-dependent %-18s rms %.4f vacil_HMS, max abs %.4f vacil_HMS\n', ...
            labels{k}, rmsDiff, maxDiff);
    end

end

%% function - load the reference results of one signal

function ref_results = il_load_reference_results(ref_path, wav_file)
%function ref_results = il_load_reference_results(ref_path, wav_file)
%
% Loads the reference results exported from the commercial software for the
% signal <wav_file>. The exported files of the two channels carry three
% columns (time and one column per channel), and the combined binaural ones
% carry two columns, as the combined binaural result is a single channel.
%
% The time-dependent specific quantity of the two channels comes in a
% single file, where the block of the second channel is stacked below the
% block of the first one and carries its own header row of band centre
% frequencies. The two blocks are split here, so that each one has the same
% layout as the combined binaural file: the first row holds the band centre
% frequencies and the first column holds the time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quantity = 'Fluctuation Strength (Hearing Model)';

ref_results.AvgSpec = readmatrix([ref_path wav_file '.Specific ' quantity '.asc'], 'FileType', 'text');
ref_results.AvgSpecCombBinaural = readmatrix([ref_path wav_file '.Specific ' quantity '_combined_binaural.asc'], 'FileType', 'text');

ref_results.TDep = readmatrix([ref_path wav_file '.' quantity ' vs. Time.asc'], 'FileType', 'text');
ref_results.TDepCombBinaural = readmatrix([ref_path wav_file '.' quantity ' vs. Time_combined_binaural.asc'], 'FileType', 'text');

Spec_TDep = readmatrix([ref_path wav_file '.Specific ' quantity ' vs. Time.asc'], 'FileType', 'text');
ref_results.Spec_TDep_combined_binaural = readmatrix([ref_path wav_file '.Specific ' quantity ' vs. Time_combined_binaural.asc'], 'FileType', 'text');

% the header rows are the ones without a time value in the first column
header_rows = find(isnan(Spec_TDep(:,1)));

if numel(header_rows) ~= 2
    error('%s: expected two channel blocks in the time-dependent specific file, found %d', ...
        mfilename, numel(header_rows));
end

ref_results.Spec_TDep_channel_1 = Spec_TDep(header_rows(1):header_rows(2)-1, :);
ref_results.Spec_TDep_channel_2 = Spec_TDep(header_rows(2):end, :);

end % end of <il_load_reference_results> subfunction

%% function - 90th percentile of a time series from the reference files

function value = il_get_90Pc(time, series, time_skip)
%function value = il_get_90Pc(time, series, time_skip)
%
% 90th percentile of the time-dependent fluctuation strength, computed on
% the samples that follow <time_skip>. The sample where the window opens is
% chosen with the rule that FluctuationStrength_ECMA418_2 applies to its own
% output, including the floor of Section 9.1.12, which discards the first 36
% values of l_50. The two single values are then computed on the same
% samples of the two time series.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, time_skip_idx] = min(abs(time - time_skip));
time_skip_idx = max(time_skip_idx, 37);

value = prctile(series(time_skip_idx:end), 90);

end % end of <il_get_90Pc> subfunction

%% function - difference between the implementation and the reference

function [rmsDiff, maxDiff] = il_difference(xRef, yRef, xSQAT, ySQAT)
%function [rmsDiff, maxDiff] = il_difference(xRef, yRef, xSQAT, ySQAT)
%
% Root-mean-square and maximum absolute difference between the
% implementation and the reference, taken on the time samples of the
% reference, where the implementation is interpolated linearly because the
% two are exported on different time steps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yInterp = interp1(xSQAT(:), ySQAT(:), xRef(:), 'linear');

valid = ~isnan(yInterp);
difference = yInterp(valid) - yRef(valid);

rmsDiff = sqrt(mean(difference.^2));
maxDiff = max(abs(difference));

end % end of <il_difference> subfunction

%% function / plot - time-dependent quantity with the absolute difference

function il_plt_tDep(xRef, yRef, xSQAT, ySQAT, label_fig, save_figs)

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
stretchY = 1.5; % stretch plot in the vertical direction
set(h,'Position', [pos(1), pos(2), pos(3), pos(4)*stretchY]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*stretchY])

cmap  = load('cmap_viridis.txt');
cmap1 = 166;
cmap2 = 34;

tiledlayout(3,2);

%% first tile - the two time series
ax = nexttile([2 2]);

plot(xRef, yRef, 'Color', cmap(cmap1, :), 'Linewidth', 1); hold all;
plot(xSQAT, ySQAT, ':', 'Color', cmap(cmap2, :), 'Linewidth', 1.5);

ylabel('Fluctuation strength (vacil_{HMS})');

legend('Reference', 'Implementation', 'Location', 'NE');
legend boxoff

xlim([min(xRef) max(xRef)]);
ylim([0 ceil(max([yRef(:); ySQAT(:)])*15)/10]);

ax.FontName = 'Times';
ax.FontSize = 16;
ax.XTickLabel = [];

%% second tile - the absolute difference on the time samples of the reference
ax2 = nexttile([1 2]);

yInterp = interp1(xSQAT(:), ySQAT(:), xRef(:), 'linear');
plot(xRef, abs(yInterp - yRef(:)), '-', 'Color', [0 0 1]);

xlabel('Time (s)');
ylabel('Abs. difference');

xlim([min(xRef) max(xRef)]);

ax2.FontName = 'Times';
ax2.FontSize = 16;

set(gcf,'color','w');

il_save_fig(save_figs, label_fig);

end % end of <il_plt_tDep> subfunction

%% function / plot - time-averaged specific quantity

function il_plt_avgSpecific(yRef, ySQAT, label_fig, save_figs)

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

barkAxis = linspace(0.5, 26.5, 53);

ax = axes;

a = plot(barkAxis, yRef, '-', 'Linewidth', 1); hold all;
b = plot(barkAxis, ySQAT, ':', 'Linewidth', 1.5);

xlabel('Critical band rate (Bark_{HMS})');
ylabel('Specific fluctuation strength (vacil_{HMS}/Bark_{HMS})');

legend([a b], {'Reference', 'Implementation'}, 'Location', 'NE');
legend boxoff

xlim([0 27]);

ax.FontName = 'Times';
ax.FontSize = 16;

set(gcf,'color','w');

il_save_fig(save_figs, label_fig);

end % end of <il_plt_avgSpecific> subfunction

%% function / plot - single values

function il_plt_singleValues(single_values, label_fig, save_figs)

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = axes;

bar(single_values.');

set(ax, 'XTickLabel', {'Channel 1', 'Channel 2', 'Combined binaural'});
ylabel('Fluctuation strength, 90th percentile (vacil_{HMS})');

legend({'Reference', 'Implementation'}, 'Location', 'NE');
legend boxoff

ax.FontName = 'Times';
ax.FontSize = 16;

set(gcf,'color','w');

il_save_fig(save_figs, label_fig);

end % end of <il_plt_singleValues> subfunction

%% function / plot - time-dependent specific quantity

function il_plt_spectrogram(xAxis, yAxis, zAxis, label_fig, save_figs)

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

ax = axes;

[xx, yy] = meshgrid(xAxis, yAxis);
pcolor(xx, yy, zAxis.');
shading interp;
colormap(load('cmap_viridis.txt'));

c = colorbar;
c.Label.String = 'Specific fluctuation strength (vacil_{HMS}/Bark_{HMS})';

set(ax, 'YScale', 'log');
yticks([63 125 250 500 1000 2000 4000 8000]);
yticklabels({'63', '125', '250', '500', '1k', '2k', '4k', '8k'});

xlabel('Time (s)');
ylabel('Band centre frequency (Hz)');

ax.FontName = 'Times';
ax.FontSize = 16;

set(gcf,'color','w');

il_save_fig(save_figs, label_fig);

end % end of <il_plt_spectrogram> subfunction

%% function - save the current figure

function il_save_fig(save_figs, label_fig)

if save_figs == 1

    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_out = [figures_dir label_fig];

    resolution = '-r600'; % Default resolution of 600 DPI

    % saveas(gcf,figname_out, 'fig');
    % print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out, '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n', mfilename, label_fig, figname_out);
end

end % end of <il_save_fig> subfunction
