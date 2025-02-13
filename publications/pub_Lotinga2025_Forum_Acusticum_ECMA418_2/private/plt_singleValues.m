function plt_singleValues(single_values1, single_values2, single_values3, save_figs,  label_fig)
% function plt_singleValues(single_values1, single_values2, single_values3, save_figs,  label_fig)
%
% Generates the following plot:
% Overall loudness, tonality, and roughness (channel 1)
%
% Plot(s) used in the following publication:
%
% Lotinga, M., Torjussen, M., and Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS
% OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED
% SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS).
% Forum Acusticum.
%
% Author: Gil Felix Greco, Braunschweig 13.02.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
movegui(h, 'center');

tiledlayout("horizontal")

%% first tile (tonality)

nexttile

str_xlabel = {'Ch1'};
br=bar(str_xlabel, single_values1);

cMap = load('cmap_inferno.txt');
colororder([cMap(166, :); cMap(34, :)])

for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, append( ' ', string(round(br(bb).YData, 2) ) ),...
        'rotation', 90, 'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', 14)
end

ylim([0 1]);
xticks([]);

ax = gca;
ax.FontName = 'Times';
ax.FontSize = 18;
ylabel( 'Tonality (tu_{HMS})', 'fontsize',18);

%% second tile (loudness)

nexttile

str_xlabel = {'Ch1'};
br=bar(str_xlabel, single_values2);

for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, append( ' ', string(round(br(bb).YData, 2) ) ),...
        'rotation', 90, 'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', 14)
end

ylim([0 10]);
xticks([]);

legend('Reference', 'Implementation', 'Location', 'northoutside', 'NumColumns', 2);

ax = gca;
ax.FontName = 'Times';
ax.FontSize = 18;
ylabel( 'Loudness (sones_{HMS})', 'fontsize',18);

%% third tile (roughness)

nexttile

str_xlabel = {'Ch1'};
br=bar(str_xlabel, single_values3);

for bb = 1:length(br)
    text(br(bb).XEndPoints, br(bb).YEndPoints, append( ' ', string(round(br(bb).YData, 2) ) ),...
        'rotation', 90, 'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', 14)
end

ylim([0 0.5]);
xticks([]);

ax = gca;
ax.FontName = 'Times';
ax.FontSize = 18;
ylabel( 'Roughness (asper_{HMS})', 'fontsize',18);

set(gcf,'color','w');

%% save fig

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];

    resolution = '-r600';
    saveas(gcf,figname_out, 'fig');
    print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end