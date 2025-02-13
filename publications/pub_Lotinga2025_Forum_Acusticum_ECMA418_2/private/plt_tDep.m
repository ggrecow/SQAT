function plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)
% function plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)
%
% Generates the following plot:
% Time-dependent quantities (tonality, loudness, and roughness) 
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

switch metric
    case 'roughness'

        xLabelTag = 'Time (s)';
        yLabelTag = 'Roughness (asper_{HMS})';
        cMap = load('cmap_inferno.txt');
        yLimits = [0 0.2];

    case 'loudness'

        xLabelTag = 'Time (s)';
        yLabelTag = 'Loudness (sone_{HMS})';
        cMap = load('cmap_viridis.txt');
        yLimits = [0 12];

    case 'tonality'

        xLabelTag = 'Time (s)';
        yLabelTag = 'Tonality (tu_{HMS})';
        cMap = load('cmap_plasma.txt');
        yLimits = [0 1.5];

end

% plot
h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

cMap1 = 166;
cMap2 = 34;

ax = nexttile;
% reference values 
plot(xRef, yRef, 'Color', cMap(cMap1, :), 'Linewidth',1); hold all;

% implementation values 
plot(xImplementation, yImplementation, ':', 'Color', cMap(cMap2, :), 'Linewidth', 1.5);  

legend('Reference', 'Implementation', 'Location', 'SW');
set(gcf,'color','w');

ylim(yLimits);

ax.FontName = 'Times';
ax.FontSize = 18;

xlabel(xLabelTag, 'fontsize',18);
ylabel(yLabelTag, 'fontsize',18);

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