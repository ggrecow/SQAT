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
% Author: Gil Felix Greco, Braunschweig 27.02.2025
% Modified: 20.03.2025 Mike Lotinga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch metric
    case 'roughness'

        xLabelTag = 'Time, s';
        yLabelTag = 'Roughness, asper_{HMS}';
        cMap = load('cmap_inferno.txt');

    case 'loudness'

        xLabelTag = 'Time, s';
        yLabelTag = 'Loudness, sone_{HMS}';
        cMap = load('cmap_viridis.txt');

    case 'tonality'

        xLabelTag = 'Time, s';
        yLabelTag = 'Tonality, tu_{HMS}';
        cMap = load('cmap_plasma.txt');
end

yLimits = [0 ceil(max(max(yRef), max(yImplementation))*15)/10];
fontSize = 22;
stretchY = 1.5; % stretch plot in the vertical direction

% plot
h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3), pos(4)*stretchY]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*stretchY])
movegui(h, 'center');

cMap1 = 166;
cMap2 = 34;

tiledlayout(3,2)

%% first tile
ax = nexttile([2 2]);

% reference values 
plot(xRef, yRef, 'Color', cMap(cMap1, :), 'Linewidth',1); hold all;
 
% implementation values 
plot(xImplementation, yImplementation, ':', 'Color', cMap(cMap2, :), 'Linewidth', 1.5);  

legend('Reference', 'Implementation', 'Location', 'NE');

xlim([0, round(xRef(end))])
ylim(yLimits);

ax.FontName = 'Times';
ax.FontSize = fontSize;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridAlpha = 0.15;

% xlabel(xLabelTag, 'fontsize', fontSize);
xticklabels(""); % Remove x axis numbering
ylabel(yLabelTag, 'fontsize', fontSize);

%% second tile

ax2 =  nexttile([1 2]);

plot( xImplementation(1:size(yRef)), abs(yImplementation(1:size(yRef))-yRef), '-', 'Color', [0 0 1] );

yticks([0 0.01]);
ax2.FontName = 'Times';
ax2.FontSize = fontSize;
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.GridLineStyle = '--';
ax2.GridAlpha = 0.15;

xlabel(xLabelTag, 'fontsize', fontSize);
ylabel('Abs. difference', 'fontsize', fontSize);

xlim([0, round(xRef(end))])
ylim([0  0.01]);

set(gcf,'color','w');

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];
    
    resolution = '-r600'; 

    % saveas(gcf,figname_out, 'fig');
    print( gcf, figname_out, '-dpdf', resolution );
    % print( gcf, figname_out,  '-dpng', resolution );
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end