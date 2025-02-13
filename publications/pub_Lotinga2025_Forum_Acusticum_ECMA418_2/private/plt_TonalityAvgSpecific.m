function plt_TonalityAvgSpecific( yRef, yImplementation, save_figs, label_fig)
% function plt_TonalityAvgSpecific( yRef, yImplementation, save_figs, label_fig)
%
% Generates the following plot:
% Time-averaged specific tonality (channel 1)
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

cmap  = load('cmap_plasma.txt');
cmap1 = 166;
cmap2 = 34;

h  =figure('Position', [200, 200, 1500, 550]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

barkAxis = linspace(0.5, 26.5, 53);

left_color = [0 0 0];
right_color = [0 0 0];
set(h,'defaultAxesColorOrder',[left_color; right_color]);

ax = nexttile;

%% left axis (bar plots)

yyaxis left;

a = bar(barkAxis,...
    yRef,...
    'EdgeColor',  cmap(cmap1, :),...
    'EdgeAlpha', 0.2, ...
    'FaceColor', cmap(cmap1, :),...
    'FaceAlpha', 0.2, ...
    'LineWidth', 1, 'LineStyle', '-',...
    'BarWidth',0.5, ...
    'DisplayName', "Reference");

hold on;

b = bar(barkAxis,...
    yImplementation,...
    'EdgeColor',  cmap(cmap2, :),...
    'FaceColor',  'none',...
    'LineWidth', 0.5, 'LineStyle', '-',...
    'BarWidth',0.5, ...
    'DisplayName', "Implementation");

xlabel( 'Critical band rate (Bark_{HMS})' );
ylabel( 'Specific tonality (tu_{HMS}/Bark_{HMS})' );

xtickangle(90);
xticks(barkAxis);

ylim([0 0.7]);

ax.XTick = 0.5:26.5;
% ax.XTick = 1:26;
ax.FontName = 'Times';
ax.FontSize = 18;

%% right axis (error plot)

yyaxis right;

c = plot( barkAxis, yImplementation-yRef, 'k*-');

% ylabel( 'Implementation-Reference (tu_{HMS}/Bark_{HMS})' );

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

    resolution = '-r600';
    saveas(gcf,figname_out, 'fig');
    print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end