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
fontSize = 20;

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3)*2, pos(4)*1.3]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)*2, pos(4)*1.3])
movegui(h, 'center');

barkAxis = linspace(0.5, 26.5, 53);

left_color = [0 0 0]; % left axis - black
right_color = [0 0 1]; % right axis - blue
set(h,'defaultAxesColorOrder',[left_color; right_color]);

ax = nexttile;

%% left axis (bar plots)

yyaxis left;

a = bar(barkAxis,...
    yRef,...
    'EdgeColor',  cmap(cmap1, :),...
    'EdgeAlpha', 1, ...
    'FaceColor', 'none',...
    'FaceAlpha', 1, ...
    'LineWidth', 1.5, 'LineStyle', '-',...
    'BarWidth',0.7, ...
    'DisplayName', "Reference");

hold on;

b = bar(barkAxis,...
    yImplementation,...
    'EdgeColor',  cmap(cmap2, :),...
    'FaceColor',  'none',...
    'LineWidth',1, 'LineStyle', '--',...
    'BarWidth',0.7, ...
    'DisplayName', "Implementation");

xtickangle(45);
xticks(barkAxis);

ylim([0 0.7]);

ax.XTick = 0.5:26.5;
ax.FontName = 'Times';
ax.FontSize = fontSize;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridAlpha = 0.15;

xlabel( 'Critical band rate, Bark_{HMS}', 'fontsize', fontSize);
ylabel( 'Specific tonality, tu_{HMS}/Bark_{HMS}', 'fontsize', fontSize);

%% right axis (error plot)

yyaxis right;

% difference = ( (yImplementation-yRef)./yRef).*100; % percentage difference
difference = (yImplementation-yRef); % relative difference

c = plot( barkAxis, difference, '*-', 'Color', right_color );

ylabel( 'Relative difference, tu_{HMS}/Bark_{HMS}', 'fontsize', fontSize);

ylim([-0.0015 0.0015]); 

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