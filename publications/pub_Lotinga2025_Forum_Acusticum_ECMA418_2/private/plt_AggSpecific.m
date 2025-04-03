function plt_AggSpecific( yRef, yImplementation, metric, save_figs, label_fig)
% function plt_AggSpecific( yRef, yImplementation, metric, save_figs, label_fig)
%
% Generates the following plot:
% Time-aggregated specific sound quality metric (channel 1)
%
% Plot(s) used in the following publication:
%
% Lotinga, M., Torjussen, M., and Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS
% OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED
% SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS).
% Forum Acusticum.
%
% Author: Gil Felix Greco, Braunschweig 27.02.2025
% Modified: 21.03.2025 Mike Lotinga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch metric
    case 'tonality'
        cmap  = load('cmap_plasma.txt');
        yLabel1 = 'tonality, tu';
        yLabel2 = 'Tonality difference, tu';
        precDig = 1;
    case 'loudness'
        cmap  = load('cmap_viridis.txt');
        yLabel1 = 'loudness, sone';
        yLabel2 = 'Loudness difference, sone';
        precDig = 1;
    case 'roughness'
        cmap  = load('cmap_inferno.txt');
        yLabel1 = 'roughness, asper';
        yLabel2 = 'Roughness difference, asper';
        precDig = 2;
end

cmap1 = 166;
cmap2 = 34;
fontSize = 22;

stretchX = 2;  % stretch plot in the horizontal direction
stretchY= 1.5; % stretch plot in the vertical direction

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3)*stretchX, pos(4)*stretchY]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)*stretchX, pos(4)*stretchY])
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
    'EdgeColor', cmap(cmap1, :),...
    'EdgeAlpha', 1, ...
    'FaceColor', 'none',...
    'FaceAlpha', 1, ...
    'LineWidth', 1.5, 'LineStyle', '-',...
    'BarWidth',0.7, ...
    'DisplayName', "Reference");

hold on;

b = bar(barkAxis,...
    yImplementation,...
    'EdgeColor', cmap(cmap2, :),...
    'FaceColor', 'none',...
    'LineWidth',1, 'LineStyle', '--',...
    'BarWidth',0.7, ...
    'DisplayName', "Implementation");

xtickangle(90);
xticks(barkAxis);

ylim([0 ceil(max(max(yRef), max(yImplementation))*10^precDig)/10^precDig]);

ax.XTick = 0.5:26.5;
ax.FontName = 'Times';
ax.FontSize = fontSize;

xlabel( 'Critical band rate, Bark_{HMS}', 'fontsize', fontSize);
ylabel( append('Specific ', yLabel1, '_{HMS}/Bark_{HMS}'), 'fontsize', fontSize);

%% right axis (error plot)

yyaxis right;
ax.YAxis(2).Exponent = 0;
ytickformat('%.3f')

% difference = ( (yImplementation-yRef)./yRef).*100; % percentage difference
difference = (yImplementation-yRef); % relative difference



c = plot( barkAxis, difference, '*-', 'Color', right_color );

ylabel( append(yLabel2 ,'_{HMS}/Bark_{HMS}'), 'fontsize', fontSize);

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

    %saveas(gcf,figname_out, 'fig');
    print( gcf, figname_out, '-dpdf', resolution );
    %print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end