function plt_tDepSpecific(x1Axis, y1Axis, z1Axis, x2Axis, y2Axis, z2Axis, metric, save_figs, label_fig)
% function plt_tDepSpecific(x1Axis, y1Axis, z1Axis, x2Axis, y2Axis, z2Axis, metric, save_figs, label_fig)
%
% Generates the following plot:
% Time-dependent specific sound quality metric (channel 1)
%
% Plot(s) used in the following publication:
%
% Lotinga, M., Torjussen, M., and Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS
% OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED
% SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS).
% Forum Acusticum.
%
% Author: Gil Felix Greco, Braunschweig 27.02.2025
% Modified: 31.03.2025 Mike Lotinga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch metric
    case 'tonality'
        cmap  = load('cmap_plasma.txt');
        yLabel1 = 'tonality, tu';
        precDig = 1;
    case 'loudness'
        cmap  = load('cmap_viridis.txt');
        yLabel1 = 'loudness, sone';
        precDig = 1;
    case 'roughness'
        cmap  = load('cmap_inferno.txt');
        yLabel1 = 'roughness, asper';
        precDig = 2;
end

stretchY = 1.5; % stretch plot in the vertical direction

h = figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'Position', [pos(1), pos(2), pos(3), pos(4)*stretchY]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*stretchY])
movegui(h, 'center');

tiledlayout(2,1);
climMax = ceil(max(max(z1Axis, [], 'all'),...
               max(z2Axis, [], 'all'))*5^precDig)/10^precDig;
fontSize = 22;

%% reference results

ax1 = nexttile;
[xx,yy] = meshgrid(x1Axis, y1Axis);
xx(:,1) = 0;  % truncating first idx of time to zero (just a visual thing because ref results may not start exactly on zero)
xx(:,end) = round(x1Axis(end)/5)*5;  % truncating last idx of time to nearest 5s (just a visual thing because ref results are almost 30s)
pcolor(xx, yy, z1Axis);
shading interp; axis tight;

ax1.FontName = 'Times';
ax1.FontSize = fontSize;
ax1.YTick = [63, 250, 1e3, 4e3, 16e3];
ax1.YTickLabel = ["63", "250",  "1k", "4k", "16k"];
ax1.YScale = 'log';
ax1.YLabel.String = 'Frequency, Hz';
ax1.YLabel.FontSize = fontSize;
ax1.XLabel.String = 'Time, s';
ax1.XLabel.FontSize = fontSize;
ax1.Title.String = 'Reference';
ax1.Title.FontWeight = 'normal';
ax1.Title.FontSize = fontSize;
colormap(ax1,cmap);
clim([0 climMax]);

%% Implementation results

ax2 = nexttile;
[xx,yy] = meshgrid(x2Axis, y2Axis);
xx(:,1) = 0;  % truncating first idx of time to zero (just a visual thing because ref results may not start exactly on zero)
xx(:,end) = round(x2Axis(end)/5)*5;  % truncating last idx of time to 30s (just a visual thing because ref results are almost 30s)
pcolor(xx, yy, z2Axis);
shading interp; axis tight;

ax2.FontName = 'Times';
ax2.FontSize = fontSize;
ax2.YTick = [63, 250, 1e3, 4e3, 16e3];
ax2.YTickLabel = ["63", "250",  "1k", "4k", "16k"];
ax2.YScale = 'log';
ax2.YLabel.String = 'Frequency, Hz';
ax2.YLabel.FontSize = fontSize;
ax2.XLabel.String = 'Time, s';
ax2.XLabel.FontSize = fontSize;
ax2.Title.String = 'Implementation';
ax2.Title.FontWeight = 'normal';
ax2.Title.FontSize = fontSize;
colormap(ax2,cmap);
clim([0 climMax]);

%% Common features

cb = colorbar;
cb.Layout.Tile = 'east';

zString = append('Specific ', yLabel1, '_{HMS}/Bark_{HMS}');
set(get(cb,'label'),'string', zString, 'fontsize', fontSize);

set(gcf,'color','w');

%% save figure

if save_figs==1
    figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = label_fig;
    figname_out = [figures_dir figname_short];

    resolution = '-r400';

    % saveas(gcf,figname_out, 'fig');
    % print( gcf, figname_out, '-dpdf', resolution );
    print( gcf, figname_out,  '-dpng', resolution );

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end