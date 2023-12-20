function plot_EPNL_comparison(ref_values, SQAT_values, tool_tag, save_figs, dir_out)

h  = figure;
set(h, 'name', ['EPNL comparison -' tool_tag ] );
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

X = categorical({'5','6','7','8','8a'});
X = reordercats(X,{'5','6','7','8','8a'});

plot(X, ref_values, 'sk', 'Linewidth', 1.2, 'Markersize', 12); hold on;
plot(X, SQAT_values, 'xk', 'Linewidth', 1.2, 'Markersize', 12);

ymin = 76; ymax = 96;
ylim([ymin ymax]);
ylabel('EPNL (EPNdB)','Interpreter','Latex');
xlabel('Case $\#$','Interpreter','Latex');

legend( [tool_tag ' (Reference)'],...
           [tool_tag ' (SQAT)'],...
           'Location','NE','Interpreter','Latex');

ax = gca;
ax.YAxis.MinorTick = 'on';
dy = 1;
ax.YAxis.MinorTickValues =  ymin:dy:ymax;

grid on;
set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    figures_dir = [dir_out 'figs' filesep];
    
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = [ 'EPNL_cross_validation_' tool_tag ];
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end