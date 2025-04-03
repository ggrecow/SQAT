function plt_Spectro(x, fs, df, olap, combine, save_figs, label_fig)
% function plt_Spectro(x, fs, df, olap, combine, save_figs, label_fig)
%
% Generates the following plot:
% Audio spectrogram
%
% Plot(s) used in the following publication:
%
% Lotinga, M., Torjussen, M., and Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS 
% OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED 
% SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS). 
% Forum Acusticum.
%
% Author: Mike Lotinga 20.03.2025
% Modified: Mike Lotinga 31.03.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Channel configuration
chansIn = size(x, 2);

if combine
    chansOut = 1;
else
    chansOut = chansIn;
end


%% STFT settings

% FFT block length
Nfft_nr = ceil(fs/df);
% if block length is odd, make it even
if mod(Nfft_nr, 2) == 1 
    Nfft = Nfft_nr - 1;
else
    Nfft = Nfft_nr;
end

% window overlap in samples
noverlap = olap*Nfft;


%% Generate spectrogram
AweightFilt = weightingFilter('A-weighting', fs);

for ii = chansIn:-1:1

        xA = AweightFilt(x(:, ii));
   
        [~, f_grm, t_grm, ps_grm(:, :, ii)] = spectrogram(xA, hann(Nfft, 'periodic'),...
                                noverlap, Nfft, fs, 'power', 'onesided');
        t_grm(1) = 0; % zero first time value
        t_grm(end) = round(t_grm(end)/5)*5; % round off last time value
end

if combine
    ps_grm = mean(ps_grm, 3);
end

%% Plot settings
% dimensions
figwidth = 24;
figheight = 16;
% octave band frequency ticks and labels
f_ticks = [31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000];
f_tick_labels = ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"];


%% Generate plot
for ii = 1:chansOut
    ps_grmdB = pow2db(ps_grm(2:end, :, ii)./(2e-5)^2);
    climMax = 5*ceil(max(ps_grmdB, [], 'all')/5);

    fig = figure;

    set(fig, 'Visible', 'on', 'units', 'centimeters',...
        'Position', [0, 0, figwidth, figheight]);
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Centimeters','PaperSize',[pos(3), pos(4)])
    movegui(fig, 'center');
    sh = surf(t_grm, f_grm(2:end), ps_grmdB);
    set(sh, 'LineStyle', 'none')
    view(0, 90)
    cbar = colorbar(gca, 'eastoutside');
    cbar_title = ylabel(cbar, "A-weighted sound pressure level,\newline           dB(A) re 2e-5 Pa");
    set(cbar_title, 'Rotation', 90, 'FontName', 'Times New Roman', 'FontSize', 18)
    xlabel("Time, s")
    ylabel("Frequency, Hz")
    set(gca, 'XLim', [0, ceil(t_grm(end))], 'YScale', 'log',...
        'YLim', [45, 17.8e3], 'YDir', 'normal',...
        'CLim', [0, climMax], 'YTick', f_ticks, 'YTickLabel', f_tick_labels,...
        'YMinorTick', 'off', 'FontName', 'Times New Roman', 'FontSize', 22);
    
    set(gcf,'color','w');
    
    %% Save plot
    if save_figs==1
        figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
        if ~exist(figures_dir,'dir')
            mkdir(figures_dir);
        end

        if ~combine
            figname_short = append(label_fig, "Ch", num2str(ii));
        else
            figname_short = label_fig;
        end

        figname_out = [figures_dir figname_short];
        
        resolution = '-r400'; 
    
        % saveas(gcf,figname_out, 'fig');
        % print( gcf, figname_out, '-dpdf', resolution , '-fillpage');
        print( gcf, figname_out,  '-dpng', resolution );
        
        fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
    end
end

end