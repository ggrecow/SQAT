function pub_Osses2023a_Forum_Acusticum_SQAT
% function pub_Osses2023a_Forum_Acusticum_SQAT
%
% Generates the tables and figures presented in the contribution to be 
%   presented at Forum Acusticum in September 2023.
%
% Reference:
% Osses, A., Felix Greco, G., and Merino-Martinez, R. (2023) Considerations
%   for the perceptual evaluation of steady-state and time vareying sounds 
%   using psychoacoustic metrics. Forum Acusticum.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

do_table2  = 0; %  Acoustic characterisation, all
do_table3  = 0; %  Psychoacoustic characterisation, all
do_fig1    = 0; % Psychoacoustic characterisation, selected (+plots)
do_fig_raw = 0; % Psychoacoustic characterisation, all (+plots)

list_tables_and_figures = {'do_table2';'do_table3';'do_fig1'; 'do_fig_raw'};
for i = 1:length(list_tables_and_figures)
    fprintf('Enter %.0f to %s\n',i,list_tables_and_figures{i});
end
bInput = input('Which table/figure you want to obtain (enter the corresponding number)?: ');
% % Setting to 1 the selected choice:
% exp2eval = sprintf('%s=1;',list_tables_and_figures{bInput});
% eval(exp2eval);
switch list_tables_and_figures{bInput}
    case 'do_table2'
        do_table2 = 1;
    case 'do_table3'
        do_table3 = 1;
    case 'do_fig1'
        do_fig1 = 1;
    case 'do_fig_raw'
        do_fig_raw = 1;
end
dir_SQAT = basepath_SQAT;

dir_curr = [cd filesep];
cd(dir_SQAT); % We change to dir_SQAT to avoid the risk of shadowing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dir_sounds = '/home/alejandro/Desktop/SQM_paper/FA2023_paper/Sounds-selection_no-commit/';
dir_sounds = [basepath_SQAT 'sound_files' filesep 'publications' filesep 'pub_Osses2023a_Forum_Acusticum_SQAT' filesep];

% file:///home/alejandro/Documents/MATLAB/MATLAB_SQAT/sound_files/publications/pub_Osses2023a_Forum_Acusticum_SQAT/
files_1 = { 'Flyover1_Airbus319-114-dBFS.wav'; 'Flyover2_Airbus319-114-dBFS.wav'; ...
            'Flyover3_Boeing737-114-dBFS.wav'; 'Flyover4_Boeing737-114-dBFS.wav'; ...
            'Flyover5_Fokker70-114-dBFS.wav';  'Flyover6_Fokker70-114-dBFS.wav'};
% file:///home/alejandro/Documents/MATLAB/MATLAB_SQAT/sound_files/publications/pub_Osses2023a_Forum_Acusticum_SQAT/2-Train-pass-by/
files_2 = { 'train_01.wav'; 'train_02.wav'; 'train_11.wav'; 'train_15.wav'; ...
            'train_18.wav'; 'train_20.wav'; 'train_21.wav'};

% file:///home/alejandro/Documents/MATLAB/MATLAB_SQAT/sound_files/publications/pub_Osses2023a_Forum_Acusticum_SQAT/3-Hummer-resonances/
files_3 = {'meas-ac-2-dist-ane.wav'; 'meas-ac-4-dist-ane-HP.wav'; 'model-ac-2-dist-ane.wav'; 'model-ac-4-dist-ane.wav'};

dir_datasets = {'1-Aircraft-fly-by'  ,114    , files_1; ...
                '2-Train-pass-by'    ,140.55 , files_2; ...
                '3-Hummer-resonances',100    , files_3};
if do_table2
    
    dBFS_4_reproduction = 120; % not too loud not too soft
    list_metrics = {'LAeq','LAFmax','LZeq','T','SEL','Delta_Leq','LAFmax_min_LAeq'};
    metrics_row = [];
    for i_dir = 1:size(dir_datasets,1)
        dBFS = dir_datasets{i_dir,2};
        dir_here = [dir_sounds dir_datasets{i_dir,1} filesep];
        files = dir_datasets{i_dir,3};
        
        N(i_dir) = length(files);
        files_data(i_dir).files = files;
        
        for i_files = 1:length(files)
            [insig,fs] = audioread([dir_here files{i_files}]);
            
            il_get_the_metrics(insig,fs,dBFS,list_metrics);
            [res,res_description,outs] = il_get_the_metrics(insig,fs,dBFS,list_metrics);
            
            gain4reproduction = dBFS-dBFS_4_reproduction;
            disp(''); close; % sound(gain4reproduction*insig(min(outs.idx):max(outs.idx),:),fs);
            metrics_row(end+1,:) = res;
        end
    end
    res_description
    il_var2latex(round(10*metrics_row)/10)
    disp('')
end

if do_fig1 || do_table3 || do_fig_raw
    % Starts by repeating the calculation as in do_table1, but here we are 
    %   interested in the exact time of the sound event 'idx'.
    
    List_psy_metrics = {'Loudness_ISO532_1', ...
                        'Sharpness_DIN45692', ...
                        'Roughness_Daniel1997', ...
                        'FluctuationStrength_Osses2016', ...
                        'Tonality_Aures1985'};
    
    list_metrics = {'LAeq','LAFmax','LZeq','T','SEL','Delta_Leq','LAFmax_min_LAeq'};
    metrics_row = [];
    
    count = 1;
    for i_dir = 1:size(dir_datasets,1)
        dBFS = dir_datasets{i_dir,2};
        dir_here = [dir_sounds dir_datasets{i_dir,1} filesep];
        
        %%%
        if do_table3 || do_fig_raw
            files = dir_datasets{i_dir,3};
        end
        if do_fig1
            % The following was a visual selection, based on the results of
            %  Table 3:
            switch dir_datasets{i_dir,1}
                case '1-Aircraft-fly-by'
                    files = {'Flyover2_Airbus319-114-dBFS.wav','Flyover3_Boeing737-114-dBFS.wav'};
                case '2-Train-pass-by'
                    files = {'train_01.wav','train_02.wav'};
                case '3-Hummer-resonances'
                    files = {'meas-ac-2-dist-ane.wav','model-ac-2-dist-ane.wav'};
            end
        end
        
        N(i_dir) = length(files);
        files_data(i_dir).files = files;
        
        for i_files = 1:length(files)
            bDo = 1; % idle variable (to manually set to 0, in debug mode)
            if bDo
                if do_fig1 || do_fig_raw || do_table3
                    h = [];
                    hname = [];
                end

                [insig,fs] = audioread([dir_here files{i_files}]);
                [res_SLM,res_description_SLM,outs] = il_get_the_metrics(insig,fs,dBFS,list_metrics);

                res = []; % where the psychoacoustic metrics will be stored
                res_description = [];
                for i_psy = 1:length(List_psy_metrics)
                    psy_metric = List_psy_metrics{i_psy};
                    params = psychoacoustic_metrics_get_defaults(psy_metric);

                    if do_fig1 || do_fig_raw
                        bPlot = 0;
                    else
                        bPlot = 1; % params.show;
                    end

                    switch psy_metric
                        case 'Loudness_ISO532_1'
                            params4input = {params.field, params.method, params.time_skip, bPlot};
                        case 'Sharpness_DIN45692'
                            params4input = {params.weight_type, params.field, params.method, params.time_skip, params.show_loudness, bPlot};
                        case 'Roughness_Daniel1997'
                            params4input = {params.time_skip, bPlot};
                        case 'FluctuationStrength_Osses2016'
                            params4input = {params.method, params.time_skip, bPlot};
                        case 'Tonality_Aures1985'
                            params4input = {params.Loudness_field, params.time_skip, bPlot};
                        otherwise
                            error('Continue!');
                    end

                    bEvent = 1;
                    bShort = ~bEvent;
                    if bShort
                        insig = insig(1:fs,1);
                        toffset = 0;
                    else
                        toffset = 0;
                    end
                    idxi = round((outs.t_min_max(1)-toffset)*fs);
                    idxf = round(outs.t_min_max(2)*fs);
                    insig_here = insig(idxi:idxf,1);
                    insig94 = 10^((dBFS-94)/20)*insig_here;
                    switch psy_metric
                        case 'Loudness_ISO532_1'
                            out = Loudness_ISO532_1(insig94,fs,params4input{:});

                            if toffset == 0
                                disp('')
                            end
                            t = out.time+(outs.t_min_max(1))-toffset;
                            N_t = out.InstantaneousLoudness;

                            [Nmax, idx_max] = max(N_t);
                            idx_perc_90(1) = find(N_t>out.N90,1,'first');
                            idx_perc_90(2) = find(N_t>out.N90,1,'last');
                            idx_perc_50(1) = find(N_t>out.N50,1,'first');
                            idx_perc_50(2) = find(N_t>out.N50,1,'last');

                            %%% Time plot:
                            figure;
                            plot(t,N_t,'b','LineWidth',2); hold on; grid on

                            ylabel('Total loudness (sone)');
                            xlabel('Time (s)');
                            title(files{i_files},'interpreter','none'); 

                            xlim(outs.t_min_max);
                            YL = get(gca,'YLim');
                            plot(t(idx_max)*[1 1],[YL(1) N_t(idx_max)],'k--');
                            plot(t(idx_perc_50(1))*[1 1],[YL(1) N_t(idx_perc_50(1))],'r-');
                            plot(t(idx_perc_50(2))*[1 1],[YL(1) N_t(idx_perc_50(2))],'r-');
                            plot(t(idx_perc_90(1))*[1 1],[YL(1) N_t(idx_perc_90(1))],'r--');
                            plot(t(idx_perc_90(2))*[1 1],[YL(1) N_t(idx_perc_90(2))],'r--');

                            h(end+1) = gcf;
                            hname{end+1} = ['c' num2str(count) '-loudness-time'];

                            %%% Frequency plot:
                            Nspec_max = out.InstantaneousSpecificLoudness(idx_max,:);
                            Nspec_50  = mean(out.InstantaneousSpecificLoudness(idx_perc_50,:));
                            Nspec_90  = mean(out.InstantaneousSpecificLoudness(idx_perc_90,:));

                            f = out.barkAxis;
                            f_Tick = .5:1:23.5;
                            % f_Hz = bark2hz(f);
                            f_Tick_Hz = round( bark2hz(f_Tick) );

                            f_Tick_Hz_txt = [];
                            for i_f = 1:length(f_Tick_Hz)
                                if mod(i_f,2) == 1
                                    if f_Tick_Hz(i_f) < 1000
                                        f_Tick_Hz_txt{i_f} = num2str(f_Tick_Hz(i_f));
                                    else
                                        f_Tick_Hz_txt{i_f} = [num2str(round(10*f_Tick_Hz(i_f)/1000)/10) ' k'];
                                    end
                                else
                                    f_Tick_Hz_txt{i_f} = '';
                                end                            
                            end

                            figure;
                            plot(f,Nspec_max,'k','LineWidth',2); hold on; grid on
                            plot(f,Nspec_50,'r-');
                            plot(f,Nspec_90,'r--');

                            set(gca,'XTick',f_Tick);
                            set(gca,'XTickLabel',f_Tick_Hz_txt);
                            ylabel('Specific loudness (sone/Bark)');
                            xlabel('Frequency (Hz)');
                            title(files{i_files},'interpreter','none'); 

                            h(end+1) = gcf;
                            hname{end+1} = ['c' num2str(count) '-loudness-freq'];

                            res(1,end+1) = out.N50;     res_description{1,end+1} = 'N50 (sone)';
                            res(1,end+1) = Nmax;        res_description{1,end+1} = 'Nmax (sone)';
                            res(1,end+1) = t(idx_max);  res_description{1,end+1} = 'tmax (s)';

                        case 'Sharpness_DIN45692'
                            out = Sharpness_DIN45692(insig94,fs,params4input{:});
                            t = out.time+(outs.t_min_max(1))-toffset;
                            S_t = out.InstantaneousSharpness;

                            [Smax, idx_max] = max(S_t);
                            idx_perc_90(1) = find(S_t>out.S90,1,'first');
                            idx_perc_90(2) = find(S_t>out.S90,1,'last');
                            idx_perc_50(1) = find(S_t>out.S50,1,'first');
                            idx_perc_50(2) = find(S_t>out.S50,1,'last');

                            %%% Time plot:
                            figure;
                            plot(t,S_t,'b','LineWidth',2); hold on; grid on

                            ylabel('Sharpness (acum)');
                            xlabel('Time (s)');
                            title(files{i_files},'interpreter','none'); 

                            xlim(outs.t_min_max);
                            YL = get(gca,'YLim');
                            plot(t(idx_max)*[1 1],[YL(1) S_t(idx_max)],'k--');
                            plot(t(idx_perc_50(1))*[1 1],[YL(1) S_t(idx_perc_50(1))],'r-');
                            plot(t(idx_perc_50(2))*[1 1],[YL(1) S_t(idx_perc_50(2))],'r-');
                            plot(t(idx_perc_90(1))*[1 1],[YL(1) S_t(idx_perc_90(1))],'r--');
                            plot(t(idx_perc_90(2))*[1 1],[YL(1) S_t(idx_perc_90(2))],'r--');

                            h(end+1) = gcf;
                            hname{end+1} = ['c' num2str(count) '-sharpness-time'];

                            res(1,end+1) = out.S50;     res_description{1,end+1} = 'S50 (sone)';
                            res(1,end+1) = Smax;        res_description{1,end+1} = 'Smax (sone)';
                            res(1,end+1) = t(idx_max);  res_description{1,end+1} = 'tmax (s)';

                    case 'Roughness_Daniel1997'
                        out = Roughness_Daniel1997(insig94,fs,params4input{:});

                        if toffset == 0
                            disp('')
                        end
                        t = out.time+(outs.t_min_max(1))-toffset;
                        R_t = out.InstantaneousRoughness;

                        [Rmax, idx_max] = max(R_t);
                        idx_perc_90(1) = find(R_t>out.R90,1,'first');
                        idx_perc_90(2) = find(R_t>out.R90,1,'last');
                        idx_perc_50(1) = find(R_t>out.R50,1,'first');
                        idx_perc_50(2) = find(R_t>out.R50,1,'last');

                        %%% Time plot:
                        figure;
                        plot(t,R_t,'b','LineWidth',2); hold on; grid on

                        ylabel('Roughness (asper)');
                        xlabel('Time (s)');
                        title(files{i_files},'interpreter','none'); 

                        xlim(outs.t_min_max);
                        YL = get(gca,'YLim');
                        plot(t(idx_max)*[1 1],[YL(1) R_t(idx_max)],'k--');
                        plot(t(idx_perc_50(1))*[1 1],[YL(1) R_t(idx_perc_50(1))],'r-');
                        plot(t(idx_perc_50(2))*[1 1],[YL(1) R_t(idx_perc_50(2))],'r-');
                        plot(t(idx_perc_90(1))*[1 1],[YL(1) R_t(idx_perc_90(1))],'r--');
                        plot(t(idx_perc_90(2))*[1 1],[YL(1) R_t(idx_perc_90(2))],'r--');

                        h(end+1) = gcf;
                        hname{end+1} = ['c' num2str(count) '-roughness-time'];

                        %%% Frequency plot:
                        Rspec_max = out.InstantaneousSpecificRoughness(:,idx_max);
                        dim_here = 2;
                        Rspec_50  = mean(out.InstantaneousSpecificRoughness(:,idx_perc_50),dim_here);
                        Rspec_90  = mean(out.InstantaneousSpecificRoughness(:,idx_perc_90),dim_here);

                        f = out.barkAxis;
                        f_Tick = .5:1:23.5;
                        % f_Hz = bark2hz(f);
                        f_Tick_Hz = round( bark2hz(f_Tick) );

                        f_Tick_Hz_txt = [];
                        for i_f = 1:length(f_Tick_Hz)
                            if mod(i_f,2) == 1
                                if f_Tick_Hz(i_f) < 1000
                                    f_Tick_Hz_txt{i_f} = num2str(f_Tick_Hz(i_f));
                                else
                                    f_Tick_Hz_txt{i_f} = [num2str(round(10*f_Tick_Hz(i_f)/1000)/10) ' k'];
                                end
                            else
                                f_Tick_Hz_txt{i_f} = '';
                            end                            
                        end

                        figure;
                        plot(f,Rspec_max,'k','LineWidth',2); hold on; grid on
                        plot(f,Rspec_50,'r-');
                        plot(f,Rspec_90,'r--');

                        set(gca,'XTick',f_Tick);
                        set(gca,'XTickLabel',f_Tick_Hz_txt);
                        ylabel('Specific roughness (asper/Bark)');
                        xlabel('Frequency (Hz)');
                        title(files{i_files},'interpreter','none'); 

                        h(end+1) = gcf;
                        hname{end+1} = ['c' num2str(count) '-roughness-freq'];

                        res(1,end+1) = out.R50;     res_description{1,end+1} = 'R50 (asper)';
                        res(1,end+1) = Rmax;        res_description{1,end+1} = 'Rmax (asper)';
                        res(1,end+1) = t(idx_max);  res_description{1,end+1} = 'tmax (s)';

                    case 'FluctuationStrength_Osses2016'
                        out = FluctuationStrength_Osses2016(insig94,fs,params4input{:});

                        t = out.time+(outs.t_min_max(1))-toffset;
                        F_t = out.InstantaneousFluctuationStrength;

                        [Fmax, idx_max] = max(F_t);
                        idx_perc_90(1) = find(F_t>out.FS90,1,'first');
                        idx_perc_90(2) = find(F_t>out.FS90,1,'last');
                        idx_perc_50(1) = find(F_t>out.FS50,1,'first');
                        idx_perc_50(2) = find(F_t>out.FS50,1,'last');

                        %%% Time plot:
                        figure;
                        plot(t,F_t,'b','LineWidth',2); hold on; grid on

                        ylabel('Fluctuation strength (vacil)');
                        xlabel('Time (s)');
                        title(files{i_files},'interpreter','none'); 

                        xlim(outs.t_min_max);
                        YL = get(gca,'YLim');
                        plot(t(idx_max)*[1 1],[YL(1) F_t(idx_max)],'k--');
                        plot(t(idx_perc_50(1))*[1 1],[YL(1) F_t(idx_perc_50(1))],'r-');
                        plot(t(idx_perc_50(2))*[1 1],[YL(1) F_t(idx_perc_50(2))],'r-');
                        plot(t(idx_perc_90(1))*[1 1],[YL(1) F_t(idx_perc_90(1))],'r--');
                        plot(t(idx_perc_90(2))*[1 1],[YL(1) F_t(idx_perc_90(2))],'r--');

                        h(end+1) = gcf;
                        hname{end+1} = ['c' num2str(count) '-fluct-time'];

                        %%% Frequency plot:
                        Fspec_max = out.InstantaneousSpecificFluctuationStrength(idx_max,:);
                        Fspec_50  = mean(out.InstantaneousSpecificFluctuationStrength(idx_perc_50,:));
                        Fspec_90  = mean(out.InstantaneousSpecificFluctuationStrength(idx_perc_90,:));

                        f = out.barkAxis;
                        f_Tick = .5:1:23.5;
                        % f_Hz = bark2hz(f);
                        f_Tick_Hz = round( bark2hz(f_Tick) );

                        f_Tick_Hz_txt = [];
                        for i_f = 1:length(f_Tick_Hz)
                            if mod(i_f,2) == 1
                                if f_Tick_Hz(i_f) < 1000
                                    f_Tick_Hz_txt{i_f} = num2str(f_Tick_Hz(i_f));
                                else
                                    f_Tick_Hz_txt{i_f} = [num2str(round(10*f_Tick_Hz(i_f)/1000)/10) ' k'];
                                end
                            else
                                f_Tick_Hz_txt{i_f} = '';
                            end                            
                        end

                        figure;
                        plot(f,Fspec_max,'k','LineWidth',2); hold on; grid on
                        plot(f,Fspec_50,'r-');
                        plot(f,Fspec_90,'r--');

                        set(gca,'XTick',f_Tick);
                        set(gca,'XTickLabel',f_Tick_Hz_txt);
                        ylabel('Specific fluctuation strength (vacil/Bark)');
                        xlabel('Frequency (Hz)');
                        title(files{i_files},'interpreter','none'); 

                        h(end+1) = gcf;
                        hname{end+1} = ['c' num2str(count) '-fluct-freq'];

                        res(1,end+1) = out.FS50;    res_description{1,end+1} = 'F50 (vacil)';
                        res(1,end+1) = Fmax;        res_description{1,end+1} = 'Fmax (vacil)';
                        res(1,end+1) = t(idx_max);  res_description{1,end+1} = 'tmax (s)';

                    case 'Tonality_Aures1985'
                        out = Tonality_Aures1985(insig94,fs,params4input{:});

                        t = out.time+(outs.t_min_max(1))-toffset;
                        T_t = out.InstantaneousTonality;

                        [Tmax, idx_max] = max(T_t);
                        idx_above_K1 = find(T_t>out.K1);

                        %%% Time plot:
                        figure;
                        plot(t,T_t,'b','LineWidth',2); hold on; grid on

                        ylabel('Tonality (t.u.)');
                        xlabel('Time (s)');
                        title(files{i_files},'interpreter','none'); 
                        xlim(outs.t_min_max);

                        for i_K = 1:length(idx_above_K1)
                            id = idx_above_K1(i_K);
                            plot(t(id),T_t(id),'md','MarkerFaceColor','m');
                        end
                        plot(t(idx_max),T_t(idx_max),'ko','MarkerFaceColor','k');

                        h(end+1) = gcf;
                        hname{end+1} = ['c' num2str(count) '-tonality-time'];

                        res(1,end+1) = out.K10;     res_description{1,end+1} = 'K10 (a.u.)';
                        res(1,end+1) = Tmax;        res_description{1,end+1} = 'Kmax (a.u.)';
                        res(1,end+1) = t(idx_max);  res_description{1,end+1} = 'tmax (s)';
                    end
                end
                metrics_row(end+1,:) = res;

                
                % dir_out = '/home/alejandro/Documents/Documenten-ENS/01-Text/05-Doc/lx2023-04-11-FA-Psycho-test-figures/Figures-all-NEW/';
                % if ~exist(dir_out,'dir')
                %     mkdir(dir_out);
                % end
                
                % figures_dir = [fileparts(mfilename('fullpath')) filesep 'figs' filesep];
                figures_dir = [dir_sounds 'figs' filesep];
                if ~exist(figures_dir,'dir')
                    mkdir(figures_dir);
                end
    
                if do_fig1 || do_fig_raw
                    for i = 1:length(h)
                        figname_short = hname{i};
                        figname_out = [figures_dir figname_short];
    
                        saveas(h(i),figname_out, 'fig' );
                        saveas(h(i),figname_out, 'epsc'); % vectorial format
                        saveas(h(i),figname_out, 'png' );
    
                        fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
                    end
                end
                close all
            end % if bDo
            count = count+1;
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [metrics_row,metrics_description,outs] = il_get_the_metrics(insig,fs,dBFS,list_metrics)

metrics_row = [];
metrics_description = [];

bCalculated = zeros([1 length(list_metrics)]);

metric = 'LAeq';
idx_dB = find(strcmp(list_metrics,'LAeq'));

weight_freq = 'A';
weight_time = 'f'; % quite independent of the time weighting

[lvl_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS);
t = (1:length(lvl_dB))/fs;

dB_below = 25; % dB. The same criterion was applied to all event types
lvl_dB_max = max(lvl_dB);

dB_below_perc = min(55,lvl_dB_max-10); % Arbitrary level (I tried before the percentile 90)
                                       % Level reduced by 10 dB is to impose a lower level if the maximum level of the sound is soft
thres = lvl_dB_max - dB_below;


[thres,i_max] = max([thres dB_below_perc]);
if i_max == 2
    fprintf('Percentile used instead of dB_below...\n');
end

idx = find(lvl_dB>=thres);
if nargout == 0
    figure; plot(t,lvl_dB,'b-'); hold on; plot(t(idx),lvl_dB(idx),'r--'); % with time
    % t_samples = 1:length(lvl_dB); figure; plot(t_samples,lvl_dB,'b-'); hold on; plot(t_samples(idx),lvl_dB(idx),'r--'); % with Samples
end
% No dt as input means: Leq over entire sound...

switch metric
    case 'LAeq'
        LAeq = Get_Leq(lvl_dB(idx),fs); % Make sure you enter only mono signals
        metrics_row(end+1) = LAeq;
        metrics_description{end+1} = metric;
        bCalculated(idx_dB) = 1;

        T = length(idx)/fs;
        if T == 0
            disp('')
        end
        idx_T = find(strcmp(list_metrics,'T'));
        if ~isempty(idx_T)
            metrics_row(end+1) = T;
            metrics_description{end+1} = 'T';
            bCalculated(idx_T) = 1;
        end

        LAFmax = max(lvl_dB);
        idx_dBmax = find(strcmp(list_metrics,'LAFmax'));
        if ~isempty(idx_dBmax)
            metrics_row(end+1) = max(lvl_dB);
            metrics_description{end+1} = 'LAFmax';
            bCalculated(idx_dBmax) = 1;
        end

        idx_dB_SEL = find(strcmp(list_metrics,'SEL'));
        if ~isempty(idx_dB_SEL)
            metrics_row(end+1) = LAeq + 10*log10(T);
            metrics_description{end+1} = 'SEL';
            bCalculated(idx_dB_SEL) = 1;
        end

        idx_dB = find(strcmp(list_metrics,'LAFmax_min_LAeq'));
        if ~isempty(idx_dB)
            metrics_row(end+1) = LAFmax - LAeq;
            metrics_description{end+1} = 'LAFmax_min_LAeq';
            bCalculated(idx_dB) = 1;
        end
end
weight_freq = 'Z';
weight_time = 'f'; % quite independent of the time weighting
[lvl_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS);
LZeq = Get_Leq(lvl_dB(idx),fs); % Make sure you enter only mono signals

idx_dB = find(strcmp(list_metrics,'LZeq'));
if ~isempty(idx_dB)
    metrics_row(end+1) = LZeq;
    metrics_description{end+1} = 'LZeq';
    bCalculated(idx_dB) = 1;
end

idx_dB = find(strcmp(list_metrics,'Delta_Leq'));
if ~isempty(idx_dB)
    metrics_row(end+1) = LZeq-LAeq;
    metrics_description{end+1} = 'Delta_Leq';
    bCalculated(idx_dB) = 1;
end

outs.idx = idx;
outs.t_idx = idx/fs;
outs.t_min_max = [min(idx/fs) max(idx/fs)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_var2latex(metrics_row)

if ~exist('var2latex.m','file')
    metrics_row
else
    var2latex(metrics_row);
end