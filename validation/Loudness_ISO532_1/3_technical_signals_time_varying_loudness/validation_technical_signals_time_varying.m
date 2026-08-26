% Script validation_technical_signals_time_varying
%
% This code computes time-varying loudness from the reference signals 
%  provided by ISO 532-1:2017 - Annex B.5. using SQAT and plot the 
%  comparison against reference values 
%
% Loudness computed using:
%   OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
%   type <help Loudness_ISO532_1> for more info
%
% In order to run this code, the user needs to download the dataset of 
%  sound files from zenodo (https://doi.org/10.5281/zenodo.7933206).
%  The obtained folder called `validation_SQAT_v1_0` has to be included in 
%  the `sound_files` folder of the toolbox. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log:
% - Author: Gil Felix Greco, Braunschweig 27.02.2023
%
% - Modifided in 07.12.2024 by Gil Felix Greco - included plot with summary 
%   of differences between reference and calculated loudness
%
% - Author: Gil Felix Greco, 21.08.2026 - new version of Loudness_ISO532_1 
%   is paired with C reference code (released in v2.0). Summary of 
%   differences now include a comparison with results from prior implementation. 
%
% - Author: Gil Felix Greco, 26.08.2026: computation of percentile values
%   is corrected (i.e. now accounts for the silence at the start/end of the 
%   signals, as given by the reference spreadsheets). Percentile values of the 
%   old loudness code version (prior to v2.0) also updated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% save figs flag

save_figs = 0;

%% input files

signal_str=[ {'Test signal 14 (propeller-driven airplane).wav'},... % free-field
             {'Test signal 15 (vehicle interior 40 kmh).wav'},...   % diffuse field 
             {'Test signal 16 (hairdryer).wav'},...                 % free-field
             {'Test signal 17 (machine gun).wav'},...               % free-field
             {'Test signal 18 (hammer).wav'},...                    % free-field
             {'Test signal 19 (door creak).wav'},...                % free-field
             {'Test signal 20 (shaking coins).wav'},...             % free-field
             {'Test signal 21 (jackhammer).wav'},...                % free-field
             {'Test signal 22 (ratchet wheel (large)).wav'},...     % free-field
             {'Test signal 23 (typewriter).wav'},...                % free-field
             {'Test signal 24 (woodpecker).wav'},...                % free-field
             {'Test signal 25 (full can rattle).wav'}];             % free-field
         
%% validation

for i=14:25 

[OUT.L{i-13},OUT.RefScalar{i-13}]=compute_and_plot(i,...     % insig_num
                                                  char(signal_str(1,i-13)),... % insig name str
                                                  save_figs,['validation_time_varying_loudness_signal_' sprintf('%g',i)]... % savefig inputs
                                                  );
end 

%% summary of differences between reference and calculated loudness

% test signals
X = categorical( {'14','15','16','17','18','19','20','21','22','23','24','25'} );
X = reordercats( X, {'14','15','16','17','18','19','20','21','22','23','24','25'} );

diff_vector_Nmax = zeros(1 , length(X)); % declare for memory allocation
diff_vector_N5 = zeros(1 , length(X)); 

% create vector with loudness differences of all test signals  
for i = 1:length(X)
    diff_vector_Nmax(i) = OUT.RefScalar{i}(1,3); % max. loudness 
    diff_vector_N5(i) = OUT.RefScalar{i}(2,3); % 5 percentile loudness
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Nmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title_fig = sprintf('Max. Loudness - summary of differences between ref. and calculated values');
h = figure('Name',title_fig);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% load differences (computed with code prior to v2.0)
% Nmax is the maximum over the complete record and is therefore not affected 
% by the time window used for the percentile, so these values are unchanged.
diff_vector_Nmax_v1 = [0.240570757820347	0.203296116740026	0.721420425934014	0.436221397800264	0.435791057533194	0.510627573506348	0.703225122162740	0.268005107297007	0.158446479058027	0.854048664706799	0.277154618957530	0.528470554768854 ];

handle_a_v1 = plot(X, diff_vector_Nmax_v1, 'xb', 'Markersize', 12); % plot results computed with v.1
hold on;
handle_a = plot(X, diff_vector_Nmax, 'sb', 'Markersize', 12); % plot results with current version

ymin = -1; ymax =1;
ylim([ymin ymax]);

variable_Nmax = '$N_{\mathrm{max,SQAT}} - N_{\mathrm{max,Ref.}}$';
ylabel( [variable_Nmax  ' (sone)'], 'Interpreter','Latex');
xlabel('Test signal','Interpreter','Latex');

legend( [handle_a_v1, handle_a], ...
    {[variable_Nmax '(v.1.x)'], [variable_Nmax '(current)']}, 'Location', 'SE')

legend box on
grid on

set(gcf,'color','w');

%%%%%%%%% save fig
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];
    
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'validation_time_varying_technical_signals_loudness_difference_Nmax';
    figname_out = [figures_dir figname_short];
    
    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%% save fig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot N5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title_fig = sprintf('N5 - summary of differences between ref. and calculated values');
h = figure('Name',title_fig);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% load differences (computed with code prior to v2.0)
% Recomputed using the same time window as the current version, so that the 
% two marker sets differ only by the change in the loudness implementation.
diff_vector_N5_v1 = [ 0.0654677427423103	0.138824097883232	0.499375217756686	0.496738810406429	0.444634442553939	0.397488740479751	0.488979193365671	0.216062319420050	0.231669567381521	0.415432100498771	0.290675684403160	0.471122066114694 ];

handle_b_v1 = plot(X, diff_vector_N5_v1, 'xk', 'Markersize', 12); % plot results computed with v.1
hold on;
handle_b = plot(X, diff_vector_N5, 'sk', 'Markersize', 12); % plot results with current version

ylim([ymin ymax]);

variable_N5 = '$N_{\mathrm{5\%,SQAT}} - N_{\mathrm{5\%,Ref.}}$';
ylabel([variable_N5 ' (sone)'], 'Interpreter', 'Latex');
xlabel('Test signal','Interpreter','Latex');

legend([handle_b_v1, handle_b], ...
    {[variable_N5 '(v.1.x)'], [variable_N5 '(current)']}, 'Location', 'SE')

grid on

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'validation_time_varying_technical_signals_loudness_difference_N5';
    figname_out = [figures_dir figname_short];
    
    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% function (compute loudness and plot comparison

function [OUT,table] = compute_and_plot(insig_num,fname_insig,save_figs,tag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this function computes the loudness using SQAT and plot the comparison
% against the reference values from the ISO 532-1:2017 - Annex B.5. 
%
% INPUTS:
%   insig_num : scalar
%       number of the reference signal to be tested
%
%   insig : string
%       name of the reference signals
%
%   save_figs : scalar
%       1 to save; <else> dont save figures 
%
%   tag : string
%       tag with the name of the figures to be saved
%
% OUTPUTS:
%   OUT : struct
%       contain all outputs from the computed loudness
%
%   table : matrix containing scalar values of Nmax (1st row) and N5 (2nd row) 
%           1st column=reference
%           2nd column=computed by SQAT
%           3rd column=difference (SQAT minus ref.)
%           4th column=relative percentage difference (SQAT minus ref.)
%
% Author: Gil Felix Greco, Braunschweig 27.02.2023
% modifided in 07.12.2024 by Gil Felix Greco - included plot with summary of differences between reference and calculated loudness
% modified in 26.08.2026 by Gil Felix Greco - N5 computed over the reference time window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signals from ISO 532-2:2017

dir_analysis_name = '3_technical_signals_time varying loudness';
dir_out = [fileparts(mfilename('fullpath')) filesep];
  
SQAT_version=1; % v1.0
dir_sounds = get_dir_validation_sounds('Loudness_ISO532_1',SQAT_version);
dir_ref_values = get_dir_reference_values('Loudness_ISO532_1',dir_analysis_name);

% calibration signal provided in the Annex C of the ISO 532-1:2017

[RefSignal,~]=audioread([dir_sounds 'calibration signal sine 1kHz 60dB.wav']);

% Test signal provided in the Annex B.5 of the ISO 532-1:2017

[signal,fs]=audioread([dir_sounds fname_insig]);

%% calibrated .wav signal

[ycal]=calibrate(signal,RefSignal,60); 

%% Loudness calculation using SQAT

if insig_num==15 % only signal to be tested considering difuse field
    
OUT = Loudness_ISO532_1( ycal,   fs,...   % input signal and sampling freq.
                                  1,...   % field; free field = 0; diffuse field = 1;
                                  2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                                  0,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                  0);     % show results, 'false' (disable, default value) or 'true' (enable)  
else % free-field
    
OUT = Loudness_ISO532_1( ycal,   fs,...   % input signal and sampling freq.
                                  0,...   % field; free field = 0; diffuse field = 1;
                                  2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                                  0,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                  0);     % show results, 'false' (disable, default value) or 'true' (enable)  
end

%% calculate difference from reference values given by ISO 532-1:2017

% time window used by ISO 532-1:2017 - Annex B.5 for the N5 statistics
%
% The header of the reference spreadsheets states that the test signals contain
% 100 ms of silence at the beginning and 500 ms at the end, and that N5 is
% obtained disregarding these intervals. The leading 100 ms is indeed skipped 
% for all signals (the reference range always starts at row 61, t = 0.100 s), 
% but the trailing 500 ms is NOT skipped for the signals 18, 22 and 24:
%   - signal 18: the reference range B61:B1093 ends exactly at the last data row;
%   - signals 22 and 24: the reference ranges (B61:B1494 and B61:B4492) extend 
%     beyond the last data row (935 and 1210), and PERCENTILE ignores the empty 
%     cells, so the effective upper limit is again the end of the record.
% These are stale ranges in the reference spreadsheet, but since the tabulated 
% reference values inherit them, the same window has to be used here in order to
% reproduce the published N5 values. Applying a uniform 0.5 s trailing trim to
% the signals 18, 22 and 24 introduces spurious deviations of about +0.42, 
% +0.14 and +0.09 sone.

trailing_trim = [0.5 0.5 0.5 0.5 0.0 0.5 0.5 0.5 0.0 0.5 0.0 0.5]; % (s), signals 14 to 25

t_start = 0.1;                                          % (s) skip leading silence
t_end   = OUT.time(end) - trailing_trim(insig_num-13);  % (s) skip trailing silence

tol = 1e-9;  % avoid dropping the boundary samples due to round-off
idx_valid = ( OUT.time >= t_start - tol ) & ( OUT.time <= t_end + tol );

N_valid = OUT.InstantaneousLoudness(idx_valid);

% overwrite N5 computed inside Loudness_ISO532_1.m 
OUT.N5 = get_exceeded_value(N_valid, 5);

% reference values provided by ISO 532-1:2017 for signals 14 to 25
reference_Nmax = [22.640 9.606 38.536 11.211 12.647 10.882 14.880 9.719 8.906 11.186 9.275 7.259];   
reference_N5 = [18.081 8.755 36.888 9.722 10.392 10.009 13.093 9.00 8.158 10.421 8.519 5.761]; 

reference_Nmax=reference_Nmax(insig_num-13);
reference_N5=reference_N5(insig_num-13);

% compute difference (SQAT minus ref.)
difference_Nmax = OUT.Nmax - reference_Nmax; % max loudness over time
difference_N5 = OUT.N5 - reference_N5; % 5% percentile loudness

% compute relative percentage difference (SQAT minus ref.)
percentage_difference_Nmax=( (OUT.Nmax-reference_Nmax)/reference_Nmax )*100; % max loudness over time
percentage_difference_N5=( (OUT.N5-reference_N5)/reference_N5 )*100; % 5% percentile loudness

table = [reference_Nmax, OUT.Nmax, difference_Nmax, percentage_difference_Nmax;
             reference_N5, OUT.N5, difference_N5, percentage_difference_N5 ];

%% plot results (total loudness over time)

h = figure('Name',['Loudness - signal ' sprintf('%g',insig_num)]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

reference = []; % to be loaded in the next line...
fname = sprintf('%sreference_values_ISO532_1_2017_signal_%g.mat', dir_ref_values, insig_num);
load(fname); % load reference vectors

% load([pwd '\reference_values\' 'reference_values_ISO532_1_2017_signal_' sprintf('%g',insig_num) '.mat']); % load reference vectors

% reference values

% plot( reference(:,1), reference(:,2),'b','Linewidth',0.5);hold on; % ref N

a=plot( reference(:,1), reference(:,3),'r:','Color',[1 0 0],'Linewidth',1);hold on; % ref N_min (5% tolerance)
plot( reference(:,1), reference(:,4),'r:','Color',[1 0 0],'Linewidth',1); % ref N_max (5% tolerance)

b=plot( reference(:,1), reference(:,5),'Color',[1 0 0],'Linewidth',0.5); % ref N_min (10% tolerance)
plot( reference(:,1), reference(:,6),'Color',[1 0 0],'Linewidth',0.5); % ref N_max (10% tolerance)

% SQAT values
c=plot( OUT.time, OUT.InstantaneousLoudness,'k','Linewidth',1); % calculated specific loudness

leg=legend([a,b,c],'5\% tolerance','10\% tolerance','SQAT','Location','Best');
set(leg,'Orientation','horizontal','Location','northoutside');

ylabel('Loudness, $N$ (sone)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 
grid off

axis([ 0 max(OUT.time) 0 max(reference(:,6)) ])

set(gcf,'color','w');

if save_figs==1
    if ~exist(dir_out,'dir')
        mkdir(dir_out);
    end
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = tag;
    figname_out = [figures_dir figname_short];
    
%     saveas(gcf,figname_out, 'fig');
%     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

end
