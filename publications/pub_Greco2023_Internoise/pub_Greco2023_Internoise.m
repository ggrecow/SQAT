function pub_Greco2023_Internoise
% function pub_Greco2023_Internoise
%
% Generates the figures presented in the contribution to be 
%   presented at Inter-noise in August 2023.
%
% This function is a wrapper to the validation codes from ../validation/ by
%   Gil Felix Greco.
%
% Reference:
% Felix Greco, G., Merino-Martinez, R., Osses, A., and Langer, S. (2023). 
%   SQAT: a MATLAB-based toolbox for quantitative sound quality analysis. 
%   Inter-noise.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

h = [];
hname = [];

do_fig2a   = 0; % Validation Loudness_ISO532_1.m
do_fig2b   = 0; % Validation Loudness_ISO532_1.m
do_fig3    = 0; % Validation Sharpness_DIN45692.m, both panels at the same time
do_fig4a   = 0; % Validation Roughness_Daniel1997.m
do_fig4b   = 0; % Validation Roughness_Daniel1997.m
do_fig5a   = 0; % Validation FluctuationStrength_Osses2016.m
do_fig5b   = 0; % Validation FluctuationStrength_Osses2016.m
% do_fig6a   = 0; % Validation Tonality_Aures1985.m (waveform)
do_fig6b   = 0; % Validation Tonality_Aures1985.m

list_figures = {'do_fig2a'; 'do_fig2b'; 'do_fig3'; 'do_fig4a'; 'do_fig4b';'do_fig5a'; 'do_fig6b'};
for i = 1:length(list_figures)
    fprintf('Enter %.0f to %s\n',i,list_figures{i});
end
bInput = input('Which figure do you want to obtain (enter the corresponding number)?: ');

switch list_figures{bInput}
    case 'do_fig2a'
        do_fig2a = 1;
    case 'do_fig2b'
        do_fig2b = 1;
    case 'do_fig3'
        do_fig3 = 1;
    case 'do_fig4a'
        do_fig4a = 1;
    case 'do_fig4b'
        do_fig4b = 1;
    case 'do_fig5a'
        do_fig5a = 1;
    case 'do_fig5b'
        do_fig5b = 1;
    case 'do_fig6b'
        do_fig6b = 1;
end
dir_validation = [basepath_SQAT 'validation' filesep];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig2a
    
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'Loudness_ISO532_1' filesep '1_synthetic_signals_stationary_loudness' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    OUT = il_get_validation_stationary_loudness_synthetic_signals;
    
    cd(curr_dir);
    
    idx = 1; % a priori knowledge, the first handle
    h(end+1) = idx;
    figure(h(end)); % sets the focus to fig. 1
    
    ref_value = OUT.RefScalar{idx}(idx);
    SQAT_value = OUT.L{idx}(idx).Loudness;
    
    fprintf('The results are:\n')
    fprintf('\t Ref: N=%.2f sone\n',ref_value);
    fprintf('\t SQAT: N=%.2f sone\n',SQAT_value);
    
    N_figs = length(OUT.L);
    idx_figs = 1:N_figs;
    idx_figs(idx) = [];
    for i = 1:length(idx_figs)
        % Closes all figures that are not used
        close(figure(idx_figs(i)));
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig2b
    
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'Loudness_ISO532_1' filesep '3_technical_signals_time_varying_loudness' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    OUT = il_get_validation_technical_signals_time_varying;
    
    cd(curr_dir);
    
    idx = 1; % a priori knowledge, the first handle
    h(end+1) = idx;
    figure(h(end)); % sets the focus to fig. 1
    
    ref_value = OUT.RefScalar{idx}(idx);
    SQAT_value = OUT.L{idx}(idx).Nmax;
    
    fprintf('The results are:\n')
    fprintf('\t Ref: Nmax=%.2f sone\n',ref_value);
    fprintf('\t SQAT: Nmax=%.2f sone\n',SQAT_value);
    
    N_figs = length(OUT.L);
    idx_figs = 1:N_figs;
    idx_figs(idx) = [];
    for i = 1:length(idx_figs)
        % Closes all figures that are not used
        close(figure(idx_figs(i)));
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig3
    
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'Sharpness_DIN45692' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    OUT = il_get_validation_Sharpness_DIN45692_NB_and_BB_signals;
    
    cd(curr_dir);
    
    idx = 1; % a priori knowledge, the first handle
    h(end+1) = idx;
    figure(h(end)); % sets the focus to fig. 1
    
    idx = 2; % a priori knowledge, the second handle
    h(end+1) = idx;
    figure(h(end)); % sets the focus to fig. 1
    
    N_figs = 4; % a priori knowledge
    idx_figs = 1:N_figs;
    idx_figs([h(1) h(2)]) = [];
    for i = 1:length(idx_figs)
        % Closes all figures that are not used
        close(figure(idx_figs(i)));
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig4a
    
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'Roughness_Daniel1997' filesep '2_AM_modulation_depth' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    
    il_get_run_verification_roughness_modulation_depth;
    
    cd(curr_dir);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig4b
    
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'Roughness_Daniel1997' filesep '1_AM_modulation_freq' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    
    OUT = il_get_run_validation_roughness_fmod;
    
    cd(curr_dir);
    
    idx = 3; % a priori knowledge, plot of 1 kHz and 8 kHz
    h(end+1) = idx;
    figure(h(end));
        
    N_figs = 4; % a priori knowledge
    idx_figs = 1:N_figs;
    idx_figs([h(1)]) = [];
    for i = 1:length(idx_figs)
        % Closes all figures that are not used
        close(figure(idx_figs(i)));
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig5a
 
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'FluctuationStrength_Osses2016' filesep '1_AM_tones_fmod' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    
    il_get_run_validation_FS_fmod;
    
    cd(curr_dir);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig5b
 
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'FluctuationStrength_Osses2016' filesep '2_AM_BBN_fmod' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    
    il_get_run_validation_FS_AM_BBN_fmod;
    
    cd(curr_dir);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_fig6b
    
    curr_dir = [cd filesep];
    dir_this_validation = [dir_validation 'Tonality_Aures1985' filesep];
    cd(dir_this_validation);
    
    close all; % to make sure that Fig. 1 is assigned to handle 1
    
    il_get_validation_signal_to_noise_ratio;
    
    cd(curr_dir);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = il_get_validation_stationary_loudness_synthetic_signals
% This inline function allows to keep the workspace separate from that of
%   the validation

OUT = []; % trick to visibly debug the variable OUT
validation_stationary_loudness_synthetic_signals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = il_get_validation_technical_signals_time_varying

OUT = []; % trick to visibly debug the variable OUT
validation_technical_signals_time_varying;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = il_get_validation_Sharpness_DIN45692_NB_and_BB_signals
% This inline function allows to keep the workspace separate from that of
%   the validation

OUT = []; % trick to visibly debug the variable OUT
validation_Sharpness_DIN45692_narrowband_and_broadband_signals;

OUT.L_narrow = L_narrow;
OUT.s_narrow = s_narrow;
OUT.s_narrow_aures = s_narrow_aures;
OUT.s_narrow_bismarck = s_narrow_bismarck;

OUT.L_broad = L_broad;
OUT.s_broad = s_broad;
OUT.s_broad_aures = s_broad_aures;
OUT.s_broad_bismarck = s_broad_bismarck;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = il_get_run_validation_roughness_fmod

OUT = [];
run_validation_roughness_fmod;

OUT.res_125hz = res_125hz;
OUT.res_250hz = res_250hz;
OUT.res_500hz = res_500hz;
OUT.res_1khz = res_1khz;
OUT.res_4khz = res_4khz;
OUT.res_8khz = res_8khz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = il_get_run_verification_roughness_modulation_depth

OUT = [];
run_verification_roughness_modulation_depth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_get_run_validation_FS_fmod

run_validation_FS_fmod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_get_run_validation_FS_AM_BBN_fmod

run_validation_FS_AM_BBN_fmod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_get_validation_signal_to_noise_ratio

validation_signal_to_noise_ratio;
