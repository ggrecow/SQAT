function dir_ref_values = get_dir_reference_values(psychoacoustic_model,dir_analysis)
% function dir_ref_value = get_dir_reference_values(psychoacoustic_model,dir_analysis)
%
% It returns the main (base) path of the toolbox.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    dir_analysis = '';
end

dir_validation = get_dir_validation(psychoacoustic_model);
switch psychoacoustic_model
    case 'Loudness_ISO532_1'
        switch dir_analysis
            case '1_synthetic_signals_stationary_loudness'
                dir_ref_values = [dir_validation dir_analysis filesep 'reference_values' filesep];
            case '2_synthetic_signals_time_varying_loudness'
                dir_ref_values = [dir_validation dir_analysis filesep 'reference_values' filesep];
            case {'3_technical_signals_time_varying_loudness','3_technical_signals_time varying loudness'}
                dir_analysis = '3_technical_signals_time_varying_loudness'; % making sure there are no spaces in the middle
                dir_ref_values = [dir_validation dir_analysis filesep 'reference_values' filesep];
        end
        
    case 'Roughness_Daniel1997'
        switch dir_analysis
            case '1_AM_modulation_freq'
                dir_ref_values = [dir_validation dir_analysis filesep 'reference_values' filesep];
            case '3_FM_modulation_depth'
                dir_ref_values = [dir_validation dir_analysis filesep 'reference_values' filesep];
            case '4_FM_level'
                dir_ref_values = [dir_validation dir_analysis filesep 'reference_values' filesep];
        end
end

if ~exist(dir_ref_values,'dir')
    dir_ref_values = '';
    fprintf('%s.m: The directory %s was not found. Maybe the folder for ref. values will still be found....\n',mfilename,dir_sound);
end