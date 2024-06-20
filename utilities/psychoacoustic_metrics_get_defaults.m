function params = psychoacoustic_metrics_get_defaults(model_name)
% function params = psychoacoustic_metrics_get_defaults_ML(model_name)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch model_name
    case {'FluctuationStrength_Osses2016','FluctuationStrength_Osses2016_from_wavfile'}
        method = 1;
        method_description = '0 = stationary method (win size=length of sound); 1 = time varying (window size = 2s)';
        start_skip = 0;
        start_skip_description = 'start_skip, in seconds for statistical calculations';
        end_skip = 0;
        end_skip_description = 'end_skip, in seconds for statistical calculations';
        show = 1;
        show_description = 'Plots the outputs (optional parameter)';
        
        params.method = method;
        params.method_description = method_description;
        params.start_skip = start_skip;
        params.start_skip_description = start_skip_description;
        params.end_skip = end_skip;
        params.end_skip_description = end_skip_description;
        params.show = show;
        params.show_description = show_description;
                                                       
    case {'Loudness_ISO532_1','Loudness_ISO532_1_from_wavfile'}
        field = 0; 
        field_description = '0 = free field; 1 = diffuse field';
        method = 2;
        method_description = '1 = stationary method; 2 time varying'; 
        start_skip = 0.5; % 
        start_skip_description = 'start_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations';
        end_skip = 0; % 
        end_skip_description = 'end_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations';
        show = 1;
        show_description = 'Plots the outputs (optional parameter)';
        
        params.field = field;
        params.field_description = field_description;
        params.method = method;
        params.method_description = method_description;
        params.start_skip = start_skip;
        params.start_skip_description = start_skip_description;
        params.end_skip = end_skip;
        params.end_skip_description = end_skip_description;
        params.show = show;
        params.show_description = show_description;

    case {'Sharpness_DIN45692','Sharpness_DIN45692_from_loudness'}
        field = 0; 
        field_description = '0 = free field; 1 = diffuse field for loudness';
        method = 2;
        method_description = '1 = stationary method; 2 time varying for loudness'; 
        start_skip = 0.5; % 
        start_skip_description = 'start_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations';
        end_skip = 0; % 
        end_skip_description = 'end_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations';
        show_loudness = 0;
        show_loudness_description = 'Plots the loudness outputs (optional parameter)';
        show = 1;
        show_description = 'Plots the outputs (optional parameter)';
        weight_type = 'DIN45692';
        weight_type_description = 'Type of sharpness models. Options: ''DIN45692'',''aures'',''bismarck''';
        
        params.field = field;
        params.field_description = field_description;
        params.method = method;
        params.method_description = method_description;
        params.start_skip = start_skip;
        params.start_skip_description = start_skip_description;
        params.end_skip = end_skip;
        params.end_skip_description = end_skip_description;
        params.show_loudness = show_loudness;
        params.show_loudness_description = show_loudness_description;
        params.show = show;
        params.show_description = show_description;
        params.weight_type = weight_type;
        params.weight_type_description = weight_type_description;
        
    case {'Roughness_Daniel1997','Roughness_Daniel1997_from_wavfile'}
        params.start_skip = 0;
        params.start_skip_description = 'start_skip';
        params.end_skip = 0;
        params.end_skip_description = 'end_skip';
        params.show = 0;
        params.show_description = 'Plots the outputs (optional parameter)';
        
    case {'Tonality_Aures1985','Tonality_Aures1985_from_wavfile'}
        params.Loudness_field = 0; 
        params.Loudness_field_description = 'Loudness: 0 = free field; 1 = diffuse field';
        params.start_skip = 0;
        params.start_skip_description = 'start_skip';
        params.end_skip = 0;
        params.end_skip_description = 'end_skip';
        params.show = 0;
        params.show_description = 'Plots the outputs (optional parameter)';
        
    case {'PsychoacousticAnnoyance_Di2016','PsychoacousticAnnoyance_More2010','PsychoacousticAnnoyance_Zwicker1999'}
        params.Loudness_field = 0; 
        params.Loudness_field_description = 'Loudness: 0 = free field; 1 = diffuse field';
        params.start_skip = 0.2;
        params.start_skip_description = 'end_skip';
        params.end_skip = 0;
        params.end_skip_description = 'end_skip';
        params.showPA = 0; % show results of PA, 'false' (disable, default value) or 'true' (enable)
        params.showPA_description = 'Plots the outputs for psychoacoustic annoyance (optional parameter)';
        params.show = 0;
        params.show_description = 'Plots the outputs (optional parameter)';
        
end