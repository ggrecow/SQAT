function params = psychoacoustic_metrics_get_defaults(model_name)
% function params = psychoacoustic_metrics_get_defaults(model_name)
%
% Author: Alejandro Osses
% Modified: Mike Lotinga, 12.06.2025 - updated to include Widmann PA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch model_name
    case {'FluctuationStrength_Osses2016','FluctuationStrength_Osses2016_from_wavfile'}
        method = 1;
        method_description = '0 = stationary method (win size=length of sound); 1 = time varying (window size = 2s)';
        time_skip = 0;
        time_skip_description = 'time_skip, in seconds for statistical calculations';
        show = 1;
        show_description = 'Plots the outputs (optional parameter)';
        
        params.method = method;
        params.method_description = method_description;
        params.time_skip = time_skip;
        params.time_skip_description = time_skip_description;
        params.show = show;
        params.show_description = show_description;
                                                       
    case {'Loudness_ISO532_1','Loudness_ISO532_1_from_wavfile'}
        field = 0; 
        field_description = '0 = free field; 1 = diffuse field';
        method = 2;
        method_description = '1 = stationary method; 2 time varying'; 
        time_skip = 0.5; % 
        time_skip_description = 'time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations';
        show = 1;
        show_description = 'Plots the outputs (optional parameter)';
        
        params.field = field;
        params.field_description = field_description;
        params.method = method;
        params.method_description = method_description;
        params.time_skip = time_skip;
        params.time_skip_description = time_skip_description;
        params.show = show;
        params.show_description = show_description;

    case {'Sharpness_DIN45692','Sharpness_DIN45692_from_loudness'}
        field = 0; 
        field_description = '0 = free field; 1 = diffuse field for loudness';
        method = 2;
        method_description = '1 = stationary method; 2 time varying for loudness'; 
        time_skip = 0.5; % 
        time_skip_description = 'time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations';
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
        params.time_skip = time_skip;
        params.time_skip_description = time_skip_description;
        params.show_loudness = show_loudness;
        params.show_loudness_description = show_loudness_description;
        params.show = show;
        params.show_description = show_description;
        params.weight_type = weight_type;
        params.weight_type_description = weight_type_description;
        
    case {'Roughness_Daniel1997','Roughness_Daniel1997_from_wavfile'}
        params.time_skip = 0;
        params.time_skip_description = 'time_skip';
        params.show = 0;
        params.show_description = 'Plots the outputs (optional parameter)';
        
    case {'Tonality_Aures1985','Tonality_Aures1985_from_wavfile'}
        params.Loudness_field = 0; 
        params.Loudness_field_description = 'Loudness: 0 = free field; 1 = diffuse field';
        params.time_skip = 0;
        params.time_skip_description = 'time_skip';
        params.show = 0;
        params.show_description = 'Plots the outputs (optional parameter)';
        
    case {'PsychoacousticAnnoyance_Di2016','PsychoacousticAnnoyance_More2010','PsychoacousticAnnoyance_Widmann1992'}
        params.Loudness_field = 0; 
        params.Loudness_field_description = 'Loudness: 0 = free field; 1 = diffuse field';
        params.time_skip = 0.2;
        params.time_skip_description = 'time_skip';
        params.showPA = 0; % show results of PA, 'false' (disable, default value) or 'true' (enable)
        params.showPA_description = 'Plots the outputs for psychoacoustic annoyance (optional parameter)';
        params.show = 0;
        params.show_description = 'Plots the outputs (optional parameter)';
        
end