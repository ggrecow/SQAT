function OUT = get_statistics(input, metric )
% function OUT = get_statistics(input, metric )
%
% This function computes several statistical indicator from
% an input vector [Nx1], [Nx2] (stereo case), or [Nx3] (stereo including comb. binaural),
%
% Last checked; Gil Felix Greco, Braunschweig 19.02.2025
%%%%%%%%%%%%%%%%%%%%%

switch metric
    case 'Loudness_ISO532_1'
        var_string = 'N';
    case 'Sharpness_DIN45692'
        var_string = 'S';
    case 'Roughness_Daniel1997'
        var_string = 'R';
    case 'FluctuationStrength_Osses2016'
        var_string = 'FS';
    case 'Tonality_Aures1985'
        var_string = 'K';
    case 'PsychoacousticAnnoyance_Di2016'
        var_string = 'PA';
    case 'PsychoacousticAnnoyance_More2010'
        var_string = 'PA';
    case 'PsychoacousticAnnoyance_Zwicker1999'
        var_string = 'PA';
    case 'Loudness_ECMA418_2'
        var_string = 'N';
    case 'Tonality_ECMA418_2'
        var_string = 'T';
    case 'Roughness_ECMA418_2'
        var_string = 'R';         
end

string_vector = { 'max'; 'min'; 'mean'; 'std'; ...
    '1'; '2'; '3'; '4'; '5'; ...
    '10'; '20'; '30'; '40'; '50'; ...
    '60'; '70'; '80'; '90'; '95';} ;

for k = 1:length(string_vector)

    if k==1
        temp_value = max(input);
    elseif k==2
        temp_value = min(input);
    elseif k==3
        temp_value = mean(input);
    elseif k==4
        temp_value = std(input);
    elseif k==5
        temp_value = get_exceeded_value(input,1);
    elseif k==6
        temp_value = get_exceeded_value(input,2);
    elseif k==7
        temp_value = get_exceeded_value(input,3);
    elseif k==8
        temp_value = get_exceeded_value(input,4);
    elseif k==9
        temp_value = get_exceeded_value(input,5);
    elseif k==10
        temp_value = get_exceeded_value(input,10);
    elseif k==11
        temp_value = get_exceeded_value(input,20);
    elseif k==12
        temp_value = get_exceeded_value(input,30);
    elseif k==13
        temp_value = get_exceeded_value(input,40);
    elseif k==14
        temp_value = median(input);
    elseif k==15
        temp_value = get_exceeded_value(input,60);
    elseif k==16
        temp_value = get_exceeded_value(input,70);
    elseif k==17
        temp_value = get_exceeded_value(input,80);
    elseif k==18
        temp_value = get_exceeded_value(input,90);
    elseif k==19
        temp_value = get_exceeded_value(input,95);
    end

    temp_varName =  strcat( var_string, char(string_vector{k}) ) ;
    OUT.( temp_varName ) = temp_value;

end
end % end function
