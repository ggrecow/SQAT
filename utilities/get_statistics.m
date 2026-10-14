function OUT = get_statistics(input, metric )
% function OUT = get_statistics(input, metric )
%
% This function computes several statistical indicator from
% an input vector [Nx1], [Nx2] (stereo case), or [Nx3] (stereo including comb. binaural),
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log
% 
% Last checked; Gil Felix Greco, Braunschweig 19.02.2025
%
% Modified: Mike Lotinga, 12.06.2025 - updated to include Widmann PA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file and code is part of the SQAT toolbox and
% is subject to the MIT license in its entirety, as detailed in the license
% text reproduced at the end of this file (see also <licenses/mit-license.txt> 
% file in the SQAT repository root.
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    case 'PsychoacousticAnnoyance_Widmann1992'
        var_string = 'PA';
    case 'Loudness_ECMA418_2'
        var_string = 'N';
    case 'Tonality_ECMA418_2'
        var_string = 'T';
    case 'Roughness_ECMA418_2'
        var_string = 'R';         
    case 'FluctuationStrength_ECMA418_2'
        var_string = 'FS';
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

%**************************************************************************
%
% copyright © 2025 Gil Felix Greco
% 
% This software is licenced under the MIT license:
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the “Software”), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, a
% nd/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%
%**************************************************************************