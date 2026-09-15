function gain = from_dB(gain_dB,divisor)
% function gain = from_dB(gain_dB,divisor)
%
% 1. Description:
%       From_dB: Convert decibels to voltage gain (if div = 20, default).
%       gain = From_dB(gain_dB)
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

if nargin < 2
    divisor = 20;
end

gain = 10 .^ (gain_dB / divisor);