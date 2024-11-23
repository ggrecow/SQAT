function gain = from_dB(gain_dB,divisor)
% function gain = from_dB(gain_dB,divisor)
%
% 1. Description:
%       From_dB: Convert decibels to voltage gain (if div = 20, default).
%       gain = From_dB(gain_dB)

if nargin < 2
    divisor = 20;
end

gain = 10 .^ (gain_dB / divisor);