function f = bark2hz(z)
% function f = bark2hz(z)
%
% 1. Description:
%       Converts from Bark to Hertz. The minimum frequency value to be converted
%       is the equivalent to 10 Hz (approximation of 9.8431 Hz) = 0.0972 Bark
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       z = 5; % Bark
%       f = bark2hz(z);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/09/2014
% Last update on: 25/09/2014 
% Last use on   : 03/11/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = 1000;
k = -20:12;
f = f0*2.^(k/3); 

zt = hz2bark(f);

f = interp1(zt,f,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
