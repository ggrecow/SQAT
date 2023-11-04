function f = bark2hz_local(z)
% function f = bark2hz_local(z)
%
% 1. Description:
%       Converts from Bark to Hertz. The minimum frequency value to be converted
%       is the equivalent to 10 Hz (approximation of 9.8431 Hz) = 0.0972 Bark.
%       Named bark2hz_local to avoid the shadowing of the 'new' Audio Toolbox 
%       function with the same name.
%
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

zt = hz2bark_local(f);

f = interp1(zt,f,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
