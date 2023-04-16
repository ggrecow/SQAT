function z = hz2bark(f)
% function z = hz2bark(f)
%
% 1. Description:
%       It converts a frequency f in Hz into a frequency z in Barks
% 
% 2. Stand-alone example:
%       f = 1000;
%       z = hz2bark(f)
%
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% 4. Reference:
%       Zwicker1961, Osses2010
%
%  See also: FREQ2AUD.m (LTFAT)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 15/08/2014
% Last update on: 15/08/2014 
% Last use on   : 03/11/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 13*atan( 0.76*(f/1000) ) + 3.5*atan( (f/(1000*7.5)).^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
