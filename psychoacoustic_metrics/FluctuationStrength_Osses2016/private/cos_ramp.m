function ramp = cos_ramp(sig_len, fs, attack, release)
% function ramp = cos_ramp(sig_len, fs, attack, release)
%
% 1. Description:
%       Generates cosine ramp with raising/falling ramp at the begninning/end 
%       (ramps have similar lengths)
%  
%   INPUTS:
%       fs     :  Sampling frequency 
%       sig_len:  Length of signal (in samples)
%       attack :  Attack time of raising envelope (in ms)
%       release:  Release time of falling envelope (optional)
%
%   OUTPUT:
%       ramp:     Resulting envelope 
% 
% 2. Stand-alone example:
%       cos_ramp;
% 
% Author        : Matthias Milczynski
% Created in    : 2008-2012
% Downloaded on : 06/09/2012
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Last update on: 18/08/2014 
% Last use on   : 19/08/2015 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   sig_len = 44100; % only for demo
end

if nargin < 2
    fs = 44100;
end

if nargin < 3
    attack = 25;
end

if nargin < 4
    release = 25;
end

attack  = round(fs*attack/1000);
ramp = ones(1, sig_len);
if attack >= sig_len
    fprintf(1,'Attack: %i, Signal: %i\n', attack, length(ramp) )
    attack = round( ( sig_len - 1 ) / 2 );
    fprintf(1,'Attack changed to: %i, %.4f (msec)\n', attack, 1000/fs*attack )
end

switch nargin
    case 4 % attack and release time are different
        for i=1:attack
            ramp(i)=ramp(i)*0.5*(1-cos(pi*i/attack));
        end
        release = round(fs*release/1000);    
        for i=(sig_len-release):sig_len
            ramp(i)=ramp(i)*0.5*(1-cos(pi*(i-sig_len)/release));
        end
    otherwise % attack and release time are the same
        for i=1:attack
            ramp(i)=ramp(i)*0.5*(1-cos(pi*i/attack));
        end
        for i=(sig_len-attack):sig_len
            ramp(i)=ramp(i)*0.5*(1-cos(pi*(i-sig_len)/attack));
        end
end

if attack == 0 | release == 0
    ramp(end) = 1;
end

if nargout == 0
    t = ( 1:sig_len )/fs;
    figure
    plot(t, ramp);
    xlabel('time [s]')
    ylabel('Amplitude')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
