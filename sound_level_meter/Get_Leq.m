function Leq = Get_Leq(levels,fs,dt,framelen_s)
% function Leq = Get_Leq(levels,fs,dt,framelen_s)
% function Leq = Get_Leq(levels)
% 
% 1. Description:
%
% 2. Stand-alone example:
%       % Obtain a one-second Leq:
%       lvls = Do_SLM(insig,fs,'A','f',100); % A-weighted, 'fast'
%       dt = 1; %s
%       Leq = Get_Leq(lvls,fs,dt); % Make sure you enter only mono signals
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also DO_SLM, Get_levels_SPL.m (inline function: il_get_Leq)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 13/07/2016
% Last update on: 13/07/2016 
% Last use on   : 24/01/2017 
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

if nargin > 2
    if nargin < 4
        framelen_s = dt;
    end
    
    H = round(dt*fs);
    N = round(framelen_s*fs);
    
    if size(levels,2) ~= 1
        error('Calculation validated only for mono inputs');
    end
    Ov = N-H; % Overlap
    levels = buffer(levels,N,Ov,'nodelay');
end

for i = 1:size(levels,2)
    idx   = find(~isnan(levels(:,i)));
    lvls  = levels(idx,i);
    lvls = (10.^(lvls/10));
    Leq(i)= 10*log10( mean(lvls) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
