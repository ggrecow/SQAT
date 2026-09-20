function gzi = Get_gzi_fluctuation(Chno)
% function gzi = Get_gzi_fluctuation(Chno)
%
% Returns gzi parameters using the specified number of channels.
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
Chstep = 0.5;
    
% Hz:   100 250   519   717 926 1084 1255 1465 1571   1972 2730 4189   15550
g0 = [0,  1,  2.5,  4.9,6.5,  8,   9,  10,  11,  11.5,  13,  15,  17.5,   24;
      1,  1,  1  ,  1  ,1  ,  1,   1,   1,   1,   1  ,   1, 0.9,   0.7, 0.5];
g0 = transpose(g0);

gzi = interp1(g0(:,1),g0(:,2),(1:Chno)*Chstep);
gzi(isnan(gzi)) = g0(end,2); % 0