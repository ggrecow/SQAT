function gzi = Get_gzi_fluctuation(Chno)
% function gzi = Get_gzi_fluctuation(Chno)
%
% Returns gzi parameters using the specified number of channels.

Chstep = 0.5;
    
% Hz:   100 250   519   717 926 1084 1255 1465 1571   1972 2730 4189   15550
g0 = [0,  1,  2.5,  4.9,6.5,  8,   9,  10,  11,  11.5,  13,  15,  17.5,   24;
      1,  1,  1  ,  1  ,1  ,  1,   1,   1,   1,   1  ,   1, 0.9,   0.7, 0.5];
g0 = transpose(g0);

gzi = interp1(g0(:,1),g0(:,2),(1:Chno)*Chstep);
gzi(isnan(gzi)) = g0(end,2); % 0