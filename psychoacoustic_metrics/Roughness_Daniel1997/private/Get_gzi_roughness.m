function gzi = Get_gzi_roughness(Chno)
% function gzi = Get_gzi_roughness(Chno)
%
% Returns the g(z) weighting of the specific roughness at the Chno
% half-Bark channel positions. The table is the revised g(z) that Dik
% Hermes derived in his 2025 implementation together with the correction
% of the upper excitation slope (see CalcGFactors.m of his routines and
% issue #47); the previous table from the 2002 implementation was tuned
% around the uncorrected slope and leaves with it.
%
% The returned gzi is the square root of the tabulated g(z): the specific
% roughness squares the product gzi*mdept*ki, so the table enters the
% specific roughness linearly, as in Hermes' formulation
% g(z)*(mdept*ki)^2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log
%
% - Author: Dik Hermes (2025), table of CalcGFactors.m
% - Author: Sergio Aguirre and Gil Felix Greco, 31.08.2026 - 
%   private function for Roughness_Daniel1997
%
%   AI disclosure: modifications performed in September 2026 were assisted 
%   by Claude Fable 5.1 and Opus 5 (Anthropic). All codes were verified by 
%   the authors.
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

gr = [ ...
     0.00    0.15
     1.00    0.26
     2.00    0.38
     3.00    0.47
     4.00    0.54
     5.00    0.65
     6.00    0.76
     7.00    0.83
     8.00    0.90
     9.00    0.98
    10.00    0.98
    11.00    0.90
    12.00    0.80
    13.00    0.70
    14.00    0.62
    15.00    0.54
    16.00    0.49
    17.00    0.43
    18.00    0.39
    19.00    0.35
    20.00    0.30
    21.00    0.30
    22.00    0.30
    23.00    0.30
    24.00    0.30];

gzi = sqrt(interp1(gr(:,1),gr(:,2),(1:Chno)/2));

end

%**************************************************************************
%
% copyright © 1999-2023 Dik Hermes
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