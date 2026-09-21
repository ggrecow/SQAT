function c = SQAT_GUI_colormap_heat(m)
% function c = SQAT_GUI_colormap_heat(m)
%
%   Colour scale used by SQAT_GUI for spectrograms and maps: black, blue,
%   magenta, red, yellow and white, linearly interpolated. The six colours
%   sit at the same positions as in a colour map of the ITA-Toolbox, so both
%   give the same scale; this file is an independent implementation.
%
% INPUT ARGUMENTS
%   m : number of colours (default 256)
%
% OUTPUTS
%   c : [m x 3] RGB values in [0, 1]
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September 2026 assisted
% by Claude Opus 5 (Anthropic). All codes were verified by
% the authors.
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
    m = 256;
end
anchors = [0 0 0      % black
           0 0 1      % blue
           1 0 1      % magenta
           1 0 0      % red
           1 1 0      % yellow
           1 1 1];    % white
positions = [1, floor((1:4) * m / 5) + 1, m];
c = interp1(positions, anchors, (1:m)');
end
