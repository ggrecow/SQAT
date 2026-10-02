function SQAT_GUI_paint(w, style)
% function SQAT_GUI_paint(w, style)
%
%   The light theme of the interface in white: the window, its grids,
%   panels and tabs, the room around the plots, the buttons and the menus
%   take white over the light grey of MATLAB; the Run button keeps its
%   green, and the tables their stripes. The dark theme gives them back the
%   colours of the theme. Called after a window or a part of it is built,
%   since a new component starts with the colour of the theme.
%
% INPUT ARGUMENTS
%   w : a uifigure, or a container in one
%   style : 'light' or 'dark'
%
% Author: Sergio Aguirre and Gil Felix Greco, October 2026
%
% AI disclosure: code development in October 2026 assisted by
% Claude Opus 5.5 (Anthropic). All codes were verified by the
% authors.
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

white = strcmp(style, 'light');
if strcmp(w.Type, 'figure')
    if white
        w.Color = [1 1 1];
    elseif isprop(w, 'ColorMode')
        w.ColorMode = 'auto';
    end
end
for h = findall(w, '-property', 'BackgroundColor')'
    if ~ismember(h.Type, {'uigridlayout', 'uipanel', 'uitab', 'axes', 'uibutton', 'uistatebutton', 'uidropdown'}) ...
            || strcmp(h.Tag, 'run')
        continue                                   % the green Run and the striped tables keep their colours
    end
    if white
        h.BackgroundColor = [1 1 1];
    elseif isprop(h, 'BackgroundColorMode')
        h.BackgroundColorMode = 'auto';
    end
end
end
