function w = SQAT_GUI_window(name, n, custom)
% function w = SQAT_GUI_window(name, n, custom)
%
%   Window of n samples for the spectrogram of the interface. The classic
%   windows are periodic (the form used for spectral analysis) and are
%   written out here, so they need no toolbox besides SQAT.
%
% INPUT ARGUMENTS
%   name : 'hann', 'hamming', 'rect' (or 'rectangular'), 'blackmanharris'
%          or 'custom'
%   n : number of samples
%   custom : vector of the window when name is 'custom'; it is resampled
%            linearly to n samples, and returned as it is when it has n
%
% OUTPUTS
%   w : [nx1] window, with its peak at 1
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

switch lower(name)
    case 'hann'
        k = (0:n-1)';
        w = 0.5 - 0.5*cos(2*pi*k/n);
    case 'hamming'
        k = (0:n-1)';
        w = 0.54 - 0.46*cos(2*pi*k/n);
    case {'rect', 'rectangular'}
        w = ones(n, 1);
    case 'blackmanharris'
        k = (0:n-1)';
        w = 0.35875 - 0.48829*cos(2*pi*k/n) + 0.14128*cos(4*pi*k/n) - 0.01168*cos(6*pi*k/n);
    case 'custom'
        if nargin < 3 || ~isnumeric(custom) || numel(custom) < 2
            error('SQAT_GUI:window', 'A custom window needs a vector of at least 2 samples.');
        end
        if numel(custom) == n
            w = custom(:);
        else
            w = interp1(linspace(0, 1, numel(custom)), custom(:), linspace(0, 1, n))';
        end
    otherwise
        error('SQAT_GUI:window', 'Unknown window %s.', name);
end
end
