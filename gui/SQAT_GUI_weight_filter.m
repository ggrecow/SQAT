function [b, a] = SQAT_GUI_weight_filter(fs, type)
% function [b, a] = SQAT_GUI_weight_filter(fs, type)
%
%   Digital filter of the A or C frequency weighting (IEC 61672-1). It is the
%   design of Gen_weighting_filters (analog zeros, poles and gain, and the
%   bilinear transformation without pre-warping) written with base MATLAB
%   only, because that function calls bilinear, of the Signal Processing
%   Toolbox. The coefficients are the same (a test compares them).
%
% INPUT ARGUMENTS
%   fs : sampling frequency (Hz)
%   type : 'A', 'C' or 'Z' (flat: b = a = 1)
%
% OUTPUTS
%   b, a : numerator and denominator of the filter, as for filter
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

w1 = 129.42731565506293;      % rad/s, the poles of IEC 61672-1 Annex E
w2 = 676.4015402329549;
w3 = 4636.125126885012;
w4 = 76618.52601685845;
switch type
    case 'A'
        K = 7.3901e9;
        z = zeros(4, 1);
        p = -[w1; w1; w2; w3; w4; w4];
    case 'C'
        K = 5.9124e9;
        z = zeros(2, 1);
        p = -[w1; w1; w4; w4];
    case 'Z'
        b = 1;
        a = 1;
        return
    otherwise
        error('SQAT_GUI:weighting', 'Unknown weighting %s: use A, C or Z.', type);
end
% bilinear transformation s = 2*fs*(z-1)/(z+1), the zeros that the analog filter has at
% infinity going to z = -1
fs2 = 2 * fs;
zd = [(fs2 + z) ./ (fs2 - z); -ones(numel(p) - numel(z), 1)];
pd = (fs2 + p) ./ (fs2 - p);
kd = real(K * prod(fs2 - z) / prod(fs2 - p));
b = real(kd * poly(zd));
a = real(poly(pd));
end
