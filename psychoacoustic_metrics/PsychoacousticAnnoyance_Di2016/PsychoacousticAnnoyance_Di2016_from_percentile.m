function OUT = PsychoacousticAnnoyance_Di2016_from_percentile(N,S,R,FS,K)
% function OUT = PsychoacousticAnnoyance_Di2016_from_percentile(N,S,R,FS,K)
%
%   This function calculates the Di's modified psychoacoustic annoyance model from scalar inputs
%   corresponding to the percentile values of loudness, sharpness, roughness, fluctuation strength and tonality
%
%   The modified psychoacoustic annoyance model is according to:
%   [1] Di et al., Improvement of Zwicker’s psychoacoustic annoyance model aiming at tonal noises, Applied Acoustics 105 (2016) 164-170
%
% - This metric combines 5 psychoacoustic metrics to quantitatively describe annoyance:
%
%    1) Loudness, N (sone)
%
%    2) Sharpness, S (acum)
%
%    3) Roughness, R (asper)
%
%    4) Fluctuation strength, FS (vacil)
%
%    5) Tonality, K (t.u.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   N: scalar
%   loudness percentile value (sone)
%
%   S: scalar
%   sharpness percentile value (acum)
%
%   R: scalar
%   roughness percentile value (asper)
%
%   FS: scalar
%   fluctuation strength percentile value (vacil)
%
%   K: scalar
%   tonality percentile value (t.u.)
%
% OUTPUTS:
%   OUT : scalar
%   modified psychoacoustic annoyance computed using the input percentile values of each metric
%
% Author: Gil Felix Greco, Braunschweig 05.04.2023
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

%% modified PA model constants (Ref. [1] pg. 168, eq (9))

alpha = 0.52;
beta = 6.41;

%% (scalar) modified psychoacoustic annoyance - computed directly from percentile values

% sharpness and loudness influence
if S > 1.75
    ws = (S-1.75)*(log10(N+10))/4; % in the Fastl&zwicker book, ln is used but it is not clear if it is natural log or log10, but most of subsequent literature uses log10
else
    ws = 0;
end

ws( isinf(ws) | isnan(ws) ) = 0;  % replace inf and NaN with zeros

% influence of roughness and fluctuation strength
wfr = ( 2.18/(N^(0.4)) ) * (0.4*FS + 0.6*R);

wfr( isinf(wfr) | isnan(wfr) ) = 0;  % replace inf and NaN with zeros

% Tonality influence
wt = (beta/(N^(alpha))) * K;

wt( isinf(wt) | isnan(wt) ) = 0;  % replace inf and NaN with zeros

% Di's modified psychoacoustic annoyance
PA_scalar = N*( 1 + sqrt( ws^2 + wfr^2 + wt^2 ) );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main output results

OUT=PA_scalar; % Annoyance calculated from the percentiles of each variable

end % end PA function
