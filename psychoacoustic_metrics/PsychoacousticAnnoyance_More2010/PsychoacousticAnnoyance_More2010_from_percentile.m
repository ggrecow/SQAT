function OUT = PsychoacousticAnnoyance_More2010_from_percentile(N,S,R,FS,K)
% function OUT = PsychoacousticAnnoyance_More2010_from_percentile(N,S,R,FS,K)
%
%   This function calculates the More's modified psychoacoustic annoyance model from scalar inputs
%   corresponding to the percentile values of loudness, sharpness, roughness, fluctuation strength and tonality
%
%   The modified psychoacoustic annoyance model is according to: (page 201)
%   [1] More, Shashikant R. Aircraft noise characteristics and metrics. PhD Thesis, Purdue University, 2010
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

%% modified PA model constants (Ref. [1] pg. 204)

gamma_0=    -0.16;
gamma_1=    11.48;
gamma_2=    0.84;
gamma_3=    1.25;
gamma_4=    0.29;
gamma_5=    5.49;

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
wt = abs( ( 1-exp(-gamma_4*N) )^2 * ( 1-exp(-gamma_5*K) )^2 );

wt( isinf(wt) | isnan(wt) ) = 0;  % replace inf and NaN with zeros

% More's modified psychoacoustic annoyance
PA_scalar = abs(N*( 1 + sqrt( gamma_0 + (gamma_1*ws^2) + (gamma_2* wfr^2) + (gamma_3*wt) ) ));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main output results

OUT=PA_scalar;               % Annoyance calculated from the percentiles of each variable

end % end PA function