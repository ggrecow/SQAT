function OUT = PsychoacousticAnnoyance_Zwicker1999_from_percentile(N,S,R,FS)
% function OUT = PsychoacousticAnnoyance_Zwicker1999_from_percentile(N,S,R,FS)
%
%   This function calculates the Zwicker's psychoacoustic annoyance model from scalar inputs
%   corresponding to the percentile values of loudness, sharpness, roughness and fluctuation strength
%
%   The psychoacoustic annoyance model is according to: (page 327) Zwicker, E. and Fastl, H. Second ed,
%   Psychoacoustics, Facts and Models, 2nd ed. M.R. Schroeder. Springer-Verlag, Berlin, 1999.
%
% - This metric combines 4 psychoacoustic metrics to quantitatively describe annoyance:
%
%    1) Loudness, N (sone)
%
%    2) Sharpness, S (acum)
%
%    3) Roughness, R (asper)
%
%    4) Fluctuation strength, FS (vacil)
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
% OUTPUTS:
%   OUT : scalar
%   Psychoacoustic Annoyance computed using the input percentile values of each metric
%
% Author: Gil Felix Greco, Braunschweig 14.03.2023
%
%  This code is part of SQAT v1.0, released 14.05.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (scalar) psychoacoustic annoyance - computed directly from percentile values

% sharpness and loudness influence
if S > 1.75
    ws = (S-1.75)*(log10(N+10))/4; % in the Fastl&zwicker book, ln is used but it is not clear if it is natural log or log10, but most of subsequent literature uses log10
else
    ws = 0;
end

ws( isinf(ws) | isnan(ws) ) = 0;  % replace inf and NaN with zeros

% influence of roughness and fluctuation strength
wfr = ( 2.18/(N^(0.4)) )*(0.4*FS + 0.6*R);

wfr( isinf(wfr) | isnan(wfr) ) = 0;  % replace inf and NaN with zeros

% psychoacoustic annoyance
PA_scalar = N*( 1 + sqrt (ws^2 + wfr^2) );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main output results
OUT=PA_scalar;               % Annoyance calculated from the percentiles of each variable


end % end PA function

%**************************************************************************
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%**************************************************************************