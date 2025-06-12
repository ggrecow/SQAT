function OUT = PsychoacousticAnnoyance_Widmann1992_from_percentile(N,S,R,FS)
% function OUT = PsychoacousticAnnoyance_Widmann1992_from_percentile(N,S,R,FS)
%
%   This function calculates Widmann's psychoacoustic annoyance model from an input acoustic signal
%
%   The psychoacoustic annoyance model is according to: (page 66) Widmann, U. (1992). Ein Modell der
%   Psychoakustischen Lästigkeit von Schallen und seine Anwendung in der Praxis der Lärmbeurteilung
%   (A model of the psychoacoustic annoyance of sounds and its application in noise assessment practice)
%   [Doctoral thesis, Technische Universität München (Technical University of Munich)].
%
%   As clarified by Lotinga, M. J. B. and A. J. Torija (2025) in
%   "Comment on "A study on calibration methods of noise annoyance data from listening tests"
%   [J. Acoust. Soc. Am. 156, 1877–1886 (2024)]." Journal of the Acoustical
%   Society of America 157(5): 3282–3285, this model is the same as that commonly
%   misattributed to (page 327) Zwicker, E. and Fastl, H. Second ed,
%   Psychoacoustics, Facts and Models, 2nd ed. Springer-Verlag, Berlin, 1999.
%
% - This metric combines 4 psychoacoustic metrics to quantitatively describe annoyance:
%
%    1) Loudness (sone) - calculated hereafter following ISO 532-1:2017
%       type <help Loudness_ISO532_1> for more info
%
%    2) Sharpness (acum) - calculated hereafter following DIN 45692:2009
%       NOTE: uses DIN 45692 weighting function by default, please change code if
%       the use of a different weighting function is desired. However, note
%       that the original PA model sharpness weighting is equal to the DIN
%       45692 weighting (i.e., Widmann weighting).
%       type <help Sharpness_DIN45692_from_loudness>
%
%    3) Roughness (asper) - calculated hereafter following Daniel & Weber model
%       type <help Roughness_Daniel1997> for more info
%
%    4) Fluctuation strength (vacil) - calculated hereafter following Osses et al. model
%       type <help FluctuationStrength_Osses2016> for more info
%
%   It should be noted however that the original model used metrics for
%   roughness and fluctuation strength that differ from those employed in
%   this implementation. The metrics employed in the original PA model
%   comprised Fastl's roughness and fluctuation strength (see Fastl &
%   Zwicker, 2007. Psychoacoustics: Facts and models.)
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
% Modified: Mike Lotinga, 12.06.2025 - created from
% PsychoacousticAnnoyance_Zwicker1999_from_percentile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (scalar) psychoacoustic annoyance - computed directly from percentile values

% sharpness and loudness influence
if S > 1.75
    % in Widmann (1992), log is used without specifying the base. In
    % Fastl&Zwicker (2007), lg is used and subsequent literature also uses log10
    ws = (S-1.75)*(log10(N+10))/4;
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
