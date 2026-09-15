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
%   Widmann defined a 1-kHz tone with 40 dB SPL as the reference signal for
%   his PA model (see page 65 in the above mentioned reference), to which 
%   he assigned an annoyance value of 1 au (annoyance unit).
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
%   Zwicker, 2007. Psychoacoustics: Facts and models.). Neverthless, calculation
%   of Widmann's psychoacoustic annoyance from the reference signal (i.e., 40 dBSPL tone at 1 kHz)
%   using this implementation yields a value of 1 (au).
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
