function OUT = PsychoacousticAnnoyance_Zwicker1999_from_wavfile(wavfilename,dBFS,LoudnessField,time_skip,showPA,show)
% function OUT = PsychoacousticAnnoyance_Zwicker1999_from_wavfile(wavfilename,dBFS,LoudnessField,time_skip,showPA,show)
%
%   This function is a compatibility-wrapper that calculates Widmann's
%   psychoacoustic annoyance model from an input acoustic signal ---
%   As clarified by Lotinga, M. J. B. and A. J. Torija (2025) in
%   "Comment on "A study on calibration methods of noise annoyance data from listening tests"
%   [J. Acoust. Soc. Am. 156, 1877–1886 (2024)]." Journal of the Acoustical
%   Society of America 157(5): 3282–3285, this model is the same as that commonly
%   misattributed to (page 327) Zwicker, E. and Fastl, H. Second ed,
%   Psychoacoustics, Facts and Models, 2nd ed. Springer-Verlag, Berlin, 1999.
%
%   The original psychoacoustic annoyance model is according to: (page 66) Widmann, U. (1992). Ein Modell der
%   Psychoakustischen Lästigkeit von Schallen und seine Anwendung in der Praxis der Lärmbeurteilung
%   (A model of the psychoacoustic annoyance of sounds and its application in noise assessment practice)
%   [Doctoral thesis, Technische Universität München (Technical University of Munich)].
%
%   Widmann defined a 1-kHz tone with 40 dB SPL as the reference signal for
%   his PA model (see page 65 in the above mentioned reference), to which 
%   he assigned an annoyance value of 1 au (annoyance unit).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   WARNING: this function is merely a wrapper of the
%   <PsychoacousticAnnoyance_Widmann1992> function, kept to maintain
%   compatibility with SQAT v1.3 and below. Nevertheless, this function,
%   will be removed in future releases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   wavfilename : char
%   wavfilename specifies the file name of a wav file to be processed
%
%   dBFS : integer
%          Full scale convention. Internally this algorithm works with 
%          a convention of full scale being equal to 94 dB SPL, or dBFS=94.
%          if the specified dBFS is different from 94 dB SPL, then a gain 
%          factor will be applied
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations
%
%   LoudnessField : integer
%   chose field for loudness calculation; free field = 0; diffuse field = 1;
%   type <help Loudness_ISO532_1> for more info
%
%   show : logical(boolean)
%   optional parameter, display results of loudness, sharpness, roughness and fluctuation strength
%   'false' (disable, default value) or 'true' (enable).
%
%   showPA : logical(boolean)
%   optional parameter, display results of psychoacoustic annoyance
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUTS:
%   OUT: struct
%      * include results from the psychoacoustic annoyance
%             ** InstantaneousPA: instantaneous quantity (unity) vs time
%             ** ScalarPA : PA (scalar value) computed using the percentile values of each metric.
%                           NOTE: if the signal's length is smaller than 2s, this is the only output as no time-varying PA is calculated
%             ** time : time vector in seconds
%             ** wfr : fluctuation strength and roughness weighting function (not squared)
%             ** ws : sharpness and loudness weighting function (not squared)
%
%             ** Statistics
%               *** PAmean : mean value of psychoacoustic annoyance (unit)
%               *** PAstd : standard deviation of instantaneous psychoacoustic annoyance (unit)
%               *** PAmax : maximum of instantaneous psychoacoustic annoyance (unit)
%               *** PAmin : minimum of instantaneous psychoacoustic annoyance (unit)
%               *** PAx : value exceeded x percent of the time
%
%      * include structs with the results from the other metrics computed
%        **  L : struct with Loudness results, type <help Loudness_ISO532_1> for more info
%        **  S : struct with Sharpness, type <help Sharpness_DIN45692_from_loudness>
%        **  R : strcut with roughness results, type <help Roughness_Daniel1997> for more info
%        ** FS : struct with fluctuation strength results, type <help FluctuationStrength_Osses2016> for more info
%
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_Loudness_ISO532_1.wav'];
%   dBFS = 94; % default for SQAT
%   PsychoacousticAnnoyance_Zwicker1999_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses
% Modified: Mike Lotinga, 12.06.2025 - moved content to
% PsychoacousticAnnoyance_Widmann1992_from_wavfile and made function a wrapper
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

OUT = PsychoacousticAnnoyance_Widmann1992_from_wavfile(wavfilename,dBFS,LoudnessField,time_skip,showPA,show);

