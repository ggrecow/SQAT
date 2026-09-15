function OUT = PsychoacousticAnnoyance_More2010_from_wavfile(wavfilename,dBFS,LoudnessField,time_skip,showPA,show)
% function OUT = PsychoacousticAnnoyance_More2010_from_wavfile(wavfilename,dBFS,LoudnessField,time_skip,showPA,show)
%
%   This function calculates the More's modified psychoacoustic annoyance 
%   model from an input acoustic signal
%
%   The modified psychoacoustic annoyance model is according to: (page 201)
%   [1] More, Shashikant. Aircraft noise characteristics and metrics. 
%       PhD Thesis, Purdue University, 2010
%
% - This metric combines 5 psychoacoustic metrics to quantitatively describe annoyance:
%
%    1) Loudness, N (sone) - calculated hereafter following ISO 532-1:2017
%       type <help Loudness_ISO532_1> for more info
%
%    2) Sharpness, S (acum) - calculated hereafter following DIN 45692:2009
%       NOTE: uses DIN 45692 weighting function by default, please change code if
%       the use of a different withgitng function is desired).
%       type <help Sharpness_DIN45692_from_loudness>
%
%    3) Roughness, R (asper) - calculated hereafter following Daniel & Weber model
%       type <help Roughness_Daniel1997> for more info
%
%    4) Fluctuation strength, FS (vacil) - calculated hereafter following 
%       Osses et al. model, type <help FluctuationStrength_Osses2016> for more info
%
%    5) Tonality, K (t.u.) - calculated hereafter following Aures' model
%       type <help Tonality_Aures1985> for more info
%
%  This script, PsychoacousticAnnoyance_More2010_from_wavfile, calls 
%    internally the main algorithm, PsychoacousticAnnoyance_More2010. The 
%    only difference is that PsychoacousticAnnoyance_More2010_from_wavfile 
%    requires a file name as first input argument and the dBFS convention 
%    value as the second input argument.
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
%   chose field for loudness calculation; free field = 0; diffuse field = 1; (used in the loudness and tonality codes)
%   type <help Loudness_ISO532_1> for more info
%
%   show : logical(boolean)
%   optional parameter, display results of loudness, sharpness, roughness, fluctuation strength and tonality
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
%             ** wt : tonality and loudness weighting function (not squared)
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
%        **  K : struct with tonality results, type <help Tonality_Aures1985> for more info
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_Loudness_ISO532_1.wav'];
%   dBFS = 94; % default for SQAT
%   PsychoacousticAnnoyance_More2010_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses
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
if nargin == 0
    help PsychoacousticAnnoyance_More2010_from_wavfile;
    return;
end
if nargin < 6
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end
if nargin < 5
    if nargout == 0
        showPA = 1; 
    else
        showPA = 0;
    end
end
if nargin <4
    pars = psychoacoustic_metrics_get_defaults('PsychoacousticAnnoyance_More2010');
    time_skip = pars.time_skip;
    fprintf('%s.m: Default time_skip value = %.0f is being used\n',mfilename,pars.time_skip);
end
if nargin <3
    pars = psychoacoustic_metrics_get_defaults('PsychoacousticAnnoyance_More2010');
    LoudnessField = pars.Loudness_field;
    fprintf('\n%s.m: Default Loudness_field value = %.0f is being used\n',mfilename,pars.Loudness_field);
end

[insig,fs] = audioread(wavfilename);
if nargin < 2 || isempty(dBFS)
    dBFS = 94; % dB
    fprintf('\n%s.m: Assuming the default full scale convention, with dBFS = %.0f\n',mfilename,dBFS);
end
gain_factor = 10^((dBFS-94)/20);
insig = gain_factor*insig;

OUT = PsychoacousticAnnoyance_More2010(insig,fs,LoudnessField,time_skip,showPA,show);

end % end of file