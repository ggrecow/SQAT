function OUT = Sharpness_DIN45692_from_wavfile(wavfilename, dBFS, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness)
% function OUT = Sharpness_DIN45692_from_wavfile(wavfilename, dBFS, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness)
%
%  Stationary and time-varying sharpness calculation according to DIN 45692
%    (2009) from a waveform (wav file). The loudness calculation, required 
%    as pre-processing for sharpness, is included in this code.
%
%  This script, Sharpness_DIN45692_from_wavfile, calls internally the main
%    sharpness algorithm, Sharpness_DIN45692. The only difference is that 
%    Sharpness_DIN45692_from_wavfile requires a file name as first input
%    argument and the dBFS convention value as the second input argument.
%
%  Loudness calculation is conducted according to ISO 532:1-2017
%  (type <help Loudness_ISO532_1> for more info)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT ARGUMENTS
%   wavfilename : file name containing the input signal 'insig' to be 
%       processed and the sampling frequency 'fs' of insig. Wavfilename 
%       should be visible to MATLAB.
%
%   dBFS: dB sound pressure level equivalent to the full scale value. If 
%       dBFS is 94 dB SPL (default), the amplitudes of 'insig' are interpreted
%       as expressed in pressure units (Pa).
%
%   weight_type : string
%   sharpness calculation using weighting function according to:
%       - 'DIN45692'
%       - 'bismarck'
%       - 'aures' (dependent on the specific loudness level)
%
%   LoudnessField : integer
%   type of field used for loudness calculation; free field = 0; diffuse field = 1;
%
%   LoudnessMethod : integer
%   method used for loudness calculation - method used for loudness 
%        calculation: stationary (from input 1/3 octave unweighted SPL)=0 (not 
%        accepted in this context); stationary = 1; time varying = 2;
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations (method=1 (time-varying) only)
%
%   show_loudness : logical(boolean)
%   optional parameter to display loudness results (only method=1)
%   'false' (disable, default value) or 'true' (enable).
%
%   show_sharpness : logical(boolean)
%   optional parameter to display sharpness results (only method=1)
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUTS (method==0; stationary)
%   OUT : struct containing the following fields
%
%       * Sharpness: sharpness (acum)
%
% OUTPUTS (method==1; time-varying)
%   OUT : struct containing the following fields
%
%       * loudness: output struct from loudness calculation (type 
%         <help loudness_ISO532_1> for more info)
%       * InstantaneousSharpness: instantaneous sharpness (acum) vs time
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousSharpness (acum)
%         ** Smean : mean value of InstantaneousSharpness (acum)
%         ** Sstd : standard deviation of InstantaneousSharpness (acum)
%         ** Smax : maximum of InstantaneousSharpness (acum)
%         ** Smin : minimum of InstantaneousSharpness (acum)
%         ** Sx : sharpness value exceeded during x percent of the time (acum)
%
%           *** HINT: time-varying loudness calculation takes some time to
%                     have a steady-response (thus sharpness too!). 
%                     Therefore, it is a good practice to consider a 
%                     time_skip to compute the statistics
%
% % Stand-alone example:
%     fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_Sharpness_DIN45692.wav'];
%     dBFS = 90; % We know this information in advance
%     Sharpness_DIN45692_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses, based on the code by Gil Felix Greco
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
    help Sharpness_DIN45692_from_wavfile;
    return;
end
if nargin < 8
    if nargout == 0
        show_loudness = 1;
    else
        show_loudness = 0;
    end
end
if nargin < 7
    if nargout == 0
        show_sharpness = 1;
    else
        show_sharpness = 0;
    end
end
if nargin < 6
    pars = psychoacoustic_metrics_get_defaults('Sharpness_DIN45692');
    time_skip = pars.time_skip;
    fprintf('\n%s.m: Default time_skip value = %.0f is being used\n',mfilename,pars.time_skip);
end
if nargin < 5
    pars = psychoacoustic_metrics_get_defaults('Sharpness_DIN45692');
    LoudnessMethod = pars.method;
    fprintf('\n%s.m: Default method value = %.0f is being used\n',mfilename,pars.method);
end
if nargin < 4
    pars = psychoacoustic_metrics_get_defaults('Sharpness_DIN45692');
    LoudnessField = pars.field;
    fprintf('\n%s.m: Default field value = %.0f is being used\n',mfilename,pars.field);
end
if nargin < 3
    pars = psychoacoustic_metrics_get_defaults('Sharpness_DIN45692');
    weight_type = pars.weight_type;
    fprintf('\n%s.m: Default weight_type value = %s is being used\n',mfilename,pars.weight_type);
end
[insig,fs] = audioread(wavfilename);
if nargin < 2 || isempty(dBFS)
    dBFS = 94; % dB
    fprintf('\n%s.m: Assuming the default full scale convention, with dBFS = %.0f\n',mfilename,dBFS);
end
gain_factor = 10^((dBFS-94)/20);
insig = gain_factor*insig;

OUT = Sharpness_DIN45692(insig, fs, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness);

end % End of file
