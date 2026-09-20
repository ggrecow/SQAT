function OUT = FluctuationStrength_Osses2016_from_wavfile(wavfilename,dBFS,method,time_skip,show)
% function OUT = FluctuationStrength_Osses2016_from_wavfile(wavfilename,dBFS,method,time_skip,show)
%
%  This function calculates the fluctuation strength using the model 
%    developed by: [1] Osses, A., Garcia A., and Kohlrausch, A.. 
%    "Modelling the sensation of fluctuation strength." Proceedings of 
%    Meetings on Acoustics 22 ICA. Vol. 28, 050005. doi:10.1121/2.0000410
%
%  This script, FluctuationStrength_Osses2016_from_wavfile, calls internally
%    the main algorithm, FluctuationStrength_Osses2016. The only difference
%    is that FluctuationStrength_Osses2016_from_wavfile requires a file 
%    name as first input argument and the dBFS convention value as the 
%    second input argument.
%
%  Reference signal: 60 dBSPL 1 kHz tone 100% modulated at 4 Hz should yield 1 vacil.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   wavfilename : char
%   wavfilename specifies the file name of a wav file to be processed
%
%   dBFS : integer
%          Full scale convention. Internally the this algorithm works with 
%          a convention of full scale being equal to 94 dB SPL, or dBFS=94.
%          if the specified dBFS is different from 94 dB SPL, then a gain 
%          factor will be applied
%   method : integer
%   method=0, stationary analysis - window size=length(insig) (s) kind of 
%             an rms value
%   method=1, time_varying analysis - window size=2 (s)
%             NOTE: if the signal's length is smaller than 2s, the analysis
%             is automatically changed to method=0
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUT:
%   OUT : struct containing the following fields
%
%       * InstantaneousFluctuationStrength: instantaneous fluctuation 
%           strength (vacil) vs time
%       * InstantaneousSpecificFluctuationStrength: specific fluctuation 
%           strength (vacil/Bark) vs time and Bark scale
%       * TimeAveragedSpecificFluctuationStrength: time-averaged specific 
%           fluctuation strength (vacil/Bark) vs Bark scale
%       * barkAxis : vector of Bark band numbers used for specific 
%           fluctuation strength computation
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousFS
%         ** FSmean : mean value of InstantaneousFS (vacil)
%         ** FSstd : standard deviation of InstantaneousFluctuationStrength
%              (vacil)
%         ** FSmax : maximum of InstantaneousFluctuationStrength (vacil)
%         ** FSmin : minimum of InstantaneousFluctuationStrength (vacil)
%         ** FSx : fluctuation strength exceeded during x percent of the time (vacil)
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_FluctuationStrength_Osses2016.wav'];
%   dBFS = 94; % default for SQAT
%   FluctuationStrength_Osses2016_from_wavfile(fname,dBFS);
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
    help FluctuationStrength_Osses2016_from_wavfile;
    return;
end

if nargin < 5
    % Default for show
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

if nargin <4
    pars = psychoacoustic_metrics_get_defaults('FluctuationStrength_Osses2016');
    time_skip = pars.time_skip;
    fprintf('\n%s.m: Default time_skip value = %.0f is being used\n',mfilename,pars.time_skip);
end
if nargin <3
    pars = psychoacoustic_metrics_get_defaults('FluctuationStrength_Osses2016');
    method = pars.method;
    fprintf('\n%s.m: Default method = %.0f is being used\n',mfilename,pars.method);
end

[insig,fs] = audioread(wavfilename);
if nargin < 2 || isempty(dBFS)
    dBFS = 94; % dB
    fprintf('\n%s.m: Assuming the default full scale convention, with dBFS = %.0f\n',mfilename,dBFS);
end
gain_factor = 10^((dBFS-94)/20); %
insig = gain_factor*insig;

OUT = FluctuationStrength_Osses2016(insig, fs, method, time_skip, show);
