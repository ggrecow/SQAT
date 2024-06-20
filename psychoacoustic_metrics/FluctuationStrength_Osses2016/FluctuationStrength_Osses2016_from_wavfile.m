function OUT = FluctuationStrength_Osses2016_from_wavfile(wavfilename, dBFS, method, start_skip, end_skip, show)
% function OUT = FluctuationStrength_Osses2016_from_wavfile(wavfilename, dBFS, method, start_skip, end_skip, show)
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
%   start_skip : number
%   end_skip : number
%   skip start/end of the signal in seconds for statistics calculations
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
%         ** FSx : percentile fluctuation strength exceeded during x percent
%              of the signal (vacil)
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_FluctuationStrength_Osses2016.wav'];
%   dBFS = 94; % default for SQAT
%   FluctuationStrength_Osses2016_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help FluctuationStrength_Osses2016_from_wavfile;
    return;
end

if nargin < 6
    % Default for show
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

if nargin < 5
    pars = psychoacoustic_metrics_get_defaults('FluctuationStrength_Osses2016');
    end_skip = pars.end_skip;
    fprintf('\n%s.m: Default end_skip value = %.0f is being used\n',...
        mfilename, pars.end_skip);
end

if nargin < 4
    pars = psychoacoustic_metrics_get_defaults('FluctuationStrength_Osses2016');
    start_skip = pars.start_skip;
    fprintf('\n%s.m: Default start_skip value = %.0f is being used\n',...
        mfilename, pars.start_skip);
end

if nargin < 3
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

OUT = FluctuationStrength_Osses2016(insig, fs, method, start_skip, end_skip, show);

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