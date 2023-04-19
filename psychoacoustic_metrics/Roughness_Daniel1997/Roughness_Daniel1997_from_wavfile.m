function OUT = Roughness_Daniel1997_from_wavfile(wavfilename,dBFS,time_skip,show)
% function OUT = Roughness_Daniel1997_from_wavfile(wavfilename,dBFS,time_skip,show)
%
%   This function calculates time-varying roughness and time-averaged specific
%     roughness using the roughness model by Daniel & Weber:
%     Daniel, P., & Weber, R. (1997). Psychoacoustical roughness: implementation
%     of an optimized model. Acustica(83), 113-123.
%
%   Reference signal: 60 dB 1 kHz tone 100% modulated at 70 Hz should yield 1 asper.
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
%       * InstantaneousRoughness: instantaneous roughness (asper) as a 
%         function of time
%       * InstantaneousSpecificRoughness: specific roughness(asper/Bark) as
%         a function of time and frequency (Bark scale)
%       * TimeAveragedSpecificRoughness: time-averaged specific roughness 
%         (asper/Bark) as a function of frequency (Bark scale)
%       * barkAxis : vector of Bark band numbers used for the computation
%         of specific roughness computation
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousRoughness
%         ** Rmean : mean value of instantaneous roughness (asper)
%         ** Rstd : standard deviation of instantaneous roughness (asper)
%         ** Rmax : maximum of instantaneous roughness (asper)
%         ** Rmin : minimum of instantaneous roughness (asper)
%         ** Rx : percentile roughness exceeded during x percent of the signal (asper)
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'Roughness_Daniel1997' filesep 'RefSignal_Roughness_1asper_48kHz_32bit.wav'];
%   dBFS = 94; % default for SQAT
%   Roughness_Daniel1997_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses
% See also: Roughness_Daniel1997.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help Roughness_Daniel1997_from_wavfile;
    return;
end

if nargin < 4
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end
if nargin <3
    pars = psychoacoustic_metrics_get_defaults('Roughness_Daniel1997');
    time_skip = pars.time_skip;
    fprintf('%s.m: Default time_skip value = %.0f is being used\n',mfilename,pars.time_skip);
end

[insig,fs] = audioread(wavfilename);
if nargin < 2 || isempty(dBFS)
    dBFS = 94; % dB
    fprintf('%s.m: Assuming the default full scale convention, with dBFS = %.0f\n',mfilename,dBFS);
end
gain_factor = 10^((dBFS-94)/20);
insig = gain_factor*insig;

OUT = Roughness_Daniel1997(insig,fs,time_skip,show);

end % end of file
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