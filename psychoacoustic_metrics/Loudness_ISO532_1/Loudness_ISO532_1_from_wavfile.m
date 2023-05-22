function OUT = Loudness_ISO532_1_from_wavfile(wavfilename, dBFS, field, method, time_skip, show)
% function OUT = Loudness_ISO532_1_from_wavfile(wavfilename, dBFS, field, method, time_skip, show)
%
%  Zwicker Loudness model according to ISO 532-1:2017 for stationary 
%    signals (Method A) and arbitrary signals (Method B).
%
%  This script, Loudness_ISO532_1_from_wavfile, calls internally the main
%    loudness algorithm, Loudness_ISO532_1. The only difference is that 
%    Loudness_ISO532_1_from_wavfile requires a file name as first input
%    argument and the dBFS convention value as the second input argument.
%
%  Reference signal: 40 dBSPL 1 kHz tone yields 1 sone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   field : integer 
%   free field = 0; diffuse field = 1;
%
%   method : integer 
%   0 = stationary (from input 1/3 octave unweighted SPL)
%   1 = stationary (from audio file) 
%   2 = time varying (from audio file)  
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUTS (method==0 and method==1; stationary method)
%   OUT : struct containing the following fields
%
%       * time_insig - time vector of the audio input, in seconds
%       * barkAxis - bark vector
%       * SpecificLoudness - time-averaged specific loudness (sone/Bark) 
%       * Loudness - loudness (sone) 
%       * LoudnessLevel - loudness level (phon)
%       * TimeAveragedSPL - time-averaged overall SPL (1/3 octave bands, DBSPL)
% 
% OUTPUTS (method==2; time-varying method)
%   OUT : struct containing the following fields
%
%       * barkAxis - vector of Bark band numbers used for specific loudness computation
%       * time - time vector of the final loudness calculation, in seconds
%       * time_insig - time vector of insig, in seconds
%       * InstantaneousLoudness - instantaneous loudness (sone) vs time
%       * InstantaneousSpecificLoudness - specific loudness (sone/Bark) vs time
%       * InstantaneousLoudnessLevel - instantaneous loudness level (phon) vs time
%       * SpecificLoudness - time-averaged specific loudness (sone/Bark) 
%       * InstantaneousSPL - overall SPL (1/3 octave bands) for each time step, in dBSPL
%       * Several statistics based on the InstantaneousLoudness 
%         ** Nmean : mean value of InstantaneousLoudness (sone)
%         ** Nstd : standard deviation of InstantaneousLoudness (sone)
%         ** Nmax : maximum of InstantaneousLoudness (sone)
%         ** Nmin : minimum of InstantaneousLoudness (sone)
%         ** Nx : percentile loudness exceeded during x percent of the signal (sone)
%         ** N_ratio : ratio between N5/N95 ( 1.1 (stationary)> N_ratio > 1.1 (time varying) )
%           *** HINT: loudness calculation takes some time to have a steady-response
%                     therefore, it is a good practice to consider a time_skip to compute the statistics 
%                     due to transient effects in the beginning of the loudness calculations
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_Loudness_ISO532_1.wav'];
%   dBFS = 94; % default for SQAT
%   Loudness_ISO532_1_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help Loudness_ISO532_1_from_wavfile;
    return;
end

if nargin < 6
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

if nargin <5
    pars = psychoacoustic_metrics_get_defaults('Loudness_ISO532_1');
    time_skip = pars.time_skip;
    fprintf('\n%s.m: Default time_skip value = %.0f is being used\n',mfilename,pars.time_skip);
end
if nargin <4
    pars = psychoacoustic_metrics_get_defaults('Loudness_ISO532_1');
    method = pars.method;
    fprintf('\n%s.m: Default method value = %.0f is being used\n',mfilename,pars.method);
end
if nargin <3
    pars = psychoacoustic_metrics_get_defaults('Loudness_ISO532_1');
    field = pars.field;
    fprintf('\n%s.m: Default field value = %.0f is being used\n',mfilename,pars.field);
end

[insig,fs] = audioread(wavfilename);
if nargin < 2 || isempty(dBFS)
    dBFS = 94; % dB
    fprintf('\n%s.m: Assuming the default full scale convention, with dBFS = %.0f\n',mfilename,dBFS);
end
gain_factor = 10^((dBFS-94)/20);
insig = gain_factor*insig;

OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show);

end

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
%**************************************************************************