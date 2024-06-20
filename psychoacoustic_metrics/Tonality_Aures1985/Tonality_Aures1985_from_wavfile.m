function OUT = Tonality_Aures1985_from_wavfile(wavfilename, dBFS, LoudnessField, start_skip, end_skip, show)
% function OUT = Tonality_Aures1985_from_wavfile(wavfilename, dBFS, LoudnessField, start_skip, end_skip, show)
%
%   This function calculates tonality metric by:
%
%   [1] Aures, Wilhelm (1985). "Berechnungsverfahren für den sensorischen Wohlklang 
%       beliebiger Schallsignale." Acta Acustica united with Acustica 59: p. 130-141.
%
%   The Aures' tonality is based on Terhard's virtual pitch theory, given by:
%
%   [2] Terhardt, E., Stoll, G. and Seewann, M. (1982). Algorithm for 
%       extraction of pitch and pitch salience from complex tonal signals. 
%       J. Acoust. Soc. Am., 71, 679-688. doi:10.1121/1.387544
%
%  Loudness calculation is conducted according to ISO 532:1-2017
%  (type <help Loudness_ISO532_1> for more info)
%
%  This script, Tonality_Aures1985_from_wavfile, calls internally the main
%    algorithm, Tonality_Aures1985. The only difference is that 
%    Tonality_Aures1985_from_wavfile requires a file name as first input
%    argument and the dBFS convention value as the second input argument.
%
%   Reference: a pure tone with 1000 Hz and 60 dBSPL has a tonality of 1 t.u.
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
%   LoudnessField : integer
%   chose field for loudness calculation; free field = 0; diffuse field = 1;
%   type <help Loudness_ISO532_1> for more info
%
%   start_skip : number
%   end_skip : number
%   skip start/end of the signal in seconds for statistic calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUT:
%   OUT : struct containing the following fields
%
%       * InstantaneousTonality: instantaneous tonality (t.u.)  vs time
%       * TonalWeighting: tonal weighting as a function of time
%       * LoudnessWeighting: loudness weighting as a function of time
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousTonality
%         ** Kmean : mean value of InstantaneousTonality (t.u.)
%         ** Kstd : standard deviation of InstantaneousTonality (t.u.)
%         ** Kmax : maximum of InstantaneousTonality (t.u.)
%         ** Kmin : minimum of InstantaneousTonality (t.u.)
%         ** Kx : percentile InstantaneousTonality exceeded during x percent of the signal (t.u.)
%
% Stand-alone example:
%   fname = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep 'RefSignal_Tonality_Aures1985.wav'];
%   dBFS = 94; % default for SQAT
%   Tonality_Aures1985_from_wavfile(fname,dBFS);
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help Tonality_Aures1985_from_wavfile;
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
    pars = psychoacoustic_metrics_get_defaults('Tonality_Aures1985');
    end_skip = pars.end_skip;
    fprintf('\n%s.m: Default end_skip value = %.0f is being used\n', ...
            mfilename, pars.end_skip);
end

if nargin < 4
    pars = psychoacoustic_metrics_get_defaults('Tonality_Aures1985');
    start_skip = pars.start_skip;
    fprintf('\n%s.m: Default start_skip value = %.0f is being used\n', ...
            mfilename, pars.start_skip);
end

if nargin < 3
    pars = psychoacoustic_metrics_get_defaults('Tonality_Aures1985');
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

OUT = Tonality_Aures1985(insig, fs, LoudnessField, start_skip, end_skip, show);

end % End of file

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
