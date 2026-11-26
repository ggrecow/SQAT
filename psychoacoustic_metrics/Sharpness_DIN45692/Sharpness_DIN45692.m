function OUT = Sharpness_DIN45692(insig, fs, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness)
% function OUT = Sharpness_DIN45692(insig, fs, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness)
%
%  Stationary and time-varying sharpness calculation according to DIN 45692
%    (2009) from an input signal. The loudness calculation, required as pre-
%    processing for sharpness, is included in this code.
%
%  Loudness calculation is conducted according to ISO 532:1-2017
%  (type <help Loudness_ISO532_1> for more info)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT ARGUMENTS
%   insig : [Nx1] array
%   calibrated audio signal (Pa), 1 channel only
%
%   fs : integer
%   sampling frequency (Hz). For method = 3, provide a dummy scalar
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
%       calculation: stationary (from input 1/3 octave unweighted SPL)=0 (not 
%       accepted in this context); stationary = 1; time varying = 2;
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics 
%        calculations (method=1 (time-varying) only)
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
% Author: Gil Felix Greco, Braunschweig 09.03.2023
% Author: Gil Felix Greco, Braunschweig 16.02.2025 - introduced get_statistics function
% Modified: Mike Lotinga 29.05/2026 - vectorised for 95% reduction in
% compute time
%
% Modified: Gil Felix Greco, November 2026 (the sharpness calculation is
%   now done by Sharpness_DIN45692_from_loudness, so that both functions
%   share one implementation and always give the same results)
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
    help Sharpness_DIN45692;
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

% loudness calculation
L = Loudness_ISO532_1(insig, fs,... % input signal and sampling freq.
                  LoudnessField,... % free field = 0; diffuse field = 1;
                 LoudnessMethod,... % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                      time_skip,... % time_skip
                  show_loudness);   % show loudness results

% sharpness calculation (see Sharpness_DIN45692_from_loudness)
if LoudnessMethod == 1 % (stationary sharpness)
    
    OUT = Sharpness_DIN45692_from_loudness(L.SpecificLoudness, weight_type, [], time_skip, show_sharpness);
    
elseif LoudnessMethod == 2 % (time-varying sharpness)
    
    OUT = Sharpness_DIN45692_from_loudness(L.InstantaneousSpecificLoudness, weight_type, L.time, time_skip, show_sharpness);
    OUT.loudness = L; % output struct from the loudness calculation
    
end

end % End of Sharpness_DIN45692
