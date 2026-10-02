function [dBFS, label] = SQAT_GUI_calibration(method, wavfilename, level, calfilename)
% function [dBFS, label] = SQAT_GUI_calibration(method, wavfilename, level, calfilename)
%
%   The full-scale level of each channel of a .wav file, in dB SPL (the dBFS
%   that SQAT_GUI_load takes), from one of the three ways to calibrate a
%   sound file described by Greco (2026, Section 3.3):
%
%   'dbfs'       : the level of a sample of value 1 is known (Scaling of WAV
%                  files, Eq. 3.2); level is that value, for every channel.
%   'calibrator' : a recording of a sound level calibrator made with the same
%                  setup (Indirect method, Eq. 3.1); level is the level of the
%                  calibrator, e.g. 94 dB (1 Pa), 114 dB (10 Pa) or the value
%                  of an adapter. dBFS = level - 20 log10(rms of the
%                  recording): one value per channel when the recording has as
%                  many channels as the file (each ear of a head and torso
%                  simulator, for example), else the first channel for all.
%   'relative'   : no information on the level (Level adjustment, Eq. 3.3);
%                  the rms of the whole file, all channels together, is set
%                  to level, so the channels keep their difference of level.
%                  The results then only compare signals among themselves.
%
% INPUT ARGUMENTS
%   method : 'dbfs', 'calibrator' or 'relative'
%   wavfilename : path of the .wav file to calibrate
%   level : dB SPL (see above)
%   calfilename : path of the calibrator recording ('calibrator' only)
%
% OUTPUTS
%   dBFS : [1 x nch] full-scale level of each channel, in dB SPL
%   label : short text for the list of signals, e.g. 'calib. 94 dB'
%
% Reference: Greco, G. F. (2026). Sound quality analysis of virtual aircraft
%   prototypes: framework development and application. Dissertation, TU
%   Braunschweig, Schriften des Instituts fuer Akustik und Dynamik, Band 5,
%   Section 3.3, Calibration and level adjustment of sound signals.
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September 2026 assisted
% by Claude Opus 5.5 (Anthropic). All codes were verified by
% the authors.
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

nch = audioinfo(wavfilename).NumChannels;
switch method
    case 'dbfs'
        dBFS = level * ones(1, nch);
        label = sprintf('%g dBFS', level);
    case 'calibrator'
        c = audioread(calfilename);
        r = sqrt(mean(c.^2, 1));                % rms of each channel of the recording
        if any(r == 0)
            error('SQAT_GUI:calibration', 'The calibrator recording %s is silent.', calfilename);
        end
        if size(c, 2) ~= nch
            r = r(1) * ones(1, nch);
        end
        dBFS = level - 20*log10(r);
        label = sprintf('calib. %g dB', level);
    case 'relative'
        y = audioread(wavfilename);
        r = sqrt(mean(y(:).^2));                % one rms over all channels
        if r == 0
            error('SQAT_GUI:calibration', 'The file %s is silent: no level to adjust.', wavfilename);
        end
        dBFS = (level - 20*log10(r)) * ones(1, nch);
        label = sprintf('rel. %g dB', level);
    otherwise
        error('SQAT_GUI:calibration', 'Unknown calibration method %s.', method);
end
end
