function[CalibratedSignal,CalFactor,dBFS]=calibrate(InputSignal,RefSignal,ReferenceLevel)
% function[CalibratedSignal,CalFactor,dBFS]=calibrate(InputSignal,RefSignal,ReferenceLevel)
%
% Converts a digital signal (numbers without units, e.g. read from a .wav
%   file) into a sound pressure signal in Pa. The calibration factor is
%   chosen so that 'RefSignal', obtained with the same setup as
%   'InputSignal', has an rms level equal to 'ReferenceLevel'. This is the
%   calibration approach of the software provided with ISO 532-1:2017.
%
%   Typical uses:
%   - calibration: 'RefSignal' is a recording of a sound level calibrator
%     and 'ReferenceLevel' is its level (e.g. 94 dB SPL)
%
%   - level adjustment: 'RefSignal' is 'InputSignal' itself (or a part of
%     it) and 'ReferenceLevel' is the desired level
%   
%   The rms value is computed over the whole 'RefSignal': trim it before
%   (e.g. remove silence at the start of a calibrator recording).
%
% USAGE:
%
%   [CalibratedSignal,CalFactor,dBFS]=calibrate(InputSignal,RefSignal,ReferenceLevel)
%
%   INPUTS: InputSignal      - signal to be calibrated (no unit)
%           RefSignal        - reference signal (no unit)
%           ReferenceLevel   - rms level of RefSignal (dB SPL)
%
%   OUTPUT: CalibratedSignal - calibrated signal (Pa), CalFactor*InputSignal
%           CalFactor        - calibration factor (Pa per unit)
%           dBFS             - full-scale level (dB SPL), i.e. the SPL of a
%                              signal with an rms value of 1
%
% Examples using SQAT:
%   See the codes in the examples/signal_calibration folder
%
% Author: Gil Greco, Braunschweig, 31/01/2020
% Author: Alejandro Osses, code with extra comments and extra outputs
% Author: Gil Felix Greco, 26.11.2026 - header revised (description,
%   inputs, outputs, examples and copyright statement); the code is unchanged
%
% AI disclosure: header revision in November 2026 assisted by Claude
% Opus 5.5 (Anthropic). All codes were verified by the authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help calibrate;
    return;
end
% Note by Alejandro: This file is indeed used by the ISO 532-1:2017

CalFactor=sqrt( 10.^(ReferenceLevel./10).* (4e-10) ./ mean(RefSignal.^2));
															 
CalibratedSignal=CalFactor.*InputSignal;	

if nargout >=3
    dBFS = ReferenceLevel-20*log10(rms(RefSignal));
end
end
