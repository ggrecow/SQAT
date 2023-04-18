function[CalibratedSignal,CalFactor,dBFS]=calibrate(InputSignal,RefSignal,ReferenceLevel)
% function[CalibratedSignal,CalFactor,dBFS]=calibrate(InputSignal,RefSignal,ReferenceLevel)
%
% This function adjusts the level of 'InputSignal' using the full scale 
%   convention given by 'RefSignal', which has a level equal to 
%   'ReferenceLevel'. 
%   This calibration function is based on the level adjustment method 
%   indicated in ISO 532-1:2017.
%
% USAGE: 
%   function[CalibratedSignal]=calibrate(InputSignal,RefSignal,ReferenceLevel) 
%
%   INPUTS: RefSignal        - reference calibration .wav (16-bit)
%           ReferenceLevel   - calibration level (dBSPL rms)
%           InputSignal      - .wav to be calibrated (16-bit)
%
%   OUTPUT: CalibratedSignal - calibrated .wav signal
%           CalFactor        - calibration factor used to scale InputSignal
%           dBFS             - value of the full scale value, expressed in
%                              dB SPL.
% Example using SQAT:
%   See validation_Sharpness_DIN45692_narrowband_and_broadband_signals.m
%
% Author: Gil Greco, Braunschweig, 31/01/2020
% Author: Alejandro Osses, code with extra comments and extra outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
