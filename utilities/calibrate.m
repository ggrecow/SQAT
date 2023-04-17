%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   calibration function based on ISO 532-1:2017  
%
% USAGE: 
%   function[CalibratedSignal]=calibrate(InputSignal,RefSignal,ReferenceLevel) 
%
%   INPUTS: RefSignal        - reference calibration .wav (16-bit)
%           ReferenceLevel   - calibration level (dBSPL rms)
%           InputSignal      - .wav to be calibrated (16-bit)
%
%   OUTPUT: CalibratedSignal - calibrated .wav signal
%
%   Gil Greco, Braunschweig, 31/01/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function[CalibratedSignal]=calibrate(InputSignal,RefSignal,ReferenceLevel) 

CalFactor=sqrt( 10.^(ReferenceLevel./10).* (4e-10) ./ mean(RefSignal.^2));
															 
CalibratedSignal=CalFactor.*InputSignal;	

return 
 
