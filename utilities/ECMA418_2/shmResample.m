function [resampledSignal, resampledRate]...
          = shmResample(signal, sampleRatein)
% [resampledSignal, resampledRate]...
%         = shmResample_(signal, sampleRatein)
%
% Returns signal resampled to 48 kHz, according to ECMA-418-2:2024
% (the Sottek Hearing Model) for an input signal.
%
% Inputs
% ------
% signal : vector or 2D matrix
%          the input signal
%
% sampleRatein : integer
%                the sample rate (frequency) of the input signal(s)
% 
% Returns
% -------
% For each channel in the input signal:
%
% resampledSignal : number or vector
%               average (overall) tonality value
% 
% resampledRate : integer
%                 the resampled signal sample rate, ie, 48 kHz
%
% Assumptions
% -----------
% The input signal is oriented with time on axis 1 (and channel # on axis
% 2), ie, the resample operation is applied along axis 1.
%
% Requirements
% ------------
% Signal Processing Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 26/09/2023
% Date last modified: 19/03/2025
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        signal double {mustBeReal}
        sampleRatein (1, 1) double {mustBePositive, mustBeInteger}
    end

%% Define constants

% Section 5.1.1 ECMA-418-2:2024
resampledRate = 48e3;  % Signal sample rate prescribed to be 48 kHz

%% Signal processing

% Input pre-processing
% --------------------
if sampleRatein ~= resampledRate  % Resample signal
    up = resampledRate/gcd(resampledRate, sampleRatein);  % upsampling factor
    down = sampleRatein/gcd(resampledRate, sampleRatein);  % downsampling factor
    resampledSignal = resample(signal, up, down);  % apply resampling
else  % don't resample
    resampledSignal = signal;
end

% end of function
