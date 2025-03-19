function signalOut = shmPreProc(signal, blockSize, hopSize, padStart, padEnd)
% signalFadePad = shmPreProc(signal, blockSize, hopSize, padStart, padEnd)
%
% Returns signal with fade-in and zero-padding pre-processing according to
% ECMA-418-2:2024 (the Sottek Hearing Model) for an input signal.
%
% Inputs
% ------
% signal : vector or 2D matrix
%          the input signal/s
%
% blockSize : integer
%             the maximum signal segmentation block size
%
% hopSize : integer
%           the maximum signal segmentation hop size
%           = (1 - overlap)*blockSize
%
% padStart : Boolean
%            flag to indicate whether to pad the start of the signal
%
% padEnd : Boolean
%          flag to indicate whether to pad the end of the signal
% 
% Returns
% -------
% signalFadePad : vector or 2D matrix
%                 the output faded, padded signal
%
% Assumptions
% -----------
% The input signal is oriented with time on axis 1 (and channel # on axis
% 2), ie, the fade and padding operation is applied along axis 1.
% The input signal must be sampled at 48 kHz.
%
% Requirements
% ------------
% None
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
        signal (:, :) double {mustBeReal}
        blockSize (1, 1) {mustBeInteger}
        hopSize (1, 1) {mustBeInteger}
        padStart {mustBeNumericOrLogical} = true
        padEnd {mustBeNumericOrLogical} = true
    end


%% Signal processing

% Input pre-processing
% --------------------
%
% Fade in weighting function Section 5.1.2 ECMA-418-2:2024
fadeWeight = repmat(transpose(0.5 - 0.5*cos(pi*(0:239)/240)), 1, size(signal, 2));
% Apply fade in
signalFade = [fadeWeight.*signal(1:240, :);
              signal(241:end, :)];

% Zero-padding Section 5.1.2 ECMA-418-2:2024
if padStart
    n_zeross = blockSize;  % start zero-padding
else
    n_zeross = 0;
end  % end of if branch to zero pad start

if padEnd
    n_samples = size(signal, 1);
    n_new = hopSize*(ceil((n_samples + hopSize + n_zeross)/hopSize) - 1);
    n_zerose = n_new - n_samples;  % end zero-padding
else
    n_zerose = 0;
end  % end of if branch to zero pad end

% Apply zero-padding
signalOut = [zeros(n_zeross, size(signalFade, 2));
                 signalFade;
                 zeros(n_zerose, size(signalFade, 2))];

% end of function
