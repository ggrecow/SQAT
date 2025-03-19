function specRoughness = shmRoughLowPass(specRoughEstTform, sampleRate, ...
                                         riseTime, fallTime)
% specRoughness = shmRoughLowPass(specRoughEstTform, sampleRatein, riseTime, fallTime)
%
% Returns specific roughness low pass filtered for smoothing according to
% ECMA-418-2:2024 (the Sottek Hearing Model) for an input transformed
% estimate of the specific roughnesss.
%
% Inputs
% ------
% specRoughEstTform : 2D matrix
%                         the input specific roughness estimate (from
%                         Equation 104)
%
% sampleRate : double
%              the sample rate (frequency) of the input specific
%              roughness (NB: this is not the original signal sample
%              rate; currently it should be set to 50 Hz)
% 
% Returns
% -------
% specRoughness : 2D matrix
%                 the filtered specific roughness
%
% Assumptions
% -----------
% The input specific roughness estimate is orientated with time on axis 1,
% and critical bands on axis 2.
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 22/09/2023
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
% Arguments validation
    arguments (Input)
        specRoughEstTform (:, :) double {mustBeReal}
        sampleRate (1, 1) double {mustBePositive}
        riseTime (1, 1) double {mustBePositive}
        fallTime (1, 1) double {mustBePositive}
    end

riseExponent = exp(-1/(sampleRate*riseTime))*ones([1, size(specRoughEstTform, 2)]);
fallExponent = exp(-1/(sampleRate*fallTime))*ones([1, size(specRoughEstTform, 2)]);

specRoughness = specRoughEstTform;

for llBlock = 2:size(specRoughEstTform, 1)

    riseMask = specRoughEstTform(llBlock, :) >= specRoughness(llBlock - 1, :);
    fallMask = ~riseMask;


    if ~isempty(specRoughEstTform(llBlock, riseMask))
        specRoughness(llBlock, riseMask) = specRoughEstTform(llBlock, riseMask).*(1 - riseExponent(riseMask))...
                                           + specRoughness(llBlock - 1, riseMask).*riseExponent(riseMask);
    end
    if ~isempty(specRoughEstTform(llBlock, fallMask))
        specRoughness(llBlock, fallMask) = specRoughEstTform(llBlock, fallMask).*(1 - fallExponent(fallMask))...
                                           + specRoughness(llBlock - 1, fallMask).*fallExponent(fallMask);
    end

    

end

