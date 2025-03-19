function roughWeight = shmRoughWeight(modRate, modfreqMaxWeight, roughWeightParams)
% roughWeight = shmRoughWeight(modRate, modfreqMaxWeight, roughWeightParams)
%
% Returns roughness weighting for high- and low-frequency (modulation
% rates) according to ECMA-418-2:2024 (the Sottek Hearing Model) for a set
% of modulation rates and parameters.
%
%
% Inputs
% ------
%
% modRate : matrix
%           the estimated modulation rates used to determine the weighting
%           factors
%
% modfreqMaxWeight : the modulation rate at which the weighting reaches its
%                    maximum value (one)
%
% roughWeightParams : the parameters for the each of the weightings (high
%                     or low)
%
% Returns
% -------
%
%
% Assumptions
% -----------
% Inputs are in compatible parallelised forms
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 10/07/2024
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
% Equation 85 [G_l,z,i(f_p,i(l,z))]
roughWeight = 1./...
                (1 +...
                 ((modRate./modfreqMaxWeight...
                   - modfreqMaxWeight./modRate)...
                  .*roughWeightParams(1, :, :)).^2).^roughWeightParams(2, :, :);
