function fluctWeight = shmFluctWeight(modRate, bandCentreFreq)
% fluctWeight = shmFluctWeight(modRate, bandCentreFreq)
%
% Returns fluctuation strength weighting for low and high modulation rates
% (band-pass characteristic of the fluctuation strength) according to
% ECMA-418-2:2025 (the Sottek Hearing Model), Section 9.1.6, Equation 148,
% for a set of modulation rates in a critical band.
%
% The carrier frequency term applies to the branch of the modulation rates
% above f_max only, and the weighting is zero for a modulation rate of zero.
%
% Inputs
% ------
%
% modRate : vector
%   the modulation rates (Hz) used to determine the weighting factors
%
% bandCentreFreq : number
%   the centre frequency (Hz) of the critical band
%
% Returns
% -------
% fluctWeight : vector
%   the weighting values for the input modulation rates
%
% Assumptions
% -----------
% Modulation rates are non-negative
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Sergio Aguirre &
%          Gil Felix Greco
%
% Date created: 17.09.2026
% Date last modified: 17.09.2026
% MATLAB version: 2026a
%
% Copyright statement: This file is part of the SQAT toolbox and is
% subject to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt>
% in the SQAT repository root.
%
% As per the licensing information, this file is provided "as is",
% WITHOUT WARRANTY OF ANY KIND, express or implied, including but not
% limited to the warranties of MERCHANTABILITY and FITNESS FOR A
% PARTICULAR PURPOSE.
%
% Checked by:
% Date last checked:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation 148 parameters [q_1l, q_2l, q_1h, q_2h, f_max]
q1Low = 0.33048;
q2Low = 0.85902;
q1High = 0.21792;
q2High = 4.6728;
modfreqMaxWeight = 4.8659;  % Hz

% carrier frequency term of the higher modulation rates
freqTerm = 1/(1 + 0.092623*abs(log2(bandCentreFreq/1000))^1.24);

% Equation 148 [w_lh(f_c,i(l,z))]
fluctWeight = zeros(size(modRate));
lowRate = modRate > 0 & modRate <= modfreqMaxWeight;
highRate = modRate > modfreqMaxWeight;
fluctWeight(lowRate) = (1 + ((modRate(lowRate)/modfreqMaxWeight...
                              - modfreqMaxWeight./modRate(lowRate))...
                             *q1Low).^2).^(-q2Low);
fluctWeight(highRate) = freqTerm*(1 + ((modRate(highRate)/modfreqMaxWeight...
                                        - modfreqMaxWeight./modRate(highRate))...
                                       *q1High).^2).^(-q2High);

end
