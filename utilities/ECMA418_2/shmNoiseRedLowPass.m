function signalFiltered = shmNoiseRedLowPass(signal, sampleRateIn)
% signalFiltered = shmNoiseRedLowPass(signal, sampleRateIn)
%
% Returns signal low pass filtered for noise reduction according to
% ECMA-418-2:2025 (the Sottek Hearing Model) for an input signal.
%
% Inputs
% ------
% signal : vector or 2D matrix
%   the input signal as single mono or stereo audio (sound pressure) signals
%
% sampleRateIn : double
%   the sample rate (frequency) of the input signal(s)
% 
% Returns
% -------
% signalFiltered : vector or 2D matrix
%   the filtered signal
%
% Assumptions
% -----------
% The input signal is oriented with time on axis 1 (and channel # on axis
% 2), ie, the filtering operation is applied along axis 1.
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (m.c.torjussen@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 22/09/2023
% Date last modified: 15/11/2025
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
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by:
% Date last checked:
%
% Arguments validation
    arguments (Input)
        signal (:, :) double {mustBeReal}
        sampleRateIn (1, 1) double {mustBePositive}
    end

k = 3; % Footnote 21 ECMA-418-2:2025
e_i = [0, 1, 1]; % Footnote 21 ECMA-418-2:2025

% Footnote 20 ECMA-418-2:2025
tau = 1/32*6/7;

d = exp(-1/(sampleRateIn*tau)); % Section 5.1.4.2 ECMA-418-2:2025

% Feed-backward coefficients, Equation 14 ECMA-418-2:2025
m_a = 1:k;
a = [1, ((-d).^m_a).*arrayfun(@(m_) nchoosek(k, m_), m_a)];

% Feed-forward coefficients, Equation 15 ECMA-418-2:2025
m_b = 0:k-1;
i = 1:k-1;
b = (((1 - d)^k)./sum(e_i(i + 1).*(d.^i))).*(d.^m_b).*e_i;

% Recursive filter Equation 13 ECMA-418-2:2025
signalFiltered = filter(b, a, signal);

end

