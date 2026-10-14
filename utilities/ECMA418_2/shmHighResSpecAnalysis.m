function [specError, constPart, lineAmp] = shmHighResSpecAnalysis(envSpectrum, nZerosBegin, nZerosEnd, modRates)
% [specError, constPart, lineAmp] = shmHighResSpecAnalysis(envSpectrum, nZerosBegin, nZerosEnd, modRates)
%
% Returns the High-resolution Spectral Analysis (HSA) of a windowed,
% downsampled envelope block according to ECMA-418-2:2025 (the Sottek
% Hearing Model), Section 9.1.4: the constant part and one spectral line
% pair per candidate modulation rate, fitted to the DFT of the envelope
% by solving the linear system of Equations 130 to 134, and the error
% function of Equation 135 (Equation 142 for a single line pair).
%
% Line amplitude convention: the line of rate f_c is x_2m + j*x_2m+1, the
% coefficient of the component at +f_c (physical amplitude 2*|line|).
% With this convention the argument of the nonlinearity in Equation 160,
% sqrt(0.5*(p0^2 + 2*sum(|line|^2))), equals the RMS value of the
% band-pass signal. Equation 123 prints the elements of x as half the real
% and imaginary parts of the line; that reading gives a line of
% 2*(x_2m + j*x_2m+1).
%
% The phase term of the window kernel of Equation 127 refers to the centre
% of the analysis window, (s_b - n_ze + n_zb - 1)/2, so that the kernel
% equals the DFT of the rectangular analysis window applied to a complex
% exponential at f_c.
%
% Inputs
% ------
%
% envSpectrum : vector
%   DFT of the windowed envelope block (s_b = 2048 at 1500 Hz), at least
%   for k = 0 to 48
%
% nZerosBegin : integer
%   number of zeros at the beginning of the analysis window [n_zb]
%
% nZerosEnd : integer
%   number of zeros at the end of the analysis window [n_ze]
%
% modRates : vector
%   candidate modulation rates (Hz) of the spectral line pairs [f_c]
%
% Returns
% -------
% specError : number
%   error function value [E_l,z(f_c)]
%
% constPart : number
%   constant part of the envelope [p_0]
%
% lineAmp : column vector
%   complex amplitudes of the spectral lines, one per modulation rate
%
% Assumptions
% -----------
% The envelope is downsampled to 1500 Hz and segmented in blocks of 2048
% samples, as for the fluctuation strength of Section 9
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
blockSize = 2048;  % [s_b_tilde]
sampleRate = 1500;  % [r_s_tilde]
resDFT = sampleRate/blockSize;  % [deltaf]

nRates = numel(modRates);

% Equation 125 [K_L]
nLines = min(max(17, round(max(modRates)/resDFT) + 8), 49);
kIdx = 0:nLines - 1;

% Equation 127 [W_E,l,z,f_c,m(k)] for the rates 0, +f_c and -f_c, with the
% normalised frequency f_n(k) = k/s_b - f_c/r_s + eps
nOnes = blockSize - nZerosEnd - nZerosBegin;
kernelRates = [0; modRates(:); -modRates(:)];
kernels = zeros(numel(kernelRates), nLines);
for iRate = 1:numel(kernelRates)
    normFreq = kIdx/blockSize - kernelRates(iRate)/sampleRate + eps;
    kernels(iRate, :) = exp(-2i*pi*normFreq*(blockSize - nZerosEnd + nZerosBegin - 1)/2)...
                        .*sin(pi*normFreq*nOnes)./sin(pi*normFreq);
end

% Equations 126, 128 and 129: columns [W_0, W_+, j*W_-] of the linear model
cols = zeros(nLines, 2*nRates + 1);
cols(:, 1) = kernels(1, :).';
for mRate = 1:nRates
    kernelPos = kernels(1 + mRate, :);
    kernelNeg = kernels(1 + nRates + mRate, :);
    cols(:, 2*mRate) = (kernelPos + kernelNeg).';
    cols(:, 2*mRate + 1) = 1i*(kernelPos - kernelNeg).';
end

envSpecK = envSpectrum(kIdx + 1);
envSpecK = envSpecK(:);

% Equations 130 to 134: the real Gram matrix of the columns implements the
% index structure of Equations 131 to 133
matA = real(cols'*cols);
vecB = real(cols'*envSpecK);
vecX = matA\vecB;

% Equation 135 [E_l,z(f_c)]
specError = sum(abs(envSpecK).^2) + (vecX.'*matA*vecX - 2*(vecB.'*vecX));

constPart = vecX(1);
lineAmp = vecX(2:2:end) + 1i*vecX(3:2:end);

end
