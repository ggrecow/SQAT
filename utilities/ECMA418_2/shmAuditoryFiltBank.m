function signalFiltered = shmAuditoryFiltBank(signal, outplot)
% signalFiltered = shmAuditoryFiltBank(signal, outplot)
%
% Returns a set of signals, bandpass filtered for the inner ear response
% in each half-Bark critical band rate scale width, according to
% ECMA-418-2:2024 (the Sottek Hearing Model) for an input calibrated
% single (sound pressure) time-series signal.
%
% Inputs
% ------
% signal : vector
%          the input signal as single audio (sound pressure) signal
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure a frequency and phase
%           response figure for the filter bank
% 
% Returns
% -------
% 
% signalFiltered : 2D matrix
%                  the filtered signals 
%
% Assumptions
% -----------
% The input signal is oriented with time on axis 1, ie, the filter
% operation is applied along axis 1.
% The input signal must be sampled at 48 kHz.
%
% Requirements
% ------------
% Signal Processing Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 27/09/2023
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
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        signal (:, 1) double {mustBeReal}
        outplot {mustBeNumericOrLogical} = false
    end

%% Define constants

sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2024
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2024
c = 0.1618;  % Half-Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2024

halfBark = 0.5:0.5:26.5;  % half-critical band rate scale
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2024
dfz = sqrt(deltaFreq0^2 + (c*bandCentreFreqs).^2);  % Section 5.1.4.1 Equation 10 ECMA-418-2:2024


%% Signal processing

% Apply auditory filter bank
% --------------------------
% Filter equalised signal using 53 1/2Bark ERB filters according to 
% Section 5.1.4.2 ECMA-418-2:2024

k = 5;  % filter order = 5, footnote 5 ECMA-418-2:2024
e_i = [0, 1, 11, 11, 1];  % filter coefficients for Section 5.1.4.2 Equation 15 ECMA-418-2:2024

for zBand = 53:-1:1
    % Section 5.1.4.1 Equation 8 ECMA-418-2:2024
    tau = (1/(2^(2*k - 1))).*nchoosek(2*k - 2, k - 1).*(1./dfz(zBand));
    
    d = exp(-1./(sampleRate48k.*tau)); % Section 5.1.4.1 ECMA-418-2:2024
    
    % Band-pass modifier Section 5.1.4.2 Equation 16/17 ECMA-418-2:2024
    bp = exp((1i.*2.*pi.*bandCentreFreqs(zBand).*(0:k+1))./sampleRate48k);
    
    % Feed-backward coefficients, Section 5.1.4.2 Equation 14 ECMA-418-2:2024
    m = 1:k;
    a_m = ([1, ((-d).^m).*arrayfun(@(m_) nchoosek(k, m_), m)]).*bp(1:k+1);

    % Feed-forward coefficients, Section 5.1.4.2 Equation 15 ECMA-418-2:2024
    m = 0:k-1;
    i = 1:k-1;
    b_m = ((((1-d).^k)./sum(e_i(i+1).*(d.^i))).*(d.^m).*e_i).*bp(1:k);

    % Recursive filter Section 5.1.4.2 Equation 13 ECMA-418-2:2024
    % Note, the results are complex so 2x the real-valued band-pass signal
    % is required.
    signalFiltered(:, zBand) = 2*real(filter(b_m, a_m, signal));

    %% Plot figures

    if outplot
        [H, f] = freqz(b_m, a_m, 10e3, 'whole', 48e3);
        phir = angle(H);
        phirUnwrap = unwrap(phir,[], 1);
        phiUnwrap = phirUnwrap/pi*180;
        % Plot frequency and phase response for filter
        if zBand == 53
            fig = figure;
            tiledlayout(fig, 2, 1);
            movegui(fig, 'center');
            ax1 = nexttile(1);
            ax2 = nexttile(2);
            hold(ax1, 'on')
            hold(ax2, 'on')
        end
        semilogx(ax1, f, 20*log10(abs(H)))
        semilogx(ax2, f, phiUnwrap);

        if zBand == 1
            ax1.XLim = [20, 20e3];
            ax1.XTick = [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
            ax1.XMinorTick = "off";
            ax1.XTickLabel = ["31.5", "63", "125", "250", "500", "1k", "2k", "4k",...
                              "8k", "16k"];
            ax1.XLabel.String = "Frequency, Hz";
            ax1.YLabel.String = "\it{H}\rm, dB";
            ax1.FontName =  'Arial';
            ax1.FontSize = 12;
            ax1.XGrid = 'on';
            ax1.YGrid = 'on';
    
            ax2.XLim = [20, 20e3];
            ax2.XTick = [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3];
            ax2.XMinorTick = "off";
            ax2.XTickLabel = ["31.5", "63", "125", "250", "500", "1k", "2k", "4k",...
                              "8k", "16k"]; 
            ax2.XLabel.String = "Frequency, Hz";
            ax2.YLabel.String = "Phase angle \circ";
            ax2.XGrid = 'on';
            ax2.YGrid = 'on';
            ax2.FontName = 'Arial';
            ax2.FontSize = 12;
        end
    end


end

% end of function
