function signalFiltered  = shmOutMidEarFilter(signal, fieldtype, outplot)
% signalFiltered  = shmOutMidEarFilter(signal, fieldtype, outplot)
%
% Returns signal filtered for outer and middle ear response according to
% ECMA-418-2:2024 (the Hearing Model of Sottek) for an input calibrated
% audio (sound pressure) time-series signal. Optional plot output shows
% frequency and phase response of the filter.
%
% Inputs
% ------
% signal : vector or 2D matrix
%          the input signal (sound pressure)
%
% fieldtype : keyword string (default: 'free-frontal')
%             determines whether the 'free-frontal' or 'diffuse' field 
%             stages are applied in the outer-middle ear filter
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a frequency and phase
%           response figure for the filter
% 
% Returns
% -------
% signalFiltered : vector or 2D matrix
%                  the output filtered signal
%
% Assumptions
% -----------
% The input signal is oriented with time on axis 1 (and channel # on axis
% 2), ie, the filtering operation is applied along axis 1.
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
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        signal (:, :) double {mustBeReal}
        fieldtype (1, 1) string {mustBeMember(fieldtype,...
                                              {'free-frontal'...
                                               'diffuse'})} = 'free-frontal'
        outplot {mustBeNumericOrLogical} = false
    end

%% Signal processing

% Apply outer & middle ear filter bank
% ------------------------------------
%
% Filter coefficients from Section 5.1.3.2 Table 1 ECMA-418-2:2024
% b_0k = [1.015896, 0.958943, 0.961372, 2.225804, 0.471735, 0.115267, 0.988029,...
%         1.952238];
% b_1k = [-1.925299, -1.806088, -1.763632, -1.43465, -0.366092, 0.0, -1.912434,...
%         0.16232]; 
% b_2k = [0.922118, 0.876439, 0.821788, -0.498204, 0.244145, -0.115267,...
%         0.926132, -0.667994];
% a_0k = ones(size(b_0k));
% a_1k = [-1.925299, -1.806088, -1.763632, -1.43465, -0.366092, -1.796003,...
%         -1.912434, 0.16232];
% a_2k = [0.938014, 0.835382, 0.78316, 0.727599, -0.28412, 0.805838, 0.914161,...
%         0.284244]; 
%
% Accurate coefficient values
b_0k = [1.015896020255593, 0.958943219304445, 0.961371976333197,...
        2.225803503609735, 0.471735128494163, 0.115267139824401,...
        0.988029297230954, 1.952237687301361];
b_1k = [-1.925298877776079, -1.806088011849494, -1.763632154338248,...
        -1.434650484792157, -0.366091796830044, 0.0,...
        -1.91243380293387, 0.162319983017519];
b_2k = [0.922118060364679, 0.876438777856084, 0.821787991845146,...
        -0.498204282194628, 0.244144703885020, -0.115267139824401,...
        0.926131550180785, -0.667994113035186];
a_0k = ones(size(b_0k));
a_1k = [-1.925298877776079, -1.806088011849494, -1.763632154338248,...
        -1.434650484792157, -0.366091796830044, -1.796002566692014,...
        -1.912433802933871, 0.162319983017519];
a_2k = [0.938014080620272, 0.835381997160530, 0.783159968178343,...
        0.727599221415107, -0.284120167620817, 0.805837815618546,...
        0.914160847411739, 0.284243574266175]; 

if fieldtype == "free-frontal"
    sos = [b_0k.', b_1k.', b_2k.', a_0k.', a_1k.', a_2k.'];
elseif fieldtype == "diffuse"
    % omit free field filter stages
    sos = [b_0k(3:end).', b_1k(3:end).', b_2k(3:end).',...
           a_0k(3:end).', a_1k(3:end).', a_2k(3:end).'];
end

% Section 5.1.3.2 ECMA-418-2:2024 Outer and middle/inner ear signal filtering
signalFiltered = sosfilt(sos, signal, 1);

%% Plot figures

if outplot
    [H, f] = freqz(sos, 10e3, 'whole', 48e3);
    phir = angle(H);
    phirUnwrap = unwrap(phir,[], 1);
    phiUnwrap = phirUnwrap/pi*180;
    % Plot frequency and phase response for filter
    fig = figure;
    tiledlayout(fig, 2, 1);
    movegui(fig, 'center');
    ax1 = nexttile(1);
    semilogx(f, 20*log10(abs(H)), 'color', [0.0, 0.2, 0.8])
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

    ax2 = nexttile(2);
    semilogx(f, phiUnwrap, 'color', [0.8, 0.1, 0.8]);
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

    title(ax1, fieldtype)
end


% end of function
