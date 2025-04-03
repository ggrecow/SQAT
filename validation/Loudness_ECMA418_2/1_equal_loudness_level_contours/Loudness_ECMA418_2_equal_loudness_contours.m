% Script Loudness_ECMA418_2_equal_loudness_contours
%
% Recreate Fig. A.1 from ECMA-418-2:2024. Compute Loudness (ECMA 418-2:2024) using the
% <Loudness_ECMA418_2.m> implementation in SQAT and
% compare with the equal-loudness-level contours from ISO 226:2003
%
%   1) compute reference equal-loudness-level contours using formulation from ISO 226:2003
%
%   2) compute overall loudness of a 1-kHz tone, for different phon curves,
%       using the model from ECMA418-2:2024 
%
%   3) perform an optimization procedure to find the SPL of tones (with
%       different frquencies) with the same loudness as the 1-kHz tone.
%       Here, the loudness is calculated using the model from ECMA418-2:2024 
%
%   Loudness using the model from ECMA418-2:2024 is computed using the
%   following function:
%   OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Loudness_ECMA418_2> for more info
%
%   All signals are generated within this code.
%
% Author: Gil Felix Greco, Braunschweig 26.01.2025
% Updated results: Gil Felix Greco, Cala Rajada 21.03.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings

dir_analysis_name = '1_equal_loudness_level_contours';
dir_out = [fileparts(mfilename('fullpath')) filesep];

% Figure where the figures (and the results) will be stored:
figures_dir = [dir_out 'figs' filesep];
if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

fname_res = 'spl_ecma_vec.mat'; % name for the result file
fname_res_full = [figures_dir fname_res];

bCalculation = ~exist(fname_res_full,'file');
if bCalculation == 0
    fprintf('%s.m: Results file found on disk, confirm that you want to load those prestored\n');
    fprintf('\t results or whether you want to re-run the calculations.\n');
    bCalculation = input('Enter your choice (1=re-run; 0=read stored results): ');
end
bLoad = ~bCalculation;

%% Compute equal-loudness-contours according to ISO 226:2003 (reference curves)

% loudness level (phon)
Ln = [80 60 40 20];

% equal-loudness-level contours from ISO 226:2003 for each Ln
reference_curve(1,:) = il_equal_loudness_contours( Ln (1:length(Ln) ) );

%% Compute equal-loudness-contours using loudness model from ECMA

if bCalculation

    freq_vec=reference_curve.frequency; % Center frequency (Hz)

    spl_ecma_vec = zeros(size(freq_vec,1), length(Ln));     % define ECMA vectors
    spl_ecma_vec(18, 1:length(Ln)) = Ln; % set spl of 1-kHz tones

    for p = 1:length(Ln) % Loop through each phon level

        fprintf('%d Phon equal loudness contour\n', Ln(p));

        % Compute ECMA loudness of a 1 kHz tone
        N_1khz = il_comp_loudness_1khz( Ln(p) ); % loudnessPowAvg (sone), of 1-kHz tone 

        % define range of optimization variable
        range = [Ln(p)-10 Ln(p)+40];

        % Iterate through frequencies
        for i = 7:length(freq_vec) % start at 80 Hz
            if i ~= 18 % fc = 1 kHz

                fprintf('\tFreq: %.2f Hz\n',  freq_vec(i));

                % Find SPL of the tone with same loudness as 1 kHz tone
                spl_ecma_vec(i,p) = fzero(@(x) il_comp_loudness(freq_vec(i), x, N_1khz), range);

            end
        end

    end % end Loop through each phon level

    save(fname_res_full, 'spl_ecma_vec');

end % end bCalculation

%% if available and desired, load equal-loudness-contours computed using loudness model from ECMA

if bLoad % load saved results
    load(fname_res_full);
end

%% Plot

color_80phon=[0.00,0.45,0.749]; % blue
color_60phon=[0.85,0.33,0.10]; % red 
color_40phon=[0.93,0.69,0.13]; % yellow
color_20phon=[0.49,0.18,0.56];  %purple
line_ISO = 1.5; % linewidth of reference results

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% 80 phon
a = semilogx( reference_curve.frequency, spl_ecma_vec(:,1), '-', 'Color', color_80phon, 'LineWidth', line_ISO ); hold all; % SQAT (loudness model from ECMA-418-2)
b = semilogx( reference_curve.frequency, reference_curve.Lp(:,1), 'k:', 'LineWidth', line_ISO); % ISO 226:2003 (reference curves)

% 60 phon
semilogx( reference_curve.frequency, spl_ecma_vec(:,2), '-', 'Color', color_60phon, 'LineWidth', line_ISO ); % SQAT (loudness model from ECMA-418-2)
semilogx( reference_curve.frequency, reference_curve.Lp(:,2), 'k:', 'LineWidth', line_ISO); % ISO 226:2003 (reference curves)

% 40 phon
semilogx( reference_curve.frequency, spl_ecma_vec(:,3), '-', 'Color', color_40phon, 'LineWidth', line_ISO ); % SQAT (loudness model from ECMA-418-2)
semilogx( reference_curve.frequency, reference_curve.Lp(:,3), 'k:', 'LineWidth', line_ISO); % ISO 226:2003 (reference curves)

% 20 phon
semilogx( reference_curve.frequency, spl_ecma_vec(:,4), '-', 'Color', color_20phon, 'LineWidth', line_ISO ); % SQAT (loudness model from ECMA-418-2)
semilogx( reference_curve.frequency, reference_curve.Lp(:,4), 'k:', 'LineWidth', line_ISO); % ISO 226:2003 (reference curves)

xlabel( 'Frequency (Hz)','Interpreter','Latex' );
ylabel( 'Sound pressure level (dB SPL)','Interpreter','Latex' );

ylim( [-15 105] );
xlim( [80 11000] );

legend([b,a], {'ISO 226:2003', 'ECMA-418-2:2024 (SQAT)'}, 'Location','southwest' );

% Create textbox
annotation(h,'textbox',...
    [0.408142857142853 0.438095238095241 0.182928571428575 0.0214285714285715],...
    'Color', color_20phon,...
    'String','20 phon',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(h,'textbox',...
    [0.408142857142853 0.492857142857143 0.182928571428576 0.100000000000006],...
    'Color', color_40phon,...
    'String','40 phon',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(h,'textbox',...
    [0.408142857142853 0.697619047619056 0.182928571428577 0.0214285714285717],...
    'Color', color_60phon,...
    'String','60 phon',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(h,'textbox',...
    [0.408142857142853 0.823809523809534 0.182928571428577 0.0214285714285717],...
    'Color', color_80phon,...
    'String','80 phon',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

legend boxoff
set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    figname_short = 'Loudness_ECMA418_2_equal_loudness_contours';
    figname_out = [figures_dir figname_short];

    % saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function - get   equal_loudness_contours according to  ISO 226:2003 (Table B1)

function OUT = il_equal_loudness_contours(Ln)
% Function OUT = il_get_equal_loudness_contours(Ln)
%
% Given a desired phon level, Ln, compute equal-loudness-levels contours using
% the formulation from ISO 226:2003, Section 4.1. The equal-loudness-levels contours
% are given in sound pressure levels, Lp (dB SPL), as a function of the
% frequency, in Hz (29 values from 20 Hz - 12.5 kHz)
%
% INPUTS:
%           Ln : scalar
%                  Desired phon level
%
% OUTPUT:
%           OUT : struct
%                      contains the following fields
%
%           Lp : [Nx1] column vector
%                  sound pressure level, in dB SPL
%
%            frequency : [Nx1] column vector
%                             frequency, in Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency vector
frequency = [ 20; 25; 31.5; 40; 50; 63; 80; 100; 125; 160; ...
    200; 250; 315; 400; 500; 630; 800; 1000; 1250; 1600; ...
    2000; 2500; 3150; 4000; 5000; 6300; 8000; 10000; 12500 ];

% Coefficients from Table 1
af = [ 0.532; 0.506; 0.480; 0.455; 0.432; 0.409; 0.387; 0.367; 0.349; 0.330; ...
    0.315; 0.301; 0.288; 0.276; 0.267; 0.259; 0.253; 0.250; 0.246; 0.244;
    0.243; 0.243; 0.243; 0.242; 0.242; 0.245; 0.254; 0.271; 0.301 ];

Lu = [ -31.6; -27.2; -23.0; -19.1; -15.9; -13.0; -10.3; -8.1; -6.2; -4.5; -3.1; ...
    -2.0; -1.1; -0.4; 0.0; 0.3; 0.5; 0.0; -2.7; -4.1; -1.0; ...
    1.7; 2.5; 1.2; -2.1; -7.1; -11.2; -10.7; -3.1];

Tf = [ 78.5; 68.7; 59.5; 51.1; 44.0; 37.5; 31.5; 26.5; 22.1; 17.9; ...
    14.4; 11.4; 8.6; 6.2; 4.4; 3.0; 2.2; 2.4; 3.5; 1.7; -1.3; ...
    -4.2; -6.0; -5.4; -1.5; 6.0; 12.6; 13.9; 12.3 ];

% ISO 226:2003 - 4.1. Derive sound pressure level (Lp = dB SPL) from loudness level (Ln = phone)
Af = 4.47e-3 .* ( (10.^(0.025.*Ln) ) - 1.15) + ( (0.4 .* ( 10.^( ( (Tf + Lu) ./ 10 ) - 9 ) ) ).^af ) ;

Lp = ( (10 ./ af) .* log10( Af ) ) - Lu + 94; % Eq. (1)

OUT.frequency = frequency;
OUT.Lp = Lp;

end % end of <il_equal_loudness_contours> subfunction


%% function - get pure tone signal

function OUT = il_generate_pure_tone(fc, Lp, length, fs)
%function OUT = il_generate_pure_tone(fc, Lp, length, fs)
%
% Generates pure tone signal with center frequency, fc (Hz), and a
% (RMS) sound pressure level, Lp (dBSPL).
%
%   INPUTS:
%   fc : Center frequency (Hz)
%   Lp : RMS SPL (dBSPL)
%   Length : signal length, in (seconds)
%   fs : Sampling frequency (Hz)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify signal characteristics

dt=1/fs;    % Time step
time = 0:dt:length; % Time vector
time =  time(:);

%% Make signal

pref = 20e-6; % reference sound pressure for air
A=pref.*10.^(Lp/20).*sqrt(2); % Amplitude

signal=A.*sin(2.*pi.*fc.*time);     % Generate signal
signal = signal(:);
spl_signal=20.*log10(rms(signal)./pref); % Verify dBSPL of the generated signal

%% output assignment

OUT.signal = signal;
OUT.spl_signal = spl_signal;
OUT.fc= fc;
OUT.fs= fs;
OUT.length=length;
OUT.time= time; % time vector

end % end of <il_generate_pure_tone> subfunction


%% function - get loudness of (reference) 1-kHz tone using ECMA loudness model

function out = il_comp_loudness_1khz(spl)
%function out = il_comp_loudness_1khz(spl)
%
% Compute loudness of a 1-kHz tone with given RMS <spl> (dB SPL)
% using the ECMA-418-2:2024 model
%
% FUNCTION:
%   OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Loudness_ECMA418_2> for more info
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal generation parameters
duration = 0.5;
fs = 48000;
f_1khz = 1000;

% Generate 1-kHz tone
signal_1kHz = il_generate_pure_tone(f_1khz, spl, duration, fs);

% Compute loudness
fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3;% time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 0; % show results, 'false' (disable, default value) or 'true' (enable)

N_1kHz = Loudness_ECMA418_2(signal_1kHz.signal, fs, fieldtype, time_skip, show);

out = N_1kHz.loudnessPowAvg;
end % end of <il_comp_loudness_1khz> subfunction

%% function - get loudness of test tones using ECMA loudness model
% and output the loudness difference between test tones and 1-kHz tone

function n_diff = il_comp_loudness(freq, spl,  N_1kHz)
%function n_diff = il_comp_loudness(freq, spl,  N_1kHz)
%
% Return the difference between the loudness of a tone at <freq> (Hz)
% and <spl> (dB SPL) and the loudness  <N_1kHz> (sone) of reference tone at 1 (KHz).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal generation parameters
duration = 0.5;
fs = 48000;

% Generate test signal
signal = il_generate_pure_tone(freq, spl, duration, fs);

% Compute loudness
fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3;% time_skip, in seconds for statistical calculations (default: 304ms - avoids transient responses of the digital filters)
show = 0; % show results, 'false' (disable, default value) or 'true' (enable)

N = Loudness_ECMA418_2(signal.signal, fs, fieldtype, time_skip, show);

% Compute difference
n_diff = N.loudnessPowAvg - N_1kHz;

end % end of <il_comp_loudness> subfunction

