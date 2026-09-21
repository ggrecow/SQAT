% Script run_validation_tonality_tone_in_noise
%
%  Verification of the tonality implementation according to ECMA-418-2:2025
%  Verification case: 1 kHz sinusoidal tone mixed with pink noise, which is
%  the case used to evaluate the psychoacoustic tonality in Annex B.2 of the
%  standard. This script reproduces Figure B.3.
%
% - Inputs, ref. data: 'Model results' curve and listening test results (mean
%   values and 95% confidence intervals) digitised from Figure B.3 of
%   ECMA-418-2:2025. Annex B.2 describes the listening tests: 16 test subjects
%   rated the tonality of each sound on a 13-point categorical scale, and the
%   mean ratings were mapped to tonality units through a linear scaling factor
%   derived by minimising the root-mean-square error over the five experiments.
%   The digitised values are stored in
%   <reference_values/ECMA418_2_FigB3.mat>, with the reading resolution
%   reported in the <description> field of that file.
%
% - Inputs, signals: mixtures of a 1 kHz sinusoidal tone at
%   Lp,tone = [55 60 65 70 75] dB SPL with pink noise at
%   Lp,noise = 40:5:80 dB SPL, i.e. 45 signals covering the five experiments
%   of Annex B.2. All signals are generated within this code (refer to the
%   subfunctions at the end of this script to see how). Annex B.2 specifies
%   the frequency of the tone and the two level ranges. The remaining
%   properties of the stimuli are left open by the standard, and the choices
%   made here are: 2 s duration, 48 kHz sampling frequency, pink noise
%   band-limited between 20 Hz and 20 kHz, one single noise realisation
%   (fixed random seed) shared by the 45 mixtures, and signals gated on and
%   off without a ramp.
%
%   The prominence ratio curve of Figure B.3 is omitted here, because the
%   prominence ratio is specified in ECMA-418-1 and is outside the scope of
%   the present implementation. The error measures reported in Annex B.2 for
%   the five experiments are 0,21 for the psychoacoustic tonality, 0,70 for
%   the prominence ratio and 0,74 for the tone-to-noise ratio, all related to
%   the 13-point categorical scale.
%
% Tonality computed using:
%   OUT = Tonality_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%   type <help Tonality_ECMA418_2> for more info
%
% Authors: Sergio Aguirre & Gil Felix Greco, 18.09.2026
%
% AI disclosure: modifications performed in September 2026 were assisted 
% by Claude Fable 5.1 and Opus 5 (Anthropic). All codes were verified by 
% the authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];

% Folder where the figures (and the results) will be stored:
figures_dir = [dir_out 'figs' filesep];
if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

fname_res = 'tonality_results.mat'; % name for the result file
fname_res_full = [figures_dir fname_res];

bCalculation = ~exist(fname_res_full,'file');
if bCalculation == 0
    fprintf('%s.m: Results file found on disk!\n', mfilename);
    fprintf('Do you want to load those prestored results or re-run the calculations?\n');
    bCalculation = input('Enter your choice (1=re-run; 0=read stored results): ');
end
bLoad = ~bCalculation;

%% reference data (ECMA-418-2:2025, Figure B.3)

ref_file = [dir_out 'reference_values' filesep 'ECMA418_2_FigB3.mat'];
load(ref_file); % loads <ref_FigB3>

Ltone = ref_FigB3.Ltone;    % level of the 1 kHz tone (dB SPL)
Lnoise = ref_FigB3.Lnoise;  % level of the pink noise (dB SPL)

%% signal and analysis parameters

fc = 1000;              % frequency of the sinusoidal tone (Hz)
fs = 48000;             % sampling frequency (Hz)
duration = 2;           % signal duration (s)
noise_band = [20 20000];% band limits of the pink noise (Hz)
seed = 1;               % seed of the pink noise realisation

fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
time_skip = 304e-3; % time_skip, in seconds for statistical calculations (default: 304 ms - avoids transient responses of the digital filters)

%% compute tonality of the tone in noise mixtures

if bCalculation

    % unit-RMS waveforms, generated once and scaled to each level below
    tone = il_generate_pure_tone(fc, duration, fs);
    noise = il_generate_pink_noise(noise_band, duration, fs, seed);

    tonality_SQAT = zeros(length(Ltone), length(Lnoise));

    tic
    for i = 1:length(Ltone)
        for j = 1:length(Lnoise)

            insig = il_set_level(tone, Ltone(i)) + il_set_level(noise, Lnoise(j));

            OUT = Tonality_ECMA418_2(insig, fs, fieldtype, time_skip);

            tonality_SQAT(i,j) = OUT.tonalityAvg;

            fprintf('%s.m: Lp,tone = %d dB SPL, Lp,noise = %d dB SPL, T = %.3f tu_HMS\n', ...
                mfilename, Ltone(i), Lnoise(j), tonality_SQAT(i,j));
        end
    end
    t_calculation = toc/60; % time to compute the 45 signals, in minutes

    %% saving results so the tonality calculation does not need to be run again

    save(fname_res_full, 'tonality_SQAT', 't_calculation', 'fc', 'fs', ...
         'duration', 'noise_band', 'seed', 'fieldtype', 'time_skip');

end
if bLoad
    load(fname_res_full); % the stored parameters replace the ones set above
    fprintf('%s.m: stored results computed with fs = %d Hz, duration = %g s, pink noise band = [%g %g] Hz, seed = %d, field type = %s\n', ...
        mfilename, fs, duration, noise_band(1), noise_band(2), seed, fieldtype);
end

%% comparison with the reference data

% deviation from the model results curve of Figure B.3
dev = tonality_SQAT - ref_FigB3.tonality_model;

% agreement with the listening test results, using the error measure defined
% in Annex B.2: the error of a point lying inside the 95% confidence interval
% is zero, and the error of a point lying outside it is the distance to the
% closest bound of the interval. Here the measure is expressed in tu_HMS,
% because the scaling factor between the 13-point categorical scale and the
% tonality units is not reported in the standard
err_SQAT = il_confidence_interval_error(tonality_SQAT, ...
    ref_FigB3.tonality_test_ci_lower, ref_FigB3.tonality_test_ci_upper);
err_ECMA = il_confidence_interval_error(ref_FigB3.tonality_model, ...
    ref_FigB3.tonality_test_ci_lower, ref_FigB3.tonality_test_ci_upper);

inside = tonality_SQAT >= ref_FigB3.tonality_test_ci_lower & ...
         tonality_SQAT <= ref_FigB3.tonality_test_ci_upper;

fprintf('\n%s.m: deviation from the model results curve of Figure B.3\n', mfilename);
fprintf('\tLp,tone (dB SPL) : rms (tu_HMS) : max abs (tu_HMS)\n');
for i = 1:length(Ltone)
    fprintf('\t%15d : %11.3f : %15.3f\n', Ltone(i), ...
        sqrt(mean(dev(i,:).^2)), max(abs(dev(i,:))));
end
fprintf('\t%15s : %11.3f : %15.3f\n', 'all', sqrt(mean(dev(:).^2)), max(abs(dev(:))));

fprintf('\n%s.m: comparison with the listening test results\n', mfilename);
fprintf('\tpoints inside the 95%% confidence interval: %d of %d (SQAT), %d of %d (Figure B.3)\n', ...
    sum(inside(:)), numel(inside), ...
    sum(sum(ref_FigB3.tonality_model >= ref_FigB3.tonality_test_ci_lower & ...
            ref_FigB3.tonality_model <= ref_FigB3.tonality_test_ci_upper)), numel(inside));
fprintf('\terror measure of Annex B.2: %.3f tu_HMS (SQAT), %.3f tu_HMS (Figure B.3)\n', ...
    err_SQAT, err_ECMA);

%% plot results, one figure per level of the sinusoidal tone

color_ref = [0.00 0.45 0.74]; % blue, as the model results curve of Figure B.3

for i = 1:length(Ltone)

    h = figure;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    % listening test results (mean and 95% confidence interval)
    a = errorbar(Lnoise, ref_FigB3.tonality_test_mean(i,:), ...
        ref_FigB3.tonality_test_mean(i,:) - ref_FigB3.tonality_test_ci_lower(i,:), ...
        ref_FigB3.tonality_test_ci_upper(i,:) - ref_FigB3.tonality_test_mean(i,:), ...
        'k-*', 'MarkerSize', 6, 'Linewidth', 0.5); hold all;

    % model results curve of Figure B.3
    b = plot(Lnoise, ref_FigB3.tonality_model(i,:), '--', ...
        'Color', color_ref, 'Linewidth', 1.5);

    % tonality computed with the implementation in SQAT
    c = plot(Lnoise, tonality_SQAT(i,:), 'ko:', 'MarkerSize', 8);

    legend([a b c], {'Listening test (mean $\pm$ 95 \% CI)', ...
        'ECMA-418-2:2025 (Fig. B.3)', 'SQAT'}, ...
        'Location', 'Best', 'Interpreter', 'Latex');
    legend boxoff

    axis([38 82 0 4]);

    ax = gca;
    set(ax,'XTick', Lnoise);
    set(ax,'YTick', 0:1:4);
    ax.YAxis.MinorTick = 'on';
    ax.YAxis.MinorTickValues = 0:0.25:4;

    title(sprintf('$f_{\\mathrm{c}}=%g$~kHz, $L_{\\mathrm{p,tone}}=%d$~dB~SPL', fc/1000, Ltone(i)), ...
        'Interpreter', 'Latex');
    ylabel('Tonality, $T$ (tu$_{\mathrm{HMS}}$)','Interpreter','Latex');
    xlabel('Level of the pink noise, $L_{\mathrm{p,noise}}$ (dB SPL)','Interpreter','Latex');

    set(gcf,'color','w');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if save_figs==1
        figname_short = sprintf('validation_tonality_tone_in_noise_%ddB', Ltone(i));
        figname_out = [figures_dir figname_short];

        % saveas(gcf,figname_out, 'fig');
        % saveas(gcf,figname_out, 'pdf');
        saveas(gcf,figname_out, 'png');

        fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% function - get pure tone signal with unit RMS

function outsig = il_generate_pure_tone(fc, duration, fs)
%function outsig = il_generate_pure_tone(fc, duration, fs)
%
% Generates a sinusoidal signal with frequency, fc (Hz), and unit RMS value.
%
%   INPUTS:
%   fc : frequency of the tone (Hz)
%   duration : signal length (s)
%   fs : sampling frequency (Hz)
%
%   OUTPUT:
%   outsig : [Nx1] column vector, unit RMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = (0:round(duration*fs)-1).'/fs; % time vector

outsig = sqrt(2).*sin(2.*pi.*fc.*time);

end % end of <il_generate_pure_tone> subfunction

%% function - get pink noise signal with unit RMS

function outsig = il_generate_pink_noise(band, duration, fs, seed)
%function outsig = il_generate_pink_noise(band, duration, fs, seed)
%
% Generates a pink noise signal with unit RMS value, band-limited between
% band(1) and band(2) (Hz). The signal is synthesised in the frequency
% domain, where the magnitude of the spectral components inside the band is
% proportional to 1/sqrt(f), which gives a power spectral density
% proportional to 1/f, and the phases are drawn from a uniform distribution.
% The phases come from a private random stream initialised with <seed>, so
% that the noise realisation is reproducible and the random number generator
% of the caller is left untouched.
%
%   INPUTS:
%   band : [1x2] vector, lower and upper band limits (Hz)
%   duration : signal length (s)
%   fs : sampling frequency (Hz)
%   seed : seed of the random number generator
%
%   OUTPUT:
%   outsig : [Nx1] column vector, unit RMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = round(duration*fs);
N = N + mod(N,2); % even number of samples

nHalf = N/2 + 1;                  % number of non-negative frequencies
freq = (0:nHalf-1).'*fs/N;        % frequency vector (Hz)

magnitude = zeros(nHalf,1);
idx = (freq >= band(1)) & (freq <= band(2));
magnitude(idx) = 1./sqrt(freq(idx));

stream = RandStream('mt19937ar','Seed',seed);
phase = 2*pi*rand(stream,nHalf,1);

spectrum = magnitude.*exp(1i*phase);
spectrum(1) = 0;                  % no dc component
spectrum(end) = real(spectrum(end)); % real-valued Nyquist component

outsig = real(ifft([spectrum; conj(flipud(spectrum(2:end-1)))], N, 1));
outsig = outsig(1:round(duration*fs));
outsig = outsig./rms(outsig);

end % end of <il_generate_pink_noise> subfunction

%% function - scale a unit-RMS signal to a given sound pressure level

function outsig = il_set_level(insig, Lp)
%function outsig = il_set_level(insig, Lp)
%
% Scales a unit-RMS signal to the RMS sound pressure level Lp (dB SPL),
% returning the signal in Pascal.
%
%   INPUTS:
%   insig : [Nx1] column vector, unit RMS
%   Lp : RMS sound pressure level (dB SPL)
%
%   OUTPUT:
%   outsig : [Nx1] column vector, sound pressure (Pa)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pref = 20e-6; % reference sound pressure for air (Pa)

outsig = insig.*pref.*10.^(Lp/20);

end % end of <il_set_level> subfunction

%% function - error measure with respect to a confidence interval

function err = il_confidence_interval_error(values, ci_lower, ci_upper)
%function err = il_confidence_interval_error(values, ci_lower, ci_upper)
%
% Root-mean-square of the distance between <values> and the confidence
% interval bounded by <ci_lower> and <ci_upper>, as defined in Annex B.2 of
% ECMA-418-2:2025. A value lying inside the interval contributes zero.
%
%   INPUTS:
%   values : matrix of values to be assessed
%   ci_lower : matrix with the lower bounds of the confidence interval
%   ci_upper : matrix with the upper bounds of the confidence interval
%
%   OUTPUT:
%   err : scalar, root-mean-square error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance = max(max(ci_lower - values, values - ci_upper), 0);

err = sqrt(mean(distance(:).^2));

end % end of <il_confidence_interval_error> subfunction
