% Script validation_frequency_weighting
%
% - This routine compares the weighting w2 of the tonality model of Aures [1],
%   eq. (9), with Fig. 5 of [1]: the relative tonality of sine tones, of
%   bandpass noise 30 Hz and 1 kHz wide, and of highpass noise above 2 kHz,
%   as a function of the critical band rate, all at 14 sone. Fig. 5 reads
%   1.00 for the tone at 8.6 Bark, and [1] normalises its measured tonality
%   to the 1 kHz tone at 14 sone (p. 137), so the implementation is
%   normalised by its own value for that tone. The reference curves are
%   analytic, eq. (9) as w2(f)/w2(1 kHz) for the sine tones and its product
%   with eq. (7) for the 30 Hz bands, so nothing is read off a figure to draw
%   the verdict; the points of Fig. 5 are shown as read off the figure, with
%   an uncertainty of about 0.1 Bark and 0.02 in the ordinate, to place the
%   model beside the measurement.
%
% - Test signals: sine tones from 1 to 18 Bark in steps of 1 Bark, which
%   traces eq. (9); following [1], at the four positions of Fig. 5 below the
%   5 kHz upper limit of the implementation, 5.1, 8.6, 13 and 17 Bark, a sine
%   tone, ideal bandpass noise 30 Hz wide and ideal bandpass noise 1 kHz wide;
%   and ideal highpass noise above 2 kHz, which Fig. 5 places at 24 Bark.
%   Every signal is brought to 14 sone as in [1]. The spectra of the noise
%   are random, so three seeds are run and the spread is shown.
%
% - The fifth position of Fig. 5, 21 Bark, is about 7.6 kHz, above the 5 kHz
%   upper limit of the implementation, so it has no counterpart here. The
%   critical band rate is computed from frequency with the formula of Zwicker,
%   as in validation_bandwidth_weighting.m.
%
% - Tonality computed using
%   OUT = Tonality_Aures1985(insig,fs,LoudnessField,time_skip,show)
%   type <help Tonality_Aures1985> for more info
%
% - Reference:
%   [1] Aures, W. (1985). Berechnungsverfahren fuer den sensorischen Wohlklang
%       beliebiger Schallsignale. Acustica 59(2), 130-141. Fig. 5, eq. (7)
%       and eq. (9).
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September 2026 assisted 
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
clear all; close all; clc;  

%% save settings
save_figs = 0;
dir_out = [fileparts(mfilename('fullpath')) filesep];

fs  = 48000;
N   = 2*fs;               % two seconds per signal
N_target = 14;            % loudness of the test signals, sone, as in [1]
seeds = 1:3;              % noise realisations per stimulus

bark  = @(f) 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2);  % Zwicker
hertz = @(z) fzero(@(f) bark(f) - z, [20 20000]);           % inverse of bark
w1    = @(dz) 0.13./(dz + 0.13);                            % eq. (7) of [1]
w2    = @(f) (1./sqrt(1 + 0.2*(f/700 + 700./f).^2)).^0.29;  % eq. (9) of [1]
w2rel = @(f) w2(f)/w2(1000);                                % relative to 1 kHz

% Fig. 5 of [1] read off the figure: critical band rate in Bark, relative tonality
fig5_tone = [5.1 1.00; 8.6 1.00; 13.0 0.871; 17.0 0.750; 21.0 0.569];
fig5_b30  = [5.1 0.304; 8.7 0.465; 13.1 0.440; 17.0 0.440; 21.0 0.461];
fig5_b1k  = [5.1 0.028; 8.6 0.032; 13.1 0.032; 17.0 0.044; 21.0 0.085];
fig5_hp   = [24 0.02];

f_pos  = arrayfun(hertz, fig5_tone(1:4,1).');      % positions below 5 kHz
f_tone = unique([arrayfun(hertz, 1:18) f_pos]);    % sine tones

f  = (0:N-1).'*fs/N; f(f > fs/2) = f(f > fs/2) - fs;
fa = abs(f);
t  = (0:N-1).'/fs;

%% the reference, the 1 kHz tone at 14 sone
K_1k = il_tonality_at_loudness(sin(2*pi*1000*t), fs, N_target);

%% sine tones
rel_tone = zeros(size(f_tone));
for i = 1:numel(f_tone)
    rel_tone(i) = il_tonality_at_loudness(sin(2*pi*f_tone(i)*t), fs, N_target)/K_1k;
end

%% bandpass noise 30 Hz and 1 kHz wide at the positions of Fig. 5
rel_b30 = zeros(numel(seeds), numel(f_pos));
rel_b1k = zeros(numel(seeds), numel(f_pos));
for i = 1:numel(f_pos)
    for s = 1:numel(seeds)
        x = il_bandpass_noise(f_pos(i), 30, fa, N, seeds(s));
        rel_b30(s,i) = il_tonality_at_loudness(x, fs, N_target)/K_1k;
        x = il_bandpass_noise(f_pos(i), 1000, fa, N, seeds(s));
        rel_b1k(s,i) = il_tonality_at_loudness(x, fs, N_target)/K_1k;
    end
end

%% highpass noise above 2 kHz, a band from 2 kHz to fs/2
rel_hp = zeros(numel(seeds), 1);
for s = 1:numel(seeds)
    x = il_bandpass_noise((2000 + fs/2)/2, fs/2 - 2000, fa, N, seeds(s));
    rel_hp(s) = il_tonality_at_loudness(x, fs, N_target)/K_1k;
end

%% report
fprintf('Tonality at %d sone relative to the 1 kHz sine tone at %d sone, which reads %.4f t.u.\n\n', ...
        N_target, N_target, K_1k);
fprintf('Sine tones\n');
fprintf('%8s %9s %10s %10s\n', 'f (Hz)', 'z (Bark)', 'SQAT', 'eq. (9)');
for i = 1:numel(f_tone)
    fprintf('%8.0f %9.2f %10.3f %10.3f\n', f_tone(i), bark(f_tone(i)), rel_tone(i), w2rel(f_tone(i)));
end

dz_30   = bark(f_pos + 15) - bark(f_pos - 15);
ref_30  = w1(dz_30).*w2rel(f_pos);
mean_30 = mean(rel_b30, 1);
fprintf('\nBandpass noise 30 Hz wide\n');
fprintf('%8s %9s %10s %10s %10s %13s %8s\n', 'fc (Hz)', 'z (Bark)', 'dz (Bark)', 'SQAT mean', 'SQAT span', 'eq. (7)*(9)', 'Fig. 5');
for i = 1:numel(f_pos)
    fprintf('%8.0f %9.2f %10.3f %10.3f %10.3f %13.3f %8.3f\n', f_pos(i), bark(f_pos(i)), dz_30(i), ...
            mean_30(i), max(rel_b30(:,i)) - min(rel_b30(:,i)), ref_30(i), fig5_b30(i,2));
end

mean_1k = mean(rel_b1k, 1);
fprintf('\nBandpass noise 1 kHz wide\n');
fprintf('%8s %9s %10s %10s %8s\n', 'fc (Hz)', 'z (Bark)', 'SQAT mean', 'SQAT span', 'Fig. 5');
for i = 1:numel(f_pos)
    fprintf('%8.0f %9.2f %10.3f %10.3f %8.3f\n', f_pos(i), bark(f_pos(i)), ...
            mean_1k(i), max(rel_b1k(:,i)) - min(rel_b1k(:,i)), fig5_b1k(i,2));
end
fprintf('\nHighpass noise above 2 kHz: SQAT mean %.3f, span %.3f, Fig. 5 %.2f at 24 Bark\n', ...
        mean(rel_hp), max(rel_hp) - min(rel_hp), fig5_hp(2));

%% verdict
tol_w2 = 0.02;            % allowance on the distance to eq. (9), relative units
dev_tone = rel_tone - w2rel(f_tone);
[~, iworst] = max(abs(dev_tone));
if max(abs(dev_tone)) <= tol_w2
    fprintf('\nCHECK 1 PASSED: over the sine tones the model stays within %.2f of eq. (9)\n', tol_w2);
else
    fprintf('\nCHECK 1 FAILED: over the sine tones the model departs from eq. (9) by more than %.2f\n', tol_w2);
end
fprintf('         largest distance %+.3f at %.0f Hz, root mean square %.3f\n', ...
        dev_tone(iworst), f_tone(iworst), sqrt(mean(dev_tone.^2)));

tol_w1w2 = 0.10;          % allowance on the distance to eq. (7) times eq. (9), relative units
dev_30 = mean_30 - ref_30;
[~, iworst] = max(abs(dev_30));
if max(abs(dev_30)) <= tol_w1w2
    fprintf('CHECK 2 PASSED: over the 30 Hz bands the model stays within %.2f of eq. (7) times eq. (9)\n', tol_w1w2);
else
    fprintf('CHECK 2 FAILED: over the 30 Hz bands the model departs from eq. (7) times eq. (9) by more than %.2f\n', tol_w1w2);
end
fprintf('         largest distance %+.3f at %.1f Bark, root mean square %.3f\n', ...
        dev_30(iworst), bark(f_pos(iworst)), sqrt(mean(dev_30.^2)));

[~, ipos] = ismember(f_pos, f_tone);
fprintf('CHECK 3 REPORT: against the points of Fig. 5 at 5.1, 8.6, 13 and 17 Bark\n');
fprintf('                sine tones,     model %s, Fig. 5 %s\n', ...
        mat2str(round(rel_tone(ipos),3)), mat2str(fig5_tone(1:4,2).'));
fprintf('                30 Hz bands,    model %s, Fig. 5 %s\n', ...
        mat2str(round(mean_30,3)), mat2str(fig5_b30(1:4,2).'));
fprintf('                1 kHz bands,    model %s, Fig. 5 %s\n', ...
        mat2str(round(mean_1k,3)), mat2str(fig5_b1k(1:4,2).'));

%% plot
figure('color','w');
zc = linspace(1, 21, 200);
fc = arrayfun(hertz, zc);
plot(zc, w2rel(fc), 'k--'); hold on;
plot(zc, w1(bark(fc + 15) - bark(fc - 15)).*w2rel(fc), 'k:');
grey = [0.5 0.5 0.5];
plot(fig5_tone(:,1), fig5_tone(:,2), 'o', 'MarkerSize', 8, 'Color', grey, 'MarkerFaceColor', grey);
plot(fig5_b30(:,1),  fig5_b30(:,2),  'd', 'MarkerSize', 8, 'Color', grey, 'MarkerFaceColor', grey);
plot(fig5_b1k(:,1),  fig5_b1k(:,2),  'v', 'MarkerSize', 8, 'Color', grey, 'MarkerFaceColor', grey);
plot(fig5_hp(1),     fig5_hp(2),     '^', 'MarkerSize', 8, 'Color', grey, 'MarkerFaceColor', grey);
for i = 1:numel(f_pos)
    plot(bark(f_pos(i))*[1 1], [min(rel_b30(:,i)) max(rel_b30(:,i))], 'k-', 'HandleVisibility', 'off');
    plot(bark(f_pos(i))*[1 1], [min(rel_b1k(:,i)) max(rel_b1k(:,i))], 'k-', 'HandleVisibility', 'off');
end
plot(24*[1 1], [min(rel_hp) max(rel_hp)], 'k-', 'HandleVisibility', 'off');
plot(bark(f_tone), rel_tone, 'ko', 'MarkerSize', 6);
plot(bark(f_pos), mean_30, 'kd', 'MarkerSize', 8);
plot(bark(f_pos), mean_1k, 'kv', 'MarkerSize', 8);
plot(24, mean(rel_hp), 'k^', 'MarkerSize', 8);
xlim([0 25]); ylim([0 1.1]); grid off;
xlabel('Critical band rate, $z$ (Bark)', 'Interpreter', 'Latex');
ylabel('Tonality relative to the 1 kHz sine tone', 'Interpreter', 'Latex');
legend({'Aures, eq. (9), sine tones', 'Aures, eq. (7) times eq. (9), 30 Hz bands', ...
        'Aures, Fig. 5, sine tones', 'Aures, Fig. 5, 30 Hz bands', 'Aures, Fig. 5, 1 kHz bands', ...
        'Aures, Fig. 5, highpass 2 kHz', 'SQAT, sine tones', 'SQAT, 30 Hz bands', ...
        'SQAT, 1 kHz bands', 'SQAT, highpass 2 kHz'}, ...
       'Location', 'EastOutside', 'Interpreter', 'Latex');
legend boxoff;

if save_figs==1
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = 'tonality_validation_frequency_weighting_Aures_fig5';
    figname_out = [figures_dir figname_short];
    saveas(gcf, figname_out, 'png');
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end

%% local functions
function x = il_bandpass_noise(fc, bw, fa, N, seed)
% ideal bandpass noise of bandwidth bw at fc, random phase, unit rms
rng(seed);
mag = zeros(N,1);
mag(fa >= fc-bw/2 & fa <= fc+bw/2) = 1;
ph  = exp(1j*2*pi*rand(N,1)); ph(1) = 1;
X   = mag.*ph; X(N/2+2:end) = conj(flipud(X(2:N/2)));
x   = real(ifft(X));
x   = x/rms(x);
end

function K = il_tonality_at_loudness(x, fs, N_target)
% brings x to the target loudness by bisection on the gain, then returns the
% time averaged tonality
x  = x/rms(x);
lo = 20e-6*10^(30/20); hi = 20e-6*10^(110/20);
for it = 1:14
    g = sqrt(lo*hi);
    if Loudness_ISO532_1(x*g, fs, 0, 1, 0, 0).Loudness < N_target, lo = g; else, hi = g; end
end
K = Tonality_Aures1985(x*sqrt(lo*hi), fs, 0, 0, false).Kmean;
end
