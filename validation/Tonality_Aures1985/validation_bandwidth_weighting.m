% Script validation_bandwidth_weighting
%
% - This routine compares the weighting w1 of the tonality model of Aures [1],
%   eq. (7), with the data it was fitted to, Fig. 6 of [1]: the tonality of
%   bandpass noise relative to that of a sine tone of the same frequency, as a
%   function of the bandwidth of the noise expressed in Bark. The reference
%   curve is analytic, so nothing is read off a figure to draw the verdict; the
%   points of Fig. 6 are shown as read off the figure, with an uncertainty of
%   about 0.01 Bark and 0.02 in the ordinate, to place the model beside the
%   measurement.
%
% - Test signals, following [1]: ideal bandpass noise 30 Hz wide at ten centre
%   frequencies from 150 Hz to 4.5 kHz, which spans a bandwidth from 0.29 down
%   to 0.04 Bark because the critical bandwidth grows with frequency, and
%   bandpass noise 1 kHz wide at 4.2 kHz, 1.37 Bark. Every signal is brought to
%   14 sone as in [1], and the reference for each is a sine tone at the centre
%   frequency at the same loudness. The spectra are random, so three seeds are
%   run and the spread is shown.
%
% - Of the two points of [1] for 1 kHz bands, at 0.69 and 1.44 Bark, the first
%   needs a centre frequency near 6.5 kHz, above the 5 kHz upper limit of the
%   implementation, so only the second has a counterpart here.
%
% - Tonality computed using
%   OUT = Tonality_Aures1985(insig,fs,LoudnessField,time_skip,show)
%   type <help Tonality_Aures1985> for more info
%
% - Reference:
%   [1] Aures, W. (1985). Berechnungsverfahren fuer den sensorischen Wohlklang
%       beliebiger Schallsignale. Acustica 59(2), 130-141. Fig. 6 and eq. (7).
%
% Author: Sergio Aguirre, September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; %#ok<CLALL>

%% save settings
save_figs = 0;
dir_out = [fileparts(mfilename('fullpath')) filesep];

fs  = 48000;
N   = 2*fs;               % two seconds per signal
N_target = 14;            % loudness of the test signals, sone, as in [1]
seeds = 1:3;              % noise realisations per stimulus

fc_30  = [4500 3700 2500 1900 1300 1000 700 400 250 150];  % 30 Hz bands
fc_1k  = 4200;                                             % 1 kHz band
bark = @(f) 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2);  % Zwicker
w1   = @(dz) 0.13./(dz + 0.13);                            % eq. (7) of [1]

% Fig. 6 of [1] read off the figure: bandwidth in Bark, relative tonality
fig6_30 = [0.025 0.83; 0.052 0.58; 0.108 0.605; 0.186 0.465; 0.273 0.30];
fig6_1k = [0.69 0.16; 1.44 0.065];

f  = (0:N-1).'*fs/N; f(f > fs/2) = f(f > fs/2) - fs;
fa = abs(f);
t  = (0:N-1).'/fs;

%% the 30 Hz bands
dz_30  = bark(fc_30 + 15) - bark(fc_30 - 15);
rel_30 = zeros(numel(seeds), numel(fc_30));
for i = 1:numel(fc_30)
    K_tone = il_tonality_at_loudness(sin(2*pi*fc_30(i)*t), fs, N_target);
    for s = 1:numel(seeds)
        x = il_bandpass_noise(fc_30(i), 30, fa, N, seeds(s));
        rel_30(s,i) = il_tonality_at_loudness(x, fs, N_target)/K_tone;
    end
end

%% the 1 kHz band
dz_1k  = bark(fc_1k + 500) - bark(fc_1k - 500);
K_tone = il_tonality_at_loudness(sin(2*pi*fc_1k*t), fs, N_target);
rel_1k = zeros(numel(seeds),1);
for s = 1:numel(seeds)
    x = il_bandpass_noise(fc_1k, 1000, fa, N, seeds(s));
    rel_1k(s) = il_tonality_at_loudness(x, fs, N_target)/K_tone;
end

%% report
mean_30 = mean(rel_30, 1);
fprintf('Bandpass noise at %d sone, tonality relative to the sine tone of the same frequency\n\n', N_target);
fprintf('%8s %8s %10s %10s %10s %10s\n', 'fc (Hz)', 'BW (Hz)', 'dz (Bark)', 'SQAT mean', 'SQAT span', 'eq. (7)');
for i = 1:numel(fc_30)
    fprintf('%8d %8d %10.3f %10.3f %10.3f %10.3f\n', fc_30(i), 30, dz_30(i), ...
            mean_30(i), max(rel_30(:,i)) - min(rel_30(:,i)), w1(dz_30(i)));
end
fprintf('%8d %8d %10.3f %10.3f %10.3f %10.3f\n', fc_1k, 1000, dz_1k, ...
        mean(rel_1k), max(rel_1k) - min(rel_1k), w1(dz_1k));

%% verdict
tol_mono = 0.02;          % allowance on the monotonic decrease, relative units
[~, order] = sort(dz_30);
if all(diff(mean_30(order)) <= tol_mono)
    fprintf('\nCHECK 1 PASSED: the relative tonality falls as the band widens in Bark\n');
else
    fprintf('\nCHECK 1 FAILED: the relative tonality rises as the band widens, %s\n', ...
            mat2str(round(mean_30(order),3)));
end

tol_w1 = 0.10;            % allowance on the distance to eq. (7), relative units
dev = mean_30 - w1(dz_30);
[~, iworst] = max(abs(dev));
if max(abs(dev)) <= tol_w1
    fprintf('CHECK 2 PASSED: over the 30 Hz bands the model stays within %.2f of eq. (7)\n', tol_w1);
else
    fprintf('CHECK 2 FAILED: over the 30 Hz bands the model departs from eq. (7) by more than %.2f\n', tol_w1);
end
fprintf('         largest distance %+.3f at %.3f Bark, root mean square %.3f\n', ...
        dev(iworst), dz_30(iworst), sqrt(mean(dev.^2)));

fprintf('CHECK 3 REPORT: 1 kHz band at %.2f Bark, model %.3f, eq. (7) %.3f, Fig. 6 reads %.3f at %.2f Bark\n', ...
        dz_1k, mean(rel_1k), w1(dz_1k), fig6_1k(2,2), fig6_1k(2,1));
fprintf('                the extraction step of [1] counts a component wider than a critical band as noise,\n');
fprintf('                so the model reads zero there while the measurement keeps a little tonality\n');

%% plot
figure('color','w');
dzc = linspace(0, 1.6, 400);
plot(dzc, w1(dzc), 'k--'); hold on;
plot(fig6_30(:,1), fig6_30(:,2), 'd', 'MarkerSize', 8, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
plot(fig6_1k(:,1), fig6_1k(:,2), 'v', 'MarkerSize', 8, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
for i = 1:numel(fc_30)
    plot(dz_30(i)*[1 1], [min(rel_30(:,i)) max(rel_30(:,i))], 'k-', 'HandleVisibility', 'off');
end
plot(dz_1k*[1 1], [min(rel_1k) max(rel_1k)], 'k-', 'HandleVisibility', 'off');
plot(dz_30, mean_30, 'ko', 'MarkerSize', 8);
plot(dz_1k, mean(rel_1k), 'kv', 'MarkerSize', 8);
xlim([0 1.6]); ylim([0 1.1]); grid off;
xlabel('Bandwidth of the noise, $\Delta z$ (Bark)', 'Interpreter', 'Latex');
ylabel('Tonality relative to the sine tone', 'Interpreter', 'Latex');
legend({'Aures, eq. (7)', 'Aures, Fig. 6, 30 Hz bands', 'Aures, Fig. 6, 1 kHz bands', ...
        'SQAT, 30 Hz bands', 'SQAT, 1 kHz band'}, 'Location', 'NorthEast', 'Interpreter', 'Latex');
legend boxoff;

if save_figs==1
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = 'tonality_validation_bandwidth_weighting_Aures_fig6';
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
