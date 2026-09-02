% Script run_verification_FS_anchor_regression
%
% Self-contained verification gate for the fluctuation strength model
% from Osses et al. [1] as implemented in SQAT
% (see <FluctuationStrength_Osses2016>). Four layers of checks:
%
% 1) Definitional anchor: a 1 kHz tone at 60 dB SPL, 100% amplitude
%    modulated at fmod = 4 Hz, yields a fluctuation strength of 1 vacil
%    [1,2]. The check applies the fluctuation strength JND of 10% [2] as
%    tolerance.
%
% 2) Regression pins: the fluctuation strength of the six AM tones of the
%    1_AM_tones_fmod verification (fc = 1 kHz, md = 1, 70 dB SPL,
%    fmod = 1 to 32 Hz), of the anchor, and of one very loud signal
%    (125 dB SPL, fmod = 4 Hz) is pinned to the values produced by the
%    current implementation within 1e-4 vacil. The very loud signal
%    exercises the corrected masked assignment of the vectorised Terhardt
%    filterbank: before that correction the computation crashed above a
%    component level of about 121 dB SPL.
%
% 3) Conformance statistics: the six-tone grid is compared against the
%    reference values of Ref. [1] (transcribed in the script of
%    1_AM_tones_fmod). The achieved root-mean-square error is pinned, so
%    the agreement with the reference can only improve. The achieved
%    state matches the published figure of 1_AM_tones_fmod, including the
%    overshoot at fmod = 2 Hz.
%
% 4) Warning behaviour: the very loud signal drives Terhardt's upper slope
%    into the regime where the model clamps it to zero (component level
%    above 120 + 1150/f dB, about 121 dB at 1 kHz), and the implementation
%    must raise SQAT:FluctuationStrength:TerhardtSlopeClamped there, since
%    the value returned is an extrapolation outside the range over which
%    the metric was validated. The same warning must stay silent at 60 and
%    70 dB SPL.
%
%  HOW TO RUN THIS CODE: this code requires the
%  <FluctuationStrength_Osses2016> function implemented in SQAT. All test
%  signals are generated inside this script, so no download of sound
%  files is necessary. Run startup_SQAT first.
%
% Author: Sergio Aguirre, Saarbruecken 31.08.2026
% Modified: Sergio Aguirre, 02.09.2026 (check 4, warning behaviour of the
% Terhardt filterbank at very high levels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% signal parameters

fs  = 44100;   % sampling frequency (Hz)
dur = 4;       % signal duration (s)
fc  = 1000;    % carrier frequency (Hz)
md  = 1;       % modulation depth
p0  = 2e-5;    % reference pressure (Pa)

fmod = [1 2 4 8 16 32];  % modulation frequencies of the grid (Hz)
spl_grid = 70;           % overall level of the grid signals (dB SPL)
spl_anchor = 60;         % overall level of the anchor (dB SPL)
spl_high = 125;          % overall level of the very loud signal (dB SPL)

%% reference and pinned values

% reference values of Ref. [1], as transcribed in 1_AM_tones_fmod
FS_ref = [0.39 0.84 1.25 1.30 0.36 0.06]; % vacil

% pinned values of the current implementation (vacil)
FS_pinned = [0.377017 1.125175 1.339299 1.306099 0.323247 0.004718];
FS_anchor_pinned = 1.019484;
FS_high_pinned   = 2.568469;

% achieved RMSE against the reference values (vacil)
RMSE_pinned = 0.1251;

tol_regression = 1e-4;   % vacil, absolute (regression sensitivity)
tol_rmse_margin = 5e-4;  % vacil, margin on the pinned RMSE (platform noise)
tol_jnd = 0.10;          % fluctuation strength JND [2], for the anchor

%% compute the fluctuation strength of every signal (method = 0, stationary)

t = (0:round(dur*fs)-1)/fs;
FSg = zeros(1,numel(fmod));

warn_id_clamp = 'SQAT:FluctuationStrength:TerhardtSlopeClamped';

lastwarn(''); % check 4 inspects the last warning raised at moderate levels
for k = 1:numel(fmod)
    insig = (1 + md*sin(2*pi*fmod(k)*t)).*sin(2*pi*fc*t);
    insig = (p0*10^(spl_grid/20))/rms(insig)*insig;  % overall level = spl_grid
    OUT = FluctuationStrength_Osses2016(insig', fs, 0, 0, 0);
    FSg(k) = OUT.FSmean;
end

insig = (1 + md*sin(2*pi*4*t)).*sin(2*pi*fc*t);
insig_anchor = (p0*10^(spl_anchor/20))/rms(insig)*insig;
OUT = FluctuationStrength_Osses2016(insig_anchor', fs, 0, 0, 0);
FS_anchor = OUT.FSmean;
[~, warn_id_moderate] = lastwarn;

insig_high = (p0*10^(spl_high/20))/rms(insig)*insig;
lastwarn('');
OUT = FluctuationStrength_Osses2016(insig_high', fs, 0, 0, 0); % raises the clamp warning once (method = 0 uses a single frame)
FS_high = OUT.FSmean;
[~, warn_id_high] = lastwarn;

%% Checks

check_pass = true;

% Check 1: definitional anchor, 1 kHz, 60 dB SPL, fmod 4 Hz, md 1 -> 1 vacil
check1 = abs(FS_anchor - 1) <= tol_jnd;
if ~check1
    fprintf('CHECK 1 FAILED: anchor FS = %.6f vacil, expected 1 vacil within the %.0f%% JND\n', FS_anchor, 100*tol_jnd);
else
    fprintf('CHECK 1 PASSED: anchor FS = %.6f vacil within the %.0f%% JND of 1 vacil (1 kHz, 60 dB SPL, fmod 4 Hz, md 1)\n', FS_anchor, 100*tol_jnd);
end
check_pass = check_pass && check1;

% Check 2: regression pins (grid, anchor, very loud signal)
allv = [FSg FS_anchor FS_high];
allp = [FS_pinned FS_anchor_pinned FS_high_pinned];
dev = max(abs(allv - allp));
check2 = dev <= tol_regression;
if ~check2
    fprintf('CHECK 2 FAILED: largest deviation from the pinned values is %.2e vacil (tolerance %.0e). The implementation changed.\n', dev, tol_regression);
else
    fprintf('CHECK 2 PASSED: all 8 pinned values match within %.0e vacil (largest deviation %.2e); the 125 dB SPL signal completes and returns %.4f vacil\n', tol_regression, dev, FS_high);
end
check_pass = check_pass && check2;

% Check 3: achieved agreement with the reference values
rmse = sqrt(mean((FSg - FS_ref).^2));
reldev = abs(FSg - FS_ref)./FS_ref;
check3 = rmse <= RMSE_pinned + tol_rmse_margin;
fprintf('%14s', 'fmod (Hz)'); fprintf('%9d', fmod); fprintf('\n');
fprintf('%14s', 'FS SQAT'); fprintf('%9.4f', FSg); fprintf('\n');
fprintf('%14s', 'FS reference'); fprintf('%9.4f', FS_ref); fprintf('\n');
fprintf('%14s', 'rel dev (%)'); fprintf('%9.1f', 100*reldev); fprintf('\n');
if ~check3
    fprintf('CHECK 3 FAILED: RMSE %.4f vacil exceeds the pinned %.4f vacil\n', rmse, RMSE_pinned);
else
    fprintf('CHECK 3 PASSED: RMSE %.4f vacil is at the pinned level (%.4f vacil)\n', rmse, RMSE_pinned);
end
check_pass = check_pass && check3;

% Check 4: warning behaviour of the Terhardt filterbank at very high levels
check4a = ~strcmp(warn_id_moderate, warn_id_clamp);
check4b = strcmp(warn_id_high, warn_id_clamp);
if ~check4a
    fprintf('CHECK 4 FAILED: %s was raised at 60 or 70 dB SPL, inside the range over which the metric was validated\n', warn_id_clamp);
elseif ~check4b
    fprintf('CHECK 4 FAILED: the 125 dB SPL signal did not raise %s (last warning identifier: "%s")\n', warn_id_clamp, warn_id_high);
else
    fprintf('CHECK 4 PASSED: no clamp warning at 60 and 70 dB SPL; the 125 dB SPL signal raised %s\n', warn_id_clamp);
end
check_pass = check_pass && check4a && check4b;

if check_pass
    fprintf('\nAll checks PASSED: the anchor holds, the implementation matches the pinned state, the agreement with Ref. [1] is at the pinned level, and the very loud signal is flagged as an extrapolation\n');
else
    fprintf('\nAt least one check FAILED\n');
end

%% plot results

figure('name','Fluctuation strength of AM tones - anchor and regression state',...
    'units','normalized','outerposition',[0 0 0.7 1]);

subplot(1,2,1);
errorbar(fmod, FS_ref, 0.1*FS_ref, 'k-*', 'MarkerSize', 6, 'Linewidth', 0.5); hold on;
semilogx(fmod, FSg, 'ro:', 'MarkerSize', 8, 'Linewidth', 1);
set(gca,'XScale','log');
legend({'Reference $\pm\:10\:\%\:(\mathrm{JND})$','SQAT'},'Location','NorthWest','Interpreter','Latex'); legend boxoff;
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');
ylabel('Fluctuation strength, $\mathrm{FS}$ (vacil)','Interpreter','Latex'); grid on;
title('AM tones, $f_{\mathrm{c}}=1$ kHz, $m_{\mathrm{d}}=1$, 70 dB SPL','Interpreter','Latex');

subplot(1,2,2);
bar(1:numel(fmod), 100*reldev); hold on;
plot([0.4 numel(fmod)+0.6], [10 10], 'k--');
set(gca,'XTickLabel',arrayfun(@num2str,fmod,'UniformOutput',false));
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');
ylabel('Relative deviation from the reference (\%)','Interpreter','Latex'); grid on;
title('Achieved agreement (dashed line: 10\% JND)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1

    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'FS_verification_anchor_regression';
    figname_out = [figures_dir figname_short];

    saveas(gcf,figname_out, 'png');

    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fail loudly when called in batch mode, so this script can gate changes
if ~check_pass
    error('run_verification_FS_anchor_regression: at least one check FAILED');
end

%% References
% [1] Osses, A., Garcia, R., & Kohlrausch, A. (2016). Modelling the
%     sensation of fluctuation strength. Proceedings of Meetings on
%     Acoustics, 28(1), 050005.
% [2] Fastl, H., & Zwicker, E. (2007). Psychoacoustics: facts and models,
%     Third edition. Springer-Verlag.
