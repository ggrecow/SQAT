% Script validation_reference_tone_1tu
%
% Verification of the tonality model from Aures against its definitional
% anchor: a pure tone with center frequency 1 kHz and sound pressure level
% 60 dB SPL yields a tonality of 1 tonality unit (t.u.) [1]. The reference
% signal RefSignal_Tonality_Aures1985.wav shipped with SQAT contains
% exactly this signal.
%
% Additionally, the tone bandwidth estimation procedure implemented in the
% model is pinned by a tight regression check on the tonal weighting
% w_tonal: the estimated bandwidth of each tone feeds the tonal weighting
% directly, so any change in the bandwidth semantics (e.g. the band-edge
% index mapping discussed in issue #4) moves w_tonal. The semantics under
% test: the band edges are the spectral bins adjacent to the tone peak that
% first fall below the half-power target hafmax = 0.707*SPL(tone).
%
%  HOW TO RUN THIS CODE: this code requires the <Tonality_Aures1985>
%  function implemented in SQAT and the reference signal
%  RefSignal_Tonality_Aures1985.wav shipped with SQAT (sound_files/reference_signals).
%  Run startup_SQAT first. Apart from that, no additional steps are
%  necessary to run this code.
%
% Author: Sergio Aguirre, Saarbruecken 31.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% load the reference signal: pure tone, 1 kHz, 60 dB SPL (1 t.u. anchor)

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

[RefSignal, fs] = audioread([dir_ref_sounds 'RefSignal_Tonality_Aures1985.wav']);

%% Compute the Aures tonality with SQAT

OUT = Tonality_Aures1985(RefSignal, fs, 0, 0, 0);

%% Reference values

Kmean_ref_definitional = 1.0;        % definitional anchor from Ref. [1] (t.u.)
Kmean_ref_regression = 0.99969754;   % value pinned by this verification (t.u.)
wtonal_ref_regression = 0.66551671;  % value pinned by this verification (-)

tol_definitional = 0.01; % t.u. (spectral windowing and leakage of the reference recording)
tol_regression = 1e-5;   % t.u. / relative (-) (regression sensitivity)

%% Checks

check_pass = true;

% Check 1: definitional anchor, 1 kHz tone at 60 dB SPL = 1 t.u.
check1 = abs(OUT.Kmean - Kmean_ref_definitional) <= tol_definitional;
if ~check1
    fprintf('CHECK 1 FAILED: Kmean = %.6f t.u., expected %.2f t.u. within +-%.2f\n', OUT.Kmean, Kmean_ref_definitional, tol_definitional);
else
    fprintf('CHECK 1 PASSED: Kmean = %.6f t.u. within +-%.2f of the 1 t.u. anchor (1 kHz, 60 dB SPL pure tone)\n', OUT.Kmean, tol_definitional);
end
check_pass = check_pass && check1;

% Check 2: regression pin on the time-averaged tonality
check2 = abs(OUT.Kmean - Kmean_ref_regression) <= tol_regression;
if ~check2
    fprintf('CHECK 2 FAILED: Kmean = %.8f t.u., pinned value %.8f t.u. (the tone bandwidth path may have changed)\n', OUT.Kmean, Kmean_ref_regression);
else
    fprintf('CHECK 2 PASSED: Kmean = %.8f t.u. matches the pinned value\n', OUT.Kmean);
end
check_pass = check_pass && check2;

% Check 3: regression pin on the tonal weighting (tone bandwidth path)
wtonal_sqat = max(OUT.TonalWeighting);
check3 = abs(wtonal_sqat - wtonal_ref_regression) <= tol_regression;
if ~check3
    fprintf('CHECK 3 FAILED: max(w_tonal) = %.8f, pinned value %.8f (the tone bandwidth path may have changed)\n', wtonal_sqat, wtonal_ref_regression);
else
    fprintf('CHECK 3 PASSED: max(w_tonal) = %.8f matches the pinned value\n', wtonal_sqat);
end
check_pass = check_pass && check3;

if check_pass
    fprintf('\nAll checks PASSED: the Aures tonality of the reference signal conforms to its 1 t.u. anchor and to the pinned bandwidth-path values\n');
else
    fprintf('\nAt least one check FAILED\n');
end

%% plot results

figure('name','Aures tonality - verification with the 1 t.u. reference tone',...
    'units','normalized','outerposition',[0 0 1 1]);

plot(OUT.time, OUT.InstantaneousTonality, 'b-', 'Linewidth', 1.2); hold on;
plot([OUT.time(1) OUT.time(end)], OUT.Kmean*[1 1], 'k--', 'Linewidth', 0.5);
plot([OUT.time(1) OUT.time(end)], [1 1], 'r:', 'Linewidth', 1.2);
legend({'Instantaneous tonality','$K_{\mathrm{mean}}$','1 t.u. reference'},'Location','NE','Interpreter','Latex'); %legend boxoff;
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('Aures tonality, $K$ (t.u.)','Interpreter','Latex'); grid on;
axis([OUT.time(1) OUT.time(end) 0 1.3]);
title('Aures tonality of the reference signal (1 kHz, 60 dB SPL pure tone)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1

    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'tonality_validation_reference_tone_1tu';
    figname_out = [figures_dir figname_short];

    saveas(gcf,figname_out, 'png');

    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% References
% [1] Aures, W. (1985). Berechnungsverfahren fuer den sensorischen Wohlklang
%     beliebiger Schallsignale. Acta Acustica united with Acustica, 59(2), 130-141.
% [2] Terhardt, E., Stoll, G., & Seewann, M. (1982). Algorithm for extraction
%     of pitch and pitch salience from complex tonal signals. Journal of the
%     Acoustical Society of America, 71(3), 679-688.
