% Script run_verification_roughness_anchor_regression
%
% Self-contained verification gate for the roughness model from Daniel &
% Weber [1] as implemented in SQAT (see <Roughness_Daniel1997>). Three
% layers of checks:
%
% 1) Definitional anchor: a 1 kHz tone at 60 dB SPL, 100% amplitude
%    modulated at fmod = 70 Hz, yields a roughness of 1 asper [1]. The
%    check applies the roughness JND of 17% [2] as tolerance.
%
% 2) Regression pins: the time-averaged roughness of 105 AM tones (seven
%    carrier frequencies, fifteen modulation frequencies, md = 1,
%    60 dB SPL) is pinned to the values produced by the current
%    implementation. Any change to the implementation shows up here first,
%    and an intentional model correction updates the pinned table with the
%    justification in the commit.
%
% 3) Conformance statistics: the same grid is compared against the
%    jury-test data of Fig. 3 of Ref. [1], stored in
%    1_AM_modulation_freq/reference_values. The achieved deviation per
%    carrier is pinned, so the agreement with the reference can only
%    improve. The pinned state of the current implementation: RMSE per
%    carrier between 0.027 and 0.100 asper, median relative deviation
%    13.0% where the reference exceeds 0.1 asper, and 39 of 94 points
%    beyond the 17% JND band, concentrated at the tails of the curves
%    where the reference values are small.
%
%  HOW TO RUN THIS CODE: this code requires the <Roughness_Daniel1997>
%  function implemented in SQAT and the reference values stored in
%  1_AM_modulation_freq/reference_values (part of the repository). All
%  test signals are generated inside this script, so no download of sound
%  files is necessary. Run startup_SQAT first.
%
% Author: Sergio Aguirre, Saarbruecken 31.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% signal parameters

fs   = 48000;        % sampling frequency (Hz), as in 1_AM_modulation_freq
dur  = 1.5;          % signal duration (s)
spl  = 60;           % overall sound pressure level of each AM tone (dB SPL)
md   = 1;            % modulation depth
tskip = 0.3;         % time skip for the statistics (s)
p0   = 2e-5;         % reference pressure (Pa)

fc   = [125 250 500 1000 2000 4000 8000];  % carrier frequencies (Hz)
fmod = 20:10:160;                          % modulation frequencies (Hz)

%% pinned values (current implementation, SQAT v1.3 lineage)

% time-averaged roughness Rmean (asper); rows: fmod, columns: fc
R_pinned = [ ...
      0.281628   0.216262   0.253985   0.207029   0.158813   0.117246   0.089455
      0.385186   0.376741   0.504634   0.434087   0.298608   0.227250   0.160049
      0.367971   0.451434   0.692970   0.678976   0.473514   0.363996   0.251984
      0.309150   0.437356   0.737680   0.863379   0.642099   0.500348   0.347603
      0.230890   0.356759   0.645222   1.008702   0.744373   0.590709   0.405828
      0.168472   0.267449   0.522054   1.000132   0.810875   0.620211   0.422330
      0.103308   0.199279   0.407744   0.931130   0.787071   0.615172   0.407797
      0.060111   0.151218   0.318873   0.796925   0.715429   0.580719   0.378098
      0.024879   0.115307   0.250895   0.671700   0.647732   0.545631   0.342571
      0.016330   0.084602   0.193943   0.564661   0.566279   0.496757   0.309553
      0.010085   0.063924   0.150598   0.475367   0.493579   0.451872   0.277381
      0.006637   0.048684   0.113781   0.398640   0.432638   0.414055   0.250295
      0.005621   0.037060   0.082995   0.340938   0.367798   0.368020   0.224840
      0.004670   0.027315   0.060346   0.288954   0.305539   0.321112   0.205004
      0.006172   0.021414   0.044317   0.244001   0.248104   0.278311   0.186557
    ];

% achieved RMSE per carrier against the Fig. 3 reference data (asper)
RMSE_pinned = [0.0435 0.0476 0.0998 0.0721 0.0415 0.0613 0.0268];

tol_regression = 1e-4;   % asper, absolute (regression sensitivity)
tol_rmse_margin = 5e-4;  % asper, margin on the pinned RMSE (platform noise)
tol_jnd = 0.17;          % roughness JND [2], for the definitional anchor

%% compute the roughness grid

R = zeros(numel(fmod), numel(fc));
t = (0:round(dur*fs)-1)/fs;

for j = 1:numel(fc)
    for k = 1:numel(fmod)
        insig = (1 + md*sin(2*pi*fmod(k)*t)) .* sin(2*pi*fc(j)*t);
        insig = (p0*10^(spl/20))/rms(insig) * insig;  % overall level = spl
        OUT = Roughness_Daniel1997(insig', fs, tskip, 0);
        R(k,j) = OUT.Rmean;
    end
end

%% load the Fig. 3 reference data (in the repository, per-carrier lengths differ)

dir_ref_values = get_dir_reference_values('Roughness_Daniel1997','1_AM_modulation_freq');
refnames = {'fmod_125hz','fmod_250hz','fmod_500hz','fmod_1khz','fmod_2khz','fmod_4khz','fmod_8khz'};
Rref = cell(1,numel(fc));
for j = 1:numel(fc)
    d = load([dir_ref_values refnames{j} '.mat']);
    v = struct2cell(d);
    Rref{j} = v{1}(:);  % fmod axis: 20:10:(10+10*numel)
end

%% Checks

check_pass = true;

% Check 1: definitional anchor, 1 kHz, 60 dB SPL, fmod 70 Hz, md 1 -> 1 asper
R_anchor = R(fmod==70, fc==1000);
check1 = abs(R_anchor - 1) <= tol_jnd;
if ~check1
    fprintf('CHECK 1 FAILED: anchor R = %.6f asper, expected 1 asper within the %.0f%% JND\n', R_anchor, 100*tol_jnd);
else
    fprintf('CHECK 1 PASSED: anchor R = %.6f asper within the %.0f%% JND of 1 asper (1 kHz, 60 dB SPL, fmod 70 Hz, md 1)\n', R_anchor, 100*tol_jnd);
end
check_pass = check_pass && check1;

% Check 2: regression pins over the whole grid
dev = max(abs(R(:) - R_pinned(:)));
check2 = dev <= tol_regression;
if ~check2
    fprintf('CHECK 2 FAILED: largest deviation from the pinned grid is %.2e asper (tolerance %.0e). The implementation changed.\n', dev, tol_regression);
else
    fprintf('CHECK 2 PASSED: all 105 grid values match the pinned table within %.0e asper (largest deviation %.2e)\n', tol_regression, dev);
end
check_pass = check_pass && check2;

% Check 3: achieved agreement with the Fig. 3 reference data, per carrier
rmse = zeros(1,numel(fc));
alldev = [];
for j = 1:numel(fc)
    n = numel(Rref{j});
    rmse(j) = sqrt(mean((R(1:n,j) - Rref{j}).^2));
    p = find(Rref{j} > 0.1);
    alldev = [alldev; abs(R(p,j) - Rref{j}(p))./Rref{j}(p)]; %#ok<AGROW>
end
check3 = all(rmse <= RMSE_pinned + tol_rmse_margin);
fprintf('%14s', 'carrier (Hz)'); fprintf('%9d', fc); fprintf('\n');
fprintf('%14s', 'RMSE (asper)'); fprintf('%9.4f', rmse); fprintf('\n');
if ~check3
    fprintf('CHECK 3 FAILED: at least one carrier exceeds its pinned RMSE against the Fig. 3 reference\n');
else
    fprintf('CHECK 3 PASSED: every carrier is at or below its pinned RMSE against the Fig. 3 reference\n');
end
fprintf('  relative deviation where the reference exceeds 0.1 asper: median %.1f%%, max %.1f%%, beyond the %.0f%% JND: %d of %d points\n', ...
    100*median(alldev), 100*max(alldev), 100*tol_jnd, sum(alldev > tol_jnd), numel(alldev));
check_pass = check_pass && check3;

if check_pass
    fprintf('\nAll checks PASSED: the anchor holds, the implementation matches the pinned state, and the agreement with Ref. [1] is at the pinned level\n');
else
    fprintf('\nAt least one check FAILED\n');
end

%% plot results

figure('name','Roughness of AM tones - anchor and regression state',...
    'units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1);
co = get(gca,'ColorOrder');
leg = cell(1,numel(fc));
for j = 1:numel(fc)
    c = co(mod(j-1,7)+1,:);
    n = numel(Rref{j});
    fm_ref = 20:10:(10+10*n);
    semilogx(fm_ref, Rref{j}, 'o', 'Color', c, 'MarkerFaceColor', c, 'HandleVisibility','off'); hold on;
    semilogx(fmod, R(:,j), '-', 'Color', c, 'Linewidth', 1.2);
    leg{j} = sprintf('$f_{\\mathrm{c}}=%g$ Hz', fc(j));
end
legend(leg,'Location','NorthWest','Interpreter','Latex');
xlabel('Modulation frequency, $f_{\mathrm{mod}}$ (Hz)','Interpreter','Latex');
ylabel('Roughness, $R$ (asper)','Interpreter','Latex'); grid on;
xlim([18 200]); ylim([0 1.1]);
title('SQAT (lines) and jury-test data of Ref. [1] (markers)','Interpreter','Latex');

subplot(1,2,2);
bar(1:numel(fc), [rmse; RMSE_pinned]');
set(gca,'XTickLabel',arrayfun(@num2str,fc,'UniformOutput',false));
legend({'computed','pinned'},'Location','NorthEast','Interpreter','Latex');
xlabel('Carrier frequency, $f_{\mathrm{c}}$ (Hz)','Interpreter','Latex');
ylabel('RMSE vs. Ref. [1] (asper)','Interpreter','Latex'); grid on;
title('Achieved agreement per carrier','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1

    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'roughness_verification_anchor_regression';
    figname_out = [figures_dir figname_short];

    saveas(gcf,figname_out, 'png');

    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fail loudly when called in batch mode, so this script can gate changes
if ~check_pass
    error('run_verification_roughness_anchor_regression: at least one check FAILED');
end

%% References
% [1] Daniel, P., & Weber, R. (1997). Psychoacoustical Roughness:
%     Implementation of an Optimized Model. Acta Acustica united with
%     Acustica, 83(1), 113-123.
% [2] Fastl, H., & Zwicker, E. (2007). Psychoacoustics: facts and models,
%     Third edition. Springer-Verlag.
