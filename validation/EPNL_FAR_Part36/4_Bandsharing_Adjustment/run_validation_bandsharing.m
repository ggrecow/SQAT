% Script run_validation_bandsharing
%
% Verification of the bandsharing adjustment applied to PNLTM
% in the EPNL code
%
%   - The Effective Perceived Noise Level (EPNL) is a metric used mainly in
%     the context of environmental aircraft noise, mainly for aircraft noise
%     certification. The maximum tone-corrected perceived noise level PNLTM
%     receives a bandsharing adjustment whenever the tone correction at the
%     time of the maximum is smaller than the average tone correction of the
%     five neighbouring spectra. The adjustment shall be computed as [1,2]:
%
%         Cavg = [ C(kM-2) + C(kM-1) + C(kM) + C(kM+1) + C(kM+2) ] / 5
%         if Cavg > C(kM):  DeltaB = Cavg - C(kM)   and   PNLTM = PNLT(kM) + DeltaB
%
%     where kM is the time index of the maximum of PNLT.
%
%   - The tone correction factors C(k) themselves are verified separately in
%     the [2_Tone_Correction_Factor] folder (against Table 3.7 of Ref. [3]).
%     Here the C(k) series is taken as given and the bandsharing arithmetic
%     applied on top of it is what is verified, against the analytic
%     formulation above (tolerance 1e-9 TPNdB).
%
%   - The test uses synthetic one-third-octave spectra designed so that the
%     bandsharing branch is triggered: the PNLT maximum occurs at a time step
%     with a smooth spectrum (small C at the peak) surrounded by time steps
%     containing strong tones (large C around the peak), so that
%     Cavg > C(kM). No external data is required.
%
%  HOW TO RUN THIS CODE: this code requires the <EPNL_FAR_Part36> function
%  and its <get_PNLT> helper available in SQAT (run startup_SQAT first).
%  Apart from that, no additional steps are necessary to run this code.
%
%  NOTE: this verification requires the bandsharing arithmetic as amended by
%  PR #56 (DeltaB = Cavg - Cmax(kM), subtraction). On versions prior to that
%  amendment the DeltaB check below fails by design.
%
% Author: Sergio Aguirre, Saarbruecken 31.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

tol = 1e-9; % tolerance for the analytic identity (TPNdB)

%% path settings

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% Synthetic one-third-octave spectra: 21 time steps x 24 bands (50 Hz to 10 kHz)

fc = [50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 ...
      2500 3150 4000 5000 6300 8000 10000];

nT = 21;         % number of time steps (0.5 s each)
kPeak = 11;      % time step of the PNLT maximum
iTone = 15;      % band index of the tones around the peak (1250 Hz)

SPL = 45 - 0.2*(1:numel(fc));          % smooth broadband base spectrum (dB per band)
M = repmat(SPL, nT, 1);                % [nT x nF] spectra matrix
M(kPeak, :) = M(kPeak, :) + 12;        % loud smooth peak at kPeak (small C at the peak)
M(kPeak, 8) = M(kPeak, 8) + 2;         % small spectral irregularity at the peak
nb = setdiff(8:14, kPeak);
M(nb, iTone) = M(nb, iTone) + 18;      % strong tones one step either side of the peak

%% Compute EPNL and decompose the bandsharing adjustment with SQAT

dt = 0.5; threshold = 10;

OUT = EPNL_FAR_Part36(M, [], 0, dt, threshold, 0);

[PNLT, PNLTM, idx, G] = get_PNLT(M, fc, OUT.PNL);

Cmax = max(G.C, [], 1); % largest tone correction factor per time step, C(k)
Cmax = Cmax(:);

win = linspace(idx-2, idx+2, 5); % the five consecutive spectra around kM
valid = win(win > 0 & win <= numel(PNLT));
Cavg = mean(Cmax(valid));        % bandsharing average per Ref. [1,2]

%% Analytic reference values from the formulation of Ref. [1,2]

DeltaB_ref = Cavg - Cmax(idx);
PNLTM_ref = PNLT(idx) + DeltaB_ref;

%% Checks

check_pass = true;

% Check 1: the synthetic case must trigger the bandsharing branch
check1 = Cavg > Cmax(idx);
if ~check1
    fprintf('CHECK 1 FAILED: the synthetic spectra do not trigger the bandsharing branch (Cavg <= C(kM))\n');
else
    fprintf('CHECK 1 PASSED: bandsharing branch triggered (Cavg = %.4f > C(kM) = %.4f)\n', Cavg, Cmax(idx));
end
check_pass = check_pass && check1;

% Check 2: PNLTM from get_PNLT equals PNLT(kM) + (Cavg - C(kM))
check2 = abs(PNLTM - PNLTM_ref) <= tol;
if ~check2
    fprintf('CHECK 2 FAILED: PNLTM = %.6f, expected PNLT(kM) + (Cavg - C(kM)) = %.6f (difference %.6f TPNdB)\n', ...
        PNLTM, PNLTM_ref, abs(PNLTM - PNLTM_ref));
else
    fprintf('CHECK 2 PASSED: PNLTM = %.6f = PNLT(kM) + DeltaB, DeltaB = Cavg - C(kM) = %.4f TPNdB\n', ...
        PNLTM, DeltaB_ref);
end
check_pass = check_pass && check2;

% Check 3: the adjusted PNLTM is the one reported by EPNL_FAR_Part36
check3 = abs(OUT.PNLTM - PNLTM_ref) <= tol;
if ~check3
    fprintf('CHECK 3 FAILED: EPNL_FAR_Part36 reports PNLTM = %.6f, expected %.6f\n', OUT.PNLTM, PNLTM_ref);
else
    fprintf('CHECK 3 PASSED: EPNL_FAR_Part36 reports the bandsharing-adjusted PNLTM = %.6f TPNdB\n', OUT.PNLTM);
end
check_pass = check_pass && check3;

% Check 4: single-spectrum inputs skip the bandsharing adjustment (guard)
[PNLT1, PNLTM1] = get_PNLT(M(1,:), fc, OUT.PNL(1));
check4 = (PNLTM1 == PNLT1(1));
if ~check4
    fprintf('CHECK 4 FAILED: single-spectrum PNLTM = %.6f differs from PNLT(1) = %.6f\n', PNLTM1, PNLT1(1));
else
    fprintf('CHECK 4 PASSED: single-spectrum input skips the adjustment (PNLTM = PNLT(1) = %.6f TPNdB)\n', PNLTM1);
end
check_pass = check_pass && check4;

if check_pass
    fprintf('\nAll checks PASSED: the bandsharing adjustment conforms to Refs. [1,2]\n');
else
    fprintf('\nAt least one check FAILED: the bandsharing adjustment does not conform to Refs. [1,2]\n');
end

%% plot results

figure('name','EPNL calculation - verification of the bandsharing adjustment',...
    'units','normalized','outerposition',[0 0 1 1]);

% PNLT curve with the bandsharing-adjusted PNLTM
subplot(2,1,1)
plot((1:nT)*dt, PNLT, 'b-', 'Linewidth', 1.2); hold on;
yl = ylim;
plot(idx*dt, PNLT(idx), 'bo', 'MarkerFaceColor', 'b');
plot(idx*dt, PNLTM_ref, 'rs', 'MarkerFaceColor', 'r');
plot([0 nT*dt], PNLTM_ref*[1 1], 'r--');
legend({'PNLT(k)','PNLT(k_M)','PNLTM (bandsharing-adjusted)'},'Location','NE','Interpreter','Latex'); %legend boxoff;
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('PNLT (TPNdB)','Interpreter','Latex'); grid on;
title('Bandsharing adjustment to PNLTM','Interpreter','Latex');

% tone correction factors with the averaging window
subplot(2,1,2)
plot((1:nT)*dt, Cmax, 'b-', 'Linewidth', 1.2); hold on;
plot(valid*dt, Cmax(valid), 'go', 'MarkerFaceColor', 'g');
plot([0 nT*dt], Cavg*[1 1], 'g--');
plot(idx*dt, Cmax(idx), 'rs', 'MarkerFaceColor', 'r');
legend({'C(k) = \mathrm{max}(C(.,k))','averaging window','C_{avg}','C(k_M)'},'Location','NE','Interpreter','Latex'); %legend boxoff;
xlabel('Time, $t$ (s)','Interpreter','Latex');
ylabel('$C$ (dB)','Interpreter','Latex'); grid on;
title('Tone correction factors and bandsharing average','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1

    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'validation_Do_bandsharing';
    figname_out = [figures_dir figname_short];

    saveas(gcf,figname_out, 'png');

    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% References
% [1] EASA, Environmental Protection Technical Specifications applicable to
%     VTOL-capable aircraft powered by tilting rotors, Doc. No. TE.CERT.00075-002,
%     section NVTOL-TILT.1105(d)(2), 12 Dec 2023.
%     https://www.easa.europa.eu/en/downloads/139024/en
% [2] Annex 16 to the Convention on International Civil Aviation,
%     Environmental Protection, Volume I - Aircraft Noise, Appendix 2, section 4.4.
% [3] International Civil Aviation Organization (2015) Doc 9501, Environmental
%     Technical Manual Volume I, Table 3.7.
