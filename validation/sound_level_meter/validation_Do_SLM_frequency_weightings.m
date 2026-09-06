% Script validation_Do_SLM_frequency_weightings
%
% - Verification of the frequency weightings A, C and Z used by the sound
%   level meter of SQAT against IEC 61672-1:2013.
%
% - Filters obtained from:
%   [b, a] = Gen_weighting_filters(fs, weight_freq)
%   type <help Gen_weighting_filters> for more info
%
% - Reference values come from the standard. Table 3 gives the design goal
%   of each frequency weighting at 34 nominal frequencies from 10 Hz to
%   20 kHz, rounded to a tenth of a decibel, together with the acceptance
%   limits for performance classes 1 and 2. Both are transcribed below.
%   The design goals themselves are recomputed by this script from the
%   analytical expressions of Annex E (normative), so the transcription is
%   checked against the expressions the table was built from, and the
%   deviations are then taken against the unrounded values.
%
% - As stated in the note to Table 3, the weightings are defined at the
%   exact frequencies f = 1000*10^((n-30)/10) Hz with n an integer from 10
%   to 43, so the filter response is evaluated there and not at the rounded
%   nominal frequencies.
%
% - Two criteria are applied at every frequency:
%
%     a) conformance (5.5.6): the deviation of the digital filter from the
%        design goal shall stay inside the acceptance limits of Table 3;
%
%     b) regression: from 10 Hz to 4 kHz the deviation shall also stay
%        below 0.05 dB. This is much tighter than the acceptance limits and
%        it is what detects a change in the pole frequencies or in the
%        normalization of the filters. It is not applied above 4 kHz, where
%        the bilinear transform used by Gen_weighting_filters warps the
%        frequency axis enough to produce a roll-off of its own. That
%        roll-off is a property of the sampling rate and is reported here.
%
% Author: Sergio Aguirre (& Claude code), 28.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Parameters

fs_list  = [44100 48000]; % sampling frequencies under test, Hz
fs_figs  = 48000;         % sampling frequency shown in the figures, Hz
f_tight  = 4000;          % upper frequency of the tight tolerance, Hz
tol      = 0.05;          % tight tolerance below f_tight, dB
tol_tab  = 0.05;          % rounding tolerance of the Table 3 transcription, dB

weight_freq = {'A','C','Z'};

%% IEC 61672-1:2013, Table 3

% nominal frequencies, Hz
f_nominal = [10 12.5 16 20 25 31.5 40 50 63 80 100 125 160 200 250 315 ...
             400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 ...
             6300 8000 10000 12500 16000 20000].';

% exact frequencies at which the weightings are defined, Hz
n_index = (10:43).';
f_exact = 1000*10.^((n_index-30)/10);

% design goals, dB, rounded to a tenth of a decibel
designGoal.A = [-70.4 -63.4 -56.7 -50.5 -44.7 -39.4 -34.6 -30.2 -26.2 ...
                -22.5 -19.1 -16.1 -13.4 -10.9 -8.6 -6.6 -4.8 -3.2 -1.9 ...
                -0.8 0.0 0.6 1.0 1.2 1.3 1.2 1.0 0.5 -0.1 -1.1 -2.5 ...
                -4.3 -6.6 -9.3].';
designGoal.C = [-14.3 -11.2 -8.5 -6.2 -4.4 -3.0 -2.0 -1.3 -0.8 -0.5 ...
                -0.3 -0.2 -0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ...
                -0.1 -0.2 -0.3 -0.5 -0.8 -1.3 -2.0 -3.0 -4.4 -6.2 ...
                -8.5 -11.2].';
designGoal.Z = zeros(size(f_exact));

% acceptance limits, dB, -Inf meaning that no lower limit is specified
class1_upper = [3.0 2.5 2.0 2.0 2.0 1.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 ...
                1.0 1.0 1.0 1.0 1.0 1.0 0.7 1.0 1.0 1.0 1.0 1.0 1.0 1.5 ...
                1.5 1.5 2.0 2.0 2.5 3.0].';
class1_lower = [-Inf -Inf -4.0 -2.0 -1.5 -1.5 -1.0 -1.0 -1.0 -1.0 -1.0 ...
                -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -0.7 -1.0 ...
                -1.0 -1.0 -1.0 -1.0 -1.0 -1.5 -2.0 -2.5 -3.0 -5.0 ...
                -16.0 -Inf].';
class2_upper = [5.0 5.0 5.0 3.0 3.0 3.0 2.0 2.0 2.0 2.0 1.5 1.5 1.5 1.5 ...
                1.5 1.5 1.5 1.5 1.5 1.5 1.0 1.5 2.0 2.0 2.5 2.5 3.0 3.5 ...
                4.5 5.0 5.0 5.0 5.0 5.0].';
class2_lower = [-Inf -Inf -Inf -3.0 -3.0 -3.0 -2.0 -2.0 -2.0 -2.0 -1.5 ...
                -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.5 -1.0 -1.5 ...
                -2.0 -2.0 -2.5 -2.5 -3.0 -3.5 -4.5 -5.0 -Inf -Inf ...
                -Inf -Inf].';

%% IEC 61672-1:2013, Annex E: analytical design goals

% E.2.3 and E.2.4: the two low frequency poles f1 and the two high frequency
% poles f4 of the C weighting solve a bi-quadratic equation set by the -3 dB
% points f_L and f_H of its power response, relative to the reference
% frequency f_r
f_r = 1000;     % Hz
f_L = 10^1.5;   % Hz
f_H = 10^3.9;   % Hz
D   = sqrt(1/2);

b_E = ( f_r^2 + (f_L^2*f_H^2)/f_r^2 - D*(f_L^2 + f_H^2) )/(1-D); % Eq. (E.4)
c_E = f_L^2 * f_H^2;                                             % Eq. (E.5)

f1 = sqrt( (-b_E - sqrt(b_E^2 - 4*c_E))/2 ); % Eq. (E.2),    20.5990 Hz
f4 = sqrt( (-b_E + sqrt(b_E^2 - 4*c_E))/2 ); % Eq. (E.3), 12194.2171 Hz

% E.3.3 and E.3.4: the A weighting adds two coupled first order high-pass
% filters to the C weighting, with poles at f2 and f3
f_Acut = 10^2.45; % Hz
f2 = ((3 - sqrt(5))/2)*f_Acut; % Eq. (E.7),  107.6526 Hz
f3 = ((3 + sqrt(5))/2)*f_Acut; % Eq. (E.8),  737.8622 Hz

% Eq. (E.1) and Eq. (E.6) without the normalization constants C_1000 and
% A_1000, which are applied below by forcing 0 dB at the reference frequency
% (C_1000 = -0.062 dB and A_1000 = -2.000 dB, see E.4.2)
C_unnorm = @(f) 20*log10( (f4^2*f.^2) ...
                          ./ ( (f.^2 + f1^2).*(f.^2 + f4^2) ) );
A_unnorm = @(f) 20*log10( (f4^2*f.^4) ...
                          ./ ( (f.^2 + f1^2).*sqrt(f.^2 + f2^2) ...
                               .*sqrt(f.^2 + f3^2).*(f.^2 + f4^2) ) );

annexE.C = C_unnorm(f_exact) - C_unnorm(f_r);
annexE.A = A_unnorm(f_exact) - A_unnorm(f_r);
annexE.Z = zeros(size(f_exact));   % Eq. (E.9)

%% Check of the transcribed design goals against Annex E

fprintf('\n%s\n', repmat('=',1,78));
fprintf('IEC 61672-1 Table 3 design goals against the expressions of Annex E\n');
fprintf('%s\n', repmat('=',1,78));
ok_table = true;
for iw = 1:length(weight_freq)
    w = weight_freq{iw};
    d = abs(annexE.(w) - designGoal.(w));
    fprintf('%s weighting: largest difference from the tabulated value %.3f dB\n', ...
            w, max(d));
    ok_table = ok_table && all(d <= tol_tab);
end
if ok_table
    fprintf('PASS: every tabulated design goal is the Annex E value rounded to 0.1 dB.\n');
else
    fprintf(2,'FAIL: a tabulated design goal is not the Annex E value rounded to 0.1 dB.\n');
end

%% Response of the filters generated by SQAT

ok_class1 = true;
ok_class2 = true;
ok_tight  = true;
verdict   = {'FAIL','PASS'};

for ifs = 1:length(fs_list)

    fs = fs_list(ifs);

    for iw = 1:length(weight_freq)
        w = weight_freq{iw};
        [b,a] = Gen_weighting_filters(fs, w);
        H = freqz(b, a, f_exact, fs);
        resp.(w) = 20*log10(abs(H));
        dev.(w)  = resp.(w) - annexE.(w);
    end

    idx_tight = f_exact <= f_tight;

    fprintf('\n%s\n', repmat('=',1,78));
    fprintf('Frequency weightings at fs = %g Hz: deviation from the design goal\n', fs);
    fprintf('%s\n', repmat('=',1,78));
    fprintf('%10s %10s %8s %8s %8s %16s %14s\n', ...
            'nominal', 'exact', 'dev A', 'dev C', 'dev Z', ...
            'class 1 limits', 'class 2 limits');
    fprintf('%10s %10s %8s %8s %8s %16s %14s\n', ...
            'Hz', 'Hz', 'dB', 'dB', 'dB', 'dB', 'dB');
    fprintf('%s\n', repmat('-',1,78));
    for i = 1:length(f_exact)
        fprintf('%10g %10.2f %8.3f %8.3f %8.3f %7.1f %8.1f %6.1f %7.1f\n', ...
                f_nominal(i), f_exact(i), dev.A(i), dev.C(i), dev.Z(i), ...
                class1_upper(i), class1_lower(i), ...
                class2_upper(i), class2_lower(i));
    end
    fprintf('%s\n', repmat('-',1,78));

    for iw = 1:length(weight_freq)
        w = weight_freq{iw};
        d = dev.(w);
        in1 = all(d <= class1_upper) && all(d >= class1_lower);
        in2 = all(d <= class2_upper) && all(d >= class2_lower);
        in_tight = all(abs(d(idx_tight)) < tol);

        [~, imax]       = max(abs(d));
        [~, imax_tight] = max(abs(d(idx_tight)));
        d_tight         = d(idx_tight);

        fprintf(['%s weighting: largest deviation %+.3f dB at %g Hz, ' ...
                 '%+.3f dB up to %g Hz. class 1 %s, class 2 %s.\n'], ...
                w, d(imax), f_nominal(imax), d_tight(imax_tight), f_tight, ...
                verdict{in1+1}, verdict{in2+1});

        ok_class1 = ok_class1 && in1;
        ok_class2 = ok_class2 && in2;
        ok_tight  = ok_tight  && in_tight;
    end

    if fs == fs_figs
        dev_figs = dev; % kept for the figures
    end
end

%% Verdict

fprintf('\n%s\n', repmat('=',1,78));
ok = ok_table && ok_class1 && ok_class2 && ok_tight;
if ok
    fprintf(['PASS: the A, C and Z weightings meet the acceptance limits of ' ...
             'classes 1 and 2\n      at every sampling frequency tested, and ' ...
             'stay within %.2f dB of the\n      design goal from 10 Hz to %g Hz.\n'], ...
            tol, f_tight);
else
    fprintf(2, 'FAIL: see the tables above.\n');
end
fprintf('%s\n\n', repmat('=',1,78));

%% Plot: weighting curves against the design goals of Table 3

title_fig1 = 'Do SLM - frequency weightings against IEC 61672-1 Table 3';
h1 = figure('Name', title_fig1);
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

f_dense = logspace(log10(10), log10(20000), 500).';

[bA,aA] = Gen_weighting_filters(fs_figs, 'A');
[bC,aC] = Gen_weighting_filters(fs_figs, 'C');
HA = freqz(bA, aA, f_dense, fs_figs);
HC = freqz(bC, aC, f_dense, fs_figs);

hA = semilogx(f_dense, 20*log10(abs(HA)), 'k-', 'Linewidth', 1); hold on;
hC = semilogx(f_dense, 20*log10(abs(HC)), 'b-', 'Linewidth', 1);
hZ = semilogx(f_dense, zeros(size(f_dense)), 'Color',[0.5 0.5 0.5], 'Linewidth', 1);
hG = semilogx(f_exact, designGoal.A, 'ro', 'Markersize', 6, 'Linewidth', 1);
semilogx(f_exact, designGoal.C, 'ro', 'Markersize', 6, 'Linewidth', 1);
semilogx(f_exact, designGoal.Z, 'ro', 'Markersize', 6, 'Linewidth', 1);
grid on; xlim([10 20000]); ylim([-75 10]);
set(gca,'XTick',[10 100 1000 10000]);
xlabel('Frequency (Hz)', 'Interpreter','Latex');
ylabel('Frequency weighting (dB)', 'Interpreter','Latex');
legend([hA hC hZ hG], ...
       {'A, SQAT','C, SQAT','Z, SQAT','design goals, Table 3'}, ...
       'Location','SW', 'Interpreter','Latex');
legend box on;
set(gcf,'color','w');

%% Plot: deviation from the design goals against the acceptance limits

title_fig2 = 'Do SLM - deviation of the frequency weightings from the design goals';
h2 = figure('Name', title_fig2);
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% the limits are drawn as they are, except where no lower limit is
% specified, which is clipped to the bottom of the axes
lim_plot = 18; % dB
c1u = class1_upper;  c1l = max(class1_lower, -lim_plot);
c2u = class2_upper;  c2l = max(class2_lower, -lim_plot);

hc2 = semilogx(f_exact, c2u, 'Color',[0.6 0.6 0.6], 'Linewidth', 1); hold on;
semilogx(f_exact, c2l, 'Color',[0.6 0.6 0.6], 'Linewidth', 1);
hc1 = semilogx(f_exact, c1u, 'r--', 'Linewidth', 1);
semilogx(f_exact, c1l, 'r--', 'Linewidth', 1);
hA = semilogx(f_exact, dev_figs.A, 'ks-', 'Markersize', 6, 'Linewidth', 1);
hC = semilogx(f_exact, dev_figs.C, 'bx-', 'Markersize', 8, 'Linewidth', 1);
semilogx([10 20000], [0 0], 'k:', 'Linewidth', 0.5);
grid on; xlim([10 20000]); ylim([-lim_plot 6]);
set(gca,'XTick',[10 100 1000 10000]);
xlabel('Frequency (Hz)', 'Interpreter','Latex');
ylabel('$L_{\mathrm{SQAT}} - L_{\mathrm{design\;goal}}$ (dB)', 'Interpreter','Latex');
title(sprintf('$f_s$ = %g Hz', fs_figs), 'Interpreter','Latex');
legend([hA hC hc1 hc2], ...
       {'A weighting','C weighting','class 1 limits','class 2 limits'}, ...
       'Location','SW', 'Interpreter','Latex');
set(legend, 'Fontsize', 8);
legend box on;
set(gcf,'color','w');

%%%%%%%%% save figs
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = {'validation_Do_SLM_frequency_weightings', ...
                     'validation_Do_SLM_frequency_weightings_deviation'};
    hfig = [h1 h2];

    for i = 1:length(hfig)
        figname_out = [figures_dir figname_short{i}];
        saveas(hfig(i), figname_out, 'png');
        fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n', ...
                mfilename, figname_short{i}, figname_out);
    end
end
%%%%%%%%% save figs

if ~ok
    error('%s: the frequency weightings did not pass every check.', mfilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
