% Script validation_stationary_level_time_skip
%
% - Validation of the band level calculation of the stationary loudness
%   method of Loudness_ISO532_1, for a non zero time_skip.
%
% - Loudness computed using:
%   OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
%   type <help Loudness_ISO532_1> for more info
%
% - Two properties are verified, both of which hold by definition and
%   therefore need no measured data. This script generates every signal it
%   uses and runs straight after cloning the toolbox.
%
%   a) ANCHOR. ISO 532-1 [1] fixes the loudness of a 1 kHz free field tone at
%      40 dB SPL at 1 sone. That single point calibrates the whole scale, so
%      it is the natural absolute reference for the stationary method.
%
%   b) INVARIANCE. time_skip discards the leading part of the signal before
%      the band levels are computed. For a stationary signal every retained
%      window carries the same mean square, so the computed level and the
%      computed loudness must not depend on how much was discarded. Any
%      dependence on time_skip is an error of the implementation.
%
% - The reference C code supplied with ISO 532-1 sums samples NumSkip to
%   NumSamples-1 and divides by the number of samples actually summed,
%   NumSamples-NumSkip (see f_square_and_smooth in ISO_532-1.c). Dividing by
%   the full signal length instead scales the mean square by
%   (NumSamples-NumSkip)/NumSamples, so the band levels come out low by
%
%       10*log10( T / (T - time_skip) )   dB
%
%   uniformly across all 28 bands, where T is the signal duration. That is
%   the deviation this script detects. It vanishes at time_skip = 0, which is
%   the value every other validation script in this folder uses, so the
%   deviation stayed invisible until a non zero time_skip was verified.
%
% Author: Sergio Aguirre (& Claude code), 23.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Parameters

fs        = 48000;  % sampling frequency, Hz
T         = 5;      % signal duration, s
L_anchor  = 40;     % level of the anchor tone, dB SPL
f_anchor  = 1000;   % frequency of the anchor tone, Hz
field     = 0;      % free field
method    = 1;      % stationary, from the audio signal

time_skip = [0 0.1 0.2 0.5 1.0 2.0]; % values verified, s

tol_level = 0.05;   % tolerance on the level invariance, dB
tol_sone  = 0.02;   % tolerance on the anchor, sone

p0 = 2e-5; % reference sound pressure, Pa

t = (0:1/fs:T-1/fs).';

%% Test signals, both stationary and both at the anchor level

insig  = {};
labels = {};

insig{1}  = sqrt(2)*p0*10^(L_anchor/20) * sin(2*pi*f_anchor*t);
labels{1} = sprintf('%g Hz tone', f_anchor);

rng(42); % fixed seed so the script is reproducible
w = randn(size(t));
insig{2}  = p0*10^(L_anchor/20) * w/rms(w);
labels{2} = 'white noise';

% both signals carry the same overall sound pressure level
L_true = zeros(1,length(insig));
for i = 1:length(insig)
    L_true(i) = 20*log10(rms(insig{i})/p0);
end

%% Computation

nSig = length(insig);
nSkip = length(time_skip);

SPL  = zeros(nSig, nSkip);
N    = zeros(nSig, nSkip);

for i = 1:nSig
    for k = 1:nSkip
        OUT = Loudness_ISO532_1(insig{i}, fs, field, method, time_skip(k), 0);
        SPL(i,k) = OUT.TimeAveragedSPL;
        N(i,k)   = OUT.Loudness;
    end
end

% The invariance is checked against the value obtained with the largest
% amount of signal retained, which is the value time_skip = 0 returns.
dev_SPL = SPL - SPL(:,1);
dev_N   = N   - N(:,1);

%% Report

fprintf('\n%s\n', repmat('=',1,86));
fprintf('Stationary level calculation: dependence on time_skip\n');
fprintf('signal duration %g s, anchor level %g dB SPL\n', T, L_anchor);
fprintf('%s\n', repmat('=',1,86));

fprintf('%-16s', 'time_skip (s)');
fprintf(' %9.2f', time_skip);
fprintf('\n%s\n', repmat('-',1,86));
for i = 1:nSig
    fprintf('%-16s', labels{i});
    fprintf(' %9.4f', dev_SPL(i,:));
    fprintf('   <- level deviation (dB)\n');
    fprintf('%-16s', '');
    fprintf(' %9.4f', N(i,:));
    fprintf('   <- loudness (sone)\n');
end
fprintf('%s\n\n', repmat('=',1,86));

ok = true;

% a) the anchor: a 1 kHz free field tone at 40 dB SPL is 1 sone
N_anchor = N(1,:);
fprintf('ISO 532-1 anchor, %g Hz at %g dB SPL, free field, expected 1 sone:\n', ...
        f_anchor, L_anchor);
fprintf('   computed %.5f to %.5f sone across the time_skip values\n', ...
        min(N_anchor), max(N_anchor));
if all(abs(N_anchor - 1) < tol_sone)
    fprintf('   PASS, within %g sone of the definition\n\n', tol_sone);
else
    fprintf(2,'   FAIL, outside %g sone of the definition\n\n', tol_sone);
    ok = false;
end

% b) the invariance
worst = max(abs(dev_SPL(:)));
fprintf('Invariance of the band levels under time_skip:\n');
fprintf('   largest deviation %.4f dB\n', worst);
if worst < tol_level
    fprintf('   PASS, below the tolerance of %g dB\n\n', tol_level);
else
    fprintf(2,'   FAIL, above the tolerance of %g dB\n\n', tol_level);
    ok = false;
end

%% Plot: level deviation against time_skip

% Deviations produced by the implementation released up to v1.3, measured with
% the same signals and the same parameters as above. They are kept here so
% that the effect of the correction stays visible.
dev_v1 = [0 -0.0853 -0.1748 -0.4551 -0.9666 -2.2160;   % tone
          0 -0.0899 -0.1759 -0.4607 -0.9792 -2.2276];  % white noise

ts_curve  = linspace(0, max(time_skip), 300);
dev_curve = -10*log10(T ./ (T - ts_curve));

title_fig = 'Stationary level calculation - dependence on time skip';
h = figure('Name', title_fig);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

handle_c = plot(ts_curve, dev_curve, 'b-', 'Linewidth', 1); hold on;
handle_a = plot(time_skip, dev_v1(1,:), 'xb', 'Markersize', 12, 'Linewidth', 1);
handle_d = plot(time_skip, dev_v1(2,:), 'ob', 'Markersize', 10, 'Linewidth', 1);
handle_b = plot(time_skip, dev_SPL(1,:), 'sk', 'Markersize', 12, 'Linewidth', 1);
plot(time_skip, dev_SPL(2,:), 'dk', 'Markersize', 10, 'Linewidth', 1);
plot([0 max(time_skip)], [0 0], 'k--', 'Linewidth', 0.5);

xlabel('Discarded lead-in, $t_{\mathrm{skip}}$ (s)', 'Interpreter','Latex');
ylabel('$L_{\mathrm{SQAT}}(t_{\mathrm{skip}}) - L_{\mathrm{SQAT}}(0)$ (dB)', 'Interpreter','Latex');
legend([handle_c, handle_a, handle_d, handle_b], ...
       {'$-10\log_{10}(T/(T-t_{\mathrm{skip}}))$', 'tone, up to v1.3', ...
        'noise, up to v1.3', 'current'}, ...
       'Location', 'SW', 'Interpreter', 'Latex');
legend box on; grid on;
set(gcf,'color','w');

%%%%%%%%% save fig
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'validation_stationary_level_time_skip';
    figname_out = [figures_dir figname_short];

    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%% save fig

if ~ok
    error('%s: the stationary level calculation depends on time_skip.', mfilename);
end
