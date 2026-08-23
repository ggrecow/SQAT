% Script validation_Do_SLM_time_weighting
%
% - Validation of the exponential time weighting applied by Do_SLM, using
%   synthetic signals of known crest factor.
%
% - Levels computed using:
%   [outsig_dB, dBFS] = Do_SLM(insig, fs, weight_freq, weight_time, dBFS)
%   Leq = Get_Leq(levels, fs, dt)
%   type <help Do_SLM> and <help Get_Leq> for more info
%
% - IEC 61672-1 defines exponential time weighting as a first order low-pass
%   applied to the SQUARED sound pressure, with the level taken as 10*log10
%   of the resulting mean square. Averaging that quantity over the whole
%   signal returns the equivalent continuous sound level, so for a stationary
%   input the result must equal 20*log10(rms(p_weighted)/p0) whatever the
%   waveform is. That identity is the reference used here.
%
% - The reference needs no measured data, so unlike the other validation
%   scripts in this folder this one is self-contained: it generates every
%   signal it uses and can be run straight after cloning the toolbox.
%
% - Three families of signals are used, all scaled to the same root mean
%   square and therefore all representing the same true sound pressure level:
%
%     a) a 1 kHz sinusoid                         crest factor of a sine
%     b) Gaussian white noise                     crest factor of noise
%     c) gated 1 kHz tone bursts, duty cycle d    crest factor set by d
%
%   Family (c) sweeps the crest factor over more than one decade. For a tone
%   gated with duty cycle d, the ratio between the mean absolute value and
%   the root mean square is (2*sqrt(2)/pi)*sqrt(d), so any implementation
%   that smooths the absolute value in place of the square reports a level
%   low by 20*log10(0.9003*sqrt(d)) dB. That is the deviation this script
%   detects, and it grows without bound as the signal becomes more impulsive.
%
% Author: Sergio Aguirre (& Claude code), 23.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Parameters

fs        = 48000;  % sampling frequency, Hz
dur       = 10;     % signal duration, s
dBFS      = 94;     % full scale convention, dB SPL for a unit rms signal
f_carrier = 1000;   % carrier frequency of the tonal signals, Hz
f_gate    = 200;    % repetition rate of the tone bursts, Hz
skip_fast = 1.5;    % lead-in discarded under fast weighting, s
skip_slow = 5.0;    % lead-in discarded under slow weighting, s (tau = 1 s)
tol       = 0.05;   % tolerance on the deviation from the reference, dB

p0 = 2e-5; % reference sound pressure, Pa

t = (0:1/fs:dur-1/fs).';

%% Build the test signals

duty = [0.5 0.25 0.1 0.05 0.02]; % duty cycle of the gated tone bursts

insig  = {};
labels = {};

insig{1}  = sin(2*pi*f_carrier*t);
labels{1} = 'sinusoid 1 kHz';

rng(42); % fixed seed so the script is reproducible
insig{2}  = randn(size(t));
labels{2} = 'white noise';

phase_gate = mod(t*f_gate, 1); % sawtooth between 0 and 1 at the gate rate
for i = 1:length(duty)
    gate = double(phase_gate < duty(i));
    insig{end+1}  = sin(2*pi*f_carrier*t) .* gate;
    labels{end+1} = sprintf('tone bursts, duty %.2f', duty(i));
end

% scale every signal to unit root mean square, so that all of them share the
% same true sound pressure level and only the crest factor differs
for i = 1:length(insig)
    insig{i} = insig{i} / rms(insig{i});
end

%% Reference and computed levels

weight_freq = {'Z','A'};
weight_time = {'f','s'};

nSig = length(insig);
crest    = zeros(nSig,1);
form     = zeros(nSig,1);
L_ref    = zeros(nSig, 4);
L_sqat   = zeros(nSig, 4);
col_name = cell(1,4);

col = 0;
for iw = 1:length(weight_freq)
    for it = 1:length(weight_time)

        col = col + 1;
        col_name{col} = sprintf('%s%s', weight_freq{iw}, upper(weight_time{it}));

        if strcmp(weight_time{it},'f')
            NumSkip = round(skip_fast*fs);
        else
            NumSkip = round(skip_slow*fs);
        end

        [b,a] = Gen_weighting_filters(fs, weight_freq{iw});

        for i = 1:nSig

            % level reported by the toolbox
            lvl = Do_SLM(insig{i}, fs, weight_freq{iw}, weight_time{it}, dBFS);
            lvl = lvl(NumSkip:end);
            L_sqat(i,col) = Get_Leq(lvl, fs, length(lvl)/fs);

            % reference: equivalent continuous level of the weighted pressure
            p  = filter(b, a, insig{i} * 10^((dBFS-94)/20));
            p  = p(NumSkip:end);
            L_ref(i,col) = 20*log10(rms(p)/p0);

            if col == 1
                crest(i) = max(abs(p))/rms(p);
                form(i)  = mean(abs(p))/rms(p);
            end
        end
    end
end

dev = L_sqat - L_ref; % deviation from the reference, dB

%% Report

fprintf('\n%s\n', repmat('=',1,87));
fprintf('Do_SLM time weighting: deviation from 20*log10(rms(p_weighted)/p0)\n');
fprintf('%s\n', repmat('=',1,87));
fprintf('%-24s %8s %8s', 'signal', 'crest', 'form');
fprintf(' %8s', col_name{:});
fprintf('\n%s\n', repmat('-',1,87));
for i = 1:nSig
    fprintf('%-24s %8.2f %8.4f', labels{i}, crest(i), form(i));
    fprintf(' %+8.3f', dev(i,:));
    fprintf('\n');
end
fprintf('%s\n', repmat('-',1,87));
fprintf('%-24s %8s %8s', 'largest absolute value', '', '');
fprintf(' %8.3f', max(abs(dev),[],1));
fprintf('\n%s\n\n', repmat('=',1,87));

ok = all(abs(dev(:)) < tol);
if ok
    fprintf('PASS: every deviation is below the tolerance of %.3f dB.\n\n', tol);
else
    [worst, idx] = max(abs(dev(:)));
    [iw_, ic_]   = ind2sub(size(dev), idx);
    fprintf(2, 'FAIL: largest deviation is %.3f dB on %s (%s), tolerance is %.3f dB.\n\n', ...
            worst, labels{iw_}, col_name{ic_}, tol);
end

%% Plot: deviation against crest factor

% Deviations produced by the implementation released up to v1.3, measured with
% the same signals and the same parameters as above, column ZF. They are kept
% here so that the effect of the correction stays visible.
dev_v1 = [0.0055; -1.0304; -3.0047; -6.0257; -9.9940; -12.8564; -17.5983];

title_fig = 'Do SLM - deviation of the time weighted level against form factor';
h = figure('Name', title_fig);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Analytical deviation of an implementation that smooths the absolute value
% of the pressure in place of its square. The level it reports is
% 20*log10(mean|p|), the correct level is 20*log10(rms(p)), and the release
% up to v1.3 added a fixed offset of 0.93 dB, so the deviation is a function
% of the form factor alone.
ff_curve  = logspace(log10(0.10), log10(0.95), 300);
dev_curve = 20*log10(ff_curve) + 0.93;

handle_c = semilogx(ff_curve, dev_curve, 'b-', 'Linewidth', 1); hold on;
handle_a = semilogx(form, dev_v1,  'xb', 'Markersize', 12, 'Linewidth', 1);
handle_b = semilogx(form, dev(:,1),'sk', 'Markersize', 12, 'Linewidth', 1);
semilogx([min(ff_curve) max(ff_curve)], [0 0], 'k--', 'Linewidth', 0.5);
xlim([min(ff_curve) max(ff_curve)]);

xlabel('Form factor of the weighted signal, $\overline{|p|}/p_{\mathrm{rms}}$', 'Interpreter','Latex');
ylabel('$L_{\mathrm{Zeq,SQAT}} - L_{\mathrm{Zeq,Ref.}}$ (dB)', 'Interpreter','Latex');
legend([handle_c, handle_a, handle_b], ...
       {'$20\log_{10}(\overline{|p|}/p_{\mathrm{rms}}) + 0.93$', 'up to v1.3', 'current'}, ...
       'Location', 'SE', 'Interpreter', 'Latex');
legend box on; grid on;
set(gcf,'color','w');

%%%%%%%%% save fig
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'validation_Do_SLM_deviation_vs_form_factor';
    figname_out = [figures_dir figname_short];

    % saveas(gcf,figname_out, 'fig');
    % saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%% save fig

if ~ok
    error('%s: the computed levels are outside the tolerance of %g dB.', mfilename, tol);
end
