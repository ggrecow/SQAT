% Script validation_Do_SLM_repeated_tonebursts
%
% - Verification of the response of the sound level meter of SQAT to
%   sequences of 4 kHz tonebursts, against IEC 61672-1:2013, clause 5.10 and
%   Equation (9).
%
% - Levels computed using:
%   [outsig_dB, dBFS] = Do_SLM(insig, fs, weight_freq, weight_time, dBFS)
%   Leq = Get_Leq(levels, fs, dt)
%   type <help Do_SLM> and <help Get_Leq> for more info
%
% - Reference values come from the standard. Clause 5.10.3 states that, for
%   a sequence of n tonebursts of duration Tb extracted from a steady 4 kHz
%   sinusoid within a total measurement duration Tm, the difference between
%   the time-averaged sound level of the sequence and the time-averaged
%   sound level of the steady sinusoid is
%
%       delta_ref = 10*log10(n*Tb/Tm)  dB                           Eq. (9)
%
%   Clause 5.10.1 requires the measured deviations from that value to stay
%   inside the acceptance limits given in Table 4 for the sound exposure
%   level toneburst response, which are the limits of the row of the
%   toneburst duration used. The specification applies to the A weighting
%   and to the C and Z weightings where provided, so the three of them are
%   tested here.
%
% - This check exercises the calibration and the energy averaging of the
%   time-weighted level rather than the time constants: the exponential
%   time weighting has unit gain at zero frequency, so its output averaged
%   over an integer number of repetition periods returns the mean square of
%   the sequence for either time weighting. The value of n*Tb/Tm reached
%   here is 0,0125, which is 19 dB below the steady level, so the sequences
%   are far more impulsive than the signals a fixed calibration offset can
%   accommodate.
%
% Author: Sergio Aguirre (& Claude code), 28.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Parameters

fs        = 48000; % sampling frequency, Hz
f_burst   = 4000;  % toneburst frequency, Hz, as specified in 5.10.1
dBFS      = 94;    % full scale convention, dB SPL for a unit rms signal
tol       = 0.2;   % tolerance on the deviation from Equation (9), dB

weight_freq = {'Z','A','C'};
weight_time = {'f','s'};

tau = [0.125 1]; % exponential time constants of 5.8.1, s, for F and S

% two sequences of tonebursts, each with its own repetition period Tp and
% number of repetitions n within the total measurement duration Tm = n*Tp
sequence(1).Tp    = 1;      % s
sequence(1).n     = 10;
sequence(1).Tb_ms = [500 200 100 50 20 10];
sequence(2).Tp    = 0.020;  % s
sequence(2).n     = 200;
sequence(2).Tb_ms = [10 5 2 1 0.5 0.25];

%% IEC 61672-1:2013, Table 4, acceptance limits

Tb_table_ms  = [1000 500 200 100 50 20 10 5 2 1 0.5 0.25];
class1_upper = [0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];
class1_lower = [-0.5 -0.5 -0.5 -1.0 -1.0 -1.0 -1.0 -1.0 -1.5 -2.0 -2.5 -3.0];
class2_upper = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.5];
class2_lower = [-1.0 -1.0 -1.0 -1.0 -1.5 -2.0 -2.0 -2.5 -2.5 -3.0 -4.0 -5.0];

%% Response to the toneburst sequences

verdict   = {'FAIL','PASS'};
ok_class1 = true;
ok_class2 = true;
ok_tol    = true;

for is = 1:length(sequence)

    Tp = sequence(is).Tp;
    n  = sequence(is).n;
    Tm = n*Tp;
    Tb = sequence(is).Tb_ms(:)/1000;
    nb = length(Tb);

    delta_ref = 10*log10(n*Tb/Tm); % Eq. (9)

    idx_lim = zeros(nb,1);
    for i = 1:nb
        idx_lim(i) = find(Tb_table_ms == sequence(is).Tb_ms(i), 1);
    end

    for it = 1:length(weight_time)

        tw     = weight_time{it};
        t_lead = 10*tau(it); % lead-in discarded, so that the time weighting
                             % is in its steady state over the whole of Tm
        dur    = t_lead + Tm;
        t      = (0:1/fs:dur-1/fs).';
        i_meas = round(t_lead*fs)+1;

        steady = sin(2*pi*f_burst*t);
        phase  = mod(t/Tp, 1); % sawtooth between 0 and 1 at the repetition rate

        for iw = 1:length(weight_freq)

            wf = weight_freq{iw};

            lvl = Do_SLM(steady, fs, wf, tw, dBFS);
            L_steady = Get_Leq(lvl(i_meas:end));

            delta = zeros(nb,1);
            for i = 1:nb
                insig = steady .* double(phase < Tb(i)/Tp);
                lvl   = Do_SLM(insig, fs, wf, tw, dBFS);
                delta(i) = Get_Leq(lvl(i_meas:end)) - L_steady;
            end

            dev  = delta - delta_ref;
            in1  = all(dev <= class1_upper(idx_lim).') && ...
                   all(dev >= class1_lower(idx_lim).');
            in2  = all(dev <= class2_upper(idx_lim).') && ...
                   all(dev >= class2_lower(idx_lim).');
            intl = all(abs(dev) < tol);

            ok_class1 = ok_class1 && in1;
            ok_class2 = ok_class2 && in2;
            ok_tol    = ok_tol    && intl;

            fprintf('\n%s\n', repmat('=',1,74));
            fprintf(['Repeated tonebursts, %s weighting, time weighting %s, ' ...
                     'n = %d, Tm = %g s\n'], wf, upper(tw), n, Tm);
            fprintf('steady 4 kHz sound level: %.2f dB\n', L_steady);
            fprintf('%s\n', repmat('=',1,74));
            fprintf('%10s %10s %10s %10s %10s %12s\n', ...
                    'Tb', 'n*Tb/Tm', 'SQAT', 'Eq. (9)', 'deviation', 'class 1 lim');
            fprintf('%10s %10s %10s %10s %10s %12s\n', ...
                    'ms', '-', 'dB', 'dB', 'dB', 'dB');
            fprintf('%s\n', repmat('-',1,74));
            for i = 1:nb
                fprintf('%10g %10.4f %10.3f %10.3f %10.3f %6.1f %5.1f\n', ...
                        1000*Tb(i), n*Tb(i)/Tm, delta(i), delta_ref(i), ...
                        dev(i), class1_upper(idx_lim(i)), class1_lower(idx_lim(i)));
            end
            fprintf('%s\n', repmat('-',1,74));
            [~, i_dev] = max(abs(dev));
            fprintf('largest deviation from Equation (9) %+.3f dB, class 1 %s, class 2 %s\n', ...
                    dev(i_dev), verdict{in1+1}, verdict{in2+1});

            % kept for the figure
            if it == 1
                store.(wf).duty{is}  = n*Tb/Tm;
                store.(wf).delta{is} = delta;
            end
        end
    end
end

%% Verdict

fprintf('\n%s\n', repmat('=',1,74));
ok = ok_class1 && ok_class2 && ok_tol;
if ok
    fprintf(['PASS: every toneburst sequence meets the acceptance limits of ' ...
             'classes 1 and 2,\n      and stays within %.2f dB of Equation ' ...
             '(9).\n'], tol);
else
    fprintf(2, 'FAIL: see the tables above.\n');
end
fprintf('%s\n\n', repmat('=',1,74));

%% Plot: time-averaged level of the sequences against Equation (9)

title_fig = 'Do SLM - response to repeated 4 kHz tonebursts';
h = figure('Name', title_fig);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

duty_curve = logspace(log10(0.01), log10(0.6), 100);

hc = semilogx(duty_curve, 10*log10(duty_curve), 'b-', 'Linewidth', 1); hold on;
for is = 1:length(sequence)
    hz = semilogx(store.Z.duty{is}, store.Z.delta{is}, 'ks', ...
                  'Markersize', 10, 'Linewidth', 1);
    ha = semilogx(store.A.duty{is}, store.A.delta{is}, 'rx', ...
                  'Markersize', 10, 'Linewidth', 1);
end
grid on;
xlim([0.01 0.6]); ylim([-21 0]);
xlabel('Duty cycle of the sequence, $nT_b/T_m$', 'Interpreter','Latex');
ylabel('$L_{\mathrm{eq,bursts}} - L_{\mathrm{eq,steady}}$ (dB)', 'Interpreter','Latex');
legend([hz ha hc], {'SQAT, Z weighting', 'SQAT, A weighting', 'Equation (9)'}, ...
       'Location','NW', 'Interpreter','Latex');
legend box on;
set(gcf,'color','w');

%%%%%%%%% save fig
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'validation_Do_SLM_repeated_tonebursts';
    figname_out = [figures_dir figname_short];

    saveas(gcf, figname_out, 'png');

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n', ...
            mfilename, figname_short, figname_out);
end
%%%%%%%%% save fig

if ~ok
    error('%s: the responses to repeated tonebursts did not pass every check.', mfilename);
end
