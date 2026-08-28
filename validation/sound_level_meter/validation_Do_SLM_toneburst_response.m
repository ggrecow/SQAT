% Script validation_Do_SLM_toneburst_response
%
% - Verification of the exponential time constants of the sound level meter
%   of SQAT against the 4 kHz toneburst responses of IEC 61672-1:2013,
%   Table 4 and Equation (7).
%
% - Levels computed using:
%   [outsig_dB, dBFS] = Do_SLM(insig, fs, weight_freq, weight_time, dBFS)
%   Leq = Get_Leq(levels, fs, dt)
%   type <help Do_SLM> and <help Get_Leq> for more info
%
% - Reference values come from the standard. Clause 5.9.2 defines the
%   toneburst response as the difference between the maximum time-weighted
%   sound level of a single 4 kHz toneburst of duration Tb, extracted from a
%   steady 4 kHz sinusoid, and the sound level of that steady sinusoid.
%   Table 4 gives the reference response and the acceptance limits for
%   performance classes 1 and 2, for twelve durations from 0,25 ms to 1 s
%   under time weighting F and for nine durations from 2 ms to 1 s under
%   time weighting S. The reference response is
%
%       delta_ref = 10*log10(1 - exp(-Tb/tau))  dB                  Eq. (7)
%
%   with tau the exponential time constant of 5.8.1, 0,125 s for F and 1 s
%   for S. Both the tabulated values and Equation (7) are used here, the
%   first for the conformance check and the second for a tight regression
%   check. According to NOTE 3 of Table 4 the reference responses hold for
%   the A, C and Z frequency weightings, so the three of them are tested.
%
% - This is the check that fixes the time constants. The response of a short
%   toneburst falls as 10*log10(Tb/tau), so a time constant in error by a
%   factor k moves every short burst by 10*log10(k) dB, which the
%   acceptance limits of Table 4 detect at once.
%
% Author: Sergio Aguirre (& Claude code), 28.08.2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Parameters

fs        = 48000; % sampling frequency, Hz
f_burst   = 4000;  % toneburst frequency, Hz, as specified in 5.9.1
dBFS      = 94;    % full scale convention, dB SPL for a unit rms signal
tol       = 0.2;   % tolerance on the deviation from Equation (7), dB

weight_freq = {'Z','A','C'};
weight_time = {'f','s'};

tau = [0.125 1]; % exponential time constants of 5.8.1, s, for F and S

%% IEC 61672-1:2013, Table 4

% toneburst durations, s
Tb.f = [1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.00025].';
Tb.s = [1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002].';

% reference toneburst response relative to the steady sound level, dB
ref.f = [0.0 -0.1 -1.0 -2.6 -4.8 -8.3 -11.1 -14.1 -18.0 -21.0 -24.0 -27.0].';
ref.s = [-2.0 -4.1 -7.4 -10.2 -13.1 -17.0 -20.0 -23.0 -27.0].';

% acceptance limits, dB
class1_upper.f = [0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0].';
class1_lower.f = [-0.5 -0.5 -0.5 -1.0 -1.0 -1.0 -1.0 -1.0 -1.5 -2.0 -2.5 -3.0].';
class2_upper.f = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.5].';
class2_lower.f = [-1.0 -1.0 -1.0 -1.0 -1.5 -2.0 -2.0 -2.5 -2.5 -3.0 -4.0 -5.0].';

class1_upper.s = [0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0].';
class1_lower.s = [-0.5 -0.5 -0.5 -1.0 -1.0 -1.5 -2.0 -2.5 -3.0].';
class2_upper.s = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0].';
class2_lower.s = [-1.0 -1.0 -1.0 -1.0 -1.5 -2.0 -3.0 -4.0 -5.0].';

%% Toneburst responses

verdict   = {'FAIL','PASS'};
ok_class1 = true;
ok_class2 = true;
ok_tol    = true;

for it = 1:length(weight_time)

    tw = weight_time{it};

    % the steady sound level is measured on a sinusoid long enough for the
    % leaky integrator to settle: the first 20*tau are discarded and the
    % level is averaged over the following 10*tau
    t_lead   = 20*tau(it);
    dur      = t_lead + 10*tau(it);
    t        = (0:1/fs:dur-1/fs).';
    steady   = sin(2*pi*f_burst*t);

    for iw = 1:length(weight_freq)

        wf = weight_freq{iw};

        lvl = Do_SLM(steady, fs, wf, tw, dBFS);
        L_steady = Get_Leq(lvl(round(t_lead*fs)+1:end));

        nb   = length(Tb.(tw));
        resp = zeros(nb,1);

        for i = 1:nb

            % single toneburst extracted from the steady sinusoid, so that
            % it starts and ends at a zero crossing
            N_burst = round(Tb.(tw)(i)*fs);
            N_total = round(0.05*fs) + N_burst + round(8*tau(it)*fs);
            insig   = zeros(N_total,1);
            idx     = round(0.05*fs) + (1:N_burst);
            insig(idx) = sin(2*pi*f_burst*(0:N_burst-1).'/fs);

            resp(i) = max(Do_SLM(insig, fs, wf, tw, dBFS)) - L_steady;
        end

        eq7 = 10*log10(1 - exp(-Tb.(tw)/tau(it)));
        dev = resp - ref.(tw);  % deviation from the reference of Table 4
        dif = resp - eq7;       % deviation from Equation (7)

        in1  = all(dev <= class1_upper.(tw)) && all(dev >= class1_lower.(tw));
        in2  = all(dev <= class2_upper.(tw)) && all(dev >= class2_lower.(tw));
        intl = all(abs(dif) < tol);

        ok_class1 = ok_class1 && in1;
        ok_class2 = ok_class2 && in2;
        ok_tol    = ok_tol    && intl;

        fprintf('\n%s\n', repmat('=',1,79));
        fprintf(['Toneburst response, %s weighting, time weighting %s ' ...
                 '(tau = %g s)\n'], wf, upper(tw), tau(it));
        fprintf('steady 4 kHz sound level: %.2f dB\n', L_steady);
        fprintf('%s\n', repmat('=',1,79));
        fprintf('%10s %10s %10s %10s %10s %12s\n', ...
                'Tb', 'SQAT', 'Table 4', 'Eq. (7)', 'SQAT-Eq7', 'class 1 lim');
        fprintf('%10s %10s %10s %10s %10s %12s\n', ...
                'ms', 'dB', 'dB', 'dB', 'dB', 'dB');
        fprintf('%s\n', repmat('-',1,79));
        for i = 1:nb
            fprintf('%10g %10.3f %10.1f %10.3f %10.3f %6.1f %5.1f\n', ...
                    1000*Tb.(tw)(i), resp(i), ref.(tw)(i), eq7(i), dif(i), ...
                    class1_upper.(tw)(i), class1_lower.(tw)(i));
        end
        fprintf('%s\n', repmat('-',1,79));
        [~, i_dev] = max(abs(dev));
        [~, i_dif] = max(abs(dif));
        fprintf(['largest deviation from Table 4 %+.3f dB, from Equation (7) ' ...
                 '%+.3f dB\n'], dev(i_dev), dif(i_dif));
        fprintf('class 1 %s, class 2 %s\n', verdict{in1+1}, verdict{in2+1});

        % kept for the figure
        response.([wf tw]) = resp;
    end
end

%% Verdict

fprintf('\n%s\n', repmat('=',1,79));
ok = ok_class1 && ok_class2 && ok_tol;
if ok
    fprintf(['PASS: every toneburst response meets the acceptance limits of ' ...
             'classes 1 and 2,\n      and stays within %.2f dB of Equation ' ...
             '(7).\n'], tol);
else
    fprintf(2, 'FAIL: see the tables above.\n');
end
fprintf('%s\n\n', repmat('=',1,79));

%% Plot: toneburst response against the reference of Table 4

title_fig = 'Do SLM - 4 kHz toneburst response';
h = figure('Name', title_fig);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

Tb_curve = logspace(log10(0.00025), log10(1), 300).';

hc = semilogx(1000*Tb_curve, 10*log10(1 - exp(-Tb_curve/tau(1))), ...
              'b-', 'Linewidth', 1); hold on;
semilogx(1000*Tb_curve, 10*log10(1 - exp(-Tb_curve/tau(2))), ...
         'b-', 'Linewidth', 1);
he = errorbar(1000*Tb.f, ref.f, -class1_lower.f, class1_upper.f, ...
              'Color', [0.6 0.6 0.6], 'Linestyle', 'none', 'Linewidth', 1);
errorbar(1000*Tb.s, ref.s, -class1_lower.s, class1_upper.s, ...
         'Color', [0.6 0.6 0.6], 'Linestyle', 'none', 'Linewidth', 1);
hz = semilogx(1000*Tb.f, response.Zf, 'ks', 'Markersize', 10, 'Linewidth', 1);
semilogx(1000*Tb.s, response.Zs, 'ks', 'Markersize', 10, 'Linewidth', 1);
ha = semilogx(1000*Tb.f, response.Af, 'rx', 'Markersize', 10, 'Linewidth', 1);
semilogx(1000*Tb.s, response.As, 'rx', 'Markersize', 10, 'Linewidth', 1);

grid on;
xlim([0.2 1200]); ylim([-30 3]);
text(250, 1.2, 'F', 'Interpreter','Latex', 'Fontsize', 12);
text(250, -4.6, 'S', 'Interpreter','Latex', 'Fontsize', 12);
xlabel('Toneburst duration, $T_b$ (ms)', 'Interpreter','Latex');
ylabel('$L_{\mathrm{max}} - L_{\mathrm{steady}}$ (dB)', 'Interpreter','Latex');
legend([hz ha hc he], {'SQAT, Z weighting', 'SQAT, A weighting', ...
                       'Equation (7)', 'Table 4, class 1 limits'}, ...
       'Location','NW', 'Interpreter','Latex');
legend box on;
set(gcf,'color','w');

%%%%%%%%% save fig
if save_figs==1

    figures_dir = [pwd filesep 'figs' filesep];

    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end

    figname_short = 'validation_Do_SLM_toneburst_response';
    figname_out = [figures_dir figname_short];

    saveas(gcf, figname_out, 'png');

    fprintf('\n%s.m: figure %s was saved on disk\n\t(full name: %s)\n', ...
            mfilename, figname_short, figname_out);
end
%%%%%%%%% save fig

if ~ok
    error('%s: the toneburst responses did not pass every check.', mfilename);
end
