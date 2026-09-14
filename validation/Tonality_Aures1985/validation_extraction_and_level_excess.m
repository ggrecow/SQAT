% Script validation_extraction_and_level_excess
%
% - This routine compares two stages of Tonality_Aures1985, the extraction of
%   sinusoidal components and the sound pressure level excess, both taken from
%   Terhardt et al. [2], with the intermediate values published by Zhang and
%   Shrestha [3]. That thesis implements the extraction and masking stages of
%   [2] from the papers, at IMM/DTU with Bruel and Kjaer, independently of the
%   lineage this code descends from. Its Appendix C prints the two power
%   spectra of its tests, 401 samples 10.77 Hz apart, and its Tables 6.2 and
%   6.3 print the components each test gave and their level excess. It is the
%   only source at hand that can tell an error of this implementation from one
%   shared with its ancestors.
%
% - The spectra are read from reference_values and injected into the two
%   stages through the local functions of the metric, which
%   Tonality_Aures1985('localfunctions') hands out. Test 1 uses 3 dB in the
%   criterion, which is what [3] used for its tables (p. 49 of the thesis: with
%   7 dB too few components were left to exercise the later stages); Test 2
%   uses the 7 dB of [2].
%
% - CHECK 1, extraction: the components found are those of Tables 6.2 and 6.3,
%   each within one sample in frequency after the refinement of eq. (3) of
%   [3], with the same level to 0.01 dB.
% - CHECK 2, excitation and threshold terms of the level excess: on the
%   components masked by the other components rather than by the noise, the
%   level excess of [3] is reproduced within 0.3 dB when the noise term is
%   left out.
% - CHECK 3, noise term: the term the metric uses equals the sum of the
%   intensities of the samples within half a Bark on each side of the
%   component, skipping the five central samples of every component in that
%   band, as [2] and [3] define it, within 0.1 dB. The sum is recomputed here
%   from the published spectrum.
% - REPORT: the level excess of [3] cannot verify the noise term. With that
%   term replaced by one constant, the formula reproduces all twelve published
%   values within 1.6 dB, while the sum the definition asks for reads 31 to
%   62 dB on those spectra. The noise term of [3] does not follow its own
%   eq. (4), so its level excesses hold no information about that term.
%
% - Tonality computed using
%   OUT = Tonality_Aures1985(insig,fs,LoudnessField,time_skip,show)
%   type <help Tonality_Aures1985> for more info
%
% - References:
%   [1] Aures, W. (1985). Berechnungsverfahren fuer den sensorischen Wohlklang
%       beliebiger Schallsignale. Acustica 59(2), 130-141.
%   [2] Terhardt, E., Stoll, G. and Seewann, M. (1982). Algorithm for
%       extraction of pitch and pitch salience from complex tonal signals.
%       J. Acoust. Soc. Am. 71(3), 679-688.
%   [3] Zhang, Z. and Shrestha, M. (2003). Sound Quality User-defined Cursor
%       Reading Control, Tonality Metric. Master thesis, IMM-Thesis-2003-22,
%       Technical University of Denmark, with Bruel and Kjaer.
%       http://www2.imm.dtu.dk/pubdb/edoc/imm2385.pdf
%
% Author: Sergio Aguirre, September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; %#ok<CLALL>

%% save settings
save_figs = 0;
dir_out = [fileparts(mfilename('fullpath')) filesep];

%% the local functions of the metric
fh    = Tonality_Aures1985('localfunctions');
names = cellfun(@func2str, fh, 'UniformOutput', false);
find_sinusoids = fh{strcmp(names, 'il_find_sinusoids')};
SPL_excess     = fh{strcmp(names, 'il_SPL_excess')};
Fq2Bark        = fh{strcmp(names, 'il_Fq2Bark')};
Threshold      = fh{strcmp(names, 'il_Threshold')};

%% the two tests of [3]: spectrum, threshold of the criterion, Tables 6.2 and 6.3
test(1).file = 'ZhangShrestha2003_AppendixC_Test1.txt';
test(1).thr  = 3;
test(1).comp = [ 387.89 58.13;  765.95 33.20;  989.05 26.51; 1528.19 27.16; ...
                1871.83 39.69; 2078.14 46.70; 2509.13 41.67; 2636.70 43.89; ...
                2863.81 37.18; 3091.05 35.01; 3348.17 40.62; 3467.18 39.10; ...
                3595.04 35.01; 3864.10 26.76];
test(1).LX   = [ 387.89 39.02;  765.95 12.86;  989.05  9.11; 1528.19  9.67; ...
                1871.83 10.35; 2078.14 16.20; 2509.13  1.52; 2636.70  3.04];
test(2).file = 'ZhangShrestha2003_AppendixC_Test2.txt';
test(2).thr  = 7;
test(2).comp = [ 387.59 87.87;  807.49 90.96; 1410.41 90.97; 2196.38 87.87];
test(2).LX   = [ 387.59 67.83;  807.49 27.40; 1410.41 21.61; 2196.38 13.66];

nfail = 0;
allL = []; allAEK = []; allEHS = []; allLXref = [];   % for the report

for k = 1:2

    A  = readmatrix([dir_out 'reference_values' filesep test(k).file], 'CommentStyle', '#');
    f  = A(:,1);  L = A(:,2);  df = f(2)-f(1);
    fprintf('Test %d of [3]: %d samples %.2f Hz apart, criterion at %d dB\n', k, numel(f), df, test(k).thr);

    %% CHECK 1, extraction
    idx  = find_sinusoids(L, 1, numel(L), df, test(k).thr);
    fref = f(idx) + 0.46*(L(idx+1) - L(idx-1));    % eq. (3) of [3], for the comparison only
    Lc   = L(idx);
    ok   = numel(idx) == size(test(k).comp,1);
    fprintf('%10s %10s %8s %10s %8s\n', 'found', 'refined', 'level', 'table', 'level');
    for i = 1:numel(idx)
        [d, m] = min(abs(test(k).comp(:,1) - fref(i)));
        hit = d < df && abs(Lc(i) - test(k).comp(m,2)) <= 0.01;
        ok  = ok && hit;
        fprintf('%10.2f %10.2f %8.2f %10.2f %8.2f %s\n', f(idx(i)), fref(i), Lc(i), test(k).comp(m,1), test(k).comp(m,2), char(repmat('OK',1,hit)));
    end
    if ok
        fprintf('CHECK 1 PASSED: the %d components of Table 6.%d are found, each within one sample and 0.01 dB\n', numel(idx), k+1);
    else
        fprintf('CHECK 1 FAILED: %d components found for the %d of Table 6.%d\n', numel(idx), size(test(k).comp,1), k+1); nfail = nfail+1;
    end

    %% the level excess of the metric, with and without its noise term
    in = struct('freq', f, 'ToneF', f(idx), 'ToneL', Lc, 'Lnoise', L);
    LX_full = SPL_excess(in);
    in.Lnoise = -300*ones(size(L));                 % no noise at all
    LX_nonoise = SPL_excess(in);

    %% the terms of eq. (4) of [3], recomputed here to sort the components
    zc = Fq2Bark(f(idx)); z = Fq2Bark(f); n = numel(idx);
    AEK2 = zeros(n,1); EHS = zeros(n,1); EGR = zeros(n,1);
    for i = 1:n
        a = 0;
        for j = 1:n
            if j == i, continue; end
            if j < i, s = -24 - 230/f(idx(j)) + 0.2*Lc(j); else, s = 27; end   % eq. (7) of [3]
            a = a + 10^((Lc(j) - s*(zc(j) - zc(i)))/20);                       % eq. (5) of [3]
        end
        AEK2(i) = a^2;
        EHS(i)  = 10^(Threshold(f(idx(i)))/10);
        band = z >= zc(i)-0.5 & z <= zc(i)+0.5;                                % the critical band of eq. (4)
        for j = 1:n
            if band(idx(j)), band(max(idx(j)-2,1):min(idx(j)+2,numel(L))) = false; end
        end
        EGR(i) = sum(10.^(L(band)/10));
    end

    %% CHECK 2, excitation and threshold terms, on the components masked by the others
    masked = 10*log10(AEK2) > 25;
    fprintf('%10s %10s %10s %10s %10s\n', 'f', 'LX table', 'LX metric', 'no noise', 'masking');
    dev = [];
    for i = 1:n
        [d, m] = min(abs(test(k).LX(:,1) - fref(i)));
        if d < df, ref = test(k).LX(m,2); else, ref = NaN; end
        fprintf('%10.2f %10.2f %10.2f %10.2f %10.1f\n', fref(i), ref, LX_full(i), LX_nonoise(i), 10*log10(AEK2(i)));
        if masked(i) && ~isnan(ref)
            dev(end+1) = LX_nonoise(i) - ref; %#ok<SAGROW>
        end
        if ~isnan(ref)
            allL(end+1) = Lc(i); allAEK(end+1) = AEK2(i); allEHS(end+1) = EHS(i); allLXref(end+1) = ref; %#ok<SAGROW>
        end
    end
    if ~isempty(dev) && max(abs(dev)) <= 0.3
        fprintf('CHECK 2 PASSED: on the %d components masked by the others, the level excess of [3] is reproduced within %.2f dB with the noise term left out\n', numel(dev), max(abs(dev)));
    else
        fprintf('CHECK 2 FAILED: largest distance %.2f dB on %d components\n', max(abs(dev)), numel(dev)); nfail = nfail+1;
    end

    %% CHECK 3, the noise term of the metric against the definition
    use = LX_full > 0;      % the metric returns zero where the excess is negative
    EGR_metric = 10.^((Lc(use) - LX_full(use))/10) - 10.^((Lc(use) - LX_nonoise(use))/10);
    d3 = 10*log10(EGR_metric) - 10*log10(EGR(use));
    if any(use) && max(abs(d3)) <= 0.1
        fprintf('CHECK 3 PASSED: on the %d components with a positive excess, the noise term of the metric equals the sum of eq. (4) within %.3f dB\n', sum(use), max(abs(d3)));
    elseif any(use)
        fprintf('CHECK 3 FAILED: the noise term of the metric departs from the sum of eq. (4) by %.2f dB\n', max(abs(d3))); nfail = nfail+1;
    else
        fprintf('CHECK 3: no component with a positive excess in this test, nothing to compare\n');
    end
    fprintf('         noise term of eq. (4) on this spectrum: %.0f to %.0f dB\n\n', min(10*log10(EGR)), max(10*log10(EGR)));

    test(k).f = f; test(k).L = L; test(k).idx = idx; test(k).fref = fref;
end

%% REPORT: what the level excess of [3] does hold
c = 0:0.5:30; emax = zeros(size(c));
for i = 1:numel(c)
    LXc = allL - 10*log10(allAEK + allEHS + 10^(c(i)/10));
    emax(i) = max(abs(LXc - allLXref));
end
[e, i] = min(emax);
fprintf('REPORT: with the noise term replaced by one constant of %.1f dB, the formula reproduces the %d level excesses of [3] within %.2f dB;\n', c(i), numel(allLXref), e);
fprintf('        the noise term of [3] does not follow its own eq. (4), and its level excesses hold no information about that term\n\n');

fprintf('%d checks failed\n', nfail);

%% plot
figure('color','w');
for k = 1:2
    subplot(2,1,k);
    plot(test(k).f, test(k).L, 'k-'); hold on;
    plot(test(k).comp(:,1), test(k).comp(:,2), 'd', 'MarkerSize', 9, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
    plot(test(k).fref, test(k).L(test(k).idx), 'ko', 'MarkerSize', 7);
    xlim([0 4400]);
    xlabel('Frequency (Hz)', 'Interpreter', 'Latex');
    ylabel('SPL (dB)', 'Interpreter', 'Latex');
    title(sprintf('Test %d of Zhang and Shrestha (2003), criterion at %d dB', k, test(k).thr), 'Interpreter', 'Latex');
    if k == 1
        legend({'Spectrum, Appendix C', 'Components, Table 6.2 and 6.3', 'SQAT'}, 'Location', 'NorthEast', 'Interpreter', 'Latex');
        legend boxoff;
    end
end

if save_figs==1
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    figname_short = 'tonality_validation_extraction_ZhangShrestha2003';
    figname_out = [figures_dir figname_short];
    saveas(gcf, figname_out, 'png');
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
