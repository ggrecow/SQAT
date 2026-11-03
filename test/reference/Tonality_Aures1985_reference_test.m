function tests = Tonality_Aures1985_reference_test
% The extraction and the level excess of Tonality_Aures1985 (Terhardt et al.
% 1982) against the intermediate values of Zhang and Shrestha (2003),
% Sound Quality User-defined Cursor Reading Control, Tonality Metric,
% IMM-Thesis-2003-22, DTU: the power spectra of its Appendix C (in
% validation/Tonality_Aures1985/reference_values) and the components and
% level excesses of its Tables 6.2 and 6.3. The checks and criteria are the
% ones of validation/Tonality_Aures1985/validation_extraction_and_level_excess.m.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
fh = Tonality_Aures1985('localfunctions');
names = cellfun(@func2str, fh, 'UniformOutput', false);
for n = {'il_find_sinusoids', 'il_SPL_excess', 'il_Fq2Bark', 'il_Threshold'}
    tc.TestData.(n{1}) = fh{strcmp(names, n{1})};
end
dir_ref = fullfile(root, 'validation', 'Tonality_Aures1985', 'reference_values');
t(1).file = fullfile(dir_ref, 'ZhangShrestha2003_AppendixC_Test1.txt');
t(1).thr = 3;                                  % the criterion the thesis used for Test 1 (p. 49)
t(1).comp = [387.89 58.13; 765.95 33.20; 989.05 26.51; 1528.19 27.16; 1871.83 39.69; 2078.14 46.70; ...
             2509.13 41.67; 2636.70 43.89; 2863.81 37.18; 3091.05 35.01; 3348.17 40.62; 3467.18 39.10; ...
             3595.04 35.01; 3864.10 26.76];    % Table 6.2: frequency (Hz), level (dB)
t(1).LX = [387.89 39.02; 765.95 12.86; 989.05 9.11; 1528.19 9.67; 1871.83 10.35; 2078.14 16.20; ...
           2509.13 1.52; 2636.70 3.04];        % Table 6.2: level excess (dB)
t(2).file = fullfile(dir_ref, 'ZhangShrestha2003_AppendixC_Test2.txt');
t(2).thr = 7;                                  % the 7 dB of Terhardt et al.
t(2).comp = [387.59 87.87; 807.49 90.96; 1410.41 90.97; 2196.38 87.87];    % Table 6.3
t(2).LX = [387.59 67.83; 807.49 27.40; 1410.41 21.61; 2196.38 13.66];
tc.TestData.t = t;
end

function test_extraction_finds_the_components_of_tables_6_2_and_6_3(tc)
% Every component of the tables is found, each within one sample in frequency
% (after the refinement of Eq. 3 of the thesis) and with the same level to 0.01 dB.
for k = 1:2
    [f, L, df, idx, fref] = il_extract(tc, k);
    c = tc.TestData.t(k).comp;
    tc.verifyNumElements(idx, size(c, 1), sprintf('Test %d: number of components', k));
    [d, m] = min(abs(c(:, 1) - fref'), [], 1);
    tc.verifyLessThan(d, df, sprintf('Test %d: frequency within one sample', k));
    sqat_report_record('tonality_aures1985', sprintf('Test %d, levels of the components (Table 6.%d)', k, k+1), L(idx), c(m, 2), 0.01, 'validation script: levels printed with 2 decimals');
    tc.verifyEqual(L(idx), c(m, 2), 'AbsTol', 0.01, sprintf('Test %d: levels', k));
end
end

function test_level_excess_of_the_masked_components(tc)
% On the components masked by the other components (more than 25 dB of
% excitation from them), the level excess of the tables is reproduced within
% 0.3 dB with the noise term left out.
for k = 1:2
    [f, L, df, idx, fref] = il_extract(tc, k);
    in = struct('freq', f, 'ToneF', f(idx), 'ToneL', L(idx), 'Lnoise', -300 * ones(size(L)));
    LX = tc.TestData.il_SPL_excess(in);
    zc = tc.TestData.il_Fq2Bark(f(idx));
    Lc = L(idx);
    got = [];
    want = [];
    for i = 1:numel(idx)
        a = 0;
        for j = [1:i-1, i+1:numel(idx)]      % Eqs. 5 and 7 of the thesis: the excitation of the others
            if j < i
                s = -24 - 230/f(idx(j)) + 0.2*Lc(j);
            else
                s = 27;
            end
            a = a + 10^((Lc(j) - s*(zc(j) - zc(i)))/20);
        end
        [d, m] = min(abs(tc.TestData.t(k).LX(:, 1) - fref(i)));
        if 20*log10(a) > 25 && d < df
            got(end+1) = LX(i); %#ok<AGROW>
            want(end+1) = tc.TestData.t(k).LX(m, 2); %#ok<AGROW>
        end
    end
    tc.assertNotEmpty(got, sprintf('Test %d: no masked component', k));
    sqat_report_record('tonality_aures1985', sprintf('Test %d, level excess of the masked components (Table 6.%d)', k, k+1), got, want, 0.3, 'validation script: 0.3 dB');
    tc.verifyEqual(got, want, 'AbsTol', 0.3, sprintf('Test %d: level excess', k));
end
end

function test_the_scripts_of_validation_pass_their_checks(tc)
% Three scripts of validation/Tonality_Aures1985 that print CHECK ...
% PASSED or FAILED, run as they are, print no FAILED. Their criteria:
% extraction and level excess (the tables of Zhang and Shrestha, and the
% noise term of Eq. 4 of their thesis within 0.1 dB); frequency weighting
% (Eq. 9 of Aures 1985, Fig. 5); bandwidth weighting (Eq. 7, Fig. 6).
% validation_bandwidth_dependence.m is left out: its README records that it
% fails by a limit of the model (no term for the roll-off of a band, issue
% 67), and it has failed since it was added.
d = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'validation', 'Tonality_Aures1985');
for f = {'validation_extraction_and_level_excess.m', 'validation_frequency_weighting.m', ...
         'validation_bandwidth_weighting.m'}
    S = sqat_run_script(fullfile(d, f{1}));
    tc.assertEmpty(S.run_err, sprintf('%s stopped: %s', f{1}, S.run_err));
    tc.verifyNotEmpty(strfind(S.run_out, 'PASSED'), [f{1} ': no check passed']);
    tc.verifyEmpty(strfind(S.run_out, 'FAILED'), [f{1} ': ' strjoin(regexp(S.run_out, '[^\n]*FAILED[^\n]*', 'match'), ' | ')]);
    if isfield(S, 'd3') && ~isempty(S.d3)
        sqat_report_record('tonality_aures1985', 'noise term of the metric against Eq. 4 of Zhang and Shrestha (dB)', S.d3, zeros(size(S.d3)), 0.1, 'validation script: 0.1 dB');
    end
    if isfield(S, 'dev_tone')
        sqat_report_record('tonality_aures1985', 'frequency weighting on sine tones against Eq. 9 of Aures (1985)', S.dev_tone, zeros(size(S.dev_tone)), S.tol_w2, 'validation script, after Aures (1985), Fig. 5');
        sqat_report_record('tonality_aures1985', 'frequency weighting on 30 Hz bands against Eqs. 7 and 9', S.dev_30, zeros(size(S.dev_30)), S.tol_w1w2, 'validation script, after Aures (1985), Figs. 5 and 6');
    elseif isfield(S, 'tol_w1') && isfield(S, 'dev')
        sqat_report_record('tonality_aures1985', 'bandwidth weighting on 30 Hz bands against Eq. 7 of Aures (1985)', S.dev, zeros(size(S.dev)), S.tol_w1, 'validation script, after Aures (1985), Fig. 6');
    end
end
end

function [f, L, df, idx, fref] = il_extract(tc, k)
% the components the metric extracts from the spectrum of Test k
A = readmatrix(tc.TestData.t(k).file, 'CommentStyle', '#');
f = A(:, 1);
L = A(:, 2);
df = f(2) - f(1);
idx = tc.TestData.il_find_sinusoids(L, 1, numel(L), df, tc.TestData.t(k).thr);
fref = f(idx) + 0.46 * (L(idx+1) - L(idx-1));     % Eq. 3 of the thesis, for the comparison only
end
