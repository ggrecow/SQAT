function tests = calculate_a0_test
% Unit tests of calculate_a0 (utilities), the transmission of the outer and
% middle ear (free field) used by the roughness and the fluctuation strength.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_fastl2007_follows_the_curve_of_the_book(tc)
% Reference: Fastl and Zwicker (2007), Fig. 8.18, as tabulated in
% calculate_a0: 0 dB at 1 kHz (8.5 Bark), the ear canal resonance of
% +7.38 dB at 3.4 kHz (16.5 Bark), 0 dB at 5.3 kHz (19 Bark), -2.59 dB at
% 7.7 kHz (21 Bark), -11.3 dB at 12 kHz (23 Bark), -40 dB at 15.5 kHz (24 Bark).
% The frequencies are table points of Get_Bark, so no interpolation enters.
[f, d] = il_a0_dB('fastl2007');
pick = [1000 3400 5300 7700 12000 15500];
want = [0 7.38 0 -2.59 -11.3 -40];
sqat_report_record('a0', 'fastl2007 at the points of Fig. 8.18', d(ismember(f, pick)), want, 1e-9, 'numerical: table points of Fastl and Zwicker (2007), Fig. 8.18');
tc.verifyEqual(d(ismember(f, pick)), want, 'AbsTol', 1e-9);
end

function test_osses2016_removes_the_ear_canal_resonance(tc)
% The simplified curve of Osses et al. (2016) is 0 dB up to 19 Bark (5.3 kHz)
% and equal to fastl2007 above it.
[f, d_osses] = il_a0_dB('fluctuationstrength_osses2016');
[~, d_fastl] = il_a0_dB('fastl2007');
tc.verifyEqual(d_osses(f <= 5300), zeros(1, nnz(f <= 5300)), 'AbsTol', 1e-9);
tc.verifyEqual(d_osses(f >= 5300), d_fastl(f >= 5300), 'AbsTol', 1e-9);
end

function test_the_frequencies_span_20_Hz_to_20_kHz(tc)
% calculate_a0 returns the bins from 20 Hz to 20 kHz, spaced fs/N.
[f, ~] = il_a0_dB('fastl2007');
tc.verifyEqual([f(1) f(end)], [20 20000], 'AbsTol', 1e-9);
tc.verifyEqual(diff(f), 10 * ones(1, numel(f) - 1), 'AbsTol', 1e-9);
end

function [f, d] = il_a0_dB(type)
% the curve in dB at fs = 48 kHz, N = 4800: 10 Hz bins, so the table points fall on bins
[~, f, a0] = calculate_a0(48000, 4800, type);
d = 20 * log10(abs(a0));
end
