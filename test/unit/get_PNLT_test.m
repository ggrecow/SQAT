function tests = get_PNLT_test
% Unit tests of get_PNLT (EPNL_FAR_Part36/helper), the tone correction of
% 14 CFR Part 36 Appendix A, Section A36.4.3 and Table A36-2.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.fb = [50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 ...
    1600 2000 2500 3150 4000 5000 6300 8000 10000];
end

function test_a_band_above_a_flat_background_gets_table_A36_2(tc)
% A flat background of 60 dB with one band raised by F dB: the ten steps find
% the difference F, and Table A36-2 gives C = F/3 (500 to 5000 Hz,
% 3 <= F < 20), 20/3 (F >= 20), F/6 (below 500 Hz or above 5 kHz, 3 <= F < 20)
% and 10/3 (F >= 20). PNLT = PNL + C for a single time step.
cases = [14 12 4; 14 25 20/3; 10 12 2; 10 25 10/3; 22 12 2];   % band, F, C
got = zeros(1, size(cases, 1));
for k = 1:size(cases, 1)
    got(k) = il_correction(tc, cases(k, 1), cases(k, 2));
end
sqat_report_record('epnl', 'tone correction of Table A36-2 (5 cases)', got, cases(:, 3)', 1e-9);
tc.verifyEqual(got, cases(:, 3)', 'AbsTol', 1e-9, 'numerical: Table A36-2 of 14 CFR Part 36');
end

function test_a_smooth_spectrum_gets_no_correction(tc)
% No change of slope above 5 dB (step 2): a flat or a straight sloped
% spectrum, and a band only 2 dB above a flat background, get C = 0.
tc.verifyEqual(il_correction(tc, 14, 0), 0);
tc.verifyEqual(il_correction(tc, 14, 2), 0);
L = 80 - 0.5 * (1:24);
[~, PNL] = get_PNL(L);
tc.verifyEqual(get_PNLT(L, tc.TestData.fb, PNL), PNL);
end

function C = il_correction(tc, band, F)
% PNLT - PNL of a 60 dB flat spectrum with one band F dB above it
L = 60 * ones(1, 24);
L(band) = 60 + F;
[~, PNL] = get_PNL(L);
C = get_PNLT(L, tc.TestData.fb, PNL) - PNL;
end
