function tests = EPNL_FAR_Part36_reference_test
% The tone correction of EPNL_FAR_Part36 (get_PNLT) against the worked
% example of ICAO Doc 9501, Environmental Technical Manual, Volume I (2015),
% Table 3.7 (a turbofan engine): the value of every step, band by band, as
% transcribed in validation/EPNL_FAR_Part36/2_Tone_Correction_Factor.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.T = load(fullfile(root, 'validation', 'EPNL_FAR_Part36', '2_Tone_Correction_Factor', ...
    'private', 'TONE_CORRECTION_TABLE_REF_VALUES.m'));
end

function test_every_step_of_table_3_7(tc)
% Columns of the table: band, f, SPL, S (step 1), |dS| (2), SPL' (4), S' (5),
% S-bar (6), SPL'' (7), F (8), C (9). The table prints 2 to 6 decimals, so
% steps 1 to 8 are compared within 0.01 dB where the table gives a finite
% value. Step 9 is compared with Table A36-2 of 14 CFR Part 36 applied to the
% F of the table: the C printed by the example differs from that table at
% 160 Hz (0.29 for F/3 - 1/2 = 0.278) and at 250 Hz (0.61 for F/6 = 0.667).
T = tc.TestData.T;
[~, ~, ~, O] = get_PNLT(T(:, 3)', T(:, 2)', 0);
steps = {'S', 4, 'step 1, S'; 'SPLP', 6, 'step 4, SPL prime'; 'SP', 7, 'step 5, S prime'; ...
         'SB', 8, 'step 6, S bar'; 'SPLPP', 9, 'step 7, SPL double prime'; 'F', 10, 'step 8, F'};
for k = 1:size(steps, 1)
    got = O.(steps{k, 1})(1:24);
    want = T(:, steps{k, 2});
    ok = isfinite(want);
    sqat_report_record('epnl', ['ICAO Doc 9501 Table 3.7, ' steps{k, 3}], got(ok), want(ok), 0.01, 'ICAO Doc 9501 prints 2 to 6 decimals');
    tc.verifyEqual(got(ok), want(ok), 'AbsTol', 0.01, steps{k, 3});
end
F = T(:, 10);
C = zeros(24, 1);                                   % Table A36-2
low = T(:, 2) < 500 | T(:, 2) > 5000;
C(F >= 1.5 & F < 3) = il_if(low(F >= 1.5 & F < 3), F(F >= 1.5 & F < 3)/3 - 1/2, 2*F(F >= 1.5 & F < 3)/3 - 1);
C(F >= 3 & F < 20) = il_if(low(F >= 3 & F < 20), F(F >= 3 & F < 20)/6, F(F >= 3 & F < 20)/3);
C(F >= 20) = il_if(low(F >= 20), 10/3, 20/3);
sqat_report_record('epnl', 'step 9, C, against Table A36-2 on the F of Table 3.7', O.C(1:24), C, 0.01, 'numerical: Table A36-2 of 14 CFR Part 36');
tc.verifyEqual(O.C(1:24), C, 'AbsTol', 0.01, 'step 9, C');
end

function v = il_if(cond, a, b)
% element-wise choice between a and b (b may be a scalar)
if isscalar(b)
    b = b * ones(size(cond));
end
if isscalar(a)
    a = a * ones(size(cond));
end
v = b;
v(cond) = a(cond);
end
