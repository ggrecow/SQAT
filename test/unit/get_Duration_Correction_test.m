function tests = get_Duration_Correction_test
% Unit tests of get_Duration_Correction (EPNL_FAR_Part36/helper), the
% duration correction of 14 CFR Part 36 Appendix A, Section A36.4.5:
% D = 10 log10(sum(10^(PNLT/10)) dt / T0) - PNLTM, with T0 = 10 s.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_a_plateau_of_10_s_gives_no_correction(tc)
% PNLT at 90 TPNdB for 10 s (20 steps of 0.5 s) between 70 TPNdB: D = 0,
% plus 0.002 dB from the first step below the 10 dB down point, which the
% sum includes; tolerance 0.01 dB.
[D, t1, t2] = il_plateau(20);
sqat_report_record('epnl', 'duration correction of a 10 s plateau', D, 0, 0.01, 'set by us: the step below the down point adds 0.002 dB');
tc.verifyEqual(D, 0, 'AbsTol', 0.01);
tc.verifyEqual([t1 t2], [11 31]);
end

function test_twice_the_duration_adds_3_dB(tc)
% Doubling the time above the down points doubles the energy:
% D rises by 10 log10(2) = 3.01 dB; tolerance 0.005 dB.
D = il_plateau(40) - il_plateau(20);
sqat_report_record('epnl', 'duration correction, 20 s against 10 s', D, 10*log10(2), 0.005, 'set by us (measured 0.001)');
tc.verifyEqual(D, 10*log10(2), 'AbsTol', 0.005);
end

function test_no_decay_warns_and_uses_the_end(tc)
% When PNLT never falls 10 dB below PNLTM after the peak, the function warns
% and integrates up to the last step.
PNLT = [70 * ones(10, 1); 90 * ones(20, 1)];
[m, i] = max(PNLT);
lastwarn('');                                             % the warning has no identifier
[~, ~, t2] = get_Duration_Correction(PNLT, m, i, 0.5, 10);
tc.verifySubstring(lastwarn, 'does not decay by more than the threshold');
tc.verifyEqual(t2, numel(PNLT));
end

function [D, t1, t2] = il_plateau(n)
% n steps of 0.5 s at 90 TPNdB, with 10 steps at 70 TPNdB on each side
PNLT = [70 * ones(10, 1); 90 * ones(n, 1); 70 * ones(10, 1)];
[m, i] = max(PNLT);
[D, t1, t2] = get_Duration_Correction(PNLT, m, i, 0.5, 10);
end
