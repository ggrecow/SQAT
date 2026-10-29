function tests = shmRoughLowPass_test
% Unit tests of shmRoughLowPass (utilities/ECMA418_2), the smoothing of the
% specific roughness, ECMA-418-2:2025 clause 7.1.7, Formulas 109 and 110.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_it_rises_fast_and_falls_slowly(tc)
% A step from 0 to 1 and back to 0 at the rate of 50 Hz of the roughness.
% Formula 109 is a first-order filter, so the output rises as
% 1 - exp(-l/(50 tau_rise)) with tau_rise = 0.0625 s and falls as
% exp(-l/(50 tau_fall)) with tau_fall = 0.5 s (Formula 110).
rs = 50;
x = [0; ones(50, 1); zeros(100, 1)];
y = shmRoughLowPass(x, rs, 0.0625, 0.5);
l = (1:50)';
rise = 1 - exp(-l / (rs * 0.0625));
fall = rise(end) * exp(-(1:100)' / (rs * 0.5));
sqat_report_record('shm', 'roughness smoothing, rise and fall against Formula 109', y(2:end), [rise; fall], 1e-12, 'numerical: exact Formula 109 of ECMA-418-2');
tc.verifyEqual(y(2:51), rise, 'AbsTol', 1e-12);
tc.verifyEqual(y(52:end), fall, 'AbsTol', 1e-12);
end

function test_each_band_is_filtered_on_its_own(tc)
% Every column (critical band) runs its own filter.
x = [zeros(1, 3); ones(20, 3) .* [1 2 3]];
y = shmRoughLowPass(x, 50, 0.0625, 0.5);
tc.verifyEqual(y(:, 2), 2 * y(:, 1), 'AbsTol', 1e-12);
tc.verifyEqual(y(:, 3), 3 * y(:, 1), 'AbsTol', 1e-12);
end
