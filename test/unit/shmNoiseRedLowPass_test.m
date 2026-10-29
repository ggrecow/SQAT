function tests = shmNoiseRedLowPass_test
% Unit tests of shmNoiseRedLowPass (utilities/ECMA418_2), the low-pass filter
% of the noise reduction of the tonality, ECMA-418-2:2025 clause 6.2.7.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_impulse_response_is_formula_11(tc)
% Formula 11 with order k = 3, e = [0 1 1] (footnote 21) and
% tau = (1/32)(6/7) s = 0.0268 s (footnote 20):
% h(n) = (1 - d)^3 / (d + d^2) n^2 d^n, d = exp(-1/(rs tau)), at the block
% rate of the tonality, rs = 48000/256 = 187.5 Hz.
rs = 187.5;
d = exp(-1 / (rs * (1/32) * (6/7)));
n = (0:199)';
h = (1 - d)^3 / (d + d^2) * n.^2 .* d.^n;
got = shmNoiseRedLowPass([1; zeros(199, 1)], rs);
sqat_report_record('shm', 'noise reduction low-pass, impulse response against Formula 11', got, h, 1e-12, 'numerical: exact Formula 11 of ECMA-418-2');
tc.verifyEqual(got, h, 'AbsTol', 1e-12);
end

function test_a_constant_passes_with_unit_gain(tc)
% The factor of Formula 11 makes the sum of h equal to 1: a constant input
% settles at the same constant.
y = shmNoiseRedLowPass(ones(400, 1), 187.5);
tc.verifyEqual(y(end), 1, 'AbsTol', 1e-9);
end
