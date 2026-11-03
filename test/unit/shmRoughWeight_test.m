function tests = shmRoughWeight_test
% Unit tests of shmRoughWeight (utilities/ECMA418_2), the weighting of the
% modulation rates of the roughness, ECMA-418-2:2025 Formula 85.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_weight_is_one_at_fmax_and_follows_formula_85(tc)
% G = 1 / (1 + ((f/fmax - fmax/f) q1)^2)^q2 reaches its maximum of one at
% f = fmax. At 1 kHz, fmax = 72.6937 (1 - 1.1739 exp(-5.4583)) Hz
% (Formula 86), q1 = 1.2822 and q2 = 0.2471 + 0.0129 (0 + 3.4253)^2
% (Formula 87).
fmax = 72.6937 * (1 - 1.1739 * exp(-5.4583));
q = [1.2822; 0.2471 + 0.0129 * 3.4253^2];
f = [fmax/4 fmax/2 fmax 2*fmax 4*fmax];
want = 1 ./ (1 + ((f/fmax - fmax./f) * q(1)).^2).^q(2);
got = shmRoughWeight(f, fmax, q);
sqat_report_record('shm', 'roughness weighting at 1 kHz against Formula 85', got, want, 1e-12, 'numerical: exact Formula 85 of ECMA-418-2');
tc.verifyEqual(got, want, 'AbsTol', 1e-12);
tc.verifyEqual(got(3), 1);
end

function test_the_weight_is_symmetric_on_a_log_scale(tc)
% Formula 85 depends on f/fmax - fmax/f, so a rate a times fmax and a rate
% fmax/a get the same weight.
fmax = 70;
q = [1.2822; 0.3];
a = [1.5 2 3 5];
tc.verifyEqual(shmRoughWeight(a*fmax, fmax, q), shmRoughWeight(fmax./a, fmax, q), 'AbsTol', 1e-12);
end
