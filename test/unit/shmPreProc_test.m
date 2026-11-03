function tests = shmPreProc_test
% Unit tests of shmPreProc (utilities/ECMA418_2), the fade-in and the
% zero-padding of ECMA-418-2:2025, clause 5.1.2.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_fade_in_and_the_padding_of_the_standard(tc)
% The first 240 samples are weighted by 0.5 - 0.5 cos(pi n / 240) (Formula 1).
% For tonality and loudness (5.1.2.1) the start gets s_b,max zeros and the
% end n_new - n_samples zeros, with
% n_new = s_h,max (ceil((n_samples + s_h,max + s_b,max)/s_h,max) - 1) (Formulas 2, 3).
n = 10000; sb = 8192; sh = 2048;
y = shmPreProc(ones(n, 1), sb, sh, true, true);
n_new = sh * (ceil((n + sh + sb)/sh) - 1);
tc.verifyEqual(numel(y), sb + n + (n_new - n));
tc.verifyEqual(y(1:sb), zeros(sb, 1));
tc.verifyEqual(y(sb+1:sb+240), 0.5 - 0.5 * cos(pi * (0:239)' / 240), 'AbsTol', 1e-15);
tc.verifyEqual(y(sb+241:sb+n), ones(n - 240, 1));
tc.verifyEqual(y(sb+n+1:end), zeros(n_new - n, 1));
end

function test_no_padding_when_asked(tc)
% Roughness and fluctuation strength pad only the start (5.1.2.2); with both
% pads off only the fade-in is applied.
y = shmPreProc(ones(1000, 1), 16384, 4096, false, false);
tc.verifySize(y, [1000 1]);
end
