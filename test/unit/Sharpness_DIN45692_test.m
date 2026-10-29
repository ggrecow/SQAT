function tests = Sharpness_DIN45692_test
% Unit tests of Sharpness_DIN45692, the sharpness of DIN 45692:2009.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_reference_noise_is_1_acum(tc)
% DIN 45692: narrow-band noise one critical band wide (920 to 1080 Hz) at
% 1 kHz and 60 dB SPL has a sharpness of 1 acum. The standard lets k vary
% from 0.105 to 0.115 around the 0.11 of the function, about 4.5 %, which is
% the tolerance here. The noise is random, so the seed is fixed.
rng(1);
n = 4*48000;
X = fft(randn(n, 1));
f = (0:n-1)' * 48000 / n;
X(~((f >= 920 & f <= 1080) | (f >= 48000-1080 & f <= 48000-920))) = 0;
x = real(ifft(X));
x = x / rms(x) * 2e-5 * 10^(60/20);
[~, O] = evalc('Sharpness_DIN45692(x, 48000, ''DIN45692'', 0, 2, 0.5, false, false)');
sqat_report_record('sharpness_din45692', 'sharpness of the reference noise (1 acum)', O.Smean, 1, 0.045, 'DIN 45692:2009, range of k from 0.105 to 0.115');
tc.verifyEqual(O.Smean, 1, 'AbsTol', 0.045);
end

function test_silence_gives_nan(tc)
% With no loudness the sharpness is 0/0: the function returns NaN (the
% SQAT Rust port records the same behaviour).
[~, O] = evalc('Sharpness_DIN45692(zeros(4*48000, 1), 48000, ''DIN45692'', 0, 2, 0.5, false, false)');
tc.verifyTrue(isnan(O.Smean));
end
