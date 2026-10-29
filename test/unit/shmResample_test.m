function tests = shmResample_test
% Unit tests of shmResample (utilities/ECMA418_2): the hearing model runs at
% 48 kHz, and any other rate is resampled to it (ECMA-418-2:2025, footnote 2).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_a_44_1_kHz_signal_goes_to_48_kHz(tc)
% One second at 44.1 kHz becomes 48000 samples at 48 kHz, and the rms of a
% sinusoid of amplitude 1 stays 1/sqrt(2) away from the edges (within 0.1 %).
[y, rate] = shmResample(sin(2*pi*1000*(0:44099)'/44100), 44100);
tc.verifyEqual(rate, 48000);
tc.verifySize(y, [48000 1]);
sqat_report_record('shm', 'rms of a resampled sinusoid', rms(y(1000:end-1000)), 1/sqrt(2), 1e-3/sqrt(2), 'set by us (measured 0.05 %)');
tc.verifyEqual(rms(y(1000:end-1000)), 1/sqrt(2), 'RelTol', 1e-3);
end

function test_48_kHz_passes_unchanged(tc)
% A signal already at 48 kHz, the rate of the ECMA-418-2 hearing model, comes
% back unchanged, with the rate 48000.
x = randn(4800, 1);
[y, rate] = shmResample(x, 48000);
tc.verifyEqual([y; rate], [x; 48000]);
end
