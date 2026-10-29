function tests = PsychoacousticAnnoyance_Di2016_test
% Unit tests of PsychoacousticAnnoyance_Di2016. Its help states no reference
% signal, so the tests check that its two entry points agree and that
% silence has no annoyance.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_scalar_follows_the_percentiles(tc)
% ScalarPA is the formula of the model on N5, S5, R5, FS5 and K5, the same
% formula PsychoacousticAnnoyance_Di2016_from_percentile applies: on a rough
% 1 kHz tone of 70 dB both give the same value.
t = (0:4*48000-1)' / 48000;
x = sqrt(2) * 2e-5 * 10^(70/20) * (1 + sin(2*pi*70*t)) .* sin(2*pi*1000*t) / sqrt(1.5);
[~, O] = evalc('PsychoacousticAnnoyance_Di2016(x, 48000, 0, 0.5, false, false)');
[~, P] = evalc('PsychoacousticAnnoyance_Di2016_from_percentile(O.L.N5, O.S.S5, O.R.R5, O.FS.FS5, O.K.K5)');
tc.verifyGreaterThan(O.ScalarPA, 0);
tc.verifyEqual(P, O.ScalarPA, 'RelTol', 1e-12);
end

function test_silence_has_no_annoyance(tc)
% Four seconds of silence give a mean psychoacoustic annoyance of 0.
[~, O] = evalc('PsychoacousticAnnoyance_Di2016(zeros(4*48000, 1), 48000, 0, 0.5, false, false)');
tc.verifyEqual(O.PAmean, 0);
end
