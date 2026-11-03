function tests = PsychoacousticAnnoyance_Widmann1992_test
% Unit tests of PsychoacousticAnnoyance_Widmann1992 and of
% PsychoacousticAnnoyance_Zwicker1999, the same model under another name.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_40_dB_at_1_kHz_is_1_au(tc)
% Widmann defined a 1 kHz tone of 40 dB SPL as the reference signal of his
% model, with an annoyance of 1 au (help of the function). Tolerance 0.01 au,
% set by us (measured 1.003 on 5 s).
[~, O] = evalc('PsychoacousticAnnoyance_Widmann1992(il_tone(40), 48000, 0, 0.5, false, false)');
sqat_report_record('annoyance', 'annoyance of 1 kHz, 40 dB (1 au)', O.PAmean, 1, 0.01, 'set by us; no published tolerance (measured 0.003)');
tc.verifyEqual(O.PAmean, 1, 'AbsTol', 0.01);
end

function test_zwicker1999_is_the_same_model(tc)
% The help: the model of Zwicker and Fastl (1999) is the one of Widmann
% (Lotinga and Torija 2025), so both give the same result.
x = il_tone(60);
[~, W] = evalc('PsychoacousticAnnoyance_Widmann1992(x, 48000, 0, 0.5, false, false)');
[~, Z] = evalc('PsychoacousticAnnoyance_Zwicker1999(x, 48000, 0, 0.5, false, false)');
tc.verifyEqual(Z.PAmean, W.PAmean, 'RelTol', 1e-12);
end

function test_silence_has_no_annoyance(tc)
% Four seconds of silence give a mean psychoacoustic annoyance of 0.
[~, O] = evalc('PsychoacousticAnnoyance_Widmann1992(zeros(4*48000, 1), 48000, 0, 0.5, false, false)');
tc.verifyEqual(O.PAmean, 0);
end

function x = il_tone(L)
% 4 s of a 1 kHz sinusoid of L dB SPL (rms), fs = 48 kHz
t = (0:4*48000-1)' / 48000;
x = sqrt(2) * 2e-5 * 10^(L/20) * sin(2*pi*1000*t);
end
