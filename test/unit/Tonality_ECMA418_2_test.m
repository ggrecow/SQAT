function tests = Tonality_ECMA418_2_test
% Unit tests of Tonality_ECMA418_2, the tonality of ECMA-418-2:2025 clause 6.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_40_dB_at_1_kHz_is_1_tu_HMS(tc)
% Reference signal of the help: a 1 kHz tone of 40 dB SPL yields 1 tu_HMS.
% Tolerance 0.0025, as the 0.25 % that footnote 9 allows the loudness.
[~, O] = evalc('Tonality_ECMA418_2(il_tone(40, 1000), 48000, ''free-frontal'', 0.304, false)');
sqat_report_record('tonality_ecma418_2', 'tonality of 1 kHz, 40 dB (1 tu_HMS)', O.tonalityAvg, 1, 0.0025, 'set by us, as the 0.25 % of ECMA-418-2 footnote 9');
tc.verifyEqual(O.tonalityAvg, 1, 'AbsTol', 0.0025);
end

function test_silence_has_no_tonality(tc)
% Four seconds of silence give an average tonality of 0.
[~, O] = evalc('Tonality_ECMA418_2(zeros(4*48000, 1), 48000, ''free-frontal'', 0.304, false)');
tc.verifyEqual(O.tonalityAvg, 0);
end

function x = il_tone(L, f)
% 4 s of a sinusoid of L dB SPL (rms) at f Hz, fs = 48 kHz
t = (0:4*48000-1)' / 48000;
x = sqrt(2) * 2e-5 * 10^(L/20) * sin(2*pi*f*t);
end

function x = il_am(L, fm)
% 4 s of a 1 kHz tone 100 % amplitude-modulated at fm, L dB SPL as the rms of the modulated signal
t = (0:4*48000-1)' / 48000;
x = sqrt(2) * 2e-5 * 10^(L/20) * (1 + sin(2*pi*fm*t)) .* sin(2*pi*1000*t) / sqrt(1.5);
end
