function tests = Loudness_ECMA418_2_test
% Unit tests of Loudness_ECMA418_2, the loudness of ECMA-418-2:2025 clause 8.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_40_dB_at_1_kHz_is_1_sone_HMS(tc)
% 5.1.8: cN is set so that a 1 kHz sinusoid of 40 dB has a total loudness of
% 1 sone_HMS; footnote 9 allows cN a tolerance of 0.25 %.
[~, O] = evalc('Loudness_ECMA418_2(il_tone(40, 1000), 48000, ''free-frontal'', 0.304, false)');
sqat_report_record('loudness_ecma418_2', 'loudness of 1 kHz, 40 dB (1 sone_HMS)', O.loudnessPowAvg, 1, 0.0025, 'ECMA-418-2:2025, footnote 9 (0.25 % on cN)');
tc.verifyEqual(O.loudnessPowAvg, 1, 'AbsTol', 0.0025);
end

function test_silence_has_no_loudness(tc)
% Four seconds of silence give an overall loudness of 0.
[~, O] = evalc('Loudness_ECMA418_2(zeros(4*48000, 1), 48000, ''free-frontal'', 0.304, false)');
tc.verifyEqual(O.loudnessPowAvg, 0);
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
