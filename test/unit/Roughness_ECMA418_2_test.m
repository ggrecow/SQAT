function tests = Roughness_ECMA418_2_test
% Unit tests of Roughness_ECMA418_2, the roughness of ECMA-418-2:2025 clause 7.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_reference_signal_is_1_asper(tc)
% Reference signal of the help: a 1 kHz tone of 60 dB SPL, 100 % modulated
% at 70 Hz, yields 1 asper (R90). Tolerance 0.0025, as the loudness.
[~, O] = evalc('Roughness_ECMA418_2(il_am(60, 70), 48000, ''free-frontal'', 0.304, false)');
sqat_report_record('roughness_ecma418_2', 'roughness of the reference signal (1 asper)', O.roughness90Pc, 1, 0.0025, 'set by us, as the 0.25 % of ECMA-418-2 footnote 9');
tc.verifyEqual(O.roughness90Pc, 1, 'AbsTol', 0.0025);
end

function test_roughness_peaks_at_70_Hz(tc)
% Roughness of a modulated 1 kHz tone is largest near 70 Hz of modulation
% (Fastl and Zwicker 2007, Fig. 11.1): 70 Hz beats 20, 40, 100 and 150 Hz.
fm = [20 40 70 100 150];
R = zeros(size(fm));
for k = 1:numel(fm)
    [~, O] = evalc('Roughness_ECMA418_2(il_am(60, fm(k)), 48000, ''free-frontal'', 0.304, false)');
    R(k) = O.roughness90Pc;
end
[~, i] = max(R);
tc.verifyEqual(fm(i), 70);
end

function test_silence_has_no_roughness(tc)
% Four seconds of silence give a 90th percentile roughness of 0.
[~, O] = evalc('Roughness_ECMA418_2(zeros(4*48000, 1), 48000, ''free-frontal'', 0.304, false)');
tc.verifyEqual(O.roughness90Pc, 0);
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
