function tests = Roughness_Daniel1997_test
% Unit tests of Roughness_Daniel1997, the roughness of Daniel and Weber (1997).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_reference_signal_is_1_asper(tc)
% Reference signal of the help: a 1 kHz tone of 60 dB SPL, 100 % modulated
% at 70 Hz, yields 1 asper. Tolerance 0.01 asper, set by us.
[~, O] = evalc('Roughness_Daniel1997(il_am(60, 70), 48000, 0.5, false)');
sqat_report_record('roughness_daniel1997', 'roughness of the reference signal (1 asper)', O.Rmean, 1, 0.01, 'set by us; no published tolerance (measured 0.0005)');
tc.verifyEqual(O.Rmean, 1, 'AbsTol', 0.01);
end

function test_roughness_peaks_at_70_Hz(tc)
% Fastl and Zwicker (2007), Fig. 11.1: 70 Hz beats 20, 40, 100 and 150 Hz.
fm = [20 40 70 100 150];
R = zeros(size(fm));
for k = 1:numel(fm)
    [~, O] = evalc('Roughness_Daniel1997(il_am(60, fm(k)), 48000, 0.5, false)');
    R(k) = O.Rmean;
end
[~, i] = max(R);
tc.verifyEqual(fm(i), 70);
end

function test_silence_has_no_roughness(tc)
% Four seconds of silence give a mean roughness of 0.
[~, O] = evalc('Roughness_Daniel1997(zeros(4*48000, 1), 48000, 0.5, false)');
tc.verifyEqual(O.Rmean, 0);
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
