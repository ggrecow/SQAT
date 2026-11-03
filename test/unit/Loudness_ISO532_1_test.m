function tests = Loudness_ISO532_1_test
% Unit tests of Loudness_ISO532_1, the loudness of ISO 532-1:2017.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_40_dB_at_1_kHz_is_1_sone(tc)
% Reference signal of the help: a 1 kHz tone of 40 dB SPL yields 1 sone
% (stationary method, free field). Tolerance: the criterion of ISO 532-1 for
% the total loudness of its test signals, 5 % or 0.1 sone, whichever is larger.
[~, O] = evalc('Loudness_ISO532_1(il_tone(40, 1000), 48000, 0, 1, 0.5, false)');
sqat_report_record('loudness_iso532_1', 'loudness of 1 kHz, 40 dB (1 sone)', O.Loudness, 1, 0.1, 'ISO 532-1:2017, 5 % or 0.1 sone');
tc.verifyEqual(O.Loudness, 1, 'AbsTol', 0.1);
end

function test_the_loudness_of_a_1kHz_tone_doubles_every_10_dB(tc)
% The phon is defined on the 1 kHz tone: at L dB SPL it has L phon, so
% 2^((L - 40)/10) sone from 40 dB up. Tolerance: the criterion of ISO 532-1
% for total loudness, 5 % or 0.1 sone, whichever is larger (measured 1.6 %
% at 80 dB).
L = [40 60 80];
N = zeros(size(L));
for k = 1:numel(L)
    [~, O] = evalc('Loudness_ISO532_1(il_tone(L(k), 1000), 48000, 0, 1, 0.5, false)');
    N(k) = O.Loudness;
end
want = 2.^((L - 40)/10);
tol = max(0.05 * want, 0.1);
sqat_report_record('loudness_iso532_1', 'loudness of 1 kHz tones, 40 to 80 dB, over the ISO tolerance', (N - want)./tol, zeros(size(L)), 1, 'ISO 532-1:2017, 5 % or 0.1 sone');
tc.verifyLessThanOrEqual(abs(N - want), tol);
end

function test_silence_has_no_loudness(tc)
% Four seconds of silence give a total loudness of 0.
[~, O] = evalc('Loudness_ISO532_1(zeros(4*48000, 1), 48000, 0, 1, 0.5, false)');
tc.verifyEqual(O.Loudness, 0);
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
