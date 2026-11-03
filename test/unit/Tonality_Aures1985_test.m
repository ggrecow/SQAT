function tests = Tonality_Aures1985_test
% Unit tests of Tonality_Aures1985, the tonality of Aures (1985).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_60_dB_at_1_kHz_is_1_tu(tc)
% Reference of the help: a pure tone of 1 kHz and 60 dB SPL has a tonality
% of 1 t.u. Tolerance 0.01 t.u., set by us.
[~, O] = evalc('Tonality_Aures1985(il_tone(60, 1000), 48000, 0, 0.5, false)');
sqat_report_record('tonality_aures1985', 'tonality of 1 kHz, 60 dB (1 t.u.)', O.Kmean, 1, 0.01, 'set by us; no published tolerance (measured 3e-5)');
tc.verifyEqual(O.Kmean, 1, 'AbsTol', 0.01);
end

function test_white_noise_has_no_tonality(tc)
% White noise holds no tonal component: the tonality stays below 0.01 t.u.
rng(1);
x = randn(4*48000, 1);
x = x / rms(x) * 2e-5 * 10^(60/20);
[~, O] = evalc('Tonality_Aures1985(x, 48000, 0, 0.5, false)');
tc.verifyLessThan(O.Kmean, 0.01);
end

function test_silence_has_no_tonality(tc)
% Four seconds of silence give a mean tonality of 0.
[~, O] = evalc('Tonality_Aures1985(zeros(4*48000, 1), 48000, 0, 0.5, false)');
tc.verifyEqual(O.Kmean, 0);
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
