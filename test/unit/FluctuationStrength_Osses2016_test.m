function tests = FluctuationStrength_Osses2016_test
% Unit tests of FluctuationStrength_Osses2016, the fluctuation strength of
% Osses et al. (2016).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_reference_signal_is_1_vacil(tc)
% Reference signal of the help: a 1 kHz tone of 60 dB SPL, 100 % modulated
% at 4 Hz, should yield 1 vacil. Tolerance 0.02 vacil, set by us (measured
% 1.006 on 5 s).
[~, O] = evalc('FluctuationStrength_Osses2016(il_am(60, 4), 48000, 1, 0.5, false)');
sqat_report_record('fluctuation_strength', 'fluctuation strength of the reference signal (1 vacil)', O.FSmean, 1, 0.02, 'set by us; no published tolerance (measured 0.006)');
tc.verifyEqual(O.FSmean, 1, 'AbsTol', 0.02);
end

function test_fluctuation_strength_peaks_at_4_Hz(tc)
% Fastl and Zwicker (2007), Fig. 10.1: largest near 4 Hz of modulation;
% 4 Hz beats 1, 2, 8 and 16 Hz (8 Hz is close, 0.99 against 1.01).
fm = [1 2 4 8 16];
F = zeros(size(fm));
for k = 1:numel(fm)
    [~, O] = evalc('FluctuationStrength_Osses2016(il_am(60, fm(k)), 48000, 1, 0.5, false)');
    F(k) = O.FSmean;
end
[~, i] = max(F);
tc.verifyEqual(fm(i), 4);
end

function test_silence_has_no_fluctuation(tc)
% Four seconds of silence give a mean fluctuation strength of 0.
[~, O] = evalc('FluctuationStrength_Osses2016(zeros(4*48000, 1), 48000, 1, 0.5, false)');
tc.verifyEqual(O.FSmean, 0);
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
