function tests = Terhardt_filterbank_test
% Unit tests of Terhardt_filterbank and Terhardt_filterbank_params
% (utilities), the critical-band filterbank of the roughness of Daniel and
% Weber (1997) and of the fluctuation strength of Osses et al. (2016).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
tc.TestData.fs = 48000;
tc.TestData.N = 4800;                                    % 10 Hz bins
tc.TestData.P = Terhardt_filterbank_params(tc.TestData.N, tc.TestData.fs);
end

function test_the_params_cover_47_channels_from_20_Hz_to_20_kHz(tc)
% 47 half-Bark channels; the audible range of the bins is 20 Hz to 20 kHz.
P = tc.TestData.P;
tc.verifyEqual(P.Chno, 47);
tc.verifyEqual([P.freqs(1) P.freqs(end)], [20 20000], 'AbsTol', 1e-9);
end

function test_a_1kHz_tone_spreads_with_the_slopes_of_terhardt(tc)
% A 1 kHz tone at 70 dB SPL, one component. Terhardt (1979): the excitation
% falls by 27 dB/Bark towards low frequencies (Eq. 4a) and by
% 24 + 230/f - 0.2 L = 10.23 dB/Bark towards high frequencies (Eq. 4b); the
% channels are half a Bark apart. The channels around the tone keep its level.
[e, info] = il_pattern(tc, 1000, 70);
tc.verifyEqual(info.n_components, 1);
top = find(abs(e - 70) < 1e-9);
tc.assertNotEmpty(top, 'no channel at the level of the tone');
low = diff(e(top(1)-4:top(1)));
high = diff(e(top(end):top(end)+4));
sqat_report_record('terhardt', 'lower slope per half Bark, Eq. 4a', low, 13.5 * ones(1, 4), 1e-9, 'numerical: exact Eq. 4a of Terhardt (1979)');
sqat_report_record('terhardt', 'upper slope per half Bark, Eq. 4b', high, -(24 + 230/1000 - 0.2*70)/2 * ones(1, 4), 1e-9, 'numerical: exact Eq. 4b of Terhardt (1979)');
tc.verifyEqual(low, 13.5 * ones(1, 4), 'AbsTol', 1e-9);
tc.verifyEqual(high, -(24 + 230/1000 - 0.2*70)/2 * ones(1, 4), 'AbsTol', 1e-9);
end

function test_silence_gives_no_component(tc)
% Nothing above the hearing threshold: no component and a silent output.
[ei, info] = Terhardt_filterbank(zeros(1, tc.TestData.N), tc.TestData.P);
tc.verifyEqual(info.n_components, 0);
tc.verifyEqual(ei, zeros(tc.TestData.P.Chno, tc.TestData.N));
end

function test_a_level_above_121_dB_is_reported_as_clamped(tc)
% Above 120 + 1150/f dB the upper slope of Eq. 4b would be positive; the
% function keeps it at zero and reports the component in info.clamp.
[~, info] = il_pattern(tc, 1000, 130);
tc.verifyEqual([info.clamp.n info.clamp.LdB info.clamp.freq], [1 130 1000]);
end

function test_a_component_above_24_Bark_runs(tc)
% A component at 16 kHz (above 24 Bark) stopped the function before PR 77.
[~, info, ei] = il_pattern(tc, 16000, 70);
tc.verifyEqual(info.n_components, 1);
tc.verifyTrue(all(isfinite(ei(:))));
end

function [e, info, ei] = il_pattern(tc, f, L)
% the level of every channel at the bin of a single component of L dB at f,
% calibrated as the filterbank expects (magnitude 10^(L/20))
fs = tc.TestData.fs;
N = tc.TestData.N;
k = round(f / (fs/N)) + 1;
spec = zeros(1, N);
spec(k) = 10^(L/20);
[ei, info, ef] = Terhardt_filterbank(spec, tc.TestData.P);
e = ef(:, k).';
end
