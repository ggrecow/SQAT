function tests = EPNL_FAR_Part36_test
% Unit tests of EPNL_FAR_Part36, the effective perceived noise level of
% 14 CFR Part 36 Appendix A.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_pnl_of_a_tone_follows_its_bands(tc)
% A 1 kHz tone of 40 dB SPL is 1 noy in the 1 kHz band, but the one-third
% octave filters leave about 20 dB in the 800 Hz and 1.25 kHz bands, above
% their SPL(d) of Table A36-3; with N = 0.85 n(max) + 0.15 sum(n)
% (A36.4.2.1) the PNL is 40.76 PNdB. The PNL of every time step is the one
% get_PNL gives on the band levels the function returns.
t = (0:5*48000-1)' / 48000;
x = sqrt(2) * 2e-5 * 10^(40/20) * sin(2*pi*1000*t);
[~, O] = evalc('EPNL_FAR_Part36(x, 48000, 1, 0.5, 10, false)');
[~, PNL] = get_PNL(O.SPL_TOB_spectra);
tc.verifyEqual(O.PNL, PNL, 'AbsTol', 1e-9);
tc.verifyGreaterThan(O.PNLM, 40);
tc.verifyLessThan(O.PNLM, 41);
end

function test_silence_gives_minus_infinity(tc)
% No band above its SPL(d): 0 noy, and 40 + 33.22 log10(0) is -Inf.
[~, O] = evalc('EPNL_FAR_Part36(zeros(5*48000, 1), 48000, 1, 0.5, 10, false)');
tc.verifyEqual(O.PNLM, -Inf);
end
