function tests = get_PNL_test
% Unit tests of get_PNL (EPNL_FAR_Part36/helper), the perceived noise level
% of 14 CFR Part 36 Appendix A (Sections A36.4.2.1 and A36.4.7).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_one_1kHz_band_reads_its_level_in_PNdB(tc)
% Table A36-3: at 1 kHz, SPL(b) = 40 dB gives 1 noy, and every 10 dB above
% doubles the noys (M(b) = 0.030103 = log10(2)/10). With one band, the total
% noisiness is 0.85 n + 0.15 n = n, so PNL = 40 + 33.22 log10(n) equals the
% band level from 40 dB up.
levels = 40:10:90;
PNL = zeros(size(levels));
for k = 1:numel(levels)
    L = -inf(1, 24);
    L(14) = levels(k);                                   % band 14 is 1 kHz
    [~, PNL(k)] = get_PNL(L);
end
sqat_report_record('epnl', 'PNL of a single 1 kHz band, 40 to 90 dB', PNL, levels, 1e-6, 'numerical: formulas of 14 CFR Part 36, A36.4');
tc.verifyEqual(PNL, levels, 'AbsTol', 1e-6);
end

function test_the_noys_of_the_other_bands_add_15_percent(tc)
% A36.4.2.1: N = 0.85 n(max) + 0.15 sum(n). Two bands of 1 noy each (1 kHz
% and 800 Hz at their SPL(b), 40 dB) give N = 0.85 + 0.3 = 1.15 noy.
L = -inf(1, 24);
L([13 14]) = 40;
PN = get_PNL(L);
tc.verifyEqual(PN, 1.15, 'AbsTol', 1e-12);
end

function test_bands_below_their_threshold_add_nothing(tc)
% Below SPL(d) of Table A36-3 a band has no noisiness: 1 kHz at 40 dB plus
% every other band at 0 dB gives 40 PNdB.
L = zeros(1, 24);
L(14) = 40;
[~, PNL] = get_PNL(L);
tc.verifyEqual(PNL, 40, 'AbsTol', 1e-9);
end
