function tests = phon2sone_local_test
% Unit tests of phon2sone_local and sone2phon_local (utilities), the
% conversion between loudness level and loudness of ISO 532-1:2017.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_40_phon_is_1_sone_and_10_phon_more_doubles_it(tc)
% The definition of the sone: 40 phon is 1 sone, and above it every 10 phon
% doubles the loudness, N = 2^((LN - 40)/10).
tc.verifyEqual(phon2sone_local([40 50 60 70 100]), [1; 2; 4; 8; 64], 'RelTol', 1e-12);
end

function test_sone2phon_uses_the_formula_of_the_standard(tc)
% ISO 532-1 gives LN = 40 + 33.22 log10(N) for N >= 1. The rounded constant
% 33.22 (10/log10(2) is 33.219) puts 64 sone at 100.0013 phon; tolerance 0.003
% phon, set by us (twice the measured error).
got = sone2phon_local([1 2 4 64]);
sqat_report_record('loudness_level', 'sone2phon_local at 1, 2, 4 and 64 sone', got, [40; 50; 60; 100], 0.003, 'set by us, twice the measured error (constant 33.22 of ISO 532-1)');
tc.verifyEqual(got, [40; 50; 60; 100], 'AbsTol', 0.003);
end

function test_the_round_trip_above_40_phon(tc)
% phon -> sone -> phon from 40 to 120 phon; the constant 33.22 leaves at most
% 1.7e-3 phon (measured); tolerance 0.004 phon, set by us (twice the measured).
p = (40:0.5:120)';
got = sone2phon_local(phon2sone_local(p));
sqat_report_record('loudness_level', 'round trip phon -> sone -> phon, 40 to 120 phon', got, p, 0.004, 'set by us, twice the measured error');
tc.verifyEqual(got, p, 'AbsTol', 0.004);
end

function test_the_output_is_a_column(tc)
% Both functions return a column, whatever the shape of the input.
tc.verifySize(phon2sone_local([40 50 60]), [3 1]);
tc.verifySize(sone2phon_local([1 2 4]), [3 1]);
end
