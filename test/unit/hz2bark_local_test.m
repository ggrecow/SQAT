function tests = hz2bark_local_test
% Unit tests of hz2bark_local and bark2hz_local (utilities), the conversion
% between frequency and critical-band rate.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_formula_follows_the_critical_band_table_of_zwicker(tc)
% hz2bark_local is the analytical expression of Zwicker and Terhardt (1980),
% J. Acoust. Soc. Am. 68, 1523-1525, which gives the critical-band rate
% "with an accuracy better than 0.2 Bark" against the table of Zwicker. The
% table, held by SQAT in Get_Bark, sets the rate at the edges of the 25
% bands; there the largest difference is 0.197 Bark (4.4 kHz). The centres
% (z + 0.5) are an interpolation of the table, where the formula reaches
% 0.24 Bark, and are left out.
[~, T] = Get_Bark(4, [], []);
sqat_report_record('bark', 'hz2bark_local at the band edges of the table of Zwicker', hz2bark_local(T(:,2)), T(:,1), 0.2, 'Zwicker and Terhardt (1980): better than 0.2 Bark');
tc.verifyEqual(hz2bark_local(T(:,2)), T(:,1), 'AbsTol', 0.2);
end

function test_the_critical_band_rate_rises_with_frequency(tc)
% A critical-band rate must grow with frequency over the audible range.
z = hz2bark_local(logspace(log10(20), log10(20000), 1000));
tc.verifyGreaterThan(diff(z), 0);
end

function test_bark2hz_inverts_hz2bark(tc)
% bark2hz_local interpolates hz2bark_local on the one-third octave grid
% 1000*2^(k/3), k = -20..12: it is exact on the grid, and between grid points
% the round trip z -> f -> z stays within 0.13 Bark (set by us, twice the measured 0.061).
g = 1000 * 2.^((-20:12)/3);
tc.verifyEqual(bark2hz_local(hz2bark_local(g)), g, 'RelTol', 1e-12);
z = linspace(hz2bark_local(g(1)), hz2bark_local(g(end)), 2000);
sqat_report_record('bark', 'round trip z -> bark2hz_local -> hz2bark_local', hz2bark_local(bark2hz_local(z)), z, 0.13, 'set by us, twice the measured error');
tc.verifyEqual(hz2bark_local(bark2hz_local(z)), z, 'AbsTol', 0.13);
end

function test_bark2hz_gives_nan_below_its_grid(tc)
% The help of bark2hz_local sets its lowest input at 0.0972 Bark (about 10 Hz).
tc.verifyTrue(isnan(bark2hz_local(0.05)));
end
