function tests = Get_Bark_test
% Unit tests of Get_Bark (utilities), the critical-band rate of the FFT bins
% taken from the critical-band table of Zwicker.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_frequencies_of_the_table_get_their_own_rate(tc)
% At the lower edge and the centre of every band, the rate is the one of the
% table (z and z + 0.5), with no interpolation error.
[~, T] = Get_Bark(4, [], []);
f = sort([T(:,2); T(:,3)])';
z = sort([T(:,1); T(:,4)])';
qb = 1:numel(f);
B = Get_Bark(2*(numel(f) - 1), qb, f);
tc.verifyEqual(B(qb), z, 'AbsTol', 1e-12);
end

function test_bins_outside_qb_stay_at_zero(tc)
% The output has N/2 + 1 bins; only the bins listed in qb get a rate.
N = 64;
qb = 5:10;
B = Get_Bark(N, qb, linspace(500, 1000, numel(qb)));
tc.verifySize(B, [1 N/2 + 1]);
tc.verifyEqual(B(setdiff(1:N/2 + 1, qb)), zeros(1, N/2 + 1 - numel(qb)));
tc.verifyGreaterThan(B(qb), 0);
end
