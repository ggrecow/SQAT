function tests = get_statistics_test
% Unit tests of get_statistics (utilities), the statistics that every metric
% returns (Nmax, Nmean, N5, ...).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_fields_carry_the_letter_of_the_metric(tc)
% One field per statistic, named after the quantity of the metric.
S = get_statistics((1:10)', 'Roughness_Daniel1997');
names = strcat('R', {'max', 'min', 'mean', 'std', '1', '2', '3', '4', '5', ...
    '10', '20', '30', '40', '50', '60', '70', '80', '90', '95'});
tc.verifyEqual(sort(fieldnames(S))', sort(names));
end

function test_the_statistics_of_a_known_series(tc)
% On 1..100: max, min, mean and std as MATLAB computes them, and every
% percentile Nx as get_exceeded_value; N50 is the median (50.5 here), where
% get_exceeded_value would give 50.
x = (1:100)';
S = get_statistics(x, 'Loudness_ISO532_1');
tc.verifyEqual([S.Nmax S.Nmin S.Nmean S.Nstd], [100 1 mean(x) std(x)]);
for p = [1 2 3 4 5 10 20 30 40 60 70 80 90 95]
    tc.verifyEqual(S.(sprintf('N%d', p)), get_exceeded_value(x, p), sprintf('N%d', p));
end
tc.verifyEqual(S.N50, median(x));
end
