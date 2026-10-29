function tests = get_exceeded_value_test
% Unit tests of get_exceeded_value (utilities), the value exceeded during a
% percentage of the time, the percentiles of every metric (N5, R90, ...).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_the_value_exceeded_during_p_percent_of_the_time(tc)
% Definition in the help: the value exceeded during P % of the time. On the
% samples 1..100, exactly P of them lie above the value exceeded during P %.
x = (1:100)';
for p = [1 5 10 50 90 95]
    v = get_exceeded_value(x, p);
    tc.verifyEqual(nnz(x > v), p, sprintf('samples above the value exceeded during %d %%', p));
end
tc.verifyEqual(get_exceeded_value(x, 100), 1, 'P100 is the minimum');
end

function test_the_percentiles_are_ordered(tc)
% A value exceeded more often cannot be larger: P1 >= P5 >= ... >= P95.
rng(3);
x = randn(1000, 1);
v = arrayfun(@(p) get_exceeded_value(x, p), [1 2 3 4 5 10 20 30 40 50 60 70 80 90 95]);
tc.verifyLessThanOrEqual(diff(v), 0);
end

function test_each_channel_is_its_own_column(tc)
% A stereo input [N x 2] gives one value per channel; a row is taken as a column.
x = (1:100)';
tc.verifyEqual(get_exceeded_value([x 2*x], 10), [90 180]);
tc.verifyEqual(get_exceeded_value(x', 10), 90);
end

function test_more_than_three_channels_is_an_error(tc)
% The help allows [N x 1], [N x 2] and [N x 3] (stereo with the binaural result).
tc.verifyError(@() get_exceeded_value(rand(4, 5), 5), ?MException);
end
