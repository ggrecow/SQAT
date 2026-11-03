function tests = shmSignalSegment_test
% Unit tests of shmSignalSegment (utilities/ECMA418_2), the segmentation
% of ECMA-418-2:2025 clause 5.1.5.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
end

function test_blocks_start_every_hop_and_the_last_takes_the_end(tc)
% Formula 18 with a fixed block size s_b and hop s_h (5.1.5.2): block l
% holds samples l s_h + 1 .. l s_h + s_b. With endShrink, the last block
% holds the last s_b samples of the signal, so no sample is lost.
x = (1:103)';
[B, iStart] = shmSignalSegment(x, 1, 8, 0.5, 1, true);
tc.verifyEqual(B(:, 1), (1:8)');
tc.verifyEqual(B(:, 2), (5:12)');
tc.verifyEqual(B(:, end), (96:103)');
tc.verifyEqual(iStart(1:3), [1 5 9]);
tc.verifyEqual(iStart(end), 96);
end

function test_the_start_index_shifts_the_blocks(tc)
% i_start (Formula 19) moves the first block, so that every band starts at
% the same time reference.
[B, ~] = shmSignalSegment((1:100)', 1, 8, 0.5, 11, false);
tc.verifyEqual(B(:, 1), (11:18)');
end

function test_a_signal_shorter_than_a_block_is_an_error(tc)
% A signal of 8 samples, cut into blocks of 8 samples with 50 % overlap, makes
% shmSignalSegment stop with an error.
tc.verifyError(@() shmSignalSegment((1:8)', 1, 8, 0.5, 1, false), ?MException);
end
