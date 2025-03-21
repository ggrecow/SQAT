function [signalSegmented, iBlocksOut] = shmSignalSegment(signal, axisn, blockSize, overlap,...
                                                          i_start, endShrink)
% signalSegmented = shmSignalSegment(signal, axisn, blockSize, overlap, i_start)
%
% Returns input signal segmented into blocks for processing.
%
% Inputs
% ------
% signal : vector or 2D matrix
%          the input signal/s
%
% axisn : integer (1 or 2, default: 1)
%         the axis along which to apply block segmentation
%
% blockSize : integer
%             the block size in samples
%
% overlap : double (>=0, < 1)
%           the proportion of overlap for each successive block
%
% i_start : integer (optional, default: 1)
%           the sample index from which to start the segmented signal
%
% endShrink : Boolean (optional, default: false)
%             option to include the end of the signal data in a block using
%             increased overlap with the preceding block
% 
% Returns
% -------
% For each channel in the input signal:
%
% signalSegmented : 2D or 3D matrix
%                   the segmented signal, arranged by channels over the
%                   last axis, samples (within each block) along the axis
%                   corresponding with axisn, and block number along the
%                   other remaining axis
%
% Also: 
%
% iBlocksOut : vector
%              the indices corresponding with each output block starting
%              index (NOTE: the indices corresponding with the input
%              indexing can be recovered by adding i_start to iBlocksOut)
%
% Assumptions
% -----------
% None
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 27/09/2023
% Date last modified: 19/03/2025
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        signal (:, :) double {mustBeReal}
        axisn (1, 1) {mustBeInteger, mustBeInRange(axisn, 1, 2)} = 1
        blockSize (1, 1) {mustBePositive, mustBeInteger} = false
        overlap (1, 1) {mustBeReal, mustBeGreaterThanOrEqual(overlap, 0),...
                        mustBeLessThan(overlap, 1)} = 0
        i_start (1, 1) {mustBePositive, mustBeInteger} = 1
        endShrink {mustBeNumericOrLogical} = false
    end

%% Signal pre-processing

% Orient input
if axisn == 2
    signal = signal.';
    axisFlip = true;
else
    axisFlip = false;
end

% Check sample index start will allow segmentation to proceed
if size(signal(i_start:end, :), 1) <= blockSize
    error("Signal is too short to apply segmentation using the selected parameters")
end

% Assign number of channels
nchans = size(signal, 2);

% Hop size
hopSize = (1 - overlap)*blockSize;

% Truncate the signal to start from i_start and to end at an index
% corresponding with the truncated signal length that will fill an
% integer number of overlapped blocks
signalTrunc = signal(i_start:end, :);
n_blocks = floor((size(signalTrunc, 1)...
                 - overlap*blockSize)/hopSize);
i_end = n_blocks*hopSize + overlap*blockSize;
signalTrunc = signalTrunc(1:i_end, :);

%% Signal segmentation

% Arrange the signal into overlapped blocks - each block reads
% along first axis, and each column is the succeeding overlapped
% block. 3 columns of zeros are appended to the left side of the
% matrix and the column shifted copies of this matrix are
% concatenated. The first 6 columns are then discarded as these all
% contain zeros from the appended zero columns.

for chan = nchans:-1:1
    signalSegmentedChan = [zeros(hopSize, 3),...
                           reshape(signalTrunc, hopSize, [])];
    
    signalSegmentedChan = cat(1, circshift(signalSegmentedChan, 3, 2),...
                              circshift(signalSegmentedChan, 2, 2),...
                              circshift(signalSegmentedChan, 1, 2),...
                              circshift(signalSegmentedChan, 0, 2));
    
    signalSegmentedChan = signalSegmentedChan(:, 7:end);

    % if branch to include block of end data with increased overlap
    if endShrink && (size(signal(i_start:end), 1) > size(signalTrunc, 1))
        signalSegmentedChanOut = [signalSegmentedChan, signal(end-blockSize + 1:end)];
        iBlocksOut = [1:hopSize:n_blocks*hopSize,...
                      size(signal(i_start:end), 1) - blockSize + 1];
    else
        signalSegmentedChanOut = signalSegmentedChan;
        iBlocksOut = 1:hopSize:n_blocks*hopSize;
    end

    signalSegmented(:, :, chan) = signalSegmentedChanOut;
end

% re-orient segmented signal to match input
if axisFlip
    signalSegmented = permute(signalSegmented, [2, 1, 3]);
end

% end of function

