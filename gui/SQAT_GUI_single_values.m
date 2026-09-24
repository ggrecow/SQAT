function T = SQAT_GUI_single_values(OUT, channel, n_channels)
% function T = SQAT_GUI_single_values(OUT, channel, n_channels)
%
%   Collects the single values of a SQAT output struct: its top-level real
%   numeric scalar fields, in field order. Nested structs (results of
%   sub-metrics), vectors, strings and the Bark step <dz> are left out.
%
%   An output of two channels (ECMA-418-2 with a binaural input) holds its
%   statistics as a row with one value per channel, [left right] or
%   [left right binaural], and the scalar values of the binaural result in
%   fields that end in Bin. With <channel> and <n_channels> the values of
%   one channel are returned, under the name of the quantity without the
%   ending Bin.
%
% INPUT ARGUMENTS
%   OUT : output struct of a SQAT metric
%   channel : (optional) 1, 2 or 'Binaural': the channel to collect
%   n_channels : (optional) number of channels of the input of the metric;
%                with 1 or no value the scalars are taken as they are
%
% OUTPUTS
%   T : table with the variables Quantity (cell array of char) and Value
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September 2026 assisted
% by Claude Opus 5 (Anthropic). All codes were verified by
% the authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 || isempty(n_channels)
    n_channels = 1;
end
if nargin < 2 || isempty(channel)
    channel = 1;
end
binaural = ischar(channel) || isstring(channel);
names = fieldnames(OUT);
labels = names;
keep = false(size(names));
values = nan(size(names));
for i = 1:numel(names)
    v = OUT.(names{i});
    if ~isnumeric(v) || ~isreal(v) || strcmp(names{i}, 'dz')
        continue
    end
    if n_channels <= 1
        if isscalar(v)
            keep(i) = true;
            values(i) = double(v);
        end
    elseif isscalar(v)
        is_bin = endsWith(names{i}, 'Bin');
        if is_bin == binaural                    % the value of the channel asked for
            if is_bin
                labels{i} = names{i}(1:end-3);
            end
            keep(i) = true;
            values(i) = double(v);
        end
    elseif isrow(v) && ismember(numel(v), [2 3]) && ~ismember(names{i}, {'time', 'timeOut'})
        c = channel;
        if binaural
            c = 3;
        end
        if c <= numel(v)
            keep(i) = true;
            values(i) = double(v(c));
        end
    end
end
T = table(labels(keep), values(keep), 'VariableNames', {'Quantity', 'Value'});
end
