function [t, y, name] = SQAT_GUI_series(OUT, metric)
% function [t, y, name] = SQAT_GUI_series(OUT, metric)
%
%   Extracts the time series that SQAT_GUI plots for a metric, as named in
%   its entry of SQAT_GUI_metrics. A stationary result has no time series,
%   and then t and y are empty.
%
% INPUT ARGUMENTS
%   OUT : output struct of the metric
%   metric : element of SQAT_GUI_metrics
%
% OUTPUTS
%   t : [Nx1] time (s)
%   y : [Nx1] values (first channel if OUT holds several)
%   name : label of the series, with its unit
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

t = []; y = []; name = metric.series_label;
time_field = metric.series{1};
value_field = metric.series{2};
if ~isfield(OUT, time_field) || ~isfield(OUT, value_field)
    return
end
t = OUT.(time_field)(:);
Y = OUT.(value_field);
if isvector(Y)
    y = Y(:);
elseif size(Y, 1) == numel(t)
    y = Y(:, 1);
else
    y = Y(1, :).';
end
end
