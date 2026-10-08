function OUT = SQAT_GUI_take(src, step)
% function OUT = SQAT_GUI_take(src, step)
%
%   Result of one metric taken from the output of another metric of the
%   same run, as planned by SQAT_GUI_share. Empty when that output does
%   not hold it.
%
% INPUT ARGUMENTS
%   src : output struct of the metric named in step.from
%   step : element of the plan of SQAT_GUI_share
%
% OUTPUTS
%   OUT : output struct of the metric of step.id, as the metric itself
%         returns it with the parameters of this run
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

OUT = [];
if isempty(src) || ~isstruct(src)
    return
end
if isempty(step.field)
    OUT = src;
elseif isfield(src, step.field)
    OUT = src.(step.field);
else
    return
end
if ~isempty(step.attach) && isfield(src, step.attach{2})
    OUT.(step.attach{1}) = src.(step.attach{2});
end
if ~isempty(step.restat)
    % the model computed this series with its own time_skip, which only
    % moves the start of the statistics
    OUT = SQAT_GUI_restat(OUT, step.id, step.restat);
end
end
