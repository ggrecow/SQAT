function OUT = SQAT_GUI_restat(OUT, id, time_skip)
% function OUT = SQAT_GUI_restat(OUT, id, time_skip)
%
%   Statistics of the output of a metric computed again for another
%   time_skip, with the same get_statistics the metric itself calls.
%
%   In the time-varying metrics that the psychoacoustic annoyance models
%   use, time_skip picks the point of the time series where the statistics
%   start (Loudness_ISO532_1 line 747 for method 2, Roughness_Daniel1997
%   line 273, FluctuationStrength_Osses2016 line 276, Tonality_Aures1985
%   line 471). The series itself does not depend on it, so the result a
%   model already holds serves any time_skip once its statistics are redone.
%
%   Loudness_ISO532_1 also states N_ratio from N5 and N95 (line 765), and
%   Sharpness_DIN45692 carries the loudness it used, which follows the same
%   time_skip.
%
% INPUT ARGUMENTS
%   OUT : output struct of the metric
%   id : id of the metric, as in SQAT_GUI_metrics
%   time_skip : time_skip (s) of this run
%
% OUTPUTS
%   OUT : the same output with the statistics of this time_skip
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

switch id
    case 'Loudness_ISO532_1'
        series_field = 'InstantaneousLoudness';
    case 'Sharpness_DIN45692'
        series_field = 'InstantaneousSharpness';
    case 'Roughness_Daniel1997'
        series_field = 'InstantaneousRoughness';
    case 'FluctuationStrength_Osses2016'
        series_field = 'InstantaneousFluctuationStrength';
    case 'Tonality_Aures1985'
        series_field = 'InstantaneousTonality';
    otherwise
        OUT = [];                       % no rule for this metric
        return
end
if ~isfield(OUT, series_field) || ~isfield(OUT, 'time')
    OUT = [];
    return
end

v = OUT.(series_field);
[~, idx] = min(abs(OUT.time - time_skip));
stat = get_statistics(v(idx:end), id);
for f = fieldnames(stat)'
    OUT.(f{1}) = stat.(f{1});
end

if strcmp(id, 'Loudness_ISO532_1')
    OUT.N_ratio = OUT.N5/OUT.N95;       % Loudness_ISO532_1 line 765
end
if strcmp(id, 'Sharpness_DIN45692') && isfield(OUT, 'loudness')
    OUT.loudness = SQAT_GUI_restat(OUT.loudness, 'Loudness_ISO532_1', time_skip);
end
end
