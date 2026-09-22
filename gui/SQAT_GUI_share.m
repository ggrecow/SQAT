function plan = SQAT_GUI_share(ids, params, n_samples, fs)
% function plan = SQAT_GUI_share(ids, params, n_samples, fs)
%
%   Order in which SQAT_GUI runs the selected metrics, and which of them
%   take their result from another metric of the same run instead of
%   computing it again.
%
%   The psychoacoustic annoyance models compute loudness, sharpness,
%   roughness, fluctuation strength and (More 2010, Di 2016) tonality on
%   the way to their own result, and return them in OUT.L, OUT.S, OUT.R,
%   OUT.FS and OUT.K. When the user also selects one of those metrics with
%   the parameters the annoyance model uses internally, the result is
%   already there. PsychoacousticAnnoyance_Zwicker1999 is a wrapper that
%   calls PsychoacousticAnnoyance_Widmann1992 with its own arguments, so
%   the two return the same output for the same parameters.
%
%   A metric that takes its result from another one is not drawn by the
%   toolbox: its figure is drawn by the graphs window when it is asked for.
%
% INPUT ARGUMENTS
%   ids : cell array of char, ids of the selected metrics
%   params : struct, params.(id) holds the parameters of each metric
%   n_samples : number of samples of the signal
%   fs : sampling frequency (Hz)
%
% OUTPUTS
%   plan : [1xN] struct array in the order to run, with fields
%          id : id of the metric
%          from : id of the metric that already holds this result, or ''
%          field : field of that output to take ('L', 'S', 'R', 'FS', 'K'),
%                  or '' for the whole output
%          attach : {name, field} to add to the result taken, so that it
%                   carries what the metric itself returns, or {}
%          restat : time_skip (s) whose statistics the result taken needs,
%                   or [] when the model already used this one
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

ids = cellstr(ids);
ids = ids(:)';
models = {'PsychoacousticAnnoyance_Widmann1992', 'PsychoacousticAnnoyance_Zwicker1999', ...
    'PsychoacousticAnnoyance_More2010', 'PsychoacousticAnnoyance_Di2016'};

% the annoyance models run first, in this order, because they hold the
% results the other metrics reuse
selected_models = models(ismember(models, ids));
order = [selected_models, ids(~ismember(ids, models))];

plan = struct('id', order, 'from', {''}, 'field', {''}, 'attach', {{}}, 'restat', {[]});
sources = selected_models;

for k = 1:numel(plan)
    id = plan(k).id;
    for s = 1:numel(sources)
        src = sources{s};
        if strcmp(src, id) || find(strcmp(order, src), 1) >= k
            continue   % a source has to have run before this metric is reached
        end
        [field, attach, restat] = il_shared_field(id, src, params, n_samples, fs);
        if isempty(field)
            continue
        end
        plan(k).from = src;
        plan(k).attach = attach;
        plan(k).restat = restat;
        if ~strcmp(field, 'whole')
            plan(k).field = field;
        end
        sources = sources(~strcmp(sources, id));   % a copy is no source
        break
    end
end
end

%% -------------------------------------------------------------------------
function [field, attach, restat] = il_shared_field(id, src, params, n_samples, fs)
% Which field of the output of src holds the result of id, with the
% parameters of this run. Empty when src does not hold it.
%
% In the metrics the models use, time_skip only picks where the statistics
% start on the time series (see SQAT_GUI_restat), so a different time_skip
% asks for the statistics again and nothing more.
field = '';
attach = {};
restat = [];
if isfield(params.(id), 'time_skip') && isfield(params.(src), 'time_skip') ...
        && ~isequal(params.(id).time_skip, params.(src).time_skip)
    restat = params.(id).time_skip;
end
p = params.(id);
q = params.(src);
models_with_tonality = {'PsychoacousticAnnoyance_More2010', 'PsychoacousticAnnoyance_Di2016'};

switch id
    case {'PsychoacousticAnnoyance_Widmann1992', 'PsychoacousticAnnoyance_Zwicker1999'}
        % Zwicker 1999 is a wrapper of Widmann 1992 called with the same
        % arguments (PsychoacousticAnnoyance_Zwicker1999.m, line 135)
        pair = {'PsychoacousticAnnoyance_Widmann1992', 'PsychoacousticAnnoyance_Zwicker1999'};
        if all(ismember({id, src}, pair)) && isequal(p.field, q.field) ...
                && isequal(p.time_skip, q.time_skip)
            field = 'whole';
        end

    case 'Loudness_ISO532_1'
        % Loudness_ISO532_1(insig, fs, LoudnessField, 2, time_skip, 0)
        if isequal(p.field, q.field) && isequal(p.method, 2)
            field = 'L';
        end

    case 'Sharpness_DIN45692'
        % Sharpness_DIN45692_from_loudness(specific loudness of the call
        % above, 'DIN45692', time, time_skip, 0)
        if isequal(p.field, q.field) && strcmp(p.weight_type, 'DIN45692') ...
                && isequal(p.method, 2)
            field = 'S';
            % the metric returns the loudness it used; the model holds the
            % same one in OUT.L, from the same call
            attach = {'loudness', 'L'};
        end

    case 'Roughness_Daniel1997'
        % Roughness_Daniel1997(insig, fs, time_skip, 0)
        field = 'R';

    case 'FluctuationStrength_Osses2016'
        % FluctuationStrength_Osses2016(insig, fs, method_FS, time_skip, 0),
        % with the stationary method for a signal shorter than 2 s
        method_FS = double((n_samples - 1)/fs >= 2);
        if isequal(p.method, method_FS)
            field = 'FS';
        end

    case 'Tonality_Aures1985'
        % Tonality_Aures1985(insig, fs, LoudnessField, 0, 0), in the models
        % that use tonality
        if ismember(src, models_with_tonality) && isequal(p.field, q.field)
            field = 'K';
            restat = [];
            if ~isequal(p.time_skip, 0)   % the models call it with time_skip 0
                restat = p.time_skip;
            end
        end
end
end
