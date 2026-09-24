function A = SQAT_GUI_extract(OUT, id, channel)
% function A = SQAT_GUI_extract(OUT, id, channel)
%
%   Lists the analyses that a SQAT output offers for plotting: the time
%   series, the profiles over the critical bands and the maps of time
%   against band, each one with its axes and units, in one orientation for
%   every metric (time along the first dimension of a map). An analysis is
%   returned when the output holds it with the expected shape, so the list
%   follows the method and the channels of the run (a stationary loudness
%   has a specific loudness profile and no time series).
%
% INPUT ARGUMENTS
%   OUT : output struct of the metric <id>
%   id : name of the SQAT function that returned OUT
%   channel : column of OUT to take: 1 or 2 for an output of two channels
%             (Loudness, Roughness and Tonality of ECMA-418-2), or
%             'Binaural' for the combined binaural result held in the
%             fields with the ending Bin. Default is 1.
%
% OUTPUTS
%   A : struct array, one element per analysis, with the fields
%       * id     - short name, the same for every metric that offers it
%       * label  - name shown in the interface
%       * kind   - 'series' (y against time x), 'profile' (y against band
%                  axis x) or 'map' (z, time by band, against time x and
%                  band axis y)
%       * x, y   - column vectors, the axes (y is the band axis of a map)
%       * z      - [numel(x) x numel(y)] matrix of a map, [] otherwise
%       * xlabel, ylabel, zlabel - axis labels with units
%       * bandscale - 'linear' or 'log', the scale of the axis that carries
%                  the bands or the frequency (x of a profile, y of a map)
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

if nargin < 3 || isempty(channel)
    channel = 1;
end
A = struct('id', {}, 'label', {}, 'kind', {}, 'x', {}, 'y', {}, 'z', {}, ...
           'xlabel', {}, 'ylabel', {}, 'zlabel', {}, 'bandscale', {});
defs = il_catalogue(id);
for k = 1:size(defs, 1)
    d = cell2struct(defs(k, :), {'id', 'label', 'kind', 'field', 'axis', 'unit', 'orient'}, 2);
    a = il_analysis(OUT, d, channel);
    if ~isempty(a)
        A(end+1) = a; %#ok<AGROW>
    end
end
end

%% -------------------------------------------------------------------------
function a = il_analysis(OUT, d, channel)
% one analysis of the catalogue, or [] when OUT does not hold it
a = [];
t = il_time(OUT);
switch d.kind
    case 'series'
        y = il_channel(OUT, d.field, channel);
        if isempty(t) || isempty(y) || ~isvector(y) || numel(y) ~= numel(t) || numel(t) < 2
            return
        end
        a = il_pack(d, t(:), y(:), [], 'Time (s)', d.unit, '', 'linear');
    case 'profile'
        b = il_axis(OUT, d.axis);
        y = il_channel(OUT, d.field, channel);
        if isempty(b) || isempty(y) || ~isvector(y) || numel(y) ~= numel(b) || numel(b) < 2
            return
        end
        a = il_pack(d, b(:), y(:), [], il_axis_label(d.axis), d.unit, '', il_axis_scale(d.axis));
    case 'map'
        b = il_axis(OUT, d.axis);
        Z = il_channel(OUT, d.field, channel);
        if isempty(t) || isempty(b) || isempty(Z) || ~ismatrix(Z) || numel(t) < 2 || numel(b) < 2
            return
        end
        if strcmp(d.orient, 'bt')
            Z = Z.';
        end
        if ~isequal(size(Z), [numel(t), numel(b)])
            return
        end
        a = il_pack(d, t(:), b(:), Z, 'Time (s)', il_axis_label(d.axis), d.unit, il_axis_scale(d.axis));
end
end

function a = il_pack(d, x, y, z, xlabel, ylabel, zlabel, bandscale)
a = struct('id', d.id, 'label', d.label, 'kind', d.kind, 'x', x, 'y', y, 'z', z, ...
           'xlabel', xlabel, 'ylabel', ylabel, 'zlabel', zlabel, 'bandscale', bandscale);
end

function t = il_time(OUT)
t = [];
if isfield(OUT, 'time')
    t = OUT.time;
elseif isfield(OUT, 'timeOut')
    t = OUT.timeOut;
end
end

function b = il_axis(OUT, name)
b = [];
if isfield(OUT, name)
    b = OUT.(name);
end
end

function s = il_axis_label(name)
switch name
    case 'barkAxis',        s = 'Critical band rate (Bark)';
    case 'bandCentreFreqs', s = 'Band centre frequency (Hz)';
    case 'TOB_freq',        s = 'One-third octave band (Hz)';
    otherwise,              s = name;
end
end

function s = il_axis_scale(name)
if ismember(name, {'bandCentreFreqs', 'TOB_freq'})
    s = 'log';
else
    s = 'linear';
end
end

function v = il_channel(OUT, field, channel)
% the values of one channel of a field: the column (or page) of a two
% channel output, the field with the ending Bin for the binaural result, and
% the field itself for an output of one channel
v = [];
if ischar(channel) || isstring(channel)          % 'Binaural'
    if isfield(OUT, [field 'Bin'])
        v = OUT.([field 'Bin']);
    end
    return
end
if ~isfield(OUT, field)
    return
end
v = OUT.(field);
if ~isnumeric(v)
    v = [];
    return
end
n = il_channels_of(OUT);
if n == 1
    if channel ~= 1
        v = [];
    end
    return
end
if channel > n
    v = [];
    return
end
switch ndims(v)
    case 3
        v = v(:, :, channel);
    case 2
        if size(v, 2) == n && ~isrow(v)
            v = v(:, channel);
        elseif isrow(v) && numel(v) == n
            v = v(channel);
        end
end
end

function n = il_channels_of(OUT)
% number of channels held by the fields of an ECMA-418-2 output
n = 1;
for f = {'loudnessTDep', 'roughnessTDep', 'tonalityTDep'}
    if isfield(OUT, f{1})
        n = size(OUT.(f{1}), 2);
    end
end
end

%% -------------------------------------------------------------------------
function defs = il_catalogue(id)
% {id, label, kind, field, axis, unit, orient} for each analysis of a metric
% ('orient' 'bt' when the metric stores a map as bands by time)
switch id
    case 'Loudness_ISO532_1'
        defs = {
            'loudness',               'Loudness vs time',                 'series',  'InstantaneousLoudness',              '',         'Loudness (sone)',                 'tb'
            'loudness_level',         'Loudness level vs time',           'series',  'InstantaneousLoudnessLevel',         '',         'Loudness level (phon)',           'tb'
            'specific_loudness_time', 'Specific loudness vs time',        'map',     'InstantaneousSpecificLoudness',      'barkAxis', 'Specific loudness (sone/Bark)',   'tb'
            'specific_loudness',      'Specific loudness',                'profile', 'SpecificLoudness',                   'barkAxis', 'Specific loudness (sone/Bark)',   'tb'};
    case 'Sharpness_DIN45692'
        defs = {
            'sharpness',              'Sharpness vs time',                'series',  'InstantaneousSharpness',             '',         'Sharpness (acum)',                'tb'};
    case 'Roughness_Daniel1997'
        defs = {
            'roughness',              'Roughness vs time',                'series',  'InstantaneousRoughness',             '',         'Roughness (asper)',               'tb'
            'specific_roughness_time','Specific roughness vs time',       'map',     'InstantaneousSpecificRoughness',     'barkAxis', 'Specific roughness (asper/Bark)', 'bt'
            'specific_roughness',     'Time-averaged specific roughness', 'profile', 'TimeAveragedSpecificRoughness',      'barkAxis', 'Specific roughness (asper/Bark)', 'tb'};
    case 'FluctuationStrength_Osses2016'
        defs = {
            'fs',                     'Fluctuation strength vs time',     'series',  'InstantaneousFluctuationStrength',   '',         'Fluctuation strength (vacil)',    'tb'
            'specific_fs_time',       'Specific fluctuation strength vs time', 'map', 'InstantaneousSpecificFluctuationStrength', 'barkAxis', 'Specific fluctuation strength (vacil/Bark)', 'tb'
            'specific_fs',            'Time-averaged specific fluctuation strength', 'profile', 'TimeAveragedSpecificFluctuationStrength', 'barkAxis', 'Specific fluctuation strength (vacil/Bark)', 'tb'};
    case 'Tonality_Aures1985'
        defs = {
            'tonality',               'Tonality vs time',                 'series',  'InstantaneousTonality',              '',         'Tonality (t.u.)',                 'tb'
            'tonal_weighting',        'Tonal weighting vs time',          'series',  'TonalWeighting',                     '',         'Tonal weighting (-)',             'tb'
            'loudness_weighting',     'Loudness weighting vs time',       'series',  'LoudnessWeighting',                  '',         'Loudness weighting (-)',          'tb'};
    case 'Loudness_ECMA418_2'
        defs = {
            'loudness',               'Loudness vs time',                 'series',  'loudnessTDep',                       '',         'Loudness (sone_{HMS})',           'tb'
            'specific_loudness_time', 'Specific loudness vs time',        'map',     'specLoudness',                       'bandCentreFreqs', 'Specific loudness (sone_{HMS}/Bark_{HMS})', 'tb'
            'specific_tonal_loudness_time', 'Specific tonal loudness vs time', 'map', 'specTonalLoudness',                  'bandCentreFreqs', 'Specific tonal loudness (sone_{HMS}/Bark_{HMS})', 'tb'
            'specific_noise_loudness_time', 'Specific noise loudness vs time', 'map', 'specNoiseLoudness',                  'bandCentreFreqs', 'Specific noise loudness (sone_{HMS}/Bark_{HMS})', 'tb'
            'specific_loudness',      'Time-averaged specific loudness',  'profile', 'specLoudnessPowAvg',                 'bandCentreFreqs', 'Specific loudness (sone_{HMS}/Bark_{HMS})', 'tb'};
    case 'Roughness_ECMA418_2'
        defs = {
            'roughness',              'Roughness vs time',                'series',  'roughnessTDep',                      '',         'Roughness (asper_{HMS})',         'tb'
            'specific_roughness_time','Specific roughness vs time',       'map',     'specRoughness',                      'bandCentreFreqs', 'Specific roughness (asper_{HMS}/Bark_{HMS})', 'tb'
            'specific_roughness',     'Time-averaged specific roughness', 'profile', 'specRoughnessAvg',                   'bandCentreFreqs', 'Specific roughness (asper_{HMS}/Bark_{HMS})', 'tb'};
    case 'Tonality_ECMA418_2'
        defs = {
            'tonality',               'Tonality vs time',                 'series',  'tonalityTDep',                       '',         'Tonality (tu_{HMS})',             'tb'
            'tonal_frequency',        'Frequency of the dominant tone vs time', 'series', 'tonalityTDepFreqs',             '',         'Frequency (Hz)',                  'tb'
            'specific_tonality_time', 'Specific tonality vs time',        'map',     'specTonality',                       'bandCentreFreqs', 'Specific tonality (tu_{HMS}/Bark_{HMS})', 'tb'
            'specific_tonal_loudness_time', 'Specific tonal loudness vs time', 'map', 'specTonalLoudness',                  'bandCentreFreqs', 'Specific tonal loudness (sone_{HMS}/Bark_{HMS})', 'tb'
            'specific_noise_loudness_time', 'Specific noise loudness vs time', 'map', 'specNoiseLoudness',                  'bandCentreFreqs', 'Specific noise loudness (sone_{HMS}/Bark_{HMS})', 'tb'
            'specific_tonality',      'Time-averaged specific tonality',  'profile', 'specTonalityAvg',                    'bandCentreFreqs', 'Specific tonality (tu_{HMS}/Bark_{HMS})', 'tb'};
    case {'PsychoacousticAnnoyance_Widmann1992', 'PsychoacousticAnnoyance_Zwicker1999', ...
          'PsychoacousticAnnoyance_More2010', 'PsychoacousticAnnoyance_Di2016'}
        defs = {
            'annoyance',              'Psychoacoustic annoyance vs time', 'series',  'InstantaneousPA',                    '',         'Psychoacoustic annoyance',        'tb'
            'weight_fr',              'Weighting of fluctuation strength and roughness vs time', 'series', 'wfr',          '',         'w_{FR}',                          'tb'
            'weight_s',               'Weighting of sharpness and loudness vs time', 'series', 'ws',                       '',         'w_{S}',                           'tb'};
    case 'EPNL_FAR_Part36'
        defs = {
            'pnlt',                   'PNLT vs time',                     'series',  'PNLT',                               '',         'PNLT (TPNdB)',                    'tb'
            'pnl',                    'PNL vs time',                      'series',  'PNL',                                '',         'PNL (PNdB)',                      'tb'
            'pn',                     'PN vs time',                       'series',  'PN',                                 '',         'PN (noy)',                        'tb'
            'spl',                    'SPL vs time',                      'series',  'InstantaneousSPL',                   '',         'SPL (dB)',                        'tb'
            'tob_spectra',            'One-third octave spectra vs time', 'map',     'SPL_TOB_spectra',                    'TOB_freq', 'SPL (dB)',                        'tb'};
    otherwise
        defs = cell(0, 7);
end
end
