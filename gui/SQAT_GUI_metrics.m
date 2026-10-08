function metrics = SQAT_GUI_metrics
% function metrics = SQAT_GUI_metrics
%
%   Catalogue of the metrics offered by SQAT_GUI: how each one is called,
%   its parameters with their default values, and whether it takes a binaural
%   pair in one call.
%
% OUTPUT
%   metrics : struct array, one element per metric, with the fields
%       * id           - name of the SQAT function
%       * label        - name shown in the interface
%       * params       - struct array with the fields name, label, type
%                        ('choice' or 'number'), options (cell array
%                        {label, value} for 'choice') and value (default)
%       * run          - function handle OUT = run(insig, fs, p, show),
%                        where p is a struct with one field per parameter
%       * stereo       - true for the metrics that take a binaural pair in
%                        one call and return both channels and the combined
%                        binaural result (ECMA-418-2, except that the
%                        tonality has no binaural value)
%
%   The defaults follow psychoacoustic_metrics_get_defaults where it
%   defines the metric, and the header of the function otherwise.
%
% Author: Sergio Aguirre and Gil Felix Greco, September 2026
%
% AI disclosure: code development in September and October
% 2026 assisted by Claude Opus 5, Claude Sonnet 5 and Claude
% Opus 5.5 (Anthropic). All codes were verified by the
% authors.
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

field_iso  = il_choice('field', 'Sound field', {'Free field', 0; 'Diffuse field', 1}, 0);
field_ecma = il_choice('fieldtype', 'Sound field', {'Free-frontal', 'free-frontal'; 'Diffuse', 'diffuse'}, 'free-frontal');

d_loud  = psychoacoustic_metrics_get_defaults('Loudness_ISO532_1');
d_sharp = psychoacoustic_metrics_get_defaults('Sharpness_DIN45692');
d_rough = psychoacoustic_metrics_get_defaults('Roughness_Daniel1997');
d_fs    = psychoacoustic_metrics_get_defaults('FluctuationStrength_Osses2016');
d_ton   = psychoacoustic_metrics_get_defaults('Tonality_Aures1985');
d_pa    = psychoacoustic_metrics_get_defaults('PsychoacousticAnnoyance_Widmann1992');

metrics = struct('id', {}, 'label', {}, 'params', {}, 'run', {}, 'stereo', {});

metrics(end+1) = il_entry('Do_SLM', 'Sound level (IEC 61672-1)', ...
    [il_choice('weight_freq', 'Frequency weighting', {'A', 'A'; 'C', 'C'; 'Z', 'Z'}, 'A'), ...
     il_choice('weight_time', 'Time weighting', {'Fast', 'f'; 'Slow', 's'; 'Impulse', 'i'}, 'f'), ...
     il_choice('tob_weight', 'One-third octave weighting', {'Z (bands)', 'Z'; 'A (bands)', 'A'; 'C (bands)', 'C'}, 'A')], ...
    @il_sound_level);

metrics(end+1) = il_entry('Loudness_ISO532_1', 'Loudness (ISO 532-1)', ...
    [il_set(field_iso, d_loud.field), ...
     il_choice('method', 'Method', {'Stationary', 1; 'Time-varying', 2}, d_loud.method), ...
     il_number('time_skip', 'Time skip (s)', d_loud.time_skip)], ...
    @(x, fs, p, show) Loudness_ISO532_1(x, fs, p.field, p.method, p.time_skip, show));

metrics(end+1) = il_entry('Loudness_ECMA418_2', 'Loudness (ECMA-418-2)', ...
    [field_ecma, il_number('time_skip', 'Time skip (s)', 0.304)], ...
    @(x, fs, p, show) Loudness_ECMA418_2(x, fs, p.fieldtype, p.time_skip, show));
metrics(end).stereo = true;

metrics(end+1) = il_entry('Sharpness_DIN45692', 'Sharpness (DIN 45692)', ...
    [il_choice('weight_type', 'Weighting', {'DIN 45692', 'DIN45692'; 'Aures', 'aures'; 'von Bismarck', 'bismarck'}, d_sharp.weight_type), ...
     il_set(field_iso, d_sharp.field), ...
     il_choice('method', 'Loudness method', {'Stationary', 1; 'Time-varying', 2}, d_sharp.method), ...
     il_number('time_skip', 'Time skip (s)', d_sharp.time_skip)], ...
    @(x, fs, p, show) Sharpness_DIN45692(x, fs, p.weight_type, p.field, p.method, p.time_skip, show, false));

metrics(end+1) = il_entry('Roughness_Daniel1997', 'Roughness (Daniel 1997)', ...
    il_number('time_skip', 'Time skip (s)', d_rough.time_skip), ...
    @(x, fs, p, show) Roughness_Daniel1997(x, fs, p.time_skip, show));

metrics(end+1) = il_entry('Roughness_ECMA418_2', 'Roughness (ECMA-418-2)', ...
    [field_ecma, il_number('time_skip', 'Time skip (s)', 0.32)], ...
    @(x, fs, p, show) Roughness_ECMA418_2(x, fs, p.fieldtype, p.time_skip, show));
metrics(end).stereo = true;

metrics(end+1) = il_entry('FluctuationStrength_Osses2016', 'Fluctuation strength (Osses 2016)', ...
    [il_choice('method', 'Method', {'Stationary', 0; 'Time-varying', 1}, d_fs.method), ...
     il_number('time_skip', 'Time skip (s)', d_fs.time_skip)], ...
    @(x, fs, p, show) FluctuationStrength_Osses2016(x, fs, p.method, p.time_skip, show));

metrics(end+1) = il_entry('Tonality_Aures1985', 'Tonality (Aures 1985)', ...
    [il_set(field_iso, d_ton.Loudness_field), ...
     il_number('time_skip', 'Time skip (s)', d_ton.time_skip)], ...
    @(x, fs, p, show) Tonality_Aures1985(x, fs, p.field, p.time_skip, show));

metrics(end+1) = il_entry('Tonality_ECMA418_2', 'Tonality (ECMA-418-2)', ...
    [field_ecma, il_number('time_skip', 'Time skip (s)', 0.304)], ...
    @(x, fs, p, show) Tonality_ECMA418_2(x, fs, p.fieldtype, p.time_skip, show));
metrics(end).stereo = true;

pa_params = [il_set(field_iso, d_pa.Loudness_field), il_number('time_skip', 'Time skip (s)', d_pa.time_skip)];
metrics(end+1) = il_entry('PsychoacousticAnnoyance_Widmann1992', 'Annoyance (Widmann 1992)', pa_params, ...
    @(x, fs, p, show) PsychoacousticAnnoyance_Widmann1992(x, fs, p.field, p.time_skip, show, false));
metrics(end+1) = il_entry('PsychoacousticAnnoyance_Zwicker1999', 'Annoyance (Zwicker 1999)', pa_params, ...
    @(x, fs, p, show) PsychoacousticAnnoyance_Zwicker1999(x, fs, p.field, p.time_skip, show, false));
metrics(end+1) = il_entry('PsychoacousticAnnoyance_More2010', 'Annoyance (More 2010)', pa_params, ...
    @(x, fs, p, show) PsychoacousticAnnoyance_More2010(x, fs, p.field, p.time_skip, show, false));
metrics(end+1) = il_entry('PsychoacousticAnnoyance_Di2016', 'Annoyance (Di 2016)', pa_params, ...
    @(x, fs, p, show) PsychoacousticAnnoyance_Di2016(x, fs, p.field, p.time_skip, show, false));

metrics(end+1) = il_entry('EPNL_FAR_Part36', 'EPNL (FAR Part 36)', ...
    [il_number('dt', 'Time step dt (s)', 0.5), ...
     il_number('threshold', 'Threshold (TPNdB)', 10)], ...
    @(x, fs, p, show) EPNL_FAR_Part36(x, fs, 1, p.dt, p.threshold, show));

end

function OUT = il_sound_level(x, fs, p, show)
% the sound level meter of SQAT on a signal in Pa, as ex_sound_level_meter.m
% uses it: the time-weighted level, its equivalent and maximum, the levels
% exceeded 5 and 90 % of the time, the sound exposure level, and the
% one-third octave levels of Do_OB13_ISO532_1 (each band the Leq of the
% whole signal), after the weighting chosen for them (Z: none). No minimum:
% the time weighting starts from zero, so the first samples would give it
L = Do_SLM(x, fs, p.weight_freq, p.weight_time, 94);
if show
    Do_SLM(x, fs, p.weight_freq, p.weight_time, 94);   % no output: the figure of Do_SLM
end
L = L(:);
fw = upper(p.weight_freq);
tw = upper(p.weight_time);
step = max(1, round(fs / 1000));                     % the level vs time every millisecond
OUT.time = (0:step:numel(L) - 1)' / fs;
OUT.InstantaneousSPL = L(1:step:end);
OUT.(['L' fw 'eq']) = Get_Leq(L, fs);
OUT.(['L' fw tw 'max']) = max(L);
OUT.(['L' fw tw '5']) = get_exceeded_value(L, 5);
OUT.(['L' fw tw '90']) = get_exceeded_value(L, 90);
OUT.(['L' fw 'E']) = OUT.(['L' fw 'eq']) + 10*log10(numel(L) / fs);
if ~strcmpi(il_tob_weight(p), 'Z')
    [b, a] = Gen_weighting_filters(fs, p.tob_weight);
    x = filter(b, a, x);
end
[bands, fc] = Do_OB13_ISO532_1(x, fs);
OUT.TOB_freq = fc(:);
OUT.level_unit = il_level_unit(fw);                % for the axis labels of the GUI
OUT.TOB_unit = il_level_unit(il_tob_weight(p));
OUT.TOB_level = 20*log10(rms(bands, 1)' / 2e-5);
end

function w = il_tob_weight(p)
% the weighting of the one-third octave levels (Z in a session saved before it)
w = 'Z';
if isfield(p, 'tob_weight')
    w = upper(p.tob_weight);
end
end

function u = il_level_unit(w)
% dB SPL unweighted, dBA and dBC weighted
u = 'dB SPL';
if ~strcmpi(w, 'Z')
    u = ['dB' upper(w)];
end
end

function e = il_entry(id, label, params, run)
e = struct('id', id, 'label', label, 'params', params, 'run', run, 'stereo', false);
end

function p = il_choice(name, label, options, value)
p = struct('name', name, 'label', label, 'type', 'choice', 'options', {options}, 'value', value);
end

function p = il_number(name, label, value)
p = struct('name', name, 'label', label, 'type', 'number', 'options', {{}}, 'value', value);
end

function p = il_set(p, value)
p.value = value;
end
