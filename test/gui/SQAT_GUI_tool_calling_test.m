function tests = SQAT_GUI_tool_calling_test
% Tests of the GUI functions that call the SQAT metrics: no window opens, the
% metrics of the toolbox run for real.
% How to run the three test files and how long they take: test/gui/README.md.
tests = functiontests(localfunctions);
end

%% Fixtures ----------------------------------------------------------------

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and gui/, also on a clean MATLAB path
fs = 48000;
t = (0:1/fs:3-1/fs)';
% 1 kHz tone, 60 dB SPL, amplitude-modulated at 4 Hz (roughness and FS > 0)
tc.TestData.x_mono = sqrt(2)*2e-5*10^(60/20) * (1 + 0.5*sin(2*pi*4*t)) .* sin(2*pi*1000*t);
% flyover-like noise burst for EPNL
rng(1);
tc.TestData.x_burst = sqrt(2)*2e-5*10^(80/20) * sin(pi*t/t(end)).^2 .* randn(size(t))/3;
tc.TestData.fs = fs;
end

function teardown(~)
delete(findall(groot, 'Type', 'figure'));
end

%% Running a metric --------------------------------------------------------

function test_run_equals_the_direct_call(tc)
% The interface must add nothing: its call returns exactly what the metric returns.
x = tc.TestData.x_mono; xb = tc.TestData.x_burst; fs = tc.TestData.fs;
direct = {
 'Loudness_ISO532_1',            @() Loudness_ISO532_1(x, fs, 0, 2, 0.5, false)
 'Loudness_ECMA418_2',           @() Loudness_ECMA418_2(x, fs, 'free-frontal', 0.304, false)
 'Sharpness_DIN45692',           @() Sharpness_DIN45692(x, fs, 'DIN45692', 0, 2, 0.5, false, false)
 'Roughness_Daniel1997',         @() Roughness_Daniel1997(x, fs, 0, false)
 'Roughness_ECMA418_2',          @() Roughness_ECMA418_2(x, fs, 'free-frontal', 0.32, false)
 'FluctuationStrength_Osses2016',@() FluctuationStrength_Osses2016(x, fs, 1, 0, false)
 'Tonality_Aures1985',           @() Tonality_Aures1985(x, fs, 0, 0, false)
 'Tonality_ECMA418_2',           @() Tonality_ECMA418_2(x, fs, 'free-frontal', 0.304, false)
 'PsychoacousticAnnoyance_Widmann1992', @() PsychoacousticAnnoyance_Widmann1992(x, fs, 0, 0.2, false, false)
 'PsychoacousticAnnoyance_Zwicker1999', @() PsychoacousticAnnoyance_Zwicker1999(x, fs, 0, 0.2, false, false)
 'PsychoacousticAnnoyance_More2010',    @() PsychoacousticAnnoyance_More2010(x, fs, 0, 0.2, false, false)
 'PsychoacousticAnnoyance_Di2016',      @() PsychoacousticAnnoyance_Di2016(x, fs, 0, 0.2, false, false)
 'EPNL_FAR_Part36',              @() EPNL_FAR_Part36(xb, fs, 1, 0.5, 10, false)
};
m = SQAT_GUI_metrics;
for k = 1:size(direct, 1)
    e = m(strcmp({m.id}, direct{k,1}));
    p = il_default_params(e);
    if strcmp(e.id, 'EPNL_FAR_Part36'), sig = xb; else, sig = x; end
    [~, out_gui] = evalc('e.run(sig, fs, p, false)');
    [~, out_ref] = evalc('direct{k,2}()');
    tc.verifyTrue(isequaln(out_gui, out_ref), [direct{k,1} ': output differs from the direct call']);
end
end

function test_sound_level_follows_the_sound_level_meter_example(tc)
% The sound level of the GUI calls Do_SLM and Get_Leq as
% ex_sound_level_meter.m does: a 1 kHz tone of 60 dB SPL gives LAeq 60 dB
% within 0.1 dB, LAFmax equal to the maximum of Do_SLM, LAE = LAeq + 10 log10 T,
% and the one-third octave band at 1 kHz carries the tone. The Leq of the
% example averages the Fast-weighted level, whose integrator starts from
% zero: a 2 s tone reads 59.72 dB, so the tone lasts 10 s here.
fs = 48000; t = (0:10*fs-1)'/fs;
x = sqrt(2) * 2e-5 * 10^(60/20) * sin(2*pi*1000*t);
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Do_SLM'));
[~, OUT] = evalc('e.run(x, fs, il_default_params(e), false)');
L = Do_SLM(x, fs, 'A', 'f', 94);
tc.verifyEqual(OUT.LAeq, Get_Leq(L, fs));
tc.verifyEqual(OUT.LAeq, 60, 'AbsTol', 0.1);
tc.verifyEqual(OUT.LAFmax, max(L));
tc.verifyEqual(OUT.LAE, OUT.LAeq + 10*log10(10), 'AbsTol', 1e-9);
[~, k] = max(OUT.TOB_level);
tc.verifyEqual(OUT.TOB_freq(k), 1000, 'RelTol', 0.01);
T = SQAT_GUI_single_values(OUT);
tc.verifyEqual(T.Quantity', {'LAeq', 'LAFmax', 'LAF5', 'LAF90', 'LAE'});
end

function test_sound_level_bands_take_the_chosen_weighting(tc)
% The one-third octave levels are the bands of the signal after the weighting
% chosen for them: Z leaves the 1 kHz band of a tone as it is, A takes a
% 100 Hz tone down by its A-weighting, about 19 dB (IEC 61672-1, Table 3).
fs = 48000; t = (0:2*fs-1)'/fs;
x = sqrt(2) * 2e-5 * 10^(70/20) * sin(2*pi*100*t);
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Do_SLM'));
p = il_default_params(e);
tc.verifyEqual(p.tob_weight, 'A');                    % as the level, A by default
p.tob_weight = 'Z';
[~, Z] = evalc('e.run(x, fs, p, false)');
p.tob_weight = 'A';
[~, A] = evalc('e.run(x, fs, p, false)');
k = find(abs(Z.TOB_freq - 100) < 1);
tc.verifyEqual(Z.TOB_level(k) - A.TOB_level(k), 19.1, 'AbsTol', 0.2);
% the units follow the weightings: dB SPL unweighted, dBA and dBC weighted
a = SQAT_GUI_extract(A, 'Do_SLM');
tc.verifyEqual(a(strcmp({a.id}, 'level')).ylabel, 'Sound pressure level (dBA)');
tc.verifyEqual(a(strcmp({a.id}, 'tob_level')).ylabel, 'Band level (dBA)');
a = SQAT_GUI_extract(Z, 'Do_SLM');
tc.verifyEqual(a(strcmp({a.id}, 'tob_level')).ylabel, 'Band level (dB SPL)');
end

function test_run_passes_the_chosen_parameters(tc)
% The parameters chosen in the GUI reach the metric: a stationary ISO 532-1
% loudness and an Aures sharpness in diffuse field give the same output as the
% direct calls with those parameters.
x = tc.TestData.x_mono; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Loudness_ISO532_1'));
p = il_default_params(e); p.method = 1;             % stationary
[~, out] = evalc('e.run(x, fs, p, false)');
[~, ref] = evalc('Loudness_ISO532_1(x, fs, 0, 1, 0.5, false)');
tc.verifyTrue(isequaln(out, ref));
tc.verifyTrue(isfield(out, 'Loudness') && ~isfield(out, 'InstantaneousLoudness'));
e = m(strcmp({m.id}, 'Sharpness_DIN45692'));
p = il_default_params(e); p.weight_type = 'aures'; p.field = 1;
[~, out] = evalc('e.run(x, fs, p, false)');
[~, ref] = evalc('Sharpness_DIN45692(x, fs, ''aures'', 1, 2, 0.5, false, false)');
tc.verifyTrue(isequaln(out, ref));
end

function test_run_without_show_opens_no_figure(tc)
% A metric run by the GUI opens no figure unless the figure is asked for.
x = tc.TestData.x_mono; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Roughness_Daniel1997'));
before = numel(findall(groot, 'Type', 'figure'));
evalc('e.run(x, fs, il_default_params(e), false)');
tc.verifyEqual(numel(findall(groot, 'Type', 'figure')), before);
evalc('e.run(x, fs, il_default_params(e), true)');
tc.verifyGreaterThan(numel(findall(groot, 'Type', 'figure')), before);
end

%% Extracting results ------------------------------------------------------

function test_analyses_of_every_metric_follow_its_output(tc)
% The analyses on offer are the ones each metric returns, by metric and by
% method, with the shapes the plots rely on. Ids that must exist per case:
x = tc.TestData.x_mono; x2 = [x, x*10^(-10/20)]; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
cases = {
 'Loudness_ISO532_1', struct('method', 1), x, {'specific_loudness'}
 'Loudness_ISO532_1', struct('method', 2), x, {'loudness','loudness_level','specific_loudness_time'}
 'Sharpness_DIN45692', struct('method', 1), x, {}
 'Sharpness_DIN45692', struct('method', 2), x, {'sharpness'}
 'Roughness_Daniel1997', struct(), x, {'roughness','specific_roughness_time','specific_roughness'}
 'FluctuationStrength_Osses2016', struct('method', 0), x, {'specific_fs'}
 'FluctuationStrength_Osses2016', struct('method', 1), x, {'fs','specific_fs_time','specific_fs'}
 'Tonality_Aures1985', struct(), x, {'tonality','tonal_weighting','loudness_weighting'}
 'Loudness_ECMA418_2', struct(), x, {'loudness','specific_loudness_time','specific_tonal_loudness_time','specific_noise_loudness_time','specific_loudness'}
 'Roughness_ECMA418_2', struct(), x, {'roughness','specific_roughness_time','specific_roughness'}
 'Tonality_ECMA418_2', struct(), x, {'tonality','tonal_frequency','specific_tonality_time','specific_tonal_loudness_time','specific_noise_loudness_time','specific_tonality'}
 'PsychoacousticAnnoyance_Widmann1992', struct(), x, {'annoyance','weight_fr','weight_s'}
 'EPNL_FAR_Part36', struct(), x, {'pnlt','pnl','pn','spl','tob_spectra'}};
for k = 1:size(cases, 1)
    id = cases{k, 1};
    e = m(strcmp({m.id}, id));
    p = il_default_params(e);
    for f = fieldnames(cases{k, 2})', p.(f{1}) = cases{k, 2}.(f{1}); end
    [~, OUT] = evalc('e.run(cases{k, 3}, fs, p, false)');
    A = SQAT_GUI_extract(OUT, id, 1);
    label = sprintf('%s method %s', id, mat2str(il_field_or(p, 'method', [])));
    tc.verifyEqual({A.id}, cases{k, 4}, label);
    for a = A
        tc.verifyNotEmpty(a.label, a.id);
        tc.verifyTrue(ismember(a.kind, {'series','profile','map'}), a.id);
        tc.verifyTrue(iscolumn(a.x) && numel(a.x) > 1, [label ' ' a.id]);
        switch a.kind
            case 'series'
                tc.verifyTrue(iscolumn(a.y) && numel(a.y) == numel(a.x), [label ' ' a.id]);
            case 'profile'
                tc.verifyTrue(iscolumn(a.y) && numel(a.y) == numel(a.x), [label ' ' a.id]);
            case 'map'
                tc.verifyTrue(iscolumn(a.y) && numel(a.y) > 1, [label ' ' a.id]);
                tc.verifyEqual(size(a.z), [numel(a.x), numel(a.y)], [label ' ' a.id]);
        end
        tc.verifyTrue(all(isfinite(a.x)), [label ' ' a.id]);
    end
end
end

function test_analyses_take_the_values_of_the_output(tc)
% The analyses offered in the graphs (time series, Bark profiles and maps) take
% their values from the output of the metric, whether the map is stored as
% bands by time (roughness) or time by bands (fluctuation strength).
x = tc.TestData.x_mono; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
% series and a map that the metric stores as bands x time
e = m(strcmp({m.id}, 'Roughness_Daniel1997'));
[~, OUT] = evalc('e.run(x, fs, il_default_params(e), false)');
A = SQAT_GUI_extract(OUT, 'Roughness_Daniel1997', 1);
tc.verifyEqual(A(strcmp({A.id}, 'roughness')).y, OUT.InstantaneousRoughness(:));
tc.verifyEqual(A(strcmp({A.id}, 'roughness')).x, OUT.time(:));
sp = A(strcmp({A.id}, 'specific_roughness_time'));
tc.verifyEqual(sp.z, OUT.InstantaneousSpecificRoughness.');   % time x bands
tc.verifyEqual(sp.y, OUT.barkAxis(:));
pr = A(strcmp({A.id}, 'specific_roughness'));
tc.verifyEqual(pr.y, OUT.TimeAveragedSpecificRoughness(:));
% a map stored as time x bands, and a row vector on a column time axis
e = m(strcmp({m.id}, 'FluctuationStrength_Osses2016'));
[~, OUT] = evalc('e.run(x, fs, il_default_params(e), false)');
A = SQAT_GUI_extract(OUT, 'FluctuationStrength_Osses2016', 1);
tc.verifyEqual(A(strcmp({A.id}, 'fs')).y, OUT.InstantaneousFluctuationStrength(:));
tc.verifyEqual(A(strcmp({A.id}, 'specific_fs_time')).z, OUT.InstantaneousSpecificFluctuationStrength);
end

function test_analyses_pick_the_channel_of_a_binaural_output(tc)
% For a two-channel ECMA-418-2 loudness output, the analyses of channel 1 and
% channel 2 take the series, profile and map of that channel.
x = tc.TestData.x_mono; x2 = [x, x*10^(-10/20)]; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Loudness_ECMA418_2'));
[~, OUT] = evalc('e.run(x2, fs, il_default_params(e), false)');
for c = 1:2
    A = SQAT_GUI_extract(OUT, 'Loudness_ECMA418_2', c);
    tc.verifyEqual(A(strcmp({A.id}, 'loudness')).y, OUT.loudnessTDep(:, c));
    tc.verifyEqual(A(strcmp({A.id}, 'specific_loudness')).y, OUT.specLoudnessPowAvg(:, c));
    tc.verifyEqual(A(strcmp({A.id}, 'specific_loudness_time')).z, OUT.specLoudness(:, :, c));
end
A = SQAT_GUI_extract(OUT, 'Loudness_ECMA418_2', 'Binaural');
tc.verifyEqual(A(strcmp({A.id}, 'loudness')).y, OUT.loudnessTDepBin(:));
tc.verifyEqual(A(strcmp({A.id}, 'specific_loudness_time')).z, OUT.specLoudnessBin);
% a mono output has one channel only
[~, OUT1] = evalc('e.run(x, fs, il_default_params(e), false)');
tc.verifyEmpty(SQAT_GUI_extract(OUT1, 'Loudness_ECMA418_2', 2));
tc.verifyEmpty(SQAT_GUI_extract(OUT1, 'Loudness_ECMA418_2', 'Binaural'));
% the tonality of ECMA-418-2 has no combined binaural value
e = m(strcmp({m.id}, 'Tonality_ECMA418_2'));
[~, OUT] = evalc('e.run(x2, fs, il_default_params(e), false)');
tc.verifyNotEmpty(SQAT_GUI_extract(OUT, 'Tonality_ECMA418_2', 2));
tc.verifyEmpty(SQAT_GUI_extract(OUT, 'Tonality_ECMA418_2', 'Binaural'));
end

function test_single_values_of_a_binaural_output_go_by_channel(tc)
% For a two-channel ECMA-418-2 loudness output, the single values of each
% channel (Nmean, N5, loudnessPowAvg) are the ones of that channel.
x = tc.TestData.x_mono; x2 = [x, x*10^(-10/20)]; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Loudness_ECMA418_2'));
[~, OUT] = evalc('e.run(x2, fs, il_default_params(e), false)');
for c = 1:2
    T = SQAT_GUI_single_values(OUT, c, 2);
    tc.verifyEqual(T.Value(strcmp(T.Quantity, 'Nmean')), OUT.Nmean(c));
    tc.verifyEqual(T.Value(strcmp(T.Quantity, 'N5')), OUT.N5(c));
    tc.verifyEqual(T.Value(strcmp(T.Quantity, 'loudnessPowAvg')), OUT.loudnessPowAvg(c));
end
T = SQAT_GUI_single_values(OUT, 'Binaural', 2);
tc.verifyEqual(T.Value(strcmp(T.Quantity, 'Nmean')), OUT.Nmean(3));
tc.verifyEqual(T.Value(strcmp(T.Quantity, 'loudnessPowAvg')), OUT.loudnessPowAvgBin);
% a plain output keeps the scalars it always had
[~, OUT1] = evalc('e.run(x, fs, il_default_params(e), false)');
tc.verifyEqual(SQAT_GUI_single_values(OUT1, 1, 1), SQAT_GUI_single_values(OUT1));
% tonality has left and right only
e = m(strcmp({m.id}, 'Tonality_ECMA418_2'));
[~, OUT] = evalc('e.run(x2, fs, il_default_params(e), false)');
T = SQAT_GUI_single_values(OUT, 2, 2);
tc.verifyEqual(T.Value(strcmp(T.Quantity, 'Tmean')), OUT.Tmean(2));
tc.verifyEmpty(SQAT_GUI_single_values(OUT, 'Binaural', 2));
end

%% Sharing results between metrics -----------------------------------------

function test_share_takes_the_result_the_metric_itself_returns(tc)
% Every result taken from another metric must be, field by field, what the
% metric returns when it is called with the parameters of this run. The
% defaults are used here, so the time_skip of the components (0.5, 0.5, 0,
% 0) differs from the one of the models (0.2).
x = tc.TestData.x_mono; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
params = struct();
for k = 1:numel(m)
    params.(m(k).id) = il_default_params(m(k));
end
ids = {'PsychoacousticAnnoyance_Widmann1992', 'PsychoacousticAnnoyance_Zwicker1999', ...
    'PsychoacousticAnnoyance_More2010', 'Loudness_ISO532_1', 'Sharpness_DIN45692', ...
    'Roughness_Daniel1997', 'FluctuationStrength_Osses2016', 'Tonality_Aures1985'};
plan = SQAT_GUI_share(ids, params, numel(x), fs);
tc.assertEqual(numel(plan), numel(ids));

done = struct();
n_taken = 0;
for k = 1:numel(plan)
    step = plan(k);
    e = m(strcmp({m.id}, step.id)); %#ok<NASGU>
    if isempty(step.from)
        [~, OUT] = evalc('e.run(x, fs, params.(step.id), false)');
        done.(step.id) = OUT;
        continue
    end
    taken = SQAT_GUI_take(done.(step.from), step);
    [~, ref] = evalc('e.run(x, fs, params.(step.id), false)');
    tc.verifyTrue(isequaln(taken, ref), ...
        sprintf('%s taken from %s differs from the metric itself', step.id, step.from));
    done.(step.id) = taken;
    n_taken = n_taken + 1;
end
tc.verifyEqual(n_taken, 6, 'five metrics and Zwicker 1999 are taken from a model');
end

%% Helpers -----------------------------------------------------------------

function p = il_default_params(e)
p = struct();
for q = e.params
    p.(q.name) = q.value;
end
end

function v = il_field_or(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
end
