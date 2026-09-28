function tests = tSQAT_GUI
% Tests of the SQAT graphical interface (gui/). How to run them and how long
% they take: test/README.md.
%
% Run with:
%   matlab -batch "startup_SQAT; cd test; r = runtests('tSQAT_GUI'); disp(table(r)); assert(all([r.Passed]))"
tests = functiontests(localfunctions);
end

%% Fixtures ----------------------------------------------------------------

function setupOnce(tc)
setappdata(groot, 'sqat_gui_mute', true);           % the player runs, with a silent buffer
addpath(fullfile(basepath_SQAT, 'gui'));
fs = 48000;
t = (0:1/fs:3-1/fs)';
% 1 kHz tone, 60 dB SPL, amplitude-modulated at 4 Hz (roughness and FS > 0)
x_mono = sqrt(2)*2e-5*10^(60/20) * (1 + 0.5*sin(2*pi*4*t)) .* sin(2*pi*1000*t);
% stereo: channel 2 is channel 1 attenuated by 10 dB
x_stereo = [x_mono, x_mono*10^(-10/20)];
% flyover-like noise burst for EPNL
rng(1);
x_burst = sqrt(2)*2e-5*10^(80/20) * sin(pi*t/t(end)).^2 .* randn(size(t))/3;
dir_tmp = tempname; mkdir(dir_tmp);
tc.TestData.dir_tmp = dir_tmp;
tc.TestData.fs = fs;
tc.TestData.x_mono = x_mono;
tc.TestData.x_burst = x_burst;
tc.TestData.wav_mono = fullfile(dir_tmp, 'tone_mono.wav');
tc.TestData.wav_stereo = fullfile(dir_tmp, 'tone_stereo.wav');
tc.TestData.wav_burst = fullfile(dir_tmp, 'burst.wav');
audiowrite(tc.TestData.wav_mono, x_mono, fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_stereo, x_stereo, fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_burst, x_burst, fs, 'BitsPerSample', 32);
% 1 kHz, 60 dB SPL, 100 % modulated at 70 Hz: nearly constant roughness
x_rough = sqrt(2)*2e-5*10^(60/20) * (1 + sin(2*pi*70*t)) .* sin(2*pi*1000*t) / sqrt(1.5);
tc.TestData.wav_rough = fullfile(dir_tmp, 'am70.wav');
audiowrite(tc.TestData.wav_rough, x_rough, fs, 'BitsPerSample', 32);
% two steady tones, 500 Hz and 2 kHz, for the filter of the waveform window
tc.TestData.wav_two = fullfile(dir_tmp, 'two_tones.wav');
audiowrite(tc.TestData.wav_two, 0.1*sin(2*pi*500*t) + 0.1*sin(2*pi*2000*t), fs, 'BitsPerSample', 32);
% a short tone, to reach the end of the file while a test waits
tc.TestData.wav_short = fullfile(dir_tmp, 'short.wav');
audiowrite(tc.TestData.wav_short, 0.1*sin(2*pi*440*t(1:round(0.3*fs))), fs, 'BitsPerSample', 32);
% steady 1 kHz tone at 60 dB SPL, for the level of the spectrogram
tc.TestData.wav_tone = fullfile(dir_tmp, 'tone_1k_60dB.wav');
audiowrite(tc.TestData.wav_tone, sqrt(2)*2e-5*10^(60/20)*sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
end

function teardownOnce(tc)
if isappdata(groot, 'sqat_gui_mute')
    rmappdata(groot, 'sqat_gui_mute');
end
rmdir(tc.TestData.dir_tmp, 's');
end

function teardown(~)
delete(findall(groot, 'Type', 'figure'));
end

%% Catalogue ---------------------------------------------------------------

function test_catalogue_lists_every_metric_of_sqat(tc)
m = SQAT_GUI_metrics;
expected = {'Loudness_ISO532_1','Loudness_ECMA418_2','Sharpness_DIN45692', ...
    'Roughness_Daniel1997','Roughness_ECMA418_2','FluctuationStrength_Osses2016', ...
    'Tonality_Aures1985','Tonality_ECMA418_2', ...
    'PsychoacousticAnnoyance_Widmann1992','PsychoacousticAnnoyance_Zwicker1999', ...
    'PsychoacousticAnnoyance_More2010','PsychoacousticAnnoyance_Di2016', ...
    'EPNL_FAR_Part36'};
tc.verifyEqual({m.id}, expected);
for k = 1:numel(m)
    tc.verifyEqual(exist(m(k).id, 'file'), 2, m(k).id);
    tc.verifyNotEmpty(m(k).label, m(k).id);
    tc.verifyTrue(isa(m(k).run, 'function_handle'), m(k).id);
end
end

function test_catalogue_parameters_are_well_formed(tc)
m = SQAT_GUI_metrics;
for k = 1:numel(m)
    for p = m(k).params
        tc.verifyTrue(ismember(p.type, {'choice','number'}), [m(k).id '.' p.name]);
        tc.verifyNotEmpty(p.label, [m(k).id '.' p.name]);
        if strcmp(p.type, 'choice')
            tc.verifyTrue(any(cellfun(@(o) isequal(o, p.value), p.options(:,2))), ...
                [m(k).id '.' p.name ': default not among the options']);
        else
            tc.verifyTrue(isnumeric(p.value) && isscalar(p.value), [m(k).id '.' p.name]);
        end
    end
end
end

function test_catalogue_defaults_follow_the_toolbox(tc)
m = SQAT_GUI_metrics;
v = @(id, name) m(strcmp({m.id}, id)).params(strcmp({m(strcmp({m.id}, id)).params.name}, name)).value;
d = psychoacoustic_metrics_get_defaults('Loudness_ISO532_1');
tc.verifyEqual(v('Loudness_ISO532_1','field'), d.field);
tc.verifyEqual(v('Loudness_ISO532_1','method'), d.method);
tc.verifyEqual(v('Loudness_ISO532_1','time_skip'), d.time_skip);
d = psychoacoustic_metrics_get_defaults('Sharpness_DIN45692');
tc.verifyEqual(v('Sharpness_DIN45692','weight_type'), d.weight_type);
tc.verifyEqual(v('Sharpness_DIN45692','field'), d.field);
tc.verifyEqual(v('Sharpness_DIN45692','method'), d.method);
tc.verifyEqual(v('Sharpness_DIN45692','time_skip'), d.time_skip);
d = psychoacoustic_metrics_get_defaults('Roughness_Daniel1997');
tc.verifyEqual(v('Roughness_Daniel1997','time_skip'), d.time_skip);
d = psychoacoustic_metrics_get_defaults('FluctuationStrength_Osses2016');
tc.verifyEqual(v('FluctuationStrength_Osses2016','method'), d.method);
tc.verifyEqual(v('FluctuationStrength_Osses2016','time_skip'), d.time_skip);
d = psychoacoustic_metrics_get_defaults('Tonality_Aures1985');
tc.verifyEqual(v('Tonality_Aures1985','field'), d.Loudness_field);
tc.verifyEqual(v('Tonality_Aures1985','time_skip'), d.time_skip);
d = psychoacoustic_metrics_get_defaults('PsychoacousticAnnoyance_Widmann1992');
for id = {'PsychoacousticAnnoyance_Widmann1992','PsychoacousticAnnoyance_Zwicker1999', ...
          'PsychoacousticAnnoyance_More2010','PsychoacousticAnnoyance_Di2016'}
    tc.verifyEqual(v(id{1},'field'), d.Loudness_field, id{1});
    tc.verifyEqual(v(id{1},'time_skip'), d.time_skip, id{1});
end
% ECMA-418-2 and EPNL: defaults stated in the headers of the functions
tc.verifyEqual(v('Loudness_ECMA418_2','fieldtype'), 'free-frontal');
tc.verifyEqual(v('Loudness_ECMA418_2','time_skip'), 0.304);
tc.verifyEqual(v('Roughness_ECMA418_2','fieldtype'), 'free-frontal');
tc.verifyEqual(v('Roughness_ECMA418_2','time_skip'), 0.32);
tc.verifyEqual(v('Tonality_ECMA418_2','fieldtype'), 'free-frontal');
tc.verifyEqual(v('Tonality_ECMA418_2','time_skip'), 0.304);
tc.verifyEqual(v('EPNL_FAR_Part36','dt'), 0.5);
tc.verifyEqual(v('EPNL_FAR_Part36','threshold'), 10);
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

function test_run_passes_the_chosen_parameters(tc)
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
x = tc.TestData.x_mono; fs = tc.TestData.fs;
m = SQAT_GUI_metrics;
e = m(strcmp({m.id}, 'Roughness_Daniel1997'));
before = numel(findall(groot, 'Type', 'figure'));
evalc('e.run(x, fs, il_default_params(e), false)');
tc.verifyEqual(numel(findall(groot, 'Type', 'figure')), before);
evalc('e.run(x, fs, il_default_params(e), true)');
tc.verifyGreaterThan(numel(findall(groot, 'Type', 'figure')), before);
end

%% Loading audio -----------------------------------------------------------

function test_load_calibrates_with_dbfs_94(tc)
[x, fs, nch] = SQAT_GUI_load(tc.TestData.wav_mono, 94, 1);
[ref, fs_ref] = audioread(tc.TestData.wav_mono);
tc.verifyEqual(x, ref(:,1));
tc.verifyEqual(fs, fs_ref);
tc.verifyEqual(nch, 1);
end

function test_load_applies_the_dbfs_gain(tc)
x94  = SQAT_GUI_load(tc.TestData.wav_mono, 94, 1);
x100 = SQAT_GUI_load(tc.TestData.wav_mono, 100, 1);
tc.verifyEqual(x100, x94 * 10^((100-94)/20), 'AbsTol', 1e-15);
end

function test_load_selects_the_channel(tc)
[x2, ~, nch] = SQAT_GUI_load(tc.TestData.wav_stereo, 94, 2);
ref = audioread(tc.TestData.wav_stereo);
tc.verifyEqual(nch, 2);
tc.verifyEqual(x2, ref(:,2));
tc.verifyError(@() SQAT_GUI_load(tc.TestData.wav_stereo, 94, 3), 'SQAT_GUI:channel');
end

function test_load_reads_several_channels_at_once(tc)
ref = audioread(tc.TestData.wav_stereo);
X = SQAT_GUI_load(tc.TestData.wav_stereo, 94, [1 2]);
tc.verifyEqual(X, ref);
X = SQAT_GUI_load(tc.TestData.wav_stereo, 100, [2 1]);
tc.verifyEqual(X, ref(:, [2 1]) * 10^(6/20));
tc.verifyError(@() SQAT_GUI_load(tc.TestData.wav_stereo, 94, [1 3]), 'SQAT_GUI:channel');
tc.verifyError(@() SQAT_GUI_load(tc.TestData.wav_stereo, 94, [1 1.5]), 'SQAT_GUI:channel');
end

%% Extracting results ------------------------------------------------------

function test_single_values_are_the_top_level_scalars(tc)
OUT.time = (0:9)';
OUT.InstantaneousLoudness = rand(10,1);
OUT.Nmax = 3; OUT.N5 = 2.5;
OUT.dz = 0.5;                       % step of the Bark axis, not a result
OUT.L = struct('N5', 1);            % nested result of a sub-metric
OUT.soundField = "free-frontal";
T = SQAT_GUI_single_values(OUT);
tc.verifyEqual(T.Quantity, {'Nmax'; 'N5'});
tc.verifyEqual(T.Value, [3; 2.5]);
end

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

function test_share_plan_follows_the_parameters(tc)
% A result is taken from another metric only when that metric computed the
% same thing. A different time_skip asks for the statistics again.
m = SQAT_GUI_metrics;
params = struct();
for k = 1:numel(m)
    params.(m(k).id) = il_default_params(m(k));
end
plan = SQAT_GUI_share({m.id}, params, 3*48000, 48000);
step = @(id) plan(strcmp({plan.id}, id));
tc.verifyEqual(plan(1).id, 'PsychoacousticAnnoyance_Widmann1992', 'the models run first');
% Zwicker 1999 calls Widmann 1992 with its own arguments
tc.verifyEqual(step('PsychoacousticAnnoyance_Zwicker1999').from, ...
    'PsychoacousticAnnoyance_Widmann1992');
tc.verifyEmpty(step('PsychoacousticAnnoyance_Zwicker1999').field, 'the whole output');
tc.verifyEmpty(step('PsychoacousticAnnoyance_Zwicker1999').restat);
% the components come from the first model, with the time_skip of this run
for id = {'Loudness_ISO532_1', 'Sharpness_DIN45692', 'Roughness_Daniel1997', ...
        'FluctuationStrength_Osses2016'}
    tc.verifyEqual(step(id{1}).from, 'PsychoacousticAnnoyance_Widmann1992', id{1});
    tc.verifyEqual(step(id{1}).restat, params.(id{1}).time_skip, id{1});
end
tc.verifyEqual(step('Sharpness_DIN45692').attach, {'loudness', 'L'});
% only More 2010 and Di 2016 compute tonality, and with time_skip 0, which
% is the default of the metric
tc.verifyEqual(step('Tonality_Aures1985').from, 'PsychoacousticAnnoyance_More2010');
tc.verifyEqual(step('Tonality_Aures1985').field, 'K');
tc.verifyEmpty(step('Tonality_Aures1985').restat);
% the ECMA-418-2 metrics and EPNL share nothing
for id = {'Loudness_ECMA418_2', 'Roughness_ECMA418_2', 'Tonality_ECMA418_2', 'EPNL_FAR_Part36'}
    tc.verifyEmpty(step(id{1}).from, id{1});
end
end

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

%% Spectrogram, windows, filters and weighting -----------------------------

function test_window_shapes(tc)
n = 1024;
for name = {'hann', 'hamming', 'rect', 'blackmanharris'}
    w = SQAT_GUI_window(name{1}, n);
    tc.verifyTrue(iscolumn(w) && numel(w) == n, name{1});
end
% periodic windows: the first sample is the value at phase zero, the peak sits at n/2 + 1
w = SQAT_GUI_window('hann', n);            tc.verifyEqual(w(1), 0, 'AbsTol', 1e-15);
w = SQAT_GUI_window('hamming', n);         tc.verifyEqual(w(1), 0.08, 'AbsTol', 1e-12);
w = SQAT_GUI_window('blackmanharris', n);  tc.verifyEqual(w(1), 6.0e-5, 'AbsTol', 1e-9);
for name = {'hann', 'hamming', 'blackmanharris'}
    w = SQAT_GUI_window(name{1}, n);
    tc.verifyEqual(w(n/2 + 1), 1, 'AbsTol', 1e-12, name{1});
    tc.verifyEqual(w(2:n/2), flipud(w(n/2+2:n)), 'AbsTol', 1e-12, name{1});   % symmetric about the peak
end
tc.verifyEqual(SQAT_GUI_window('rect', n), ones(n, 1));
tc.verifyEqual(SQAT_GUI_window('Rectangular', 8), ones(8, 1));
% a custom window is resampled to the length asked for, and kept as it is at its own length
v = [0 1 3 2]';
tc.verifyEqual(SQAT_GUI_window('custom', 4, v), v);
tc.verifyEqual(SQAT_GUI_window('custom', 7, v), interp1(linspace(0, 1, 4), v, linspace(0, 1, 7))', 'AbsTol', 1e-12);
tc.verifyError(@() SQAT_GUI_window('kaiser', 8), 'SQAT_GUI:window');
tc.verifyError(@() SQAT_GUI_window('custom', 8), 'SQAT_GUI:window');
end

function test_read_window_takes_text_and_mat_files(tc)
d = tc.TestData.dir_tmp;
v = [0.1 0.5 1 0.5 0.1]';
writematrix(v, fullfile(d, 'w_col.csv'));
writematrix(v', fullfile(d, 'w_row.txt'), 'Delimiter', ' ');
w = v; save(fullfile(d, 'w_var.mat'), 'w');
for f = {'w_col.csv', 'w_row.txt', 'w_var.mat'}
    got = SQAT_GUI_read_window(fullfile(d, f{1}));
    tc.verifyEqual(got, v, 'AbsTol', 1e-12, f{1});
    tc.verifyTrue(iscolumn(got), f{1});
end
% what cannot be a window is refused
writematrix([1 2; 3 4], fullfile(d, 'w_matrix.csv'));
writematrix([1 NaN 2]', fullfile(d, 'w_nan.csv'));
writematrix(1, fullfile(d, 'w_one.csv'));
copyfile(fullfile(d, 'w_col.csv'), fullfile(d, 'w_col.xyz'));
for f = {'w_matrix.csv', 'w_nan.csv', 'w_one.csv', 'w_col.xyz'}
    tc.verifyError(@() SQAT_GUI_read_window(fullfile(d, f{1})), 'SQAT_GUI:window', f{1});
end
end

function test_spectrogram_reads_the_level_of_a_tone_with_any_window(tc)
fs = 48000; t = (0:fs-1)'/fs;
x = sqrt(2)*2e-5*10^(60/20) * sin(2*pi*1500*t);      % on a bin of a 1024 point FFT (32 x 46.875 Hz)
for name = {'hann', 'hamming', 'rect', 'blackmanharris'}
    [tt, f, L, info] = SQAT_GUI_spectrogram(x, fs, name{1}, 10, 50);
    tc.verifyEqual(max(L(:)), 60, 'AbsTol', 0.05, name{1});
    [~, i] = max(max(L, [], 2));
    tc.verifyEqual(f(i), 1500, name{1});
    tc.verifyEqual(size(L), [513, numel(tt)], name{1});
    tc.verifyEqual(info.n_fft, 1024);
    tc.verifyEqual(info.overlap, 50);
    tc.verifyFalse(info.limited);
end
end

function test_spectrogram_follows_degree_and_overlap(tc)
fs = 48000; x = randn(fs, 1) * 0.01;
[~, f, L4096] = SQAT_GUI_spectrogram(x, fs, 'hann', 12, 50);
tc.verifyEqual(numel(f), 2049);
tc.verifyEqual(f(2) - f(1), fs/4096, 'AbsTol', 1e-9);
[t50, ~, ~] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 50);
[t75, ~, ~] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 75);
tc.verifyEqual(numel(t50), floor((fs - 1024)/512) + 1);
tc.verifyEqual(numel(t75), floor((fs - 1024)/256) + 1);
% a window given as a vector is used as it is at the FFT length, and resampled to it otherwise
w = SQAT_GUI_window('hann', 1024);
[~, ~, Lc] = SQAT_GUI_spectrogram(x, fs, w, 10, 50);
[~, ~, Lh] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 50);
tc.verifyEqual(Lc, Lh);
v = SQAT_GUI_window('hann', 200);
[~, ~, Lr] = SQAT_GUI_spectrogram(x, fs, v, 10, 50);
[~, ~, Le] = SQAT_GUI_spectrogram(x, fs, SQAT_GUI_window('custom', 1024, v), 10, 50);
tc.verifyEqual(Lr, Le);
% a signal shorter than the FFT is padded
[t, ~, Ls] = SQAT_GUI_spectrogram(x(1:300), fs, 'hann', 10, 50);
tc.verifyEqual(size(Ls, 2), numel(t));
tc.verifyGreaterThanOrEqual(numel(t), 1);
end

function test_spectrogram_limits_the_frames_and_says_so(tc)
fs = 8000; x = randn(200*fs, 1) * 0.01;
[t, ~, L, info] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 90);
tc.verifyTrue(info.limited);
tc.verifyLessThan(info.overlap, 90);
tc.verifyLessThanOrEqual(numel(t), 4000);
tc.verifyLessThanOrEqual(numel(L), 8e6);
% the largest FFT still fits in memory
[t, f, L, info] = SQAT_GUI_spectrogram(x, fs, 'hann', 16, 95);
tc.verifyLessThanOrEqual(numel(L), 8e6);
tc.verifyTrue(info.limited);
end

function test_spectral_filter_removes_a_box_and_only_that(tc)
fs = 48000; t = (0:3*fs-1)'/fs;
x = sin(2*pi*500*t) + sin(2*pi*2000*t);
% no box, or a box beyond the spectrum: the signal comes back
tc.verifyEqual(SQAT_GUI_spectral_filter(x, fs, zeros(0, 4)), x);
tc.verifyEqual(SQAT_GUI_spectral_filter(x, fs, [0 3 30000 40000]), x, 'AbsTol', 1e-12);
% a box around 2 kHz for the whole signal
y = SQAT_GUI_spectral_filter(x, fs, [0 3 1500 2500]);
tc.verifyTrue(iscolumn(y) && numel(y) == numel(x));
mid = fs:2*fs;                                        % away from the edges
tc.verifyLessThan(il_tone_db(y(mid), fs, 2000) - il_tone_db(x(mid), fs, 2000), -40);
tc.verifyEqual(il_tone_db(y(mid), fs, 500), il_tone_db(x(mid), fs, 500), 'AbsTol', 0.1);
% a box that lasts one second only removes the tone there
y = SQAT_GUI_spectral_filter(x, fs, [1 2 1500 2500]);
inside = round(1.2*fs):round(1.8*fs); before = 1:round(0.8*fs); after = round(2.2*fs):numel(x);
tc.verifyLessThan(il_tone_db(y(inside), fs, 2000) - il_tone_db(x(inside), fs, 2000), -40);
tc.verifyEqual(il_tone_db(y(before), fs, 2000), il_tone_db(x(before), fs, 2000), 'AbsTol', 0.1);
tc.verifyEqual(il_tone_db(y(after), fs, 2000), il_tone_db(x(after), fs, 2000), 'AbsTol', 0.1);
% the box can be kept and the rest removed
y = SQAT_GUI_spectral_filter(x, fs, [0 3 1500 2500], 'keep');
tc.verifyLessThan(il_tone_db(y(mid), fs, 500) - il_tone_db(x(mid), fs, 500), -40);
tc.verifyEqual(il_tone_db(y(mid), fs, 2000), il_tone_db(x(mid), fs, 2000), 'AbsTol', 0.1);
y = SQAT_GUI_spectral_filter(x, fs, [1 2 1500 2500], 'keep');       % outside the box's time: silence
tc.verifyLessThan(rms(y(1:round(0.8*fs))), 1e-6);
tc.verifyLessThan(rms(y(round(2.2*fs):end)), 1e-6);
tc.verifyEqual(il_tone_db(y(inside), fs, 2000), il_tone_db(x(inside), fs, 2000), 'AbsTol', 0.1);
tc.verifyLessThan(il_tone_db(y(inside), fs, 500) - il_tone_db(x(inside), fs, 500), -40);
tc.verifyEqual(SQAT_GUI_spectral_filter(x, fs, zeros(0, 4), 'keep'), x);   % no box: as it was
tc.verifyError(@() SQAT_GUI_spectral_filter(x, fs, [0 3 1500 2500], 'other'), 'SQAT_GUI:filter');
% two boxes
y = SQAT_GUI_spectral_filter(x, fs, [0 3 1500 2500; 0 3 400 600]);
tc.verifyLessThan(il_tone_db(y(mid), fs, 500) - il_tone_db(x(mid), fs, 500), -40);
tc.verifyLessThan(il_tone_db(y(mid), fs, 2000) - il_tone_db(x(mid), fs, 2000), -40);
end

function test_weighting_follows_iec_61672_at_the_reference_points(tc)
fs = 48000; t = (0:2*fs-1)'/fs;
x100 = sin(2*pi*100*t); x1k = sin(2*pi*1000*t);
tc.verifyEqual(SQAT_GUI_weight(x100, fs, 'Z'), x100);
keep = fs/2:2*fs;                                     % after the start-up of the filter
db = @(y, x) 20*log10(rms(y(keep)) / rms(x(keep)));
tc.verifyEqual(db(SQAT_GUI_weight(x100, fs, 'A'), x100), -19.1, 'AbsTol', 0.2);
tc.verifyEqual(db(SQAT_GUI_weight(x100, fs, 'C'), x100), -0.3, 'AbsTol', 0.2);
tc.verifyEqual(db(SQAT_GUI_weight(x1k, fs, 'A'), x1k), 0, 'AbsTol', 0.1);
tc.verifyEqual(db(SQAT_GUI_weight(x1k, fs, 'C'), x1k), 0, 'AbsTol', 0.1);
% the curve drawn on the spectrogram is the response of the same filter
tc.verifyEqual(SQAT_GUI_weight_curve([100 1000], fs, 'A'), [-19.1 0], 'AbsTol', 0.2);
tc.verifyEqual(SQAT_GUI_weight_curve([100 1000], fs, 'C'), [-0.3 0], 'AbsTol', 0.2);
tc.verifyEqual(SQAT_GUI_weight_curve([100 1000], fs, 'Z'), [0 0]);
tc.verifyTrue(iscolumn(SQAT_GUI_weight_curve([100; 1000], fs, 'A')));
tc.verifyError(@() SQAT_GUI_weight(x1k, fs, 'B'), 'SQAT_GUI:weighting');
end

function test_weighting_filters_equal_those_of_the_sound_level_meter(tc)
% SQAT_GUI_weight_filter rewrites the design of Gen_weighting_filters without the toolbox
% function it calls (bilinear); the coefficients must stay the same
for fs = [44100 48000 96000]
    for type = {'A', 'C'}
        [b, a] = SQAT_GUI_weight_filter(fs, type{1});
        [b0, a0] = Gen_weighting_filters(fs, type{1});
        lbl = sprintf('%s at %d Hz', type{1}, fs);
        tc.verifyEqual(b(:), b0(:), 'RelTol', 1e-9, 'AbsTol', 1e-12, lbl);
        tc.verifyEqual(a(:), a0(:), 'RelTol', 1e-9, 'AbsTol', 1e-12, lbl);
    end
end
[b, a] = SQAT_GUI_weight_filter(48000, 'Z');
tc.verifyEqual([b a], [1 1]);
tc.verifyError(@() SQAT_GUI_weight_filter(48000, 'D'), 'SQAT_GUI:weighting');
end

function test_gui_helpers_need_no_toolbox(tc)
% besides MATLAB itself: the files of the waveform window and the extraction
names = {'SQAT_GUI_window', 'SQAT_GUI_read_window', 'SQAT_GUI_spectrogram', ...
    'SQAT_GUI_spectral_filter', 'SQAT_GUI_weight', 'SQAT_GUI_weight_curve', ...
    'SQAT_GUI_weight_filter', 'SQAT_GUI_extract', 'SQAT_GUI_load', 'SQAT_GUI_single_values'};
for k = 1:numel(names)
    [~, prods] = matlab.codetools.requiredFilesAndProducts(which(names{k}));
    tc.verifyEmpty(setdiff({prods.Name}, {'MATLAB'}), [names{k} ' needs a toolbox']);
end
end

function test_enhanced_stft_settings(tc)
% the two smoothings of the function: both put a 1 kHz tone at 1 kHz, and the sharp one gives the thinner line
[x, fs] = audioread(tc.TestData.wav_tone);
w = zeros(1, 2);
for k = 1:2
    modes = {'readable', 'sharp'};
    [~, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, modes{k});
    tc.verifyEqual(info.smoothing, modes{k});
    p = max(L, [], 2);
    [~, i] = max(p);
    tc.verifyEqual(f(i), 1000, 'AbsTol', 10);
    w(k) = nnz(p >= max(p) - 3);                      % rows within 3 dB of the peak
end
tc.verifyLessThan(w(2), w(1));
tc.verifyError(@() SQAT_GUI_enhanced_stft(x, fs, 'other'), 'SQAT_GUI_enhanced_stft:smoothing');
end

function test_enhanced_stft_preview_steps_by_one_column(tc)
% the preview of a long signal: frames step by up to one output column; with 1 ms columns nothing changes
[x, fs] = audioread(tc.TestData.wav_tone);
x = x(1:round(1.5*fs));                                        % 1.5 s in 2000 columns: 1 ms per column
[t, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
[t_p, f_p, L_p, info] = SQAT_GUI_enhanced_stft(x, fs, 'readable', [], [], [], true);
tc.verifyTrue(info.preview);
tc.verifyEqual({t_p, f_p, L_p}, {t, f, L});
[t, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable', 50);    % 30 ms columns, as in a long signal
[t_p, f_p, L_p] = SQAT_GUI_enhanced_stft(x, fs, 'readable', 50, [], [], true);
tc.verifyEqual({t_p, f_p}, {t, f});
tc.verifyNotEqual(L_p, L);
[~, i] = max(max(L_p, [], 2));
tc.verifyEqual(f_p(i), 1000, 'AbsTol', 10);
end

%% Main window -------------------------------------------------------------

function test_gui_opens_with_the_expected_controls(tc)
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyClass(fig, 'matlab.ui.Figure');
for tag = {'logo','file_count','load_files','signals_list','theme', ...
           'add_analysis','analysis_list','run','stop_run','open_graphs','open_waveform','export', ...
           'show_plots','save_figures','split_figures','figures_folder', ...
           'console','results_table','status','progress'}
    tc.verifyNotEmpty(findobj(fig, 'Tag', tag{1}), ['missing control: ' tag{1}]);
end
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, 'No files loaded');
tc.verifyEqual(findobj(fig, 'Tag', 'status').Text, 'Ready');
tc.verifyEmpty(findobj(fig, 'Tag', 'sqat_figures'), 'the graphs window shows the SQAT figures');
dd = findobj(fig, 'Tag', 'analysis_metric_1');
m = SQAT_GUI_metrics;
tc.verifyEqual(dd.ItemsData, {m.id});
tc.verifyEqual(dd.Value, 'Loudness_ISO532_1');                 % one analysis at start
tc.verifyEmpty(findobj(fig, 'Tag', 'analysis_metric_2'));
tc.verifyEmpty(findobj(fig, 'Tag', 'channel'), 'the channel is set per signal');
tc.verifyEmpty(findobj(fig, 'Tag', 'dbfs'), 'the dBFS is set per signal');
end

function test_gui_called_without_output_keeps_working(tc)
% The usual call, SQAT_GUI at the prompt, returns nothing: the callbacks must still work.
evalc('SQAT_GUI({tc.TestData.wav_mono}, ''Visible'', ''off'')');
fig = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI');
tc.assertNumElements(fig, 1);
tc.addTeardown(@() delete(fig));
if il_has_theme()
    il_press(fig, 'theme');
    tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'light');
end
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
il_press(fig, 'open_graphs');
tc.verifyNumElements(il_window('SQAT_GUI_graphs'), 1);
il_press(fig, 'open_waveform');
tc.verifyNumElements(il_window('SQAT_GUI_waveform'), 1);
end

function test_gui_windows_skip_a_user_default_CreateFcn(tc)
% A startup.m may set a default CreateFcn that fails on uifigures (e.g.
% addToolbarExplorationButtons in R2024b); the GUI windows must not run it.
old = get(groot, 'defaultFigureCreateFcn');
tc.addTeardown(@() set(groot, 'defaultFigureCreateFcn', old));
setappdata(groot, 'sqat_gui_created', {});
set(groot, 'defaultFigureCreateFcn', @(f, ~) setappdata(groot, 'sqat_gui_created', ...
    [getappdata(groot, 'sqat_gui_created'), {f}]));
tc.addTeardown(@() rmappdata(groot, 'sqat_gui_created'));
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
il_press(fig, 'open_waveform');
wins = [fig, il_window('SQAT_GUI_graphs'), il_window('SQAT_GUI_waveform')];
tc.assertNumElements(wins, 3);
created = getappdata(groot, 'sqat_gui_created');
for k = 1:numel(wins)
    tc.verifyFalse(any(cellfun(@(f) isequal(f, wins(k)), created)), ...
        ['the default CreateFcn ran on ' wins(k).Tag]);
end
end

function test_gui_play_and_stop_never_break_the_interface(tc)
% Audio output may be missing (batch mode); the buttons must then only log it.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_press(w, 'play');
il_press(w, 'stop');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyTrue(contains(log, 'Playing') || contains(log, 'Audio output unavailable'));
end

function test_gui_loads_files_and_draws_the_waveform(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(il_signal_names(fig), {'tone_mono.wav', 'tone_stereo.wav'});
tc.verifyTrue(findobj(fig, 'Tag', 'signal_tick_1').Value);   % new signals come ticked
tc.verifyTrue(findobj(fig, 'Tag', 'signal_tick_2').Value);
tc.verifyNotEmpty(findobj(fig, 'Tag', 'signal_remove_2'));
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, '2 files loaded');
% each signal has its channel (All by default for a stereo file) and its calibration
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Value, '1');
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_2').Value, 'All');
c = findobj(fig, 'Tag', 'signal_cal_1');
tc.verifyEqual(c.Value, 94);
tc.verifyEqual(c.FontAngle, 'italic');                      % the default is marked
il_signal_dbfs(fig, 1, 94);
tc.verifyEqual(c.FontAngle, 'normal');                      % set by hand, even to the same value
end

function test_gui_each_signal_offers_only_its_own_channels(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Items, {'1'});
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_2').Items, {'1', '2', 'All'});
end

function test_gui_signal_list_keeps_the_x_in_view_for_long_names(tc)
% The name takes the space that is left, so a long name cannot push the x out.
long = fullfile(tc.TestData.dir_tmp, [repmat('a_very_long_signal_name_', 1, 5) '.wav']);
copyfile(tc.TestData.wav_mono, long);
fig = SQAT_GUI({long}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
g = findobj(fig, 'Tag', 'signals_list');
tc.verifyEqual(g.ColumnWidth{3}, '1x');
tc.verifyTrue(all(cellfun(@isnumeric, g.ColumnWidth([1 2 4 5 6]))));
tc.verifyLessThanOrEqual(sum([g.ColumnWidth{[1 2 4 5 6]}]), 250);   % the name keeps most of the 430 px
end

function test_gui_channels_offer_all_only_with_a_stereo_file(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(findobj(fig, 'Tag', 'signal_channel_1').Items, {'1'});
il_remove_signal(fig, 1);
tc.verifyEmpty(il_signal_names(fig));
end

function test_gui_signal_list_ticks_and_removes_signals(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_mark_signal(fig, 2, false);                       % only the first is analysed
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_mono.wav'});
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, '1 file(s)');
il_mark_signal(fig, 2, true);
il_mark_signal(fig, 1, false);
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'});
% nothing ticked: a warning and no run
il_mark_signal(fig, 2, false);
il_press(fig, 'run');
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No signals ticked');
% removing a signal takes its results away
il_mark_signal(fig, 1, true); il_mark_signal(fig, 2, true);
il_press(fig, 'run');
il_remove_signal(fig, 1);
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'});
tc.verifyEqual(il_signal_names(fig), {'tone_1k_60dB.wav'});
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, '1 file loaded');
il_remove_signal(fig, 1);
tc.verifyEqual(findobj(fig, 'Tag', 'file_count').Text, 'No files loaded');
tc.verifyEmpty(findobj(fig, 'Tag', 'results_table').Data);
il_press(fig, 'run');                                % the interface goes on
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No files loaded');
end

function test_gui_signals_carry_a_number_and_the_plots_use_it_as_a_tag(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(il_signal_numbers(fig), {'#1', '#2'});
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'loudness');
lg = findobj(g, 'Type', 'legend');
tc.assertNumElements(lg, 1);
tc.verifyEqual(lg.String, {'Signal #1, ch1', 'Signal #2, ch1'});
il_set(g, 'graph_analysis', 'stats');
tc.verifyEqual(findobj(g, 'Tag', 'graph_stats').ColumnName(:)', ...
    {'Quantity', 'Signal #1, ch1', 'Signal #2, ch1'});
% a number is not given twice: the signal that is left keeps its own
il_remove_signal(fig, 1);
tc.verifyEqual(il_signal_numbers(fig), {'#2'});
end

function test_gui_theme_switches_between_dark_and_light(tc)
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
logo = findobj(fig, 'Tag', 'logo');
b = findobj(fig, 'Tag', 'theme');
if ~il_has_theme()                                  % before R2025a: the light look, and no toggle
    tc.verifyTrue(endsWith(logo.ImageSource, 'logo.png'));
    tc.verifyTrue(isfile(logo.ImageSource));
    tc.verifyEqual(char(b.Enable), 'off');
    tc.verifySubstring(b.Tooltip, 'R2025a');
    return
end
tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'dark');         % starts dark, like pySQAT
tc.verifyTrue(endsWith(logo.ImageSource, 'logo_white.png'));
tc.verifyTrue(isfile(logo.ImageSource));
tc.verifyEqual(b.Text, char(9788));                 % a sun
tc.verifySubstring(b.Tooltip, 'Light');
il_press(fig, 'theme');
tc.verifyEqual(char(fig.Theme.BaseColorStyle), 'light');
tc.verifyTrue(endsWith(logo.ImageSource, 'logo.png'));
tc.verifyTrue(isfile(logo.ImageSource));
tc.verifyEqual(b.Text, char(9790));                 % a moon
tc.verifySubstring(b.Tooltip, 'Dark');
end

function test_gui_status_bar_follows_the_run(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
tc.verifyEqual(findobj(fig, 'Tag', 'progress').Value, 0);
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
tc.verifyEqual(findobj(fig, 'Tag', 'progress').Value, 100);
end

function test_startup_puts_the_gui_on_the_path(tc)
rmpath(fullfile(basepath_SQAT, 'gui'));
tc.addTeardown(@() addpath(fullfile(basepath_SQAT, 'gui')));
evalc('run(fullfile(basepath_SQAT, ''startup_SQAT.m''))');
tc.verifyNotEmpty(which('SQAT_GUI'));
end

%% Running the analyses ----------------------------------------------------

function test_gui_run_shows_a_progress_dialog(tc)
% the dialog needs a visible window, so this one is shown
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'on');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Roughness_Daniel1997'});
% the run gives the timer a turn between steps: record what the dialog shows then
tm = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.05, 'TimerFcn', @(~, ~) record());
tc.addTeardown(@() delete(tm));
setappdata(groot, 'sqat_gui_dialogs', {});
tc.addTeardown(@() rmappdata(groot, 'sqat_gui_dialogs'));
    function record()
        if isappdata(fig, 'sqat_progress')
            d = getappdata(fig, 'sqat_progress');
            setappdata(groot, 'sqat_gui_dialogs', [getappdata(groot, 'sqat_gui_dialogs'), ...
                {struct('value', d.Value, 'message', d.Message, 'cancelable', d.Cancelable)}]);
        end
    end
start(tm);
il_press(fig, 'run');
stop(tm);
seen = getappdata(groot, 'sqat_gui_dialogs');
tc.assertNotEmpty(seen, 'no progress dialog during the run');
tc.verifyEqual(seen{1}.cancelable, matlab.lang.OnOffSwitchState.on);
tc.verifyTrue(any(cellfun(@(r) contains(r.message, 'Roughness_Daniel1997') || contains(r.message, 'Loudness'), seen)));
tc.verifyGreaterThan(max(cellfun(@(r) r.value, seen)), 0);
tc.verifyFalse(isappdata(fig, 'sqat_progress'), 'the dialog stays after the run');
tc.verifySubstring(findobj(fig, 'Tag', 'status').Text, 'Done');
end

function test_gui_results_have_one_tab_per_signal(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tabs = @() findobj(fig, 'Tag', 'results_signal');         % in the order of the tab group
tc.verifyEqual({tabs().Title}, {'Results #1', 'Results #2'});
names = {'tone_mono.wav', 'tone_1k_60dB.wav'};
for k = 1:2                                          % each tab holds the rows of its signal
    t = tabs();
    D = findobj(t(k), 'Type', 'uitable').Data;
    tc.verifyEqual(D, T(strcmp(T.File, names{k}), {'Analysis', 'Metric', 'Channel', 'Quantity', 'Value', 'Unit', 'Parameters'}));
end
il_remove_signal(fig, 1);                            % a removed signal takes its tab away
tc.verifyEqual({tabs().Title}, {'Results #2'});
il_remove_signal(fig, 1);
tc.verifyEmpty(tabs());
end

function test_gui_results_put_the_signals_and_channels_side_by_side(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Roughness_Daniel1997'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
% the first quantity of #1: the mono file, then channels 1 and 2 of the stereo one
tc.verifyEqual(T.File(1:3), {'tone_mono.wav'; 'tone_stereo.wav'; 'tone_stereo.wav'});
tc.verifyEqual(T.Channel(1:3), {'1'; '1'; '2'});
tc.verifyEqual(T.Quantity(2:3), repmat(T.Quantity(1), 2, 1));
tc.verifyEqual(T.Signal(1:3), {'#1'; '#2'; '#2'});
tc.verifyEqual(unique(T.Unit(strcmp(T.Quantity, 'N5'))), {'sone'});
tc.verifyEqual(unique(T.Unit(strcmp(T.Metric, 'Roughness_Daniel1997'))), {'asper'});
tc.verifyNotEqual(T.Quantity{4}, T.Quantity{1});
% every analysis comes whole, #1 before #2
a = str2double(erase(T.Analysis, '#'));
tc.verifyTrue(issorted(a));
% the tab of a signal keeps its own rows, channels interleaved
tabs = findobj(fig, 'Tag', 'results_signal');
D = findobj(tabs(2), 'Type', 'uitable').Data;
tc.verifyEqual(D.Channel(1:4), {'1'; '2'; '1'; '2'});
end

function test_gui_all_channels_runs_a_binaural_pair_in_one_call(tc)
fig = SQAT_GUI({tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1', 'Loudness_ECMA418_2'});
il_signal_channel(fig, 1, 'All');
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
iso = T(strcmp(T.Metric, 'Loudness_ISO532_1'), :);
ecma = T(strcmp(T.Metric, 'Loudness_ECMA418_2'), :);
tc.verifyEqual(unique(iso.Channel, 'stable'), {'1'; '2'});
% the channels are interleaved: each quantity of channel 1, then of channel 2
tc.verifyEqual(iso.Channel(1:4), {'1'; '2'; '1'; '2'});
tc.verifyEqual(iso.Quantity(1), iso.Quantity(2));
tc.verifyEqual(iso.Quantity(3), iso.Quantity(4));
tc.verifyEqual(ecma.Channel(1:3), {'1'; '2'; 'Binaural'});
tc.verifyEqual(unique(ecma.Channel, 'stable'), {'1'; '2'; 'Binaural'});
% the values are the ones of the direct calls
ref = audioread(tc.TestData.wav_stereo); fs = tc.TestData.fs;
[~, r2] = evalc('Loudness_ISO532_1(ref(:, 2), fs, 0, 2, 0.5, false)');
row = strcmp(iso.Channel, '2') & strcmp(iso.Quantity, 'N5');
tc.verifyEqual(iso.Value(row), r2.N5);
[~, rb] = evalc('Loudness_ECMA418_2(ref, fs, ''free-frontal'', 0.304, false)');
for c = {'1', '2'}
    row = strcmp(ecma.Channel, c{1}) & strcmp(ecma.Quantity, 'Nmean');
    tc.verifyEqual(ecma.Value(row), rb.Nmean(str2double(c{1})));
end
row = strcmp(ecma.Channel, 'Binaural') & strcmp(ecma.Quantity, 'Nmean');
tc.verifyEqual(ecma.Value(row), rb.Nmean(3));
% the pair was analysed in one call
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyEqual(numel(strfind(log, 'Running Loudness_ECMA418_2 on tone_stereo.wav')), 1);
tc.verifySubstring(log, 'Running Loudness_ECMA418_2 on tone_stereo.wav, both channels');
tc.verifyEqual(numel(strfind(log, 'Running Loudness_ISO532_1 on tone_stereo.wav')), 2);
end

function test_gui_edits_the_parameters_of_a_metric(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
pw = il_open_params(fig, 1);
tc.verifyEqual(findobj(pw, 'Tag', 'unit_time_skip').Text, 's');
tc.verifyEqual(findobj(pw, 'Tag', 'param_time_skip').Limits, [0 Inf]);
c = findobj(pw, 'Tag', 'param_method');
tc.verifyEqual(c.Value, 2);
c.Value = 1; c.ValueChangedFcn(c, []);          % stationary
t = findobj(pw, 'Tag', 'param_time_skip'); t.Value = 1; t.ValueChangedFcn(t, []);
b = findobj(pw, 'Tag', 'params_reset'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(c.Value, 2);
tc.verifyEqual(t.Value, 0.5);
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'FF tv t 0.5s');
c.Value = 1; c.ValueChangedFcn(c, []);          % stationary again, for the run
il_press(fig, 'add_analysis');                  % the window stays open while the list changes
tc.verifyTrue(isvalid(pw));
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyTrue(any(strcmp(T.Quantity, 'Loudness')));
tc.verifyFalse(any(strcmp(T.Quantity, 'N5')));
end

function test_gui_runs_several_metrics_on_several_files(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1','Roughness_Daniel1997'});
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_mono.wav'; 'tone_stereo.wav'});
tc.verifyEqual(unique(T.Metric), {'Loudness_ISO532_1'; 'Roughness_Daniel1997'});
% the value shown is the value of the metric
x_file = audioread(tc.TestData.wav_mono);          % what the interface reads (float32 file)
[~, ref] = evalc('Loudness_ISO532_1(x_file, tc.TestData.fs, 0, 2, 0.5, false)');
row = strcmp(T.File, 'tone_mono.wav') & strcmp(T.Metric, 'Loudness_ISO532_1') & strcmp(T.Quantity, 'N5');
tc.verifyEqual(T.Value(row), ref.N5);
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'Done');
end

function test_gui_uses_the_dbfs_and_channel_of_each_signal(tc)
fig = SQAT_GUI({tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_signal_channel(fig, 1, '1');
il_press(fig, 'run');
n94 = il_value(fig, 'tone_stereo.wav', 'Loudness_ISO532_1', 'N5');
il_signal_dbfs(fig, 1, 100);
il_press(fig, 'run');
n100 = il_value(fig, 'tone_stereo.wav', 'Loudness_ISO532_1', 'N5');
tc.verifyGreaterThan(n100, n94);
il_signal_dbfs(fig, 1, 94);
il_signal_channel(fig, 1, '2');
il_press(fig, 'run');
n_ch2 = il_value(fig, 'tone_stereo.wav', 'Loudness_ISO532_1', 'N5');
tc.verifyLessThan(n_ch2, n94);                  % channel 2 is 10 dB lower
end

function test_gui_analyses_the_active_file_first(tc)
% The file shown in the interface is the first one analysed, so its result
% is the first one on screen.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_activate_signal(fig, 2);
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
first = strfind(log, 'Running Loudness_ISO532_1 on tone_1k_60dB.wav');
second = strfind(log, 'Running Loudness_ISO532_1 on tone_mono.wav');
tc.assertNotEmpty(first);
tc.assertNotEmpty(second);
tc.verifyLessThan(first(1), second(1), 'the active file was not analysed first');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.File), {'tone_1k_60dB.wav'; 'tone_mono.wav'});
end

function test_gui_exports_the_results_to_excel(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
xlsx = fullfile(tc.TestData.dir_tmp, 'results.xlsx');
describe = getappdata(fig, 'sqat_run_description');
SQAT_GUI_export(findobj(fig, 'Tag', 'results_table').Data, xlsx, describe());
T = findobj(fig, 'Tag', 'results_table').Data;
R = readtable(xlsx, 'Sheet', 'Results', 'TextType', 'char');
tc.verifyEqual(R.Quantity, T.Quantity);
tc.verifyEqual(R.Value, T.Value, 'AbsTol', 1e-12);
% what is needed to reproduce a value travels with it
tc.verifyEqual(unique(R.Unit), {'asper'});
tc.verifyEqual(unique(R.Cal_dB_SPL), 94);
tc.verifyEqual(unique(R.Parameters), {'Time skip: 0 s'});
tc.verifyEqual(unique(R.Path), {tc.TestData.wav_mono});
S = readtable(xlsx, 'Sheet', 'Settings', 'TextType', 'char');
tc.verifyTrue(any(startsWith(S.Item, 'Signal #1')));
v = S.Value(strcmp(S.Item, 'Signal #1'));
tc.verifySubstring(v{1}, '(default)');
v = S.Value(strcmp(S.Item, 'Analysis #1'));
tc.verifySubstring(v{1}, 'Roughness_Daniel1997');
tc.verifyTrue(any(strcmp(S.Item, 'SQAT version')));
end

function test_gui_marks_the_results_when_a_setting_changes(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tab = findobj(fig, 'Type', 'uitab', 'Title', 'Results');
tc.assertNumElements(tab, 1);
pw = il_open_params(fig, 1);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
tc.verifySubstring(tab.Title, 'settings changed');
il_press(fig, 'run');
tc.verifyEqual(tab.Title, 'Results');
il_signal_dbfs(fig, 1, 100);                         % the calibration of a signal with results
tc.verifySubstring(tab.Title, 'settings changed');
end

function test_gui_reports_a_failing_metric_and_goes_on(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1','Roughness_Daniel1997'});
pw = il_open_params(fig, 1);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
c = findobj(pw, 'Tag', 'param_time_skip'); c.Value = 10; c.ValueChangedFcn(c, []);  % longer than the signal: the metric raises an error
il_press(fig, 'run');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifySubstring(log, 'ERROR');
tc.verifySubstring(log, 'Loudness_ISO532_1');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.Metric), {'Roughness_Daniel1997'});
end

function test_gui_run_without_files_or_metrics_only_warns(tc)
fig = SQAT_GUI({}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'run');
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'No files');
fig2 = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig2));
il_select_metrics(fig2, {});
% the summary of the ECMA-418-2 parameters
il_select_metrics(fig2, {'Loudness_ECMA418_2'});
tc.verifyEqual(findobj(fig2, 'Tag', 'analysis_summary_1').Text, 'FFront t 0.304s');
il_select_metrics(fig2, {});
il_press(fig2, 'run');
tc.verifySubstring(strjoin(findobj(fig2, 'Tag', 'console').Value, newline), 'No metrics');
end

function test_gui_compares_one_metric_with_two_sets_of_parameters(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'add_analysis');                       % a copy of the first: Loudness #2
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_metric_2').Value, 'Loudness_ISO532_1');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_1').Text, '#1');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_2').Text, '#2');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'FF tv t 0.5s');
pw = il_open_params(fig, 2);
c = findobj(pw, 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);   % stationary
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_2').Text, 'FF stat t 0.5s');
il_press(fig, 'run');
T = findobj(fig, 'Tag', 'results_table').Data;
tc.verifyEqual(unique(T.Analysis, 'stable'), {'#1'; '#2'});
tc.verifyEqual(unique(T.Metric), {'Loudness_ISO532_1'});
tc.verifyTrue(any(strcmp(T.Quantity(strcmp(T.Analysis, '#1')), 'N5')));
tc.verifyTrue(any(strcmp(T.Quantity(strcmp(T.Analysis, '#2')), 'Loudness')));
il_press(fig, 'open_graphs');                        % both are on offer in the graphs
g = il_window('SQAT_GUI_graphs');
dm = findobj(g, 'Tag', 'graph_metric');
tc.verifyEqual(dm.Items, {'#1 Loudness (ISO 532-1)', '#2 Loudness (ISO 532-1)'});
il_set(g, 'graph_metric', 'Loudness_ISO532_1#2');
tc.verifyNotEmpty(findobj(g, 'Type', 'axes'));
% removing the first leaves the second with its number: #2 keeps meaning the same analysis
b = findobj(fig, 'Tag', 'analysis_remove_1'); b.ButtonPushedFcn(b, []);
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_1').Text, '#2');
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_summary_1').Text, 'FF stat t 0.5s');
tc.verifyEmpty(findobj(fig, 'Tag', 'analysis_metric_2'));
il_press(fig, 'add_analysis');                       % a number is not given twice
tc.verifyEqual(findobj(fig, 'Tag', 'analysis_number_2').Text, '#3');
end

function test_gui_saves_the_figures_of_the_metrics(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
out_dir = fullfile(tc.TestData.dir_tmp, 'figs'); mkdir(out_dir);
h = findobj(fig, 'Tag', 'save_figures');   h.Value = true;
h = findobj(fig, 'Tag', 'figures_folder'); h.Value = out_dir;
il_press(fig, 'run');
png = dir(fullfile(out_dir, '*.png'));
tc.verifyNotEmpty(png);
tc.verifyTrue(all(startsWith({png.name}, 'tone_mono_s1_Roughness_Daniel1997_a1')));
tc.verifyEmpty(il_sqat_figures());                % saved, then closed
il_press(fig, 'run');                             % a second run keeps the first files
tc.verifyNumElements(dir(fullfile(out_dir, '*.png')), 2 * numel(png));
end

function test_gui_split_figures_saves_one_file_per_axes(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
out_dir = fullfile(tc.TestData.dir_tmp, 'figs_split'); mkdir(out_dir);
h = findobj(fig, 'Tag', 'save_figures');   h.Value = true;
h = findobj(fig, 'Tag', 'split_figures');  h.Value = true;
h = findobj(fig, 'Tag', 'figures_folder'); h.Value = out_dir;
il_press(fig, 'run');
figs = il_all_sqat_figures();                      % drawn to be saved, kept hidden
tc.assertNotEmpty(figs);
n_axes = numel(findobj(figs, 'Type', 'axes'));
png = dir(fullfile(out_dir, '*.png'));
tc.verifyEqual(numel(png), n_axes);
tc.verifyGreaterThanOrEqual(n_axes, numel(figs));
for ax = findobj(figs, 'Type', 'axes')'
    tc.verifyEqual(ax.Colormap, load('cmap_inferno.txt'));
end
end

%% Graphs windows ----------------------------------------------------------

function test_gui_plots_the_series_of_the_active_file(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1','Roughness_Daniel1997'});
il_press(fig, 'run');
tc.verifyEmpty(il_window('SQAT_GUI_graphs'));      % opens on request only
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
tc.verifyNumElements(g, 1);
pm = findobj(g, 'Tag', 'graph_metric');
tc.verifyEqual(pm.ItemsData, {'Loudness_ISO532_1','Roughness_Daniel1997'});
il_set(g, 'graph_metric', 'Roughness_Daniel1997');
% the window holds the figure that Roughness_Daniel1997 draws itself
x_file = audioread(tc.TestData.wav_mono);
before = findall(groot, 'Type', 'figure');
[~, ref] = evalc('Roughness_Daniel1997(x_file, tc.TestData.fs, 0, true)');
after = findall(groot, 'Type', 'figure');
ref_fig = after(~ismember(after, before));
tc.verifyEqual(numel(findobj(g, 'Type', 'axes')), numel(findobj(ref_fig, 'Type', 'axes')));
ln = findobj(g, 'Type', 'line');
tc.verifyTrue(any(arrayfun(@(l) isequal(l.YData(:), ref.InstantaneousRoughness(:)), ln)), ...
    'no line of the window holds the instantaneous roughness');
% with the inferno colour scale of the toolbox
for ax = findobj(g, 'Type', 'axes')'
    tc.verifyEqual(ax.Colormap, load('cmap_inferno.txt'));
end
end

function test_gui_draws_the_figure_of_the_active_file_in_the_analysis(tc)
% The analysis draws the figure of the active file in the same call, so the
% metric runs once and the graphs window only copies it. The figure carries
% the settings of the analysis.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tc.verifyNumElements(il_all_sqat_figures(), 1, 'the analysis kept no figure to plot');
tc.verifyEmpty(il_sqat_figures(), 'the figure of the analysis must stay hidden');
il_signal_dbfs(fig, 1, 100);                            % changed after the run
il_press(fig, 'open_graphs');
il_press(fig, 'open_graphs');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyEqual(numel(strfind(log, 'Drawing the SQAT figure')), 0, ...
    'the graphs window ran the metric again');
x_file = audioread(tc.TestData.wav_mono);                % 94 dBFS, as in the analysis
[~, ref] = evalc('Loudness_ISO532_1(x_file, tc.TestData.fs, 0, 2, 0.5, false)');
ln = findobj(il_window('SQAT_GUI_graphs'), 'Type', 'line');
tc.verifyTrue(any(arrayfun(@(l) isequal(l.YData(:), ref.InstantaneousLoudness(:)), ln)));
end

function test_gui_draws_the_figure_of_another_file_on_request(tc)
% The analysis keeps the figure of the active file alone. The figure of any
% other file is drawn when the graphs window asks for it, once, with the
% settings of the analysis.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
tc.verifyNumElements(il_all_sqat_figures(), 1, 'the analysis drew more than the active file');
il_mark_signal(fig, 1, false);                       % the tone is the only signal ticked
il_press(fig, 'open_graphs');
il_press(fig, 'open_graphs');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.verifyEqual(numel(strfind(log, 'Drawing the SQAT figure')), 1);
tc.verifyNumElements(il_all_sqat_figures(), 2);
end

function test_gui_graphs_fall_back_when_a_sqat_figure_fails(tc)
% Roughness_ECMA418_2 cannot draw its figure for this signal (CLim error in
% the toolbox); the values stand and the window shows the time series.
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_ECMA418_2'});
il_press(fig, 'run');
tc.verifyNotEmpty(findobj(fig, 'Tag', 'results_table').Data);
il_press(fig, 'open_graphs');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeSubstring(log, 'could not be drawn', 'the toolbox now draws this figure');
ln = findobj(il_window('SQAT_GUI_graphs'), 'Type', 'line');
tc.verifyNumElements(ln, 1);
end

function test_gui_close_deletes_the_drawn_figures(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
fig.CloseRequestFcn(fig, []);
tc.verifyEmpty(il_all_sqat_figures());
tc.verifyEmpty(il_window('SQAT_GUI_graphs'));
end

function test_gui_plot_of_a_nearly_constant_series_is_readable(tc)
% A constant result must not be zoomed into its rounding noise.
fig = SQAT_GUI({tc.TestData.wav_rough}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
il_set(il_window('SQAT_GUI_graphs'), 'graph_analysis', 'roughness');
ax = findobj(il_window('SQAT_GUI_graphs'), 'Type', 'axes');
y = findobj(ax, 'Type', 'line').YData;
tc.assertLessThan(max(y) - min(y), 1e-3 * mean(y));      % the premise: nearly constant
tc.verifyGreaterThanOrEqual(diff(ax.YLim), 0.09 * mean(y));
tc.verifyTrue(ax.YLim(1) <= min(y) && ax.YLim(2) >= max(y));
end

function test_gui_show_plots_after_run_opens_the_graphs_window(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
h = findobj(fig, 'Tag', 'show_plots'); h.Value = true;
il_press(fig, 'run');
tc.verifyNumElements(il_window('SQAT_GUI_graphs'), 1);
tc.verifyEmpty(il_sqat_figures());                % no second copy of the figure
end

function test_gui_split_figures_gives_one_tab_per_panel(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997', 'Loudness_ECMA418_2'});
il_press(fig, 'run');
for id = {'Roughness_Daniel1997', 'Loudness_ECMA418_2'}     % plain axes, and a tiled layout
    il_press(fig, 'open_graphs');
    g = il_window('SQAT_GUI_graphs');
    il_set(g, 'graph_metric', id{1});
    n_axes = numel(findobj(g, 'Type', 'axes'));
    tc.verifyNumElements(findobj(g, 'Type', 'uitab'), 1, id{1});
    cb = findobj(fig, 'Tag', 'split_figures'); cb.Value = true; cb.ValueChangedFcn(cb, []);
    tabs = findobj(g, 'Type', 'uitab');
    tc.verifyNumElements(tabs, n_axes, id{1});
    for t = tabs'
        tc.verifyNumElements(findobj(t, 'Type', 'axes'), 1, id{1});
    end
    cb.Value = false; cb.ValueChangedFcn(cb, []);
    tc.verifyNumElements(findobj(g, 'Type', 'uitab'), 1, id{1});
end
end

function test_gui_graph_window_offers_the_analyses_of_the_metric(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dd = findobj(g, 'Tag', 'graph_analysis');
% one signal: the SQAT figure and all the analyses come first
tc.verifyEqual(dd.ItemsData, {'sqat','all','loudness','loudness_level','specific_loudness_time','stats'});
tc.verifyEqual(dd.Value, 'sqat');
% the analyses follow the method: a stationary loudness has no time series
c = findobj(il_open_params(fig, 1), 'Tag', 'param_method'); c.Value = 1; c.ValueChangedFcn(c, []);
il_press(fig, 'run');
dd = findobj(il_window('SQAT_GUI_graphs'), 'Tag', 'graph_analysis');
tc.verifyEqual(dd.ItemsData, {'sqat','all','specific_loudness','stats'});
end

function test_gui_graphs_overlay_the_ticked_signals(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dd = findobj(g, 'Tag', 'graph_analysis');
% two signals: no SQAT figure (it is one per signal), the analyses that overlay
tc.verifyEqual(dd.ItemsData, {'loudness','loudness_level','specific_loudness_time','stats'});
tc.verifyEqual(dd.Value, 'loudness');
ln = findobj(g, 'Type', 'line');
tc.assertNumElements(ln, 2);
fs = tc.TestData.fs;
wavs = {tc.TestData.wav_mono, tc.TestData.wav_tone};
for k = 1:2
    x_file = audioread(wavs{k});
    [~, ref] = evalc('Loudness_ISO532_1(x_file, fs, 0, 2, 0.5, false)');
    hit = arrayfun(@(l) isequal(l.YData(:), ref.InstantaneousLoudness(:)) && ...
        isequal(l.XData(:), ref.time(:)), ln);
    tc.verifyTrue(any(hit), sprintf('signal %d is not on the plot', k));
end
lg = findobj(g, 'Type', 'legend');
tc.verifyEqual(sort(lg.String(:)), {'Signal #1, ch1'; 'Signal #2, ch1'});   % tags, not file names
% unticking a signal redraws the window that follows the list
il_mark_signal(fig, 2, false);
tc.verifyNumElements(findobj(il_window('SQAT_GUI_graphs'), 'Type', 'line'), 1);
end

function test_gui_graphs_overlay_profiles_over_the_bark_axis(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_rough}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'specific_roughness');
ln = findobj(g, 'Type', 'line');
tc.assertNumElements(ln, 2);
fs = tc.TestData.fs;
x_file = audioread(tc.TestData.wav_mono);
[~, ref] = evalc('Roughness_Daniel1997(x_file, fs, 0, false)');
hit = arrayfun(@(l) isequal(l.YData(:), ref.TimeAveragedSpecificRoughness(:)) && ...
    isequal(l.XData(:), ref.barkAxis(:)), ln);
tc.verifyTrue(any(hit));
ax = findobj(g, 'Type', 'axes');
tc.verifySubstring(ax.XLabel.String, 'Bark');
end

function test_gui_graphs_put_maps_side_by_side_on_one_colour_scale(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_rough}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'specific_roughness_time');
ax = findobj(g, 'Type', 'axes');
tc.assertNumElements(ax, 2);
tc.verifyNumElements(findobj(g, 'Type', 'surface'), 2);
tc.verifyEqual(ax(1).CLim, ax(2).CLim);
tc.verifyGreaterThan(diff(ax(1).CLim), 0);
for k = 1:2
    tc.verifyEqual(ax(k).Colormap, load('cmap_inferno.txt'));
end
end

function test_gui_graphs_tabulate_the_statistics_of_the_signals(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'stats');
st = findobj(g, 'Tag', 'graph_stats');
tc.assertNumElements(st, 1);
tc.verifyEqual(st.ColumnName(:)', {'Quantity', 'Signal #1, ch1', 'Signal #2, ch1'});
q = string(st.Data(:, 1));
names = {'tone_mono.wav', 'tone_1k_60dB.wav'};
for c = 2:3
    row = find(q == "N5");
    tc.verifyEqual(st.Data{row, c}, il_value(fig, names{c - 1}, 'Loudness_ISO532_1', 'N5'));
end
end

function test_gui_graphs_all_analyses_gives_a_tab_each(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
il_set(g, 'graph_analysis', 'all');
tabs = findobj(g, 'Type', 'uitab');
tc.verifyEqual(numel(tabs), 4);                 % three analyses and the statistics
tc.verifyNumElements(findobj(g, 'Type', 'surface'), 1);
tc.verifyNumElements(findobj(g, 'Tag', 'graph_stats'), 1);
end

function test_gui_graphs_choose_the_channel_of_a_binaural_analysis(tc)
fig = SQAT_GUI({tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ECMA418_2'});
il_signal_channel(fig, 1, 'All');
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dc = findobj(g, 'Tag', 'graph_channel');
tc.verifyEqual(dc.Items, {'1','2','Binaural','All'});
il_set(g, 'graph_analysis', 'loudness');
ref = audioread(tc.TestData.wav_stereo); fs = tc.TestData.fs;
[~, rb] = evalc('Loudness_ECMA418_2(ref, fs, ''free-frontal'', 0.304, false)');
il_set(g, 'graph_channel', '2');
tc.verifyEqual(findobj(g, 'Type', 'line').YData(:), rb.loudnessTDep(:, 2));
il_set(g, 'graph_channel', 'Binaural');
tc.verifyEqual(findobj(g, 'Type', 'line').YData(:), rb.loudnessTDepBin(:));
end

function test_gui_pinned_graph_window_keeps_its_signals(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
il_press(fig, 'open_graphs');                    % a window that follows the list is reused
g1 = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g1, 1);
il_set(g1, 'graph_pin', true);
il_mark_signal(fig, 2, false);
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);       % pinned: as it was
il_set(g1, 'graph_analysis', 'loudness_level');             % drawn again with its own signals
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
il_set(g1, 'graph_analysis', 'loudness');
il_press(fig, 'open_graphs');                    % the pinned window is kept, another opens
g = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g, 2);
g2 = g(g ~= g1);
tc.verifyEqual(findobj(g2, 'Tag', 'graph_analysis').Value, 'sqat');
il_set(g2, 'graph_analysis', 'loudness');
tc.verifyNumElements(findobj(g2, 'Type', 'line'), 1);
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
% a change in the list reaches the window that follows it, not the pinned one
il_mark_signal(fig, 2, true);
tc.verifyNumElements(findobj(g2, 'Type', 'line'), 2);
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
end

function test_gui_pinned_graph_window_keeps_its_results_after_another_run(tc)
% A run replaces the results; a pinned window keeps the ones it was pinned with.
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g1 = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g1, 1);
il_set(g1, 'graph_pin', true);
il_select_metrics(fig, {'Roughness_Daniel1997'});
il_press(fig, 'run');
il_press(fig, 'open_graphs');                    % the pinned window stays, another opens
g = il_window('SQAT_GUI_graphs');
tc.assertNumElements(g, 2);
g2 = g(g ~= g1);
tc.verifyEqual(findobj(g2, 'Tag', 'graph_metric').ItemsData, {'Roughness_Daniel1997'});
% the pinned window is drawn again from its own results, not from the last run
il_set(g1, 'graph_analysis', 'loudness_level');
il_set(g1, 'graph_analysis', 'loudness');
dm = findobj(g1, 'Tag', 'graph_metric');
tc.verifyEqual(dm.ItemsData, {'Loudness_ISO532_1'});
tc.verifyEqual(dm.Value, 'Loudness_ISO532_1');
tc.verifyNumElements(findobj(g1, 'Type', 'line'), 2);
tc.verifyEqual(findobj(g2, 'Tag', 'graph_metric').Value, 'Roughness_Daniel1997');
end

function test_gui_graphs_compare_a_mono_and_a_stereo_signal(tc)
fig = SQAT_GUI({tc.TestData.wav_mono, tc.TestData.wav_stereo}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_select_metrics(fig, {'Loudness_ISO532_1'});
il_signal_channel(fig, 2, 'All');
il_press(fig, 'run');
il_press(fig, 'open_graphs');
g = il_window('SQAT_GUI_graphs');
dc = findobj(g, 'Tag', 'graph_channel');
tc.verifyEqual(dc.Items, {'1', '2', 'All'});
il_set(g, 'graph_channel', 'All');
il_set(g, 'graph_analysis', 'loudness');
tc.verifyNumElements(findobj(g, 'Type', 'line'), 3);       % the mono, and both channels of the stereo
tc.verifyEqual(findobj(g, 'Type', 'legend').String, ...
    {'Signal #1, ch1', 'Signal #2, ch1', 'Signal #2, ch2'});
il_set(g, 'graph_channel', '1');
tc.verifyNumElements(findobj(g, 'Type', 'line'), 2);       % one channel, as before
end

function test_gui_graphs_window_needs_results(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_graphs');
tc.verifyEmpty(il_window('SQAT_GUI_graphs'));
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'Run an analysis first');
end

%% Waveform window ---------------------------------------------------------

function test_colormap_is_the_inferno_scale_of_the_toolbox(tc)
c = load('cmap_inferno.txt');
tc.verifySize(c, [256 3]);
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
ax = findobj(il_window('SQAT_GUI_waveform'), 'Tag', 'spectrogram');
tc.verifyEqual(ax.Colormap, c);
end

function test_waveform_window_shows_the_signal_and_follows_playback(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
tc.verifyNumElements(w, 1);
ln = findobj(w, 'Type', 'line');
tc.verifyNumElements(ln, 1);
tc.verifyEqual(ln.YData(:), audioread(tc.TestData.wav_mono));   % calibrated at 94 dBFS
ph = findobj(w, 'Tag', 'playhead');
tc.verifyEqual(ph.Value, 0);
b = findobj(w, 'Tag', 'play');
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyEqual(b.Text, 'Pause');
t_wait = tic;                                      % the audio device may start late
while ph.Value == 0 && toc(t_wait) < 3
    pause(0.05);
end
tc.verifyGreaterThan(ph.Value, 0);
il_press(w, 'play');                               % pause
tc.verifyEqual(b.Text, 'Play');
il_press(w, 'stop');
tc.verifyEqual(ph.Value, 0);
end

function test_waveform_window_has_the_spectrogram(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
ax = findobj(w, 'Tag', 'spectrogram');
tc.assertNumElements(ax, 1);
tc.verifyEqual(ax.Colormap, load('cmap_inferno.txt'));
tc.verifyEqual(ax.YScale, 'log');
sf = findobj(ax, 'Type', 'surface');
tc.assertNumElements(sf, 1);
% level in dB SPL: the 60 dB tone peaks at 60 dB, less the scalloping of the window
tc.verifyEqual(max(sf.CData(:)), 60, 'AbsTol', 1.5);
[~, i_max] = max(max(sf.CData, [], 2));
tc.verifyEqual(sf.YData(i_max), 1000, 'AbsTol', 48000/1024);
tc.verifySubstring(ax.Colorbar.Label.String, 'dB');
tc.verifyNotEmpty(findobj(w, 'Tag', 'playhead_spectrogram'));
end

function test_waveform_space_starts_and_pauses_playback(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
b = findobj(w, 'Tag', 'play');
w.KeyPressFcn(w, struct('Key', 'a'));                      % another key does nothing
tc.verifyEqual(b.Text, 'Play');
w.KeyPressFcn(w, struct('Key', 'space'));                  % the first play, with no button pressed
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyEqual(b.Text, 'Pause');
pause(0.4);
w.KeyPressFcn(w, struct('Key', 'space'));
tc.verifyEqual(b.Text, 'Play');
pause(0.4);
w.KeyPressFcn(w, struct('Key', 'space'));                  % and plays again
tc.verifyEqual(b.Text, 'Pause');
% the same key event arriving twice at once is one press (a focused button also answers to space)
pause(0.4);
w.KeyPressFcn(w, struct('Key', 'space'));
w.KeyPressFcn(w, struct('Key', 'space'));
tc.verifyEqual(b.Text, 'Play');
il_press(w, 'stop');
end

function test_waveform_tests_play_in_silence(tc)
% the tests of the waveform window start the player, and must not be heard
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
audio = getappdata(w, 'sqat_audio');
info = getappdata(w, 'sqat_play');
tc.assertGreaterThan(max(abs(audio())), 0.05, 'the file itself is silent');
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyTrue(info().playing);
tc.verifyEqual(double(info().buffer_peak), 0, 'the buffer that goes to the player is not silent');
il_press(w, 'stop');
end

function test_waveform_click_moves_the_playhead(tc)
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
ph = findobj(w, 'Tag', 'playhead'); ph2 = findobj(w, 'Tag', 'playhead_spectrogram');
ax = findobj(w, 'Tag', 'waveform_axes'); axs = findobj(w, 'Tag', 'spectrogram');
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [1.5 0 0]));
tc.verifyEqual(ph.Value, 1.5, 'AbsTol', 1/48000);
tc.verifyEqual(ph2.Value, 1.5, 'AbsTol', 1/48000);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [2.0 1000 0]));   % the spectrogram too
tc.verifyEqual(ph.Value, 2.0, 'AbsTol', 1/48000);
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [99 0 0]));         % beyond the end: the end
tc.verifyEqual(ph.Value, 3.0, 'AbsTol', 2/48000);
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [-4 0 0]));         % before the start: the start
tc.verifyEqual(ph.Value, 0, 'AbsTol', 1/48000);
% what is drawn on the axes does not take the click
tc.verifyTrue(all(strcmp({findobj(ax, 'Type', 'line').PickableParts}, 'none')));
tc.verifyEqual(findobj(axs, 'Type', 'surface').PickableParts, 'none');
end

function test_waveform_play_starts_where_the_click_was_and_stop_goes_back(tc)
% the state of the play is checked, not the time: the audio device starts a second or two late
fig = SQAT_GUI({tc.TestData.wav_mono}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
ph = findobj(w, 'Tag', 'playhead');
ax = findobj(w, 'Tag', 'waveform_axes');
info = getappdata(w, 'sqat_play');
fs = 48000;
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [1.5 0 0]));
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyTrue(info().playing);
tc.verifyEqual(info().sample, round(1.5*fs) + 1);
il_wait_until(@() ph.Value ~= 1.5);                        % the playhead moves once the device runs
tc.verifyNotEqual(ph.Value, 1.5);
% a click while playing jumps there and goes on
ax.ButtonDownFcn(ax, struct('IntersectionPoint', [0.2 0 0]));
tc.verifyEqual(findobj(w, 'Tag', 'play').Text, 'Pause');
tc.verifyTrue(info().playing);
tc.verifyEqual(info().sample, round(0.2*fs) + 1);
il_press(w, 'stop');
tc.verifyEqual(ph.Value, 0);
tc.verifyFalse(info().playing);
il_press(w, 'play');                                       % after Stop the file starts over
tc.verifyTrue(info().playing);
tc.verifyEqual(info().sample, 1);
il_press(w, 'stop');
end

function test_waveform_loops_by_default(tc)
fig = SQAT_GUI({tc.TestData.wav_short}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
lp = findobj(w, 'Tag', 'loop');
tc.verifyTrue(lp.Value, 'the loop must be on by default');
b = findobj(w, 'Tag', 'play');
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
pause(1.2);                                                % four times the length of the file
tc.verifyEqual(b.Text, 'Pause', 'the playback ended instead of looping');
il_press(w, 'stop');
tc.verifyEqual(b.Text, 'Play');
pause(0.5);
tc.verifyEqual(b.Text, 'Play', 'the loop restarted after Stop');
% without the loop the playback ends at the end of the file
lp.Value = false;                                          % read when the file ends
il_press(w, 'play');
pause(1.2);
tc.verifyEqual(b.Text, 'Play');
tc.verifyEqual(findobj(w, 'Tag', 'playhead').Value, 0);
end

function test_waveform_box_modes(tc)
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
axs = findobj(w, 'Tag', 'spectrogram');
ph = findobj(w, 'Tag', 'playhead');
audio = getappdata(w, 'sqat_audio');
[ref, fs] = audioread(tc.TestData.wav_two);
tc.verifyEqual(audio(), ref);                              % no box: the file as it is
tc.verifyEqual(findobj(w, 'Tag', 'draw_box').Text, 'Draw filter');
tc.verifyEqual(findobj(w, 'Tag', 'clear_boxes').Text, 'Clear filters');
% a click on the spectrogram moves the playhead until the tool is armed
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [2.5 1500 0]));   % first corner
tc.verifyEqual(ph.Value, 0);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [0.5 2500 0]));   % the opposite one, in any order
box = findobj(w, 'Tag', 'box');
tc.assertNumElements(box, 1);
tc.verifyEqual([min(box.XData) max(box.XData)], [0.5 2.5], 'AbsTol', 1e-9);
tc.verifyEqual([min(box.YData) max(box.YData)], [1500 2500], 'AbsTol', 1e-9);
tc.verifyFalse(bt.Value, 'the box tool stays armed');
% by default what is inside only plays: the sound is not changed
md = findobj(w, 'Tag', 'box_mode');
tc.verifyEqual(md.ItemsData, {'loop', 'isolate', 'remove'});
tc.verifyEqual(md.Items, {'Filter: loop only', 'Filter: isolate', 'Filter: remove'});
tc.verifyEqual(md.Value, 'loop');
tc.verifyEqual(audio(), ref);
mid = round(1.2*fs):round(1.8*fs); early = 1:round(0.3*fs);
% remove: what is inside goes
il_set(w, 'box_mode', 'remove');
y = audio();
tc.verifyLessThan(il_tone_db(y(mid), fs, 2000) - il_tone_db(ref(mid), fs, 2000), -40);
tc.verifyEqual(il_tone_db(y(mid), fs, 500), il_tone_db(ref(mid), fs, 500), 'AbsTol', 0.1);
tc.verifyEqual(il_tone_db(y(early), fs, 2000), il_tone_db(ref(early), fs, 2000), 'AbsTol', 0.1);
% isolate: only what is inside stays
il_set(w, 'box_mode', 'isolate');
y = audio();
tc.verifyEqual(il_tone_db(y(mid), fs, 2000), il_tone_db(ref(mid), fs, 2000), 'AbsTol', 0.1);
tc.verifyLessThan(il_tone_db(y(mid), fs, 500) - il_tone_db(ref(mid), fs, 500), -40);
tc.verifyLessThan(rms(y(early)), 1e-6);                    % outside the box's time
il_set(w, 'box_mode', 'loop');
tc.verifyEqual(audio(), ref);
% clearing the boxes clears the filter too
il_set(w, 'box_mode', 'remove');
tc.verifyNotEqual(audio(), ref);
il_press(w, 'clear_boxes');
tc.verifyEmpty(findobj(w, 'Tag', 'box'));
tc.verifyEqual(audio(), ref);
% a click on the same axes moves the playhead again
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.0 1000 0]));
tc.verifyEqual(ph.Value, 1.0, 'AbsTol', 1/fs);
end

function test_waveform_box_can_be_dragged(tc)
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
axs = findobj(w, 'Tag', 'spectrogram');
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [0.5 1500 0]));   % press
tc.assertNotEmpty(w.WindowButtonMotionFcn, 'nothing follows the mouse');
w.WindowButtonMotionFcn(w, struct('IntersectionPoint', [1.5 2000 0]));
prev = findobj(w, 'Tag', 'box_preview');
tc.assertNumElements(prev, 1);                             % the box follows the mouse
tc.verifyEqual([min(prev.XData) max(prev.XData)], [0.5 1.5], 'AbsTol', 1e-9);
w.WindowButtonMotionFcn(w, struct('IntersectionPoint', [2.5 2500 0]));
w.WindowButtonUpFcn(w, struct('IntersectionPoint', [2.5 2500 0]));   % release far away: the box
box = findobj(w, 'Tag', 'box');
tc.assertNumElements(box, 1);
tc.verifyEqual([min(box.XData) max(box.XData)], [0.5 2.5], 'AbsTol', 1e-9);
tc.verifyEqual([min(box.YData) max(box.YData)], [1500 2500], 'AbsTol', 1e-9);
tc.verifyEmpty(findobj(w, 'Tag', 'box_preview'));
tc.verifyEmpty(w.WindowButtonMotionFcn, 'the mouse is still followed');
tc.verifyFalse(bt.Value);
% a press and a release on the same spot is the first click of two
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.0 1000 0]));
w.WindowButtonUpFcn(w, struct('IntersectionPoint', [1.0 1000 0]));
tc.verifyNumElements(findobj(w, 'Tag', 'box'), 1);         % still the first box only
tc.verifyNumElements(findobj(w, 'Tag', 'box_corner'), 1);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [2.0 3000 0]));
tc.verifyNumElements(findobj(w, 'Tag', 'box'), 2);
end

function test_waveform_loops_inside_the_box(tc)
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
axs = findobj(w, 'Tag', 'spectrogram');
ph = findobj(w, 'Tag', 'playhead');
b = findobj(w, 'Tag', 'play');
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.0 1500 0]));
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [1.4 2500 0]));
il_press(w, 'play');                                       % starts at the box, not at the start
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
il_wait_until(@() ph.Value > 0);
tc.verifyGreaterThanOrEqual(ph.Value, 1.0 - 0.02);
% it stays in the box for more than three times its length
inside = true; t0 = tic;
while toc(t0) < 1.5
    pause(0.05);
    inside = inside && ph.Value >= 1.0 - 0.02 && ph.Value <= 1.4 + 0.1;
end
tc.verifyTrue(inside, 'the playhead left the box');
tc.verifyEqual(b.Text, 'Pause');
il_press(w, 'stop');
% without the loop the box is played once
lp = findobj(w, 'Tag', 'loop');
lp.Value = false;
il_press(w, 'play');
il_wait_until(@() strcmp(b.Text, 'Play') || ph.Value > 1.0);
il_wait_until(@() strcmp(b.Text, 'Play'));
tc.verifyEqual(b.Text, 'Play', 'the box was played more than once');
end

function test_waveform_box_loop_is_left_and_entered_by_clicks(tc)
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
axs = findobj(w, 'Tag', 'spectrogram'); axw = findobj(w, 'Tag', 'waveform_axes');
ph = findobj(w, 'Tag', 'playhead');
b = findobj(w, 'Tag', 'play');
info = getappdata(w, 'sqat_play');                         % where the play comes from, whatever the audio latency
fs = 48000;
ax_click = @(ax, t) ax.ButtonDownFcn(ax, struct('IntersectionPoint', [t 1000 0]));
draw_box = @(t1, t2) draw_two_corners(w, axs, t1, t2);
% the file plays from the start, outside where the box is going to be
ax_click(axw, 0.05);
il_press(w, 'play');
log = strjoin(findobj(fig, 'Tag', 'console').Value, newline);
tc.assumeTrue(contains(log, 'Playing'), 'no audio output on this machine');
tc.verifyFalse(info().box_loop);
% a box drawn now, with the loop on, takes the playback into it
draw_box(1.0, 1.4);
tc.verifyTrue(info().playing);
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(1.0*fs) + 1);
il_wait_until(@() ph.Value > 0.5);
tc.verifyTrue(il_stays_between(ph, 0.99, 1.41, 1.5), 'the playhead left the box');
% a click outside the box leaves the loop: the file plays on from there
ax_click(axw, 2.0);
tc.verifyEqual(b.Text, 'Pause');
tc.verifyFalse(info().box_loop);
tc.verifyEqual(info().sample, round(2.0*fs) + 1);
% a click inside the box goes back into the loop
ax_click(axw, 1.2);
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(1.2*fs) + 1);
il_wait_until(@() ph.Value > 1.19 && ph.Value < 1.5);
tc.verifyTrue(il_stays_between(ph, 0.99, 1.41, 1.5), 'the playhead left the box');
% the loop tick off: the play ends where the box ends
lp = findobj(w, 'Tag', 'loop');
lp.Value = false; lp.ValueChangedFcn(lp, []);
il_wait_until(@() strcmp(b.Text, 'Play'));
tc.verifyEqual(b.Text, 'Play', 'the box was not played once');
lp.Value = true; lp.ValueChangedFcn(lp, []);
il_press(w, 'stop');
% a box drawn after a pause takes the next play to its start
il_press(w, 'clear_boxes');
ax_click(axw, 2.2);
il_press(w, 'play'); pause(0.3); il_press(w, 'play');     % paused past the box
tc.verifyEqual(b.Text, 'Play');
draw_box(0.5, 0.9);
pause(0.3);
il_press(w, 'play');
tc.verifyTrue(info().box_loop);
tc.verifyEqual(info().sample, round(0.5*fs) + 1);
il_press(w, 'stop');
end

function test_waveform_weighting_filters_the_playback_and_the_spectrogram(tc)
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
wd = findobj(w, 'Tag', 'wave_weighting');
tc.verifyEqual(wd.Items, {'Z', 'A', 'C'});
tc.verifyEqual(wd.Value, 'Z');
audio = getappdata(w, 'sqat_audio');
[ref, fs] = audioread(tc.TestData.wav_two);
sf = findobj(w, 'Tag', 'spectrogram'); sf = findobj(sf, 'Type', 'surface');
z_level = @() max(sf.CData(abs(sf.YData - 500) < 30, :), [], 'all');
level_z = z_level();
il_set(w, 'wave_weighting', 'A');
tc.verifyEqual(audio(), SQAT_GUI_weight(ref, fs, 'A'), 'AbsTol', 1e-12);
sf = findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
tc.verifyEqual(max(sf.CData(abs(sf.YData - 500) < 30, :), [], 'all') - level_z, ...
    SQAT_GUI_weight_curve(500, fs, 'A'), 'AbsTol', 0.6);      % the bin nearest 500 Hz is within 1/2 bin
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Colorbar.Label.String, 'dB(A)');
il_set(w, 'wave_weighting', 'C');
tc.verifyEqual(audio(), SQAT_GUI_weight(ref, fs, 'C'), 'AbsTol', 1e-12);
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Colorbar.Label.String, 'dB(C)');
il_set(w, 'wave_weighting', 'Z');
tc.verifyEqual(audio(), ref);
end

function test_waveform_spectrogram_options(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
wd = findobj(w, 'Tag', 'spec_window');
tc.verifyEqual(wd.ItemsData, {'hann', 'hamming', 'rect', 'blackmanharris'});
tc.verifyEqual(wd.Value, 'hann');
tc.verifyEqual(findobj(w, 'Tag', 'spec_degree').Value, 10);
tc.verifyEqual(findobj(w, 'Tag', 'spec_overlap').Value, 50);
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
c_hann = surf().CData;
[x, fs] = audioread(tc.TestData.wav_tone);
[~, ~, Lh] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 50);
tc.verifyEqual(max(c_hann(:)), max(Lh(:)), 'AbsTol', 1e-9);   % calibrated at 94 dBFS
il_set(w, 'spec_window', 'rect');
tc.verifyNotEqual(surf().CData, c_hann);
il_set(w, 'spec_window', 'hann');
il_set(w, 'spec_degree', 12);
f4096 = (0:2048)' * fs / 4096;
tc.verifyEqual(numel(surf().YData), nnz(f4096 >= 20));
il_set(w, 'spec_overlap', 75);
tc.verifyEqual(numel(surf().XData), floor((numel(x) - 4096) / 1024) + 1);
% values outside the range go back into it
d = findobj(w, 'Tag', 'spec_degree');
tc.verifyClass(d, 'matlab.ui.control.Spinner');           % arrows, no typing and Enter
tc.verifyEqual(d.Limits, [6 16]);
tc.verifyEqual(d.Step, 1);
d.Value = 16; d.ValueChangedFcn(d, []);                   % what the up arrow gives at the top
tc.verifyEqual(d.Value, 16);
tc.verifyEqual(numel(findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface').YData), nnz((0:32768)' * 48000 / 65536 >= 20));
d.Value = 10; d.ValueChangedFcn(d, []);
o = findobj(w, 'Tag', 'spec_overlap');
tc.verifyClass(o, 'matlab.ui.control.Spinner');           % arrows too
tc.verifyEqual(o.Limits, [0 95]);
tc.verifyEqual(o.Step, 5);
o.Value = 95; o.ValueChangedFcn(o, []);                   % what the up arrow gives at the top
tc.verifyEqual(o.Value, 95);
o.Value = 0; o.ValueChangedFcn(o, []);                    % and the down arrow at the bottom: no overlap
[t0, ~, ~] = SQAT_GUI_spectrogram(x, fs, 'hann', 10, 0);
tc.verifyEqual(numel(findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface').XData), numel(t0));
end

function test_waveform_is_as_wide_as_the_spectrogram(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
pause(1);                                                      % the alignment waits for the layout
p = findobj(w, 'Tag', 'spectrogram').InnerPosition;
q = findobj(w, 'Tag', 'waveform_axes').InnerPosition;
tc.verifyEqual(q([1 3]), p([1 3]), 'AbsTol', 1);
end

function test_waveform_spectrogram_takes_a_window_from_a_file(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
c_hann = surf().CData;
wd = findobj(w, 'Tag', 'spec_window');
% a Hann window written to a file gives the Hann spectrogram
file = fullfile(tc.TestData.dir_tmp, 'my_hann.csv');
writematrix(SQAT_GUI_window('hann', 1024), file);
setappdata(w, 'sqat_next_file', file);                     % stands for the file dialog
il_press(w, 'import_window');
tc.verifyEqual(wd.Value, 'custom');
tc.verifySubstring(wd.Items{end}, 'my_hann');
tc.verifyEqual(surf().CData, c_hann, 'AbsTol', 1e-2);     % the text file keeps 5 digits
w_hann = SQAT_GUI_window('hann', 1024);
file = fullfile(tc.TestData.dir_tmp, 'my_hann.mat');
save(file, 'w_hann');
setappdata(w, 'sqat_next_file', file);
il_press(w, 'import_window');
tc.verifyEqual(surf().CData, c_hann, 'AbsTol', 1e-9);
% a window of another length is resampled
file = fullfile(tc.TestData.dir_tmp, 'my_ramp.txt');
writematrix([0 1 3 1 0], file);
setappdata(w, 'sqat_next_file', file);
il_press(w, 'import_window');
tc.verifyEqual(numel(wd.Items), 5, 'a second import must replace the first one');
tc.verifySubstring(wd.Items{end}, 'my_ramp');
tc.verifyNotEqual(surf().CData, c_hann);
% a file that is no window is refused and the window stays
file = fullfile(tc.TestData.dir_tmp, 'not_a_window.csv');
writematrix([1 2; 3 4], file);
before = wd.Value;
setappdata(w, 'sqat_next_file', file);
il_press(w, 'import_window');
tc.verifyEqual(wd.Value, before);
tc.verifySubstring(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'not_a_window');
end

%% Enhanced STFT in the waveform window ------------------------------------

function test_waveform_enhanced_stft_switch(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_no_background', true);        % the map is computed at once, in the foreground
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
sw = findobj(w, 'Tag', 'spec_enhanced');
tc.verifyEqual(sw.Value, 'Off');
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
c_plain = surf().CData;
il_set(w, 'spec_enhanced', 'On');
for tag = {'spec_degree', 'spec_overlap', 'spec_window', 'import_window'}
    tc.verifyEqual(char(findobj(w, 'Tag', tag{1}).Enable), 'off', ['faded: ' tag{1}]);
end
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Title.String, 'Enhanced STFT');
md = findobj(w, 'Tag', 'spec_enhanced_mode');
tc.verifyEqual(char(md.Enable), 'on');
tc.verifyEqual(md.Value, 'readable');
c = surf().CData;
[~, y] = il_map_centres(surf());
tc.verifyEqual(surf().FaceColor, 'texturemap');
il_set(w, 'spec_enhanced_mode', 'sharp');
tc.verifySubstring(findobj(w, 'Tag', 'spectrogram').Title.String, 'sharp');
tc.verifyNotEqual(surf().CData, c);
il_set(w, 'spec_enhanced_mode', 'readable');
tc.verifyEqual(surf().CData, c);
[~, i] = max(max(c, [], 2));
tc.verifyEqual(y(i), 1000, 'AbsTol', 25);                   % the tone sits at 1 kHz on the log grid
tc.verifyEqual(max(c(:)), max(c_plain(:)), 'AbsTol', 1e-3);  % peak aligned with the level of a steady tone
il_set(w, 'spec_enhanced', 'Off');
for tag = {'spec_degree', 'spec_overlap', 'spec_window', 'import_window'}
    tc.verifyEqual(char(findobj(w, 'Tag', tag{1}).Enable), 'on', ['back: ' tag{1}]);
end
tc.verifyEqual(char(md.Enable), 'off');
tc.verifyEqual(surf().CData, c_plain);
end

function test_enhanced_stft_zoom_recomputes(tc)
% a zoom on the enhanced map recomputes the excerpt with finer columns; zooming out restores the full map
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_no_background', true);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
t_full = surf().XData;
L_full = surf().CData;
t = il_map_centres(surf());
tc.verifyGreaterThan(t(2) - t(1), 0.0011);                    % 3 s in 2000 columns: 1.5 ms per column
ax.XLim = [1 1.5];
zoom_now = getappdata(w, 'sqat_spec_zoom');
zoom_now();
[t_z, y] = il_map_centres(surf());
tc.verifyEqual(t_z(2) - t_z(1), 0.001, 'AbsTol', 1e-9);        % 1 ms columns in the excerpt
tc.verifyLessThan(t_z(1), 1);
tc.verifyGreaterThan(t_z(end), 1.5);
tc.verifyLessThan(t_z(end) - t_z(1), 1.2);                     % only the excerpt and its margins
[~, i] = max(max(surf().CData, [], 2));
tc.verifyEqual(y(i), 1000, 'AbsTol', 10);
ax.XLim = [0 3];
zoom_now();
tc.verifyEqual(surf().XData, t_full);
tc.verifyEqual(surf().CData, L_full);
end

function test_enhanced_stft_options_keep_the_zoom(tc)
% the maps are kept without weighting: weighting and mode changes keep the zoom, and the
% weighting only adds its curve to each row
fig = SQAT_GUI({tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_no_background', true);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
ax.XLim = [1 1.5];
ax.YLim = [200 5000];
zoom_now = getappdata(w, 'sqat_spec_zoom');
zoom_now();
L_z = surf().CData;
[t_z, f_z] = il_map_centres(surf());
il_set(w, 'wave_weighting', 'A');
tc.verifyEqual(ax.XLim, [1 1.5]);
tc.verifyEqual(ax.YLim, [200 5000]);
tc.verifyEqual(il_map_centres(surf()), t_z);                   % the same excerpt
d = surf().CData - L_z;
tc.verifyEqual(d, repmat(d(:, 1), 1, numel(t_z)), 'AbsTol', 1e-9);   % one offset per row
a = SQAT_GUI_weight_curve(f_z, 48000, 'A');
tc.verifyEqual(d(:, 1), a(:), 'AbsTol', 1e-6);
il_set(w, 'spec_enhanced_mode', 'sharp');
tc.verifyEqual(ax.XLim, [1 1.5]);
tc.verifyEqual(ax.YLim, [200 5000]);
tc.verifyEqual(il_map_centres(surf()), t_z);
tc.verifySubstring(ax.Title.String, 'sharp');
il_set(w, 'wave_weighting', 'Z');
il_set(w, 'spec_enhanced_mode', 'readable');
tc.verifyEqual(surf().CData, L_z);                             % back where it started
end

function test_enhanced_stft_home_shows_the_whole_file(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_no_background', true);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
x_full = surf().XData;
L_full = surf().CData;
ax.XLim = [1 1.5];
ax.YLim = [500 2000];
zoom_now = getappdata(w, 'sqat_spec_zoom');
zoom_now();
tc.verifyNotEqual(surf().XData, x_full);
il_press(w, 'spec_home');
tc.verifyEqual(ax.XLim, [0 3]);
tc.verifyEqual(ax.YLim, [20 24000]);
tc.verifyEqual(surf().XData, x_full);
tc.verifyEqual(surf().CData, L_full);
end

function test_enhanced_stft_zoom_buttons_take_turns_with_the_box(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_no_background', true);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
zin = findobj(w, 'Tag', 'spec_zoom_in');
zout = findobj(w, 'Tag', 'spec_zoom_out');
box = findobj(w, 'Tag', 'draw_box');
il_set(w, 'spec_zoom_in', true);
tc.verifyEqual(char(zoom(w).Enable), 'on');
tc.verifyEqual(char(zoom(w).Direction), 'in');
il_set(w, 'spec_zoom_out', true);
tc.verifyFalse(zin.Value);
tc.verifyEqual(char(zoom(w).Direction), 'out');
il_set(w, 'draw_box', true);                                   % the box takes the clicks back
tc.verifyFalse(zout.Value);
tc.verifyEqual(char(zoom(w).Enable), 'off');
il_set(w, 'spec_zoom_in', true);                               % and gives them to the zoom
tc.verifyFalse(box.Value);
tc.verifyEqual(char(zoom(w).Enable), 'on');
il_set(w, 'spec_zoom_in', false);
tc.verifyEqual(char(zoom(w).Enable), 'off');
end

function test_enhanced_stft_colour_floor_moves_by_5_dB(tc)
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
setappdata(fig, 'sqat_no_background', true);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
top = ax.CLim(2);
tc.verifyEqual(ax.CLim, top + [-45 0], 'AbsTol', 1e-9);
il_press(w, 'spec_black_more');
tc.verifyEqual(ax.CLim, top + [-40 0], 'AbsTol', 1e-9);
w.KeyPressFcn(w, struct('Key', 'uparrow'));
tc.verifyEqual(ax.CLim, top + [-35 0], 'AbsTol', 1e-9);
il_set(w, 'spec_enhanced_mode', 'sharp');                      % the floor stays through a redraw
tc.verifyEqual(ax.CLim(2) - ax.CLim(1), 35, 'AbsTol', 1e-9);
il_set(w, 'spec_enhanced_mode', 'readable');
w.KeyPressFcn(w, struct('Key', 'downarrow'));
il_press(w, 'spec_black_less');
tc.verifyEqual(ax.CLim, top + [-45 0], 'AbsTol', 1e-9);
for k = 1:20, il_press(w, 'spec_black_more'); end
tc.verifyEqual(ax.CLim, top + [-5 0], 'AbsTol', 1e-9);         % at least 5 dB of colour
for k = 1:30, il_press(w, 'spec_black_less'); end
tc.verifyEqual(ax.CLim, top + [-105 0], 'AbsTol', 1e-9);       % at most 60 dB below the default
end

function test_enhanced_stft_comes_from_the_background_pool(tc)
% with the pool the window waits with a note, and the map that arrives equals the direct computation
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
tc.verifyEmpty(surf());
tc.verifyNotEmpty(findobj(ax, 'Tag', 'spec_wait'));
il_wait_until(@() ~isempty(surf()), 60);
tc.assertNotEmpty(surf());
tc.verifyEmpty(findobj(ax, 'Tag', 'spec_wait'));
[x, fs] = SQAT_GUI_load(tc.TestData.wav_tone, 94, 1);
[~, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
keep = f >= 20;
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
tc.verifyEmpty(strfind(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'ERROR'));
end

function test_enhanced_stft_of_a_long_signal_shows_a_preview_first(tc)
fs = 48000;
t = (0:1/fs:61-1/fs)';
wav = fullfile(tc.TestData.dir_tmp, 'long_tone.wav');
audiowrite(wav, 0.1*sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
fig = SQAT_GUI({wav}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
seen_preview = false;
t0 = tic;
while toc(t0) < 120
    pause(0.05);
    seen_preview = seen_preview || (~isempty(surf()) && contains(ax.Title.String, 'preview'));
    if ~isempty(surf()) && ~contains(ax.Title.String, 'preview')
        break
    end
end
tc.verifyTrue(seen_preview);
tc.assertNotEmpty(surf());
tc.verifyFalse(contains(ax.Title.String, 'preview'));
[x, fs] = SQAT_GUI_load(wav, 94, 1);
[~, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
keep = f >= 20;
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
end

function test_enhanced_stft_follows_a_new_signal_while_computing(tc)
% a map still on its way for the last signal never replaces the map of the new one
fig = SQAT_GUI({tc.TestData.wav_tone, tc.TestData.wav_two}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
il_activate_signal(fig, 1);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
il_activate_signal(fig, 2);
tc.verifySubstring(w.Name, 'two_tones');
ax = findobj(w, 'Tag', 'spectrogram');
surf = @() findobj(ax, 'Type', 'surface');
il_wait_until(@() ~isempty(surf()), 60);
tc.assertNotEmpty(surf());
pause(3);                                                      % time for a late map of the tone to arrive
tc.verifySubstring(w.Name, 'two_tones');
[x, fs] = SQAT_GUI_load(tc.TestData.wav_two, 94, 1);
[~, f, L] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
keep = f >= 20;
tc.verifyEqual(surf().CData, double(single(L(keep, :))) + SQAT_GUI_weight_curve(f(keep), fs, 'Z'));
end

function test_waveform_close_while_a_map_is_computed(tc)
% closing the window with maps queued leaves no timer running and logs no error
fig = SQAT_GUI({tc.TestData.wav_tone}, 'Visible', 'off');
tc.addTeardown(@() delete(fig));
timers_before = numel(timerfindall);
il_press(fig, 'open_waveform');
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
w.CloseRequestFcn(w, []);
tc.verifyFalse(isvalid(w));
pause(5);                                                      % any map on its way arrives meanwhile
tc.verifyEqual(numel(timerfindall), timers_before);
tc.verifyEmpty(strfind(strjoin(findobj(fig, 'Tag', 'console').Value, newline), 'ERROR'));
il_press(fig, 'open_waveform');                                % and the window opens again
w = il_window('SQAT_GUI_waveform');
il_set(w, 'spec_enhanced', 'On');
surf = @() findobj(findobj(w, 'Tag', 'spectrogram'), 'Type', 'surface');
il_wait_until(@() ~isempty(surf()), 60);
tc.verifyNotEmpty(surf());
end

%% Helpers -----------------------------------------------------------------

function p = il_default_params(e)
p = struct();
for q = e.params
    p.(q.name) = q.value;
end
end

function draw_two_corners(w, axs, t1, t2)
% arms the filter tool and clicks the two corners of a box from t1 to t2 (s), 1.5 to 2.5 kHz
bt = findobj(w, 'Tag', 'draw_box');
bt.Value = true; bt.ValueChangedFcn(bt, []);
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [t1 1500 0]));
axs.ButtonDownFcn(axs, struct('IntersectionPoint', [t2 2500 0]));
end

function ok = il_stays_between(ph, lo, hi, seconds)
% true when the playhead stays between lo and hi (s) for the given time
ok = true; t0 = tic;
while toc(t0) < seconds
    pause(0.05);
    ok = ok && ph.Value >= lo && ph.Value <= hi;
end
end

function il_wait_until(cond, seconds)
% waits up to 6 s (or the given time) for a condition (the audio device may take seconds to
% start), letting the timers and the background pool callbacks run
if nargin < 2, seconds = 6; end
t0 = tic;
while ~cond() && toc(t0) < seconds
    pause(0.05);
end
end

function [t, f] = il_map_centres(srf)
% time (row) and frequency (column) of the columns and rows of an enhanced map drawn as a
% texture: uniform in time, log-spaced in frequency, between the edges of its one face
[M, F] = deal(size(srf.CData, 2), size(srf.CData, 1));
xe = [min(srf.XData(:)), max(srf.XData(:))];
ye = log([min(srf.YData(:)), max(srf.YData(:))]);
t = xe(1) + ((1:M) - 0.5) * diff(xe) / M;
f = exp(ye(1) + ((1:F)' - 0.5) * diff(ye) / F);
end

function db = il_tone_db(y, fs, f0)
% level of the component of y at f0 (dB re an arbitrary reference), from a windowed FFT
n = numel(y);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/n);
Y = fft(y(:) .* w);
k = round(f0 * n / fs) + 1;
db = 20*log10(max(abs(Y(k-2:k+2))) * 2 / sum(w) + eps);
end

function v = il_field_or(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
end

function il_activate_signal(fig, k)
b = findobj(fig, 'Tag', sprintf('signal_name_%d', k));
b.ButtonPushedFcn(b, []);
end

function il_mark_signal(fig, k, tf)
c = findobj(fig, 'Tag', sprintf('signal_tick_%d', k));
c.Value = tf;
c.ValueChangedFcn(c, []);
end

function il_remove_signal(fig, k)
b = findobj(fig, 'Tag', sprintf('signal_remove_%d', k));   % the window is hidden: no question asked
b.ButtonPushedFcn(b, []);
end

function il_signal_channel(fig, k, c)
d = findobj(fig, 'Tag', sprintf('signal_channel_%d', k));
d.Value = c;
d.ValueChangedFcn(d, []);
end

function il_signal_dbfs(fig, k, v)
e = findobj(fig, 'Tag', sprintf('signal_cal_%d', k));
e.Value = v;
e.ValueChangedFcn(e, []);
end

function names = il_signal_names(fig)
names = {};
k = 1;
while ~isempty(findobj(fig, 'Tag', sprintf('signal_name_%d', k)))
    names{end+1} = findobj(fig, 'Tag', sprintf('signal_name_%d', k)).Text; %#ok<AGROW>
    k = k + 1;
end
end

function nums = il_signal_numbers(fig)
nums = {};
for k = 1:numel(il_signal_names(fig))
    nums{end+1} = findobj(fig, 'Tag', sprintf('signal_number_%d', k)).Text; %#ok<AGROW>
end
end

function il_select_metrics(fig, ids)
set_analyses = getappdata(fig, 'sqat_set_analyses');
set_analyses(ids);
end

function w = il_open_params(fig, k)
b = findobj(fig, 'Tag', sprintf('analysis_params_%d', k));
b.ButtonPushedFcn(b, []);
w = findall(groot, 'Type', 'figure', 'Tag', 'SQAT_GUI_params');
end

function tf = il_has_theme()
% the theme function, and with it the dark theme of the GUI, came in R2025a
tf = exist('theme', 'file') > 0;
end

function il_set(win, tag, value)
c = findobj(win, 'Tag', tag);
c.Value = value;
c.ValueChangedFcn(c, []);
end

function il_press(fig, tag)
b = findobj(fig, 'Tag', tag);
b.ButtonPushedFcn(b, []);
end

function w = il_window(tag)
w = findall(groot, 'Type', 'figure', 'Tag', tag);
end

function figs = il_sqat_figures()
% visible figures opened by the metrics (the interface windows carry a SQAT_GUI tag)
figs = il_all_sqat_figures();
figs = figs(strcmp({figs.Visible}, 'on'));
end

function figs = il_all_sqat_figures()
% every figure drawn by a metric, shown or kept hidden by the interface
figs = findall(groot, 'Type', 'figure');
if isempty(figs)
    return
end
figs = figs(~startsWith({figs.Tag}, 'SQAT_GUI'));
end

function v = il_value(fig, file, metric, quantity)
T = findobj(fig, 'Tag', 'results_table').Data;
v = T.Value(strcmp(T.File, file) & strcmp(T.Metric, metric) & strcmp(T.Quantity, quantity));
end
