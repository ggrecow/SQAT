function tests = SQAT_GUI_unit_test
% Unit tests of the functions in gui/: no window opens and no SQAT metric runs.
% How to run the three test files and how long they take: test/gui/README.md.
tests = functiontests(localfunctions);
end

%% Fixtures ----------------------------------------------------------------

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and gui/, also on a clean MATLAB path
fs = 48000;
t = (0:1/fs:3-1/fs)';
x_mono = sqrt(2)*2e-5*10^(60/20) * (1 + 0.5*sin(2*pi*4*t)) .* sin(2*pi*1000*t);   % 1 kHz, 60 dB SPL, AM at 4 Hz
dir_tmp = tempname; mkdir(dir_tmp);
tc.TestData.dir_tmp = dir_tmp;
tc.TestData.wav_mono = fullfile(dir_tmp, 'tone_mono.wav');
tc.TestData.wav_stereo = fullfile(dir_tmp, 'tone_stereo.wav');     % channel 2 is channel 1 less 10 dB
tc.TestData.wav_tone = fullfile(dir_tmp, 'tone_1k_60dB.wav');      % steady 1 kHz at 60 dB SPL
tc.TestData.wav_two = fullfile(dir_tmp, 'two_tones.wav');          % 500 Hz and 2 kHz
audiowrite(tc.TestData.wav_mono, x_mono, fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_stereo, [x_mono, x_mono*10^(-10/20)], fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_tone, sqrt(2)*2e-5*10^(60/20)*sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
audiowrite(tc.TestData.wav_two, 0.1*sin(2*pi*500*t) + 0.1*sin(2*pi*2000*t), fs, 'BitsPerSample', 32);
end

function teardownOnce(tc)
rmdir(tc.TestData.dir_tmp, 's');
end

%% Catalogue ---------------------------------------------------------------

function test_catalogue_lists_every_metric_of_sqat(tc)
% The catalogue of the GUI lists the sound level meter and the 13 metrics of
% SQAT in order, each an
% existing function with a label and a run handle.
m = SQAT_GUI_metrics;
expected = {'Do_SLM','Loudness_ISO532_1','Loudness_ECMA418_2','Sharpness_DIN45692', ...
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
% Every parameter of the catalogue is a choice or a number with a label, and
% the default of a choice is among its options.
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
% The default parameters of the catalogue equal the defaults of the toolbox
% (psychoacoustic_metrics_get_defaults) for each metric.
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

%% Loading audio -----------------------------------------------------------

function test_load_calibrates_with_dbfs_94(tc)
% At 94 dBFS the loaded signal equals the samples of the file, with its
% sampling rate and channel count.
[x, fs, nch] = SQAT_GUI_load(tc.TestData.wav_mono, 94, 1);
[ref, fs_ref] = audioread(tc.TestData.wav_mono);
tc.verifyEqual(x, ref(:,1));
tc.verifyEqual(fs, fs_ref);
tc.verifyEqual(nch, 1);
end

function test_load_applies_the_dbfs_gain(tc)
% Loading at 100 dBFS scales the signal by 6 dB against 94 dBFS.
x94  = SQAT_GUI_load(tc.TestData.wav_mono, 94, 1);
x100 = SQAT_GUI_load(tc.TestData.wav_mono, 100, 1);
tc.verifyEqual(x100, x94 * 10^((100-94)/20), 'AbsTol', 1e-15);
end

function test_load_applies_one_dbfs_per_channel(tc)
% a [1 x nch] full scale gives each channel its own gain, also when one channel is read
x = SQAT_GUI_load(tc.TestData.wav_stereo, [94 100], [1 2]);
ref = audioread(tc.TestData.wav_stereo);
tc.verifyEqual(x, ref .* [1, 10^(6/20)], 'AbsTol', 1e-15);
tc.verifyEqual(SQAT_GUI_load(tc.TestData.wav_stereo, [94 100], 2), ref(:, 2) * 10^(6/20), 'AbsTol', 1e-15);
end

%% Calibration -------------------------------------------------------------

function test_calibration_from_a_full_scale_level(tc)
% 'dbfs': the level given, for every channel (Greco 2026, Eq. 3.2)
[d, label] = SQAT_GUI_calibration('dbfs', tc.TestData.wav_stereo, 104);
tc.verifyEqual(d, [104 104]);
tc.verifyEqual(label, '104 dBFS');
end

function test_calibration_from_a_calibrator_recording(tc)
% 'calibrator': dBFS = level - 20 log10(rms of the recording) (Greco 2026,
% Eq. 3.1, as calibrate.m); a mono recording serves every channel, a
% recording with as many channels as the file calibrates each of them
fs = 48000;
t = (0:fs-1)' / fs;
mono = fullfile(tc.TestData.dir_tmp, 'cal_mono.wav');
stereo = fullfile(tc.TestData.dir_tmp, 'cal_stereo.wav');
audiowrite(mono, 0.5 * sin(2*pi*1000*t), fs, 'BitsPerSample', 32);
audiowrite(stereo, [0.5 * sin(2*pi*1000*t), 0.25 * sin(2*pi*1000*t)], fs, 'BitsPerSample', 32);
[d, label] = SQAT_GUI_calibration('calibrator', tc.TestData.wav_stereo, 94, mono);
tc.verifyEqual(d, (94 - 20*log10(0.5/sqrt(2))) * [1 1], 'AbsTol', 1e-4);
tc.verifyEqual(label, 'calib. 94 dB');
[~, ~, dBFS_ref] = calibrate(1, audioread(mono), 94);      % the toolbox function of Eq. 3.1
tc.verifyEqual(d(1), dBFS_ref, 'AbsTol', 1e-9);
d = SQAT_GUI_calibration('calibrator', tc.TestData.wav_stereo, 114, stereo);
tc.verifyEqual(d, 114 - 20*log10([0.5 0.25]/sqrt(2)), 'AbsTol', 1e-4);
end

function test_calibration_to_a_relative_level_keeps_the_channels_apart(tc)
% 'relative': one gain for the whole file, so that the rms of all channels
% together is the level (Greco 2026, Eq. 3.3); channel 2 stays 10 dB below 1
d = SQAT_GUI_calibration('relative', tc.TestData.wav_stereo, 70);
x = SQAT_GUI_load(tc.TestData.wav_stereo, d, [1 2]);
tc.verifyEqual(d(1), d(2));
tc.verifyEqual(20*log10(sqrt(mean(x(:).^2))) + 94, 70, 'AbsTol', 1e-6);   % SQAT reads 94 dB as 1 Pa
tc.verifyEqual(20*log10(rms(x(:, 1)) / rms(x(:, 2))), 10, 'AbsTol', 1e-6);
end

function test_calibration_with_one_level_per_channel(tc)
% Each method takes one level per channel: dBFS gives them as they are, the
% relative level sets the rms of each channel to its own level, and a number
% of levels that is neither 1 nor the number of channels is an error.
[d, label] = SQAT_GUI_calibration('dbfs', tc.TestData.wav_stereo, [94 100]);
tc.verifyEqual(d, [94 100]);
tc.verifyEqual(label, '94/100 dBFS');
[d, label] = SQAT_GUI_calibration('relative', tc.TestData.wav_stereo, [70 65]);
x = SQAT_GUI_load(tc.TestData.wav_stereo, d, [1 2]);
tc.verifyEqual(20*log10(rms(x)) + 94, [70 65], 'AbsTol', 1e-6);   % SQAT reads 94 dB as 1 Pa
tc.verifyEqual(label, 'rel. 70/65 dB');
tc.verifyError(@() SQAT_GUI_calibration('dbfs', tc.TestData.wav_stereo, [94 100 106]), 'SQAT_GUI:calibration');
end

function test_calibration_from_one_calibrator_recording_per_channel(tc)
% A head and torso simulator calibrated one ear at a time: recording k
% calibrates channel k, mono or stereo with the calibrator in channel k
% (the other channel near silence); a level per channel is the calibrator
% level of each. Two recordings for a mono file is an error.
fs = 48000;
t = (0:fs-1)' / fs;
s = sin(2*pi*1000*t);
f = fullfile(tc.TestData.dir_tmp, {'cal_left.wav', 'cal_right.wav', 'ear_left.wav', 'ear_right.wav'});
audiowrite(f{1}, 0.5 * s, fs, 'BitsPerSample', 32);
audiowrite(f{2}, 0.25 * s, fs, 'BitsPerSample', 32);
audiowrite(f{3}, [0.5 * s, 1e-4 * s], fs, 'BitsPerSample', 32);
audiowrite(f{4}, [1e-4 * s, 0.25 * s], fs, 'BitsPerSample', 32);
expected = 94 - 20*log10([0.5 0.25]/sqrt(2));
[d, label] = SQAT_GUI_calibration('calibrator', tc.TestData.wav_stereo, 94, f(1:2));
tc.verifyEqual(d, expected, 'AbsTol', 1e-4);
tc.verifyEqual(label, 'calib. 94 dB');
tc.verifyEqual(SQAT_GUI_calibration('calibrator', tc.TestData.wav_stereo, 94, f(3:4)), expected, 'AbsTol', 1e-4);
tc.verifyEqual(SQAT_GUI_calibration('calibrator', tc.TestData.wav_stereo, [94 114], f(1:2)), ...
    expected + [0 20], 'AbsTol', 1e-4);
tc.verifyError(@() SQAT_GUI_calibration('calibrator', tc.TestData.wav_mono, 94, f(1:2)), 'SQAT_GUI:calibration');
end

function test_calibration_rejects_silence(tc)
% Calibration stops with an error on a silent file (relative level) and on a
% silent calibrator recording.
silent = fullfile(tc.TestData.dir_tmp, 'silent.wav');
audiowrite(silent, zeros(4800, 1), 48000);
tc.verifyError(@() SQAT_GUI_calibration('relative', silent, 70), 'SQAT_GUI:calibration');
tc.verifyError(@() SQAT_GUI_calibration('calibrator', tc.TestData.wav_mono, 94, silent), 'SQAT_GUI:calibration');
end

function test_load_selects_the_channel(tc)
% Loading a stereo file on channel 2 returns that channel, and a channel the
% file does not have is an error.
[x2, ~, nch] = SQAT_GUI_load(tc.TestData.wav_stereo, 94, 2);
ref = audioread(tc.TestData.wav_stereo);
tc.verifyEqual(nch, 2);
tc.verifyEqual(x2, ref(:,2));
tc.verifyError(@() SQAT_GUI_load(tc.TestData.wav_stereo, 94, 3), 'SQAT_GUI:channel');
end

function test_load_reads_several_channels_at_once(tc)
% Several channels load at once in the order asked, with the gain of the
% dBFS; a missing or non-integer channel is an error.
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
% The single values of an output are its numeric scalars at the top level:
% time series, axis steps, nested results and text are left out.
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

function test_analysis_table_of_a_series_and_a_profile(tc)
% A series or a profile is a table of two columns, its axis and its values,
% named by their labels with the units.
A = struct('id', 'loudness', 'label', 'Loudness vs time', 'kind', 'series', 'x', [0; 0.5; 1], ...
    'y', [1; 2; 3], 'z', [], 'xlabel', 'Time (s)', 'ylabel', 'Loudness (sone)', 'zlabel', '', 'bandscale', 'linear');
T = SQAT_GUI_analysis_table(A);
tc.verifyEqual(T.Properties.VariableNames, {'Time (s)', 'Loudness (sone)'});
tc.verifyEqual(T{:, :}, [A.x, A.y]);
A.kind = 'profile';
A.xlabel = 'Critical band rate (Bark)';
A.ylabel = 'Specific loudness (sone/Bark)';
T = SQAT_GUI_analysis_table(A);
tc.verifyEqual(T.Properties.VariableNames, {'Critical band rate (Bark)', 'Specific loudness (sone/Bark)'});
A.ylabel = 'Specific loudness (sone_{HMS}/Bark_{HMS})';    % the braces of TeX go
T = SQAT_GUI_analysis_table(A);
tc.verifyEqual(T.Properties.VariableNames{2}, 'Specific loudness (sone_HMS/Bark_HMS)');
end

function test_analysis_table_of_a_map_names_each_band(tc)
% A map is its time axis and one column per band, each named by the quantity,
% its unit and the band with the unit of the band axis.
A = struct('id', 'specific_loudness_time', 'label', 'Specific loudness vs time', 'kind', 'map', ...
    'x', [0; 0.002], 'y', [0.1; 0.2; 24], 'z', [1 2 3; 4 5 6], 'xlabel', 'Time (s)', ...
    'ylabel', 'Critical band rate (Bark)', 'zlabel', 'Specific loudness (sone/Bark)', 'bandscale', 'linear');
T = SQAT_GUI_analysis_table(A);
tc.verifyEqual(T.Properties.VariableNames, {'Time (s)', 'Specific loudness (sone/Bark) at 0.1 Bark', ...
    'Specific loudness (sone/Bark) at 0.2 Bark', 'Specific loudness (sone/Bark) at 24 Bark'});
tc.verifyEqual(T{:, :}, [A.x, A.z]);
A.ylabel = 'Band';                                   % an axis with no unit
T = SQAT_GUI_analysis_table(A);
tc.verifyEqual(T.Properties.VariableNames{2}, 'Specific loudness (sone/Bark) at 0.1');
A.kind = 'image';
tc.verifyError(@() SQAT_GUI_analysis_table(A), 'SQAT_GUI:analysis_table');
end

%% Export of the analyses to Excel -----------------------------------------

function test_export_data_writes_a_sheet_per_analysis_and_per_map_channel(tc)
% One workbook per signal and metric: Info first, then the single values, the
% series and the profiles with a column per channel, then the maps with a sheet
% per channel, even when a map comes first in the list asked for.
t = (0:0.1:0.4)';
b = [100; 1000; 10000];
Z = reshape(1:15, 5, 3);
ser = @(y) il_analysis('loudness', 'Loudness vs time', 'series', t, y, [], 'Time (s)', 'Loudness (sone_{HMS})', '');
pro = @(y) il_analysis('specific_loudness', 'Time-averaged specific loudness', 'profile', b, y, [], ...
    'Band centre frequency (Hz)', 'Specific loudness (sone_{HMS}/Bark_{HMS})', '');
map = @(z) il_analysis('specific_tonal_loudness_time', 'Specific tonal loudness vs time', 'map', t, b, z, ...
    'Time (s)', 'Band centre frequency (Hz)', 'Specific tonal loudness (sone_{HMS}/Bark_{HMS})');
values = @(q, v) table(q, v, 'VariableNames', {'Quantity', 'Value'});
entries = [il_entry('1', [ser(t), pro(b), map(Z)], values({'N5'; 'Nmax'}, [1; 2])), ...
           il_entry('2', [ser(2*t), pro(2*b), map(2*Z)], values({'N5'; 'Nmax'}, [3; 4])), ...
           il_entry('Binaural', [ser(3*t), map(3*Z)], values({'Nmax'}, 6))];   % no profile, one value
info = cell2table({'Signal', '#1 train.wav'}, 'VariableNames', {'Item', 'Value'});
f = fullfile(tc.TestData.dir_tmp, 'train_s1_Loudness_ECMA418_2.xlsx');
files = SQAT_GUI_export_data(entries, {'specific_tonal_loudness_time', 'values', 'loudness', 'specific_loudness'}, ...
    f, info, @(q) 'sone_HMS');
tc.verifyEqual(files, {f});
tc.verifyEqual(cellstr(sheetnames(f))', {'Info', 'Single values', 'Loudness', 'Time-averaged specific loudness', ...
    'ch1 Spec. tonal loudness', 'ch2 Spec. tonal loudness', 'binaural Spec. tonal loudness'});
c = readcell(f, 'Sheet', 'Loudness');
tc.verifyEqual(c(1, :), {'Time (s)', 'Loudness (sone_HMS), ch1', 'Loudness (sone_HMS), ch2', 'Loudness (sone_HMS), binaural'});
tc.verifyEqual(readmatrix(f, 'Sheet', 'Loudness'), [t, t, 2*t, 3*t], 'AbsTol', 1e-12);
c = readcell(f, 'Sheet', 'Time-averaged specific loudness');
tc.verifyEqual(c(1, :), {'Band centre frequency (Hz)', 'Specific loudness (sone_HMS/Bark_HMS), ch1', ...
    'Specific loudness (sone_HMS/Bark_HMS), ch2'});
c = readcell(f, 'Sheet', 'binaural Spec. tonal loudness');
tc.verifyEqual(c{1, 2}, 'Specific tonal loudness (sone_HMS/Bark_HMS) at 100 Hz');
tc.verifyEqual(readmatrix(f, 'Sheet', 'binaural Spec. tonal loudness'), [t, 3*Z], 'AbsTol', 1e-12);
c = readcell(f, 'Sheet', 'Single values');
tc.verifyEqual(c(1, :), {'Quantity', 'Unit', 'Value, ch1', 'Value, ch2', 'Value, binaural'});
tc.verifyEqual(c(2:3, 1:4), {'N5', 'sone_HMS', 1, 3; 'Nmax', 'sone_HMS', 2, 4});
tc.verifyTrue(ismissing(c{2, 5}), 'binaural has no N5');
tc.verifyEqual(c{3, 5}, 6);
c = readcell(f, 'Sheet', 'Info');
tc.verifyEqual(c(1:3, :), {'Item', 'Value'; 'Signal', '#1 train.wav'; 'Channels', 'ch1, ch2, binaural'});
tc.verifyEqual(c(strcmp(c(:, 1), 'Sheet Loudness'), 2), {'Loudness vs time, one column per channel'});
tc.verifyEqual(c(strcmp(c(:, 1), 'Sheet ch2 Spec. tonal loudness'), 2), ...
    {'Specific tonal loudness vs time, ch2, one column per band'});
end

function test_export_data_shortens_the_long_names_of_sheets(tc)
% Excel takes 31 characters in the name of a sheet: the long names of the
% analyses are shortened, the headers of the columns keep them in full.
t = (0:0.5:1)';
b = [1; 2];
A = [il_analysis('weight_fr', 'Weighting of fluctuation strength and roughness vs time', 'series', t, t, [], 'Time (s)', 'w_{FR}', ''), ...
     il_analysis('weight_s', 'Weighting of sharpness and loudness vs time', 'series', t, t, [], 'Time (s)', 'w_{S}', ''), ...
     il_analysis('specific_fs', 'Time-averaged specific fluctuation strength', 'profile', b, b, [], ...
        'Critical band rate (Bark)', 'Specific fluctuation strength (vacil/Bark)', ''), ...
     il_analysis('specific_fs_time', 'Specific fluctuation strength vs time', 'map', t, b, ones(3, 2), ...
        'Time (s)', 'Critical band rate (Bark)', 'Specific fluctuation strength (vacil/Bark)')];
f = fullfile(tc.TestData.dir_tmp, 'long_names.xlsx');
info = cell2table(cell(0, 2), 'VariableNames', {'Item', 'Value'});
SQAT_GUI_export_data(il_entry('1', A, table({}, [], 'VariableNames', {'Quantity', 'Value'})), {A.id}, f, info);
names = cellstr(sheetnames(f))';
tc.verifyEqual(names, {'Info', 'w_FR', 'w_S', 'Avg. spec. fluctuation strength', 'ch1 Spec. fluctuation strength'});
tc.verifyLessThanOrEqual(cellfun(@numel, names), 31);
c = readcell(f, 'Sheet', 'w_FR');
tc.verifyEqual(c(1, :), {'Time (s)', 'w_FR, ch1'});
end

function test_export_data_writes_a_table_too_long_for_a_sheet_to_csv(tc)
% A sheet holds 1048575 rows under its header: a longer series goes to a CSV
% file next to the workbook, and the Info sheet names it.
t = (0:1048575)' / 1000;                               % 1048576 rows, 17.5 min at 1 kHz
A = il_analysis('level', 'Sound level vs time', 'series', t, 60 + 0*t, [], 'Time (s)', 'Sound pressure level (dBZ)', '');
f = fullfile(tc.TestData.dir_tmp, 'long_s1_Do_SLM.xlsx');
info = cell2table(cell(0, 2), 'VariableNames', {'Item', 'Value'});
files = SQAT_GUI_export_data(il_entry('1', A, table({}, [], 'VariableNames', {'Quantity', 'Value'})), {'level'}, f, info);
csv = fullfile(tc.TestData.dir_tmp, 'long_s1_Do_SLM_Sound_level.csv');
tc.verifyEqual(files, {f, csv});
tc.verifyEqual(cellstr(sheetnames(f))', {'Info'});
tc.verifyEqual(size(readmatrix(csv)), [1048576 2]);
c = readcell(f, 'Sheet', 'Info');
tc.verifyEqual(c(end, :), {'File long_s1_Do_SLM_Sound_level.csv', ...
    'Sound level vs time, one column per channel (1048576 rows, more than a sheet holds)'});
end

function test_export_data_needs_one_axis_for_the_channels(tc)
% The channels of a series share a column of time; different axes stop the export.
s = @(t) il_analysis('loudness', 'Loudness vs time', 'series', t, t, [], 'Time (s)', 'Loudness (sone)', '');
none = table({}, [], 'VariableNames', {'Quantity', 'Value'});
entries = [il_entry('1', s((0:2)'), none), il_entry('2', s((0:3)'), none)];
info = cell2table(cell(0, 2), 'VariableNames', {'Item', 'Value'});
tc.verifyError(@() SQAT_GUI_export_data(entries, {'loudness'}, fullfile(tc.TestData.dir_tmp, 'axes.xlsx'), info), ...
    'SQAT_GUI:export_data');
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

%% Spectrogram, windows, filters and weighting -----------------------------

function test_window_shapes(tc)
% Hann, Hamming, rectangular and Blackman-Harris windows come as columns of
% the length asked for, with the expected shape.
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
% A window can be read from a CSV or text file, as a row or a column, and from
% a .mat file; it always comes back as a column.
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
% A 60 dB tone on a bin of the FFT reads 60 dB at its frequency with each of
% the four windows, with the size, FFT length and overlap reported.
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
% The FFT degree sets the frequency resolution and the overlap the number of
% frames; a window given as a vector is used as it is or resampled to the FFT
% length, and a signal shorter than the FFT is padded.
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
% On a long signal the spectrogram lowers the overlap to keep the map within
% 4000 frames and 8e6 values, and reports it; the largest FFT also fits.
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
% The spectral filter takes out a box of time and frequency (over 40 dB at
% the tone inside it) and leaves the rest within 0.1 dB; in keep mode it keeps
% only the box; with no box the signal comes back.
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
% The A and C weightings give the IEC 61672-1 values at 100 Hz and 1 kHz,
% within 0.2 dB, Z leaves the signal unchanged, the curve drawn on the
% spectrogram is the response of the same filter, and B is an error.
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
    tc.verifyEqual(10*log10(median(sum(10.^(L(:, 500:1500)/10), 1))), 60, 'AbsTol', 0.1);   % dB SPL: a column adds up to the tone
end
tc.verifyLessThan(w(2), w(1));
tc.verifyError(@() SQAT_GUI_enhanced_stft(x, fs, 'other'), 'SQAT_GUI_enhanced_stft:smoothing');
end

function test_enhanced_stft_of_silence_and_a_high_f_min(tc)
% a silent excerpt (a zoom over a pause) gives a flat finite map; f_min at or above fs/2 is refused by name
fs = 48000;
[~, ~, L] = SQAT_GUI_enhanced_stft(zeros(fs/2, 1), fs, {'readable', 'sharp'});
for k = 1:2
    tc.verifyTrue(all(isfinite(L{k}(:))));
    tc.verifyEqual(min(L{k}(:)), max(L{k}(:)));
end
tc.verifyError(@() SQAT_GUI_enhanced_stft(zeros(fs/2, 1), fs, 'sharp', [], fs/2), 'SQAT_GUI_enhanced_stft:f_min');
end

function test_enhanced_stft_gives_both_smoothings_in_one_call(tc)
% the reassignment is shared: both maps at once equal the two calls
[x, fs] = audioread(tc.TestData.wav_two);
[t, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, {'readable', 'sharp'});
tc.verifyEqual(info.smoothing, {'readable', 'sharp'});
[t1, f1, L1] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
[~, ~, L2] = SQAT_GUI_enhanced_stft(x, fs, 'sharp');
tc.verifyEqual({t, f}, {t1, f1});
tc.verifyEqual(L, {L1, L2});
end

function test_enhanced_stft_excerpt_reads_on_the_scale_of_the_whole_map(tc)
% a loud half (80 dB tone, 40 dB noise) and a quiet one (20 dB tone, 0 dB noise): with the
% references of the whole map, a zoom into the quiet half keeps its levels and its floor
fs = 48000;
t = (0:1/fs:1.5-1/fs)';
rng(1);
x = [sqrt(2)*2e-5*1e4*sin(2*pi*1000*t) + 2e-5*100*randn(size(t)); ...
     sqrt(2)*2e-5*10*sin(2*pi*2000*t) + 2e-5*randn(size(t))];
[t_w, f, L, info] = SQAT_GUI_enhanced_stft(x, fs, 'readable');
i1 = round(1.7*fs) + 1;
i2 = round(2.8*fs);
[t_z, ~, L_z] = SQAT_GUI_enhanced_stft(x(i1:i2), fs, 'readable', 1100, [], [], [], info.ref);
t_z = t_z + (i1 - 1) / fs;
in = @(tv) tv >= 2 & tv <= 2.5;
tone = abs(f - 2000) < 40;
tc.verifyEqual(max(L_z(tone, in(t_z)), [], 'all'), max(L(tone, in(t_w)), [], 'all'), 'AbsTol', 0.5);
tc.verifyEqual(min(L_z(:, in(t_z)), [], 'all'), min(L(:, in(t_w)), [], 'all'), 'AbsTol', 0.5);
[~, ~, L_own] = SQAT_GUI_enhanced_stft(x(i1:i2), fs, 'readable', 1100);   % on its own scale, the floor drops
tc.verifyLessThan(min(L_own, [], 'all'), min(L(:, in(t_w)), [], 'all') - 20);
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

%% Path --------------------------------------------------------------------

function test_startup_puts_the_gui_on_the_path(tc)
% startup_SQAT adds the gui folder to the MATLAB path.
rmpath(fullfile(basepath_SQAT, 'gui'));
tc.addTeardown(@() addpath(fullfile(basepath_SQAT, 'gui')));
evalc('run(fullfile(basepath_SQAT, ''startup_SQAT.m''))');
tc.verifyNotEmpty(which('SQAT_GUI'));
end

%% Helpers -----------------------------------------------------------------

function a = il_analysis(id, label, kind, x, y, z, xlabel, ylabel, zlabel)
% one analysis as SQAT_GUI_extract gives it
a = struct('id', id, 'label', label, 'kind', kind, 'x', x, 'y', y, 'z', z, ...
    'xlabel', xlabel, 'ylabel', ylabel, 'zlabel', zlabel, 'bandscale', 'linear');
end

function en = il_entry(channel, analyses, values)
% one entry of the results of SQAT_GUI: a channel of a signal and a metric
en = struct('channel', channel, 'analyses', analyses, 'values', values);
end

function p = il_default_params(e)
p = struct();
for q = e.params
    p.(q.name) = q.value;
end
end

function db = il_tone_db(y, fs, f0)
% level of the component of y at f0 (dB re an arbitrary reference), from a windowed FFT
n = numel(y);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/n);
Y = fft(y(:) .* w);
k = round(f0 * n / fs) + 1;
db = 20*log10(max(abs(Y(k-2:k+2))) * 2 / sum(w) + eps);
end
