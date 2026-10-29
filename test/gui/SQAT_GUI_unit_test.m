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

%% Extracting results ------------------------------------------------------

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

%% Spectrogram, windows, filters and weighting -----------------------------

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

%% Image PDF ---------------------------------------------------------------

function test_image_pdf_holds_the_pixels_on_a_page_of_their_size(tc)
% SQAT_GUI_image_pdf writes one page of the size of the image at the dpi given,
% with the pixels in one Flate stream, row by row, r g b, without loss.
rgb = reshape(uint8(10:10:180), [2 3 3]);              % 2 rows, 3 columns
f = fullfile(tc.TestData.dir_tmp, 'image.pdf');
SQAT_GUI_image_pdf(rgb, f, 300);
txt = fileread(f);
tc.verifyTrue(startsWith(txt, '%PDF-1.4'));
tc.verifySubstring(txt, '/MediaBox [0 0 0.7200 0.4800]');   % 3 x 2 px at 300 dpi, in points
tc.verifySubstring(txt, '/Width 3 /Height 2 /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode');
tc.verifySubstring(txt, 'q 0.7200 0 0 0.4800 0 0 cm /Im0 Do Q');
tc.verifyTrue(endsWith(strtrim(txt), '%%EOF'));
% the stream inflates to the pixels of the first row, then the second
fid = fopen(f, 'r'); bytes = fread(fid, '*uint8')'; fclose(fid);
k = strfind(bytes, uint8('/FlateDecode'));
k = k + find(bytes(k:end) == uint8(newline), 2);            % after "stream" and its line break
n = sscanf(txt(strfind(txt, '/FlateDecode /Length') + 20:end), '%d', 1);
in = java.util.zip.InflaterInputStream(java.io.ByteArrayInputStream(typecast(bytes(k(2):k(2) + n - 1), 'int8')));
out = zeros(1, 0, 'uint8');
while true
    b = in.read();                                          % a byte at a time: 18 of them
    if b < 0
        break
    end
    out(end+1) = b; %#ok<AGROW>
end
expected = reshape(permute(rgb, [3 2 1]), 1, []);
tc.verifyEqual(out, expected);
% the offsets of the cross reference table point at the objects
xref = sscanf(txt(strfind(txt, 'startxref') + 9:end), '%d', 1);
tc.verifyTrue(startsWith(txt(xref + 1:end), 'xref'));
offsets = sscanf(txt(xref + 1:end), 'xref 0 6 0000000000 65535 f %d 00000 n %d 00000 n %d 00000 n %d 00000 n %d 00000 n');
for j = 1:5
    tc.verifyTrue(startsWith(txt(offsets(j) + 1:end), sprintf('%d 0 obj', j)));
end
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
