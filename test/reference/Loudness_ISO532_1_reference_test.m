function tests = Loudness_ISO532_1_reference_test
% Loudness_ISO532_1 against the 25 test signals of ISO 532-1:2017, Annex B,
% with the criteria of the standard (clause 7 and Annex B):
%   stationary signals (B.2, B.3): total loudness within 5 % or 0.1 sone,
%     whichever is larger, and specific loudness within the envelope of
%     5 % or 0.1 sone/Bark that the standard gives;
%   time-varying signals (B.4, B.5): loudness vs time, and for B.4 the
%     specific loudness vs time, within the envelope of 5 % or 0.1 (with
%     the +/-2 ms temporal tolerance built in), and at most 1 % of the
%     samples up to the envelope of 10 % or 0.2.
% The reference values and envelopes are the ones of the standard, in
% validation/Loudness_ISO532_1/*/reference_values. The sound files are the
% ones of the standard, from https://standards.iso.org/iso/532/-1/ed-1/en/
% (ISO 532-1 - Program etc.zip); they are read from the folder named by the
% environment variable SQAT_ISO532_1_DIR, or from test/reference/data/, and
% the tests are marked incomplete when neither holds them.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
tc.applyFixture(matlab.unittest.fixtures.PathFixture(root, 'IncludingSubfolders', true));   % SQAT and test/report
iso = getenv('SQAT_ISO532_1_DIR');
if isempty(iso)
    iso = fullfile(fileparts(mfilename('fullpath')), 'data', 'ISO 532-1 - Program etc');
end
tc.TestData.iso = iso;
tc.TestData.ref = fullfile(root, 'validation', 'Loudness_ISO532_1');
cal = fullfile(iso, 'Annex C', 'calibration signal sine 1kHz 60dB.wav');
tc.TestData.ok = isfile(cal);
if tc.TestData.ok
    c = audioread(cal);                        % 1 kHz at 60 dB SPL: the level of full scale
    tc.TestData.gain = 10^((60 - 20*log10(rms(c)) - 94)/20);   % to the SQAT convention, 1 = 1 Pa
end
end

function test_stationary_signals_B2_and_B3(tc)
% The five stationary test signals of ISO 532-1:2017 (B.2 from third-octave
% levels, B.3 from the sound files): total loudness within 5 % or 0.1 sone
% of the reference, and specific loudness inside the envelope of the
% standard.
tc.assumeTrue(tc.TestData.ok, 'ISO 532-1 sound files not found (SQAT_ISO532_1_DIR)');
N_ref = [83.296 14.655 4.019 1.549 10.498];    % total loudness of signals 1 to 5 (sone)
names = {'', 'Annex B.3/Test signal 2 (250 Hz 80 dB).wav', 'Annex B.3/Test signal 3 (1 kHz 60 dB).wav', ...
         'Annex B.3/Test signal 4 (4 kHz 40 dB).wav', 'Annex B.3/Test signal 5 (pinknoise 60 dB).wav'};
for k = 1:5
    if k == 1                                   % B.2: one-third octave levels, method 0
        L = [-60 -60 78 79 89 72 80 89 75 87 85 79 86 80 71 70 72 71 72 74 69 65 67 77 68 58 45 30];
        [~, O] = evalc('Loudness_ISO532_1(L, 1, 0, 0, 0, false)');
    else
        [x, fs] = audioread(fullfile(tc.TestData.iso, names{k}));
        [~, O] = evalc('Loudness_ISO532_1(x * tc.TestData.gain, fs, 0, 1, 0, false)');
    end
    tol = max(0.05 * N_ref(k), 0.1);
    sqat_report_record('loudness_iso532_1', sprintf('signal %d, total loudness (sone)', k), O.Loudness, N_ref(k), tol, 'ISO 532-1:2017, 5 % or 0.1 sone');
    tc.verifyEqual(O.Loudness, N_ref(k), 'AbsTol', tol, sprintf('signal %d, total loudness', k));
    R = il_ref(tc, '1_synthetic_signals_stationary_loudness', k);    % [Bark, N'ref, N'min, N'max]
    s = interp1(O.barkAxis, O.SpecificLoudness, R(:, 1));
    il_envelope(tc, sprintf('signal %d, specific loudness', k), s, R(:, 2), R(:, 3), R(:, 4), [], []);
end
end

function test_time_varying_signals_B4_and_B5(tc)
% The 20 time-varying test signals of ISO 532-1:2017 (B.4 synthetic, B.5
% technical): loudness against time, and for B.4 the specific loudness
% against time at the Bark given by the standard, inside the 5 % envelope,
% with at most 1 % of the samples out to the 10 % envelope.
tc.assumeTrue(tc.TestData.ok, 'ISO 532-1 sound files not found (SQAT_ISO532_1_DIR)');
d = [dir(fullfile(tc.TestData.iso, 'Annex B.4', '*.wav')); dir(fullfile(tc.TestData.iso, 'Annex B.5', '*.wav'))];
tc.assertNumElements(d, 20);
for k = 1:numel(d)
    n = sscanf(d(k).name, 'Test signal %d');
    field = double(n == 15);                    % signal 15 (vehicle interior) is the diffuse-field one
    [x, fs] = audioread(fullfile(d(k).folder, d(k).name));
    [~, O] = evalc('Loudness_ISO532_1(x(:, 1) * tc.TestData.gain, fs, field, 2, 0, false)');
    folder = '2_synthetic_signals_time_varying_loudness';
    if n >= 14
        folder = '3_technical_signals_time_varying_loudness';
    end
    R = il_ref(tc, folder, n);                  % [t, Nref, Nmin5, Nmax5, Nmin10, Nmax10]
    N = interp1(O.time, O.InstantaneousLoudness, R(:, 1));
    il_envelope(tc, sprintf('signal %d, loudness vs time', n), N, R(:, 2), R(:, 3), R(:, 4), R(:, 5), R(:, 6));
    if n <= 13                                  % B.4 also gives the specific loudness at one Bark
        S = load(fullfile(tc.TestData.ref, folder, 'reference_values', ...
            sprintf('reference_values_ISO532_1_2017_signal_%d_specific_loudness.mat', n)));
        [~, iz] = min(abs(O.barkAxis - S.target_bark));
        s = interp1(O.time, O.InstantaneousSpecificLoudness(:, iz), R(:, 1));
        Q = S.reference_2;                      % [N'ref, N'min5, N'max5, N'min10, N'max10], on the times of R
        il_envelope(tc, sprintf('signal %d, specific loudness vs time at %.1f Bark', n, S.target_bark), ...
            s, Q(:, 1), Q(:, 2), Q(:, 3), Q(:, 4), Q(:, 5));
    end
end
end

function R = il_ref(tc, folder, n)
S = load(fullfile(tc.TestData.ref, folder, 'reference_values', sprintf('reference_values_ISO532_1_2017_signal_%d.mat', n)));
R = S.reference;
end

function il_envelope(tc, what, x, ref, lo5, hi5, lo10, hi10)
% x inside [lo5, hi5] everywhere, or, when the 10 % envelope is given, at most
% 1 % of the samples outside it and all of them inside [lo10, hi10]. The
% report gets the worst deviation as a fraction of the 5 % half-width on its
% side (1 = on the edge of the envelope).
up = x > ref;
u = zeros(size(x));
u(up) = (x(up) - ref(up)) ./ max(hi5(up) - ref(up), eps);
u(~up) = (ref(~up) - x(~up)) ./ max(ref(~up) - lo5(~up), eps);
out5 = x < lo5 - 1e-9 | x > hi5 + 1e-9;
if isempty(lo10)
    sqat_report_record('loudness_iso532_1', [what ', 5 % envelope'], u, zeros(size(u)), 1, 'ISO 532-1:2017, envelope of 5 % or 0.1 within 2 ms');
    tc.verifyEqual(nnz(out5), 0, [what ': samples outside the 5 % envelope']);
else
    in = ~out5;
    sqat_report_record('loudness_iso532_1', [what ', 5 % envelope (99 % of the samples)'], sort(u(in)), zeros(nnz(in), 1), 1, 'ISO 532-1:2017, 5 % envelope; up to 1 % of the samples to 10 %');
    tc.verifyLessThanOrEqual(nnz(out5) / numel(x), 0.01, [what ': more than 1 % outside the 5 % envelope']);
    tc.verifyEqual(nnz(x < lo10 - 1e-9 | x > hi10 + 1e-9), 0, [what ': samples outside the 10 % envelope']);
end
end
