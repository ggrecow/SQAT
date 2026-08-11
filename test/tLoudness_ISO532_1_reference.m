classdef tLoudness_ISO532_1_reference < matlab.unittest.TestCase
% Conformance tests for Loudness_ISO532_1 against the ISO 532-1:2017
% Annex A.4 reference implementation.
%
% Rather than comparing single-value statistics (Nmax, N5) - which average
% away where and why two implementations differ - these tests compare the
% loudness time series sample by sample against golden data produced by the
% unmodified reference C. See test/tools/ for how that data is regenerated.
%
% Requires the ISO test signals from Zenodo (doi 10.5281/zenodo.7933206)
% under sound_files/validation_SQAT_v1_0/Loudness_ISO532_1. Tests are
% skipped, not failed, when the signals are absent.
%
% Run with:  runtests('tLoudness_ISO532_1_reference')

    properties (Constant)
        % Largest tolerated deviation from the reference, in sone. The
        % reference implementation is deterministic and SQAT should track it
        % to well within this; it is not a perceptual tolerance.
        AbsTol = 0.02;
    end

    properties
        SoundDir
        GoldenDir
        CalSignal
    end

    properties (TestParameter)
        % Annex B.4 synthetic (6-13) and Annex B.5 technical (14-25) signals
        tvSignal = num2cell(6:25);
        % Annex B.3 stationary signals
        statSignal = num2cell(2:5);
    end

    methods (TestClassSetup)
        function locateData(tc)
            here          = fileparts(mfilename('fullpath'));
            repo          = fileparts(here);
            tc.GoldenDir  = fullfile(here, 'golden');
            tc.SoundDir   = fullfile(repo, 'sound_files', ...
                                     'validation_SQAT_v1_0', 'Loudness_ISO532_1');

            tc.assumeTrue(isfolder(tc.SoundDir), sprintf( ...
                ['ISO test signals not found at\n  %s\n' ...
                 'Download the dataset from doi 10.5281/zenodo.7933206 and ' ...
                 'place validation_SQAT_v1_0 inside sound_files.'], tc.SoundDir));

            calFile = fullfile(tc.SoundDir, 'calibration signal sine 1kHz 60dB.wav');
            tc.assumeTrue(isfile(calFile), ...
                sprintf('Calibration signal not found at\n  %s', calFile));
            tc.CalSignal = audioread(calFile);

            % Make sure the toolbox itself is reachable
            tc.assumeNotEmpty(which('Loudness_ISO532_1'), ...
                'Run startup_SQAT before the tests.');
        end
    end

    methods
        function insig = readSignal(tc, num)
            % Resolve "Test signal <num> (...).wav" and calibrate to 60 dB SPL
            d = dir(fullfile(tc.SoundDir, sprintf('Test signal %d (*.wav', num)));
            tc.assumeNotEmpty(d, sprintf('Test signal %d not found.', num));
            [x, ~] = audioread(fullfile(tc.SoundDir, d(1).name));
            insig  = calibrate(x(:,1), tc.CalSignal, 60);
        end

        function [fs] = signalRate(tc, num)
            d = dir(fullfile(tc.SoundDir, sprintf('Test signal %d (*.wav', num)));
            info = audioinfo(fullfile(tc.SoundDir, d(1).name));
            fs   = info.SampleRate;
        end

        function ref = readGolden(tc, name)
            f = fullfile(tc.GoldenDir, name);
            tc.assertTrue(isfile(f), sprintf( ...
                ['Golden data missing: %s\n' ...
                 'Regenerate with: cd test/tools && make && make golden'], f));
            m   = readmatrix(f, 'NumHeaderLines', 1);
            ref = m(:,2);
        end
    end

    methods (Test)

        function timeVaryingMatchesReference(tc, tvSignal)
            insig = tc.readSignal(tvSignal);
            fs    = tc.signalRate(tvSignal);
            ref   = tc.readGolden(sprintf('sig%02d_timevarying.csv', tvSignal));

            OUT = Loudness_ISO532_1(insig, fs, 0, 2, 0, 0);
            got = OUT.InstantaneousLoudness;

            n = min(numel(got), numel(ref));
            tc.assertGreaterThan(n, 0, 'No overlapping samples to compare.');

            d      = got(1:n) - ref(1:n);
            [mx,i] = max(abs(d));

            tc.verifyLessThanOrEqual(mx, tc.AbsTol, sprintf( ...
                ['Signal %d deviates from the ISO 532-1 reference by %.5f sone ' ...
                 '(%.2f%% of peak) at t = %.3f s.\n' ...
                 'mean|dN| = %.5f, RMS = %.5f over %d samples.'], ...
                tvSignal, mx, 100*mx/max(ref(1:n)), (i-1)*2e-3, ...
                mean(abs(d)), sqrt(mean(d.^2)), n));
        end

        function stationaryMatchesReference(tc, statSignal)
            insig = tc.readSignal(statSignal);
            fs    = tc.signalRate(statSignal);
            ref   = tc.readGolden(sprintf('sig%02d_stationary.csv', statSignal));

            OUT = Loudness_ISO532_1(insig, fs, 0, 1, 0, 0);

            tc.verifyEqual(OUT.Loudness, ref(1), 'AbsTol', tc.AbsTol, sprintf( ...
                'Stationary signal %d: got %.5f sone, reference %.5f sone.', ...
                statSignal, OUT.Loudness, ref(1)));
        end

        function levelInputMatchesReference(tc)
            levels = [-60 -60 78 79 89 72 80 89 75 87 85 79 86 80 ...
                       71  70 72 71 72 74 69 65 67 77 68 58 45 30];
            ref = tc.readGolden('sig01_stationary_levels.csv');
            OUT = Loudness_ISO532_1(levels, 1, 0, 0, 0, 0);

            tc.verifyEqual(OUT.Loudness, ref(1), 'AbsTol', tc.AbsTol, sprintf( ...
                'Level input: got %.5f sone, reference %.5f sone.', ...
                OUT.Loudness, ref(1)));
        end

    end
end
