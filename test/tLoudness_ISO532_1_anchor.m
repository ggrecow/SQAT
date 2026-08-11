classdef tLoudness_ISO532_1_anchor < matlab.unittest.TestCase
% Anchor tests for Loudness_ISO532_1.
%
% The defining calibration point of ISO 532-1 is that a 1 kHz tone at
% 40 dB SPL has a loudness of 1 sone. These tests assert exactly that, plus
% the Annex B.2 third-octave level vector, and they are deliberately
% self-contained: the signals are synthesised or hard-coded here, so nothing
% outside the toolbox is needed and the whole class runs in a few seconds.
%
% This is the guard rail. It is not the test that detects the ISO 532-1
% conformance defects - see tLoudness_ISO532_1_reference for that.
%
% Run with:  runtests('tLoudness_ISO532_1_anchor')

    properties (Constant)
        fs   = 48000;       % Hz, the rate the filter bank is designed for
        pref = 2e-5;        % Pa
    end

    methods (Static)
        function insig = tone(level_dB, freq_Hz, dur_s, fs, pref)
            % Calibrated sine of the requested SPL, as a [N x 1] Pa signal
            t     = (0:round(dur_s*fs)-1).'/fs;
            p_rms = pref * 10^(level_dB/20);
            insig = sqrt(2) * p_rms * sin(2*pi*freq_Hz*t);
        end
    end

    methods (Test)

        function referenceToneIsOneSone_stationary(tc)
            % 1 kHz @ 40 dB SPL -> 1 sone, method 1 (stationary from audio)
            insig = tc.tone(40, 1000, 1, tc.fs, tc.pref);
            OUT   = Loudness_ISO532_1(insig, tc.fs, 0, 1, 0, 0);

            tc.verifyEqual(OUT.Loudness, 1, 'RelTol', 0.01, ...
                'A 1 kHz tone at 40 dB SPL must yield 1 sone (stationary method).');
        end

        function referenceToneIsOneSone_timeVarying(tc)
            % Same tone through method 2; compare the steady-state value once
            % the temporal weighting has settled (last 100 ms).
            insig = tc.tone(40, 1000, 1, tc.fs, tc.pref);
            OUT   = Loudness_ISO532_1(insig, tc.fs, 0, 2, 0, 0);

            steady = mean(OUT.InstantaneousLoudness(end-49:end));
            tc.verifyEqual(steady, 1, 'RelTol', 0.01, ...
                'A 1 kHz tone at 40 dB SPL must yield 1 sone (time-varying method).');
        end

        function levelInputMatchesAnnexB2(tc)
            % Annex B.2 test signal 1: 28 unweighted third-octave levels,
            % 25 Hz .. 12.5 kHz. ISO gives N = 83.296 sone.
            levels = [-60 -60 78 79 89 72 80 89 75 87 85 79 86 80 ...
                       71  70 72 71 72 74 69 65 67 77 68 58 45 30];
            OUT = Loudness_ISO532_1(levels, 1, 0, 0, 0, 0);

            tc.verifyEqual(OUT.Loudness, 83.296, 'RelTol', 0.005, ...
                'Annex B.2 level vector must yield 83.296 sone.');
        end

        function loudnessGrowsMonotonicallyWithLevel(tc)
            % Sanity property: louder input, louder result. Cheap, and it
            % catches gross breakage in the core-loudness / slope stages.
            lvls = 30:10:80;
            N    = zeros(size(lvls));
            for k = 1:numel(lvls)
                insig = tc.tone(lvls(k), 1000, 0.5, tc.fs, tc.pref);
                N(k)  = Loudness_ISO532_1(insig, tc.fs, 0, 1, 0, 0).Loudness;
            end
            tc.verifyTrue(all(diff(N) > 0), ...
                sprintf('Loudness must increase with level, got [%s].', ...
                        num2str(N, '%.3f ')));
        end

        function doublingLevelBy10dBRoughlyDoublesLoudness(tc)
            % Above ~40 phon, +10 dB should roughly double loudness. Loose
            % bounds - this is a smoke test for the sone scale, not a spec.
            N60 = Loudness_ISO532_1(tc.tone(60, 1000, 0.5, tc.fs, tc.pref), ...
                                    tc.fs, 0, 1, 0, 0).Loudness;
            N70 = Loudness_ISO532_1(tc.tone(70, 1000, 0.5, tc.fs, tc.pref), ...
                                    tc.fs, 0, 1, 0, 0).Loudness;

            tc.verifyGreaterThan(N70/N60, 1.7);
            tc.verifyLessThan(   N70/N60, 2.3);
        end

    end
end
