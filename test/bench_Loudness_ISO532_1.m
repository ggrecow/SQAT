function T = bench_Loudness_ISO532_1(signals)
% function T = bench_Loudness_ISO532_1(signals)
%
% Benchmarks Loudness_ISO532_1 over the ISO test signals and prints a table.
% Reports wall-clock time only - correctness is the test suite's job
% (run_all_tests).
%
% To compare two revisions, run it on each and compare the totals:
%
%   git checkout main
%   matlab -batch "cd test; bench_Loudness_ISO532_1"
%   git checkout <branch>
%   matlab -batch "cd test; bench_Loudness_ISO532_1"
%
% If the branch you are comparing against predates test/, copy this file
% aside first, or check out just the metric:
%
%   git checkout <branch> -- psychoacoustic_metrics/
%
% INPUT
%   signals : optional vector of ISO signal numbers (default 6:25, the
%             time-varying ones - stationary signals are dominated by
%             filter-bank time and show little variation).
%
% OUTPUT
%   T : table of per-signal timings, also printed.
%
% Requires the Zenodo validation sounds in
% sound_files/validation_SQAT_v1_0/Loudness_ISO532_1 (doi 10.5281/zenodo.7933206).

    if nargin < 1
        signals = 6:25;
    end

    here = fileparts(mfilename('fullpath'));
    repo = fileparts(here);
    if isempty(which('Loudness_ISO532_1'))
        addpath(repo);
        startup_SQAT(fullfile(repo, filesep));
    end

    snd = fullfile(repo, 'sound_files', 'validation_SQAT_v1_0', 'Loudness_ISO532_1');
    if ~isfolder(snd)
        error('bench:noSounds', ...
             ['ISO test signals not found at\n  %s\n' ...
              'Download from doi 10.5281/zenodo.7933206 and place ' ...
              'validation_SQAT_v1_0 inside sound_files.'], snd);
    end
    cal = audioread(fullfile(snd, 'calibration signal sine 1kHz 60dB.wav'));

    fprintf('\n%s on %s\n', version('-release'), computer);
    fprintf('%s\n\n', repmat('=', 1, 58));
    fprintf(' sig | audio (s) | time (s) | x realtime | name\n');
    fprintf('%s\n', repmat('-', 1, 58));

    num = []; dur = []; sec = []; nm = strings(0);
    for n = signals
        d = dir(fullfile(snd, sprintf('Test signal %d (*.wav', n)));
        if isempty(d)
            fprintf(' %3d | (not found, skipped)\n', n);
            continue
        end
        [x, fs] = audioread(fullfile(snd, d(1).name));
        y = calibrate(x(:,1), cal, 60);

        Loudness_ISO532_1(y, fs, 0, 2, 0, 0);          % warm up
        t0 = tic;
        Loudness_ISO532_1(y, fs, 0, 2, 0, 0);
        el = toc(t0);

        name = regexprep(d(1).name, '^Test signal \d+ \((.*)\)\.wav$', '$1');
        fprintf(' %3d | %9.1f | %8.2f | %9.1fx | %s\n', ...
                n, numel(y)/fs, el, (numel(y)/fs)/el, name);

        num(end+1,1) = n;           %#ok<AGROW>
        dur(end+1,1) = numel(y)/fs; %#ok<AGROW>
        sec(end+1,1) = el;          %#ok<AGROW>
        nm(end+1,1)  = string(name);%#ok<AGROW>
    end

    fprintf('%s\n', repmat('-', 1, 58));
    fprintf(' TOTAL %9.1f s of audio in %.2f s  (%.1fx realtime)\n\n', ...
            sum(dur), sum(sec), sum(dur)/sum(sec));

    T = table(num, nm, dur, sec, dur./sec, ...
              'VariableNames', {'Signal','Name','AudioSeconds','Seconds','TimesRealtime'});

    if nargout == 0
        clear T;
    end
end
