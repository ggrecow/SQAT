function results = run_all_tests(varargin)
% function results = run_all_tests(...)
%
% Runs the SQAT test suite. Bootstraps the toolbox path if needed, so it can
% be called from a clean MATLAB session or headlessly:
%
%   matlab -batch "cd test; run_all_tests"
%
% Optional name/value arguments are forwarded to testsuite(), e.g.
%
%   run_all_tests('Name', '*anchor*')
%
% Tests that need the Zenodo validation sounds report as Incomplete
% (skipped) rather than Failed when those files are absent.

    here = fileparts(mfilename('fullpath'));
    repo = fileparts(here);

    if isempty(which('Loudness_ISO532_1'))
        addpath(repo);
        startup_SQAT(fullfile(repo, filesep));
    end
    addpath(here);

    suite   = testsuite(here, varargin{:});
    runner  = matlab.unittest.TestRunner.withTextOutput( ...
                  'OutputDetail', matlab.unittest.Verbosity.Terse);
    results = runner.run(suite);

    fprintf('\n%s\n', repmat('=', 1, 62));
    fprintf('  passed %d   failed %d   skipped %d   (%.1f s)\n', ...
            nnz([results.Passed]), nnz([results.Failed]), ...
            nnz([results.Incomplete]), sum([results.Duration]));
    fprintf('%s\n', repmat('=', 1, 62));

    if nargout == 0
        clear results;
    end
end
