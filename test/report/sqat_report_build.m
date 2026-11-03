function ok = sqat_report_build(varargin)
% function ok = sqat_report_build(varargin)
%
%   Runs the test suite (test/, every subfolder) with the measured-error
%   recorder on, and writes one self-contained page, test/report/out/index.html
%   (open it in any browser; no server, no network).
%
%   sqat_report_build               % the whole suite
%   sqat_report_build('unit')       % one folder of test/
%   ok = sqat_report_build(...)     % false when a test failed (the CI uses it)
%
%   The page is template.html with the data of the run injected: the result
%   of every test, its failure output, what it checks (the comment block under
%   its function line) and the error of every comparison recorded with
%   sqat_report_record. The same data goes to test/report/out/data.json, so
%   that the runs of several folders (the parallel jobs of the CI) can be
%   joined into one page by merge_reports.py.

here = fileparts(mfilename('fullpath'));
root_test = fileparts(here);
root = fileparts(root_test);
out = fullfile(here, 'out');
if ~isfolder(out)
    mkdir(out);
end
target = root_test;
if nargin > 0
    target = fullfile(root_test, varargin{1});
end

records = fullfile(out, 'records.jsonl');
fclose(fopen(records, 'w'));
setenv('SQAT_REPORT_FILE', records);
cleanup = onCleanup(@() setenv('SQAT_REPORT_FILE', ''));
addpath(here);
suite = matlab.unittest.TestSuite.fromFolder(target, 'IncludingSubfolders', true);
suite = suite(~startsWith({suite.BaseFolder}, here));   % the report folder holds no tests
% the detailed output goes to the screen as the tests run (the CI log shows
% where a run stops) and to a file, for the raw log of the report
logfile = fullfile(out, 'runner.log');
if isfile(logfile)
    delete(logfile);
end
import matlab.unittest.plugins.*
runner = matlab.unittest.TestRunner.withNoPlugins;
runner.addPlugin(TestRunProgressPlugin.withVerbosity(matlab.unittest.Verbosity.Detailed));
runner.addPlugin(TestRunProgressPlugin.withVerbosity(matlab.unittest.Verbosity.Detailed, ToFile(logfile)));
runner.addPlugin(FailureDiagnosticsPlugin);
runner.addPlugin(DiagnosticsRecordingPlugin);   % the failure output of each test
results = runner.run(suite);
log = [fileread(logfile), il_summary(results, fileread(records))];

tests = cell(1, numel(results));
docs = containers.Map();
for k = 1:numel(results)
    r = results(k);
    s = suite(k);
    file = fullfile(s.BaseFolder, [strtok(r.Name, '/') '.m']);
    if ~isKey(docs, file)
        docs(file) = il_docs(file);
    end
    d = docs(file);
    fn = extractAfter(r.Name, '/');
    status = 'ok';
    if r.Failed
        status = 'failed';
    elseif r.Incomplete
        status = 'ignored';
    end
    txt = '';
    if isfield(d, fn)
        txt = d.(fn);
    end
    [~, folder] = fileparts(s.BaseFolder);
    tests{k} = struct('name', r.Name, 'suite', folder, 'file', strtok(r.Name, '/'), ...
        'status', status, 'time', r.Duration, 'doc', txt, 'out', il_failure(r));
end

% no pager and no colour: they write terminal codes into the text
[~, git] = system(sprintf('git --no-pager -C "%s" rev-parse --short HEAD', root));
[~, branch] = system(sprintf('git --no-pager -C "%s" rev-parse --abbrev-ref HEAD', root));
git = [strtrim(git) ' ' strtrim(branch)];
[~, dirty] = system(sprintf('git --no-pager -C "%s" status --porcelain', root));
meta = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
    'git', strtrim(git), 'dirty', ~isempty(strtrim(dirty)), 'matlab', version, ...
    'machine', computer('arch'), 'cores', feature('numcores'));
data = struct('meta', meta, 'tests', {tests}, 'records', fileread(records), 'log', log);
json = strrep(jsonencode(data), '</', '<\/');

fid = fopen(fullfile(out, 'data.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', json);
fclose(fid);
page = strrep(fileread(fullfile(here, 'template.html')), '/*__DATA__*/', ['window.REPORT = ' json ';']);
fid = fopen(fullfile(out, 'index.html'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', page);
fclose(fid);

ok = ~any([results.Failed]);
fprintf('report: %s\n', fullfile(out, 'index.html'));
fprintf('tests: %d, failed %d, ignored %d | comparisons recorded: %d\n', numel(results), ...
    nnz([results.Failed]), nnz([results.Incomplete] & ~[results.Failed]), ...
    numel(regexp(fileread(records), '\n')));
end

function d = il_docs(file)
% the comment block right under each test function line: what the test checks
d = struct();
lines = splitlines(fileread(file));
for k = 1:numel(lines)
    tok = regexp(lines{k}, '^function\s+(test_\w+)\s*\(', 'tokens', 'once');
    if isempty(tok)
        continue
    end
    txt = {};
    for j = k+1:numel(lines)
        c = regexp(lines{j}, '^\s*%\s?(.*)$', 'tokens', 'once');
        if isempty(c)
            break
        end
        txt{end+1} = c{1}; %#ok<AGROW>
    end
    d.(tok{1}) = strjoin(txt, ' ');
end
end

function txt = il_summary(results, records)
% the end of the raw log: one line per test with its result and time, then
% one line per recorded comparison with its error, tolerance and utilization
status = repmat({'PASS'}, 1, numel(results));
status([results.Incomplete]) = {'INCOMPLETE'};
status([results.Failed]) = {'FAIL'};
lines = [{'', 'Results', '======='}, ...
    arrayfun(@(k) sprintf('%-10s %8.2f s  %s', status{k}, results(k).Duration, results(k).Name), ...
    1:numel(results), 'UniformOutput', false)];
recs = splitlines(strtrim(records));
recs = recs(~cellfun(@isempty, recs));
lines = [lines, {'', 'Comparisons with a reference', '============================', ...
    sprintf('%-11s %-11s %-11s %s', 'utilization', 'max error', 'tolerance', 'case (test)')}];
for k = 1:numel(recs)
    r = jsondecode(recs{k});
    lines{end+1} = sprintf('%-11.4g %-11.4g %-11.4g %s (%s)', r.util, r.max_abs, r.tol_abs, r.case, r.test); %#ok<AGROW>
end
txt = sprintf('%s\n', lines{:});
end

function txt = il_failure(r)
% the diagnostics of a failed or incomplete test, as the runner prints them
txt = '';
if isfield(r.Details, 'DiagnosticRecord')
    txt = strjoin(arrayfun(@(d) d.Report, r.Details.DiagnosticRecord, 'UniformOutput', false), newline);
end
end
