function sqat_report_record(metric, case_name, got, want, tol_abs, source)
% function sqat_report_record(metric, case_name, got, want, tol_abs, source)
%
%   Records the measured error of one comparison for the test report
%   (sqat_report_build). A test calls it before it verifies the same values.
%   It does nothing unless the environment variable SQAT_REPORT_FILE names
%   the file to append to, so a plain runtests is unchanged.
%
%   metric    : block of the report the comparison belongs to (e.g. 'slm')
%   case_name : what was compared, in a few words
%   got, want : measured and reference values (same size)
%   tol_abs   : largest absolute error the test accepts
%   source    : where the tolerance comes from (a standard, a paper, or
%               'set by us' with the reason), shown next to it in the report
%
%   Each record holds the worst absolute and relative error and the
%   utilization, the measured error divided by the allowed error (1 = at the
%   limit, 0 = identical to the reference).

file = getenv('SQAT_REPORT_FILE');
if isempty(file)
    return
end
if nargin < 6
    source = '';
end
want = want(:);
err = abs(got(:) - want);
big = abs(want) >= 1e-6 * max(abs(want));          % relative error only where the reference is not near zero
rel = err(big) ./ abs(want(big));
[suite, test] = il_caller();
r = struct('suite', suite, 'test', test, 'metric', metric, 'case', case_name, ...
    'n', numel(err), 'max_abs', max([err; 0]), 'max_rel', max([rel; 0]), ...
    'mode', 'abs', 'tol_abs', tol_abs, 'tol_rel', 0, 'util', max([err; 0]) / tol_abs, 'source', source);
fid = fopen(file, 'a', 'n', 'UTF-8');
fprintf(fid, '%s\n', jsonencode(r));
fclose(fid);
end

function [suite, test] = il_caller()
% the test function that called, as runtests names it (tFile/test_name), and its folder
suite = '';
test = '';
for s = dbstack('-completenames')'
    k = strfind(s.name, 'test_');
    if ~isempty(k) && k(1) == 1 || contains(s.name, '>test_')
        [folder, file] = fileparts(s.file);
        [~, suite] = fileparts(folder);
        test = [file '/' regexprep(s.name, '^.*>', '')];
        return
    end
end
end
