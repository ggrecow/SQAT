function S = sqat_run_script(file__, copy__)
% function S = sqat_run_script(file, copy)
%
%   Runs a validation script of validation/ in a workspace of its own (the
%   scripts start with clear all, which would clear the variables of a test)
%   and returns every variable it leaves, plus
%     S.run_err : the message of the error that stopped it, '' if none
%     S.run_out : what it printed
%   (named so that no variable of a script overwrites them)
%   The figures it opens are closed. With copy true, a copy of the script
%   runs from a temporary folder, so that the results it stored next to
%   itself (and would ask to reuse) are not found and everything is computed.

if nargin > 1 && copy__
    dir__ = tempname;
    mkdir(dir__);
    copyfile(file__, dir__);
    [~, n__, e__] = fileparts(file__);
    file__ = fullfile(dir__, [n__ e__]);
    % kept out of the workspace, which the script clears: the copy to remove
    % and the folder to come back to (a script may cd into its own folder)
    setappdata(groot, 'sqat_run_script', struct('dir', dir__, 'here', pwd));
end
try
    out__ = evalc('run(file__)');
    err__ = '';                              % set after the script, which clears all
catch e__
    out__ = '';
    err__ = e__.message;
end
names__ = who;
S = struct();
for k__ = 1:numel(names__)
    if ~endsWith(names__{k__}, '__')
        S.(names__{k__}) = eval(names__{k__});
    end
end
S.run_err = err__;
S.run_out = out__;
close all
if isappdata(groot, 'sqat_run_script')
    ctx__ = getappdata(groot, 'sqat_run_script');
    rmappdata(groot, 'sqat_run_script');
    cd(ctx__.here);
    rmdir(ctx__.dir, 's');
end
end
