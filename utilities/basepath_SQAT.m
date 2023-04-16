function bp = basepath_SQAT
% function bp = basepath_SQAT
%
% It returns the main (base) path of the toolbox.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bp = [fileparts(fileparts(mfilename('fullpath'))) filesep];