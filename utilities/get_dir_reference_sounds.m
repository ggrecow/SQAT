function dir_ref_sounds = get_dir_reference_sounds(psychoacoustic_model)
% function dir_ref_sounds = get_dir_reference_sounds(psychoacoustic_model)
%
% It returns the main (base) path of the toolbox.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    dir_analysis = '';
end

dir_ref_sounds = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep psychoacoustic_model filesep];

if ~exist(dir_ref_sounds,'dir')
    dir_ref_sounds = '';
    fprintf('%s.m: The directory %s was not found. Maybe the folder for ref. sounds will still be found....\n',mfilename,dir_ref_sounds);
end
