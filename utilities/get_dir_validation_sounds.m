function dir_sound = get_dir_validation_sounds(psychoacoustic_model)
% function dir_sound = get_dir_validation_sounds(psychoacoustic_model)
%
% It returns the main (base) path of the toolbox.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_sound = [basepath_SQAT 'sound_files' filesep 'validation' filesep psychoacoustic_model filesep];

if ~exist(dir_sound,'dir')
    dir_sound = '';
    fprintf('%s.m: The directory %s was not found. Maybe the sound files will still be found....\n',mfilename,dir_sound);
end
