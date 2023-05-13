function dir_sound = get_dir_validation_sounds(psychoacoustic_model,version)
% function dir_sound = get_dir_validation_sounds(psychoacoustic_model)
%
% It returns the main (base) path of the toolbox.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if version==1
    
dir_sound = [basepath_SQAT 'sound_files' filesep 'validation_SQAT_v1_0' filesep psychoacoustic_model filesep];
else
end

if ~exist(dir_sound,'dir')
    dir_sound = '';
    fprintf('%s.m: The directory %s was not found. Maybe the sound files will still be found....\n',mfilename,dir_sound);
end
