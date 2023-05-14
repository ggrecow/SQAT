function dir_validation = get_dir_validation(psychoacoustic_model)
% function dir_validation = get_dir_validation(psychoacoustic_model)
%
% It returns the main (base) path of the toolbox.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_validation = [basepath_SQAT 'validation' filesep psychoacoustic_model filesep];

if ~exist(dir_validation,'dir')
    warning('%s does not seem to exist',dir_validation);
end