function startup_SQAT(bp)
% function startup_SQAT(bp)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bp = [fileparts(mfilename('fullpath')) filesep]; % obtains the folder where this script is
    % 'bp' in Alejandro's computer: /home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI_sim/MATLAB/tb_Brazuca/
end

main_dirs = { 'examples', ...
              'psychoacoustic_metrics', ...
              'sound_level_meter', ...
              'utilities', ...
              'validation'};

for i_dir = 1:length(main_dirs)
    main_dir_i = [bp main_dirs{i_dir} filesep];
    if exist(main_dir_i,'dir')
        bAdd = 1;
        switch main_dirs{i_dir}
            case 'sound_level_meter'
                bAdd = ~exist('Do_SLM.m','file'); % avoids adding the folder twice
        end
        if bAdd
            addpath(main_dirs{i_dir});
            fprintf('%s.m: Main dir added to path:\n \t%s\n', mfilename, main_dir_i);
        end
    end
end

L_main        = [bp 'psychoacoustic_metrics' filesep 'Loudness_ISO532_1'             filesep];
L_validation  = [bp 'validation'             filesep 'Loudness_ISO532_1'             filesep];
L_example     = [bp 'examples'               filesep 'Loudness_ISO532_1'             filesep];

R_main        = [bp 'psychoacoustic_metrics' filesep 'Roughness_Daniel1997'          filesep];
R_validation  = [bp 'validation'             filesep 'Roughness_Daniel1997'          filesep];
R_example     = [bp 'examples'               filesep 'Roughness_Daniel1997'          filesep];

FS_main       = [bp 'psychoacoustic_metrics' filesep 'FluctuationStrength_Osses2016' filesep];
FS_validation = [bp 'validation'             filesep 'FluctuationStrength_Osses2016' filesep];
FS_example    = [bp 'examples'               filesep 'FluctuationStrength_Osses2016' filesep];

S_main        = [bp 'psychoacoustic_metrics' filesep 'Sharpness_DIN45692'            filesep];
S_validation  = [bp 'validation'             filesep 'Sharpness_DIN45692'            filesep];
S_example     = [bp 'examples'               filesep 'Sharpness_DIN45692'            filesep];

T_main        = [bp 'psychoacoustic_metrics' filesep 'Tonality_Aures1985'            filesep];
T_validation  = [bp 'validation'             filesep 'Tonality_Aures1985'            filesep];
T_example     = [bp 'examples'               filesep 'Tonality_Aures1985'            filesep];

% L_main  = [bp 'psychoacoustic_metrics' filesep 'Loudness_Chalupper2002'        filesep];
% % fs_tue_basepath = '/home/alejandro/Documents/MATLAB/fluctuation-strength-TUe/';
bAdd = ~exist('Loudness_ISO532_1.m','file');
if bAdd
    addpath(L_main);
    addpath(L_validation);
    addpath(L_example);
end

bAdd = ~exist('Roughness_Daniel1997.m','file');
if bAdd
    addpath(R_main);
    addpath(R_validation);
    addpath(R_example);
end

bAdd = ~exist('FluctuationStrength_Osses2016.m','file');
if bAdd
    addpath(FS_main);
    addpath(FS_validation);
    addpath(FS_example);
end
% bAdd = ~exist('Loudness_Chalupper2002.m','file');

bAdd = ~exist('Sharpness_DIN45692.m','file');
if bAdd
    addpath(S_main);
    addpath(S_validation);
    addpath(S_example);
end

bAdd = ~exist('Tonality_Aures1985.m','file');
if bAdd
    addpath(T_main);
    addpath(T_validation);
    addpath(T_example);
end
 
% bAdd = ~exist('rmsdb.m','file');
% if bAdd
%     addpath([bp 'utilities' filesep]);
% end
