function startup_SQAT(bp)
% function startup_SQAT(bp)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bp = [fileparts(mfilename('fullpath')) filesep]; % obtains the folder where this script is
    % 'bp' in Alejandro's computer: /home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI_sim/MATLAB/tb_Brazuca/
end

main_dirs = { 'psychoacoustic_metrics', ...
              'sound_level_meter', ...
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

FS_main       = [bp 'psychoacoustic_metrics' filesep 'FluctuationStrength_Osses2016' filesep];
FS_validation = [bp 'validation'             filesep 'FluctuationStrength_Osses2016' filesep];

% R_main  = [bp 'psychoacoustic_metrics' filesep 'Roughness_Daniel1997'          filesep];
% % L_main  = [bp 'psychoacoustic_metrics' filesep 'Loudness_Chalupper2002'        filesep];
% L_main  = [bp 'psychoacoustic_metrics' filesep 'Loudness_ISO532_1'             filesep];
% S_main  = [bp 'psychoacoustic_metrics' filesep 'Sharpness_DIN45692'            filesep];
% T_main  = [bp 'psychoacoustic_metrics' filesep 'Tonality_Aures1985'            filesep];
% % fs_tue_basepath = '/home/alejandro/Documents/MATLAB/fluctuation-strength-TUe/';
bAdd = ~exist('FluctuationStrength_Osses2016.m','file');
if bAdd
    addpath(FS_main);
    addpath(FS_validation);
end
% bAdd = ~exist('Roughness_Daniel1997.m','file');
% if bAdd
%     addpath(R_main);
% end
% % bAdd = ~exist('Loudness_Chalupper2002.m','file');
% bAdd = ~exist('Loudness_ISO532_1.m','file');
% if bAdd
%     addpath(L_main);
% end
% bAdd = ~exist('Sharpness_DIN45692.m','file');
% if bAdd
%     addpath(S_main);
% end
% bAdd = ~exist('Tonality_Aures1985.m','file');
% if bAdd
%     addpath(T_main);
% end
 
% bAdd = ~exist('rmsdb.m','file');
% if bAdd
%     addpath([bp 'utilities' filesep]);
% end
