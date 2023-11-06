function startup_SQAT(bp)
% function startup_SQAT(bp)
%
% This scripts initialises the toolbox. As recommended by the authors, the
%   added paths will be removed from the MATLAB directories when MATLAB is
%   closed. This means that startup_SQAT needs to be run once after MATLAB
%   has started.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bp = [fileparts(mfilename('fullpath')) filesep]; % obtains the folder where this script is
end

main_dirs = { 'examples', ...
              'psychoacoustic_metrics', ...
              'publications', ...
              'sound_files', ...
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
            addpath([bp main_dirs{i_dir}]);
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

PA_More_main       = [bp 'psychoacoustic_metrics' filesep 'PsychoacousticAnnoyance_More2010'    filesep];
PA_More_example    = [bp 'examples'               filesep 'PsychoacousticAnnoyance_More2010'    filesep];

PA_Di_main         = [bp 'psychoacoustic_metrics' filesep 'PsychoacousticAnnoyance_Di2016'      filesep];
PA_Di_example      = [bp 'examples'               filesep 'PsychoacousticAnnoyance_Di2016'      filesep];

PA_Zwicker_main    = [bp 'psychoacoustic_metrics' filesep 'PsychoacousticAnnoyance_Zwicker1999' filesep];
PA_Zwicker_example = [bp 'examples'               filesep 'PsychoacousticAnnoyance_Zwicker1999' filesep];

SLM_example        = [bp 'examples'               filesep 'sound_level_meter' filesep];

EPNL_main        = [bp 'psychoacoustic_metrics' filesep 'EPNL_FAR_Part36' filesep];
EPNL_helper        = [bp 'psychoacoustic_metrics' filesep 'EPNL_FAR_Part36' filesep 'helper' filesep];
EPNL_validation  = [bp 'validation'             filesep 'EPNL_FAR_Part36' filesep];
EPNL_example     = [bp 'examples'               filesep 'EPNL_FAR_Part36' filesep];

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

bAdd = ~exist('PsychoacousticAnnoyance_Zwicker1999.m','file');
if bAdd
    addpath(PA_Zwicker_main);
    addpath(PA_Zwicker_example);
end

bAdd = ~exist('PsychoacousticAnnoyance_More2010.m','file');
if bAdd
    addpath(PA_More_main);
    addpath(PA_More_example);
end

bAdd = ~exist('PsychoacousticAnnoyance_Di2016.m','file');
if bAdd
    addpath(PA_Di_main);
    addpath(PA_Di_example);
end

bAdd = ~exist('ex_sound_level_meter.m','file');
if bAdd
    addpath(SLM_example);
end

bAdd = ~exist('EPNL_FAR_Part36.m','file');
if bAdd
    addpath(EPNL_main);
    addpath(EPNL_validation);
    addpath(EPNL_example);
end

bAdd = ~exist('get_PNL.m','file');
if bAdd
    addpath(EPNL_helper);
end

%%% Adding the publications' directory (alphabetical order):
bAdd = ~exist('pub_Greco2023_Internoise.m','file');
if bAdd
    dir2add = [bp 'publications' filesep 'pub_Greco2023_Internoise' filesep];
    if exist(dir2add,'dir')
        addpath(dir2add)
    end
end

bAdd = ~exist('pub_Osses2023c_Forum_Acusticum_SQAT.m','file');
if bAdd
    dir2add = [bp 'publications' filesep 'pub_Osses2023c_Forum_Acusticum_SQAT' filesep];
    if exist(dir2add,'dir')
        addpath(dir2add)
    end
end
