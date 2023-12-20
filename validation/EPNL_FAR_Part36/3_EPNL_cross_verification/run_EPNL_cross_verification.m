% Script run_EPNL_cross_verification
%
%   In this script, a cross-verification between the EPNL results obtained using SQAT
%   and EPNL values provided in the literature is performed.
%
%   The reference results provided in the literature are obtained from Table
%   4 of the following publication:
%
%   [1] Rizzi, S. A., LeGriffon, I., Pieren, R., & Bertsch, L. (2020). A Comparison 
%       of Aircraft Flyover Auralizations by the Aircraft Noise Simulation
%       Working Group. In AIAA AVIATION 2020 FORUM
%       https://doi.org/10.2514/6.2020-2582
%
%   The (auralized) audio files are provided as complementary data from Ref. [1]
%   are used here. They can be freely donwloaded from:
%
%   [2]  https://stabserv.larc.nasa.gov/flyover/?doing_wp_cron=1703067390.0725319385528564453125
%         (Last viewed 20 December, 2023)
%
%   In this script, the EPNL values from the audio files provided in the link above,
%   are computed using SQAT. The values obtained with the SQAT EPNL implementation
%   are compared with the ones provided in Table
%
%   EPNL (and/or associated quantities) are computed using the following function:
%
%   OUT = EPNL_FAR_Part36( insig, fs, method, dt, threshold, show )
%   type <help EPNL_FAR_Part36> for more info
%
%   HOW TO RUN THIS CODE: apart from the <EPNL_FAR_Part36> function provided in SQAT,
%   you need to download the audio files from Ref. [2] in order to run this
%   code. The obtained folder called <AIAA-2020-2582_Sound_Files> has to be
%   included in the sound_files folder of the toolbox (or also loaded properly
%   in this code).
%
% Author: Gil Felix Greco, Braunschweig, 20.12.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% reference values from Table 4 of Ref. [1]

ref = [ 85.6 85.6 86.7
          91.5 91.4 92.4
          79.2 81.4 79.5
          80.5 78.4 79.0
          81.2 81.6 83.3] ;
      
%% load audio files provided by Ref. [1] in Ref. [2]

% Points to <sound_files> folder of SQAT. Thus, if you download the file from Ref. [2] and 
% unzip it in <sound_files> folder, all will work. If you want to place
% the sound files in a different directory, you need to change it here
main_path = [basepath_SQAT 'sound_files' filesep 'AIAA-2020-2582_Sound_Files' filesep];  

% string with the different names (auralization tools) to load .wav files
tool_tag = {'NAF', 'AURAFONE', 'FLAURA'};

% number of cases per auralization tool
nCases = 5;

% load .wav files (these files are already in Pa units, so no calibration is necessary)
[NAF, AURAFONE, FLAURA] = get_audio_files( main_path, nCases, tool_tag );

%% Compute EPNL using SQAT

%  Compute EPNL from loaded.wav files
% (these files are already in Pa units, so no calibration is necessary)
[NAF, AURAFONE, FLAURA] = get_EPNL( NAF, AURAFONE, FLAURA, nCases );

%% Arrange EPNL values computed with SQAT 

EPNL_SQAT = [ NAF{1, 1}.EPNL.EPNL  AURAFONE{1, 1}.EPNL.EPNL   FLAURA{1, 1}.EPNL.EPNL  
                      NAF{1, 2}.EPNL.EPNL  AURAFONE{1, 2}.EPNL.EPNL   FLAURA{1, 2}.EPNL.EPNL  
                      NAF{1, 3}.EPNL.EPNL  AURAFONE{1, 3}.EPNL.EPNL   FLAURA{1, 3}.EPNL.EPNL  
                      NAF{1, 4}.EPNL.EPNL  AURAFONE{1, 4}.EPNL.EPNL   FLAURA{1, 4}.EPNL.EPNL  
                      NAF{1, 5}.EPNL.EPNL  AURAFONE{1, 5}.EPNL.EPNL   FLAURA{1, 5}.EPNL.EPNL  ] ;
                  
%% Compute difference of EPNL results ( SQAT minus Ref [1] )

EPNL_delta = EPNL_SQAT - ref;

fprintf ( '\n\nThe differences between the EPNL values computed with SQAT and the ref. values are\n\n');
disp( round(EPNL_delta,2) );

%% Plot comparison of EPNL results obtained with SQAT with the ones provided by Ref. [1]

% dir where the figures will be stored:
dir_out = [fileparts(mfilename('fullpath')) filesep];
    
% NAF - ref. vs. SQAT results for all cases
plot_EPNL_comparison(ref(:,1), EPNL_SQAT(:,1), tool_tag{1}, save_figs, dir_out);

% AURAFONE - ref. vs. SQAT results for all cases
plot_EPNL_comparison(ref(:,2), EPNL_SQAT(:,2), tool_tag{2}, save_figs, dir_out);

% FLAURA - ref. vs. SQAT results for all cases
plot_EPNL_comparison(ref(:,3), EPNL_SQAT(:,3), tool_tag{3}, save_figs, dir_out);
 