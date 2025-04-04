function OUT = get_reference_results(wav_file)
% function OUT = get_reference_results
%
% Get reference results obtained from a commercial software. The ref.
% results are stored in the validation folders of each individual sound
% quality metric.
%
% Loudness, ref. results folder:
% validation/Loudness_ECMA418_2/2_Loudness_ECMA418_2_software_comparison/reference_results
%
% Roughness, ref. results folder:
% validation/Roughness_ECMA418_2/2_Roughness_ECMA418_2_software_comparison/reference_results
%
% Tonality, ref. results folder:
% validation/Tonality_ECMA418_2/Tonality_ECMA418_2_software_comparison/reference_results
%
% Output
%   OUT : structure
%            single strcuture containing all reference results
%
% Author: Gil Felix Greco, Braunschweig 13.02.2025
% Modified: Mike Lotinga 03/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define directories
switch wav_file
    case  'TrainStation7-0100-0130'
    R_dir = [basepath_SQAT filesep 'validation' filesep 'Roughness_ECMA418_2' filesep '2_Roughness_ECMA418_2_software_comparison' filesep 'reference_results' filesep];
    L_dir = [basepath_SQAT filesep 'validation' filesep 'Loudness_ECMA418_2' filesep '2_Loudness_ECMA418_2_software_comparison' filesep 'reference_results' filesep];
    T_dir = [basepath_SQAT filesep 'validation' filesep 'Tonality_ECMA418_2' filesep 'Tonality_ECMA418_2_software_comparison' filesep 'reference_results' filesep];
    case 'Park3-0002-0027_UAS'
    R_dir = [basepath_SQAT filesep 'publications' filesep 'pub_Lotinga2025_Forum_Acusticum_ECMA418_2' filesep 'data' filesep 'Roughness' filesep];
    L_dir = [basepath_SQAT filesep 'publications' filesep 'pub_Lotinga2025_Forum_Acusticum_ECMA418_2' filesep 'data' filesep 'Loudness' filesep];
    T_dir = [basepath_SQAT filesep 'publications' filesep 'pub_Lotinga2025_Forum_Acusticum_ECMA418_2' filesep 'data' filesep 'Tonality' filesep];
end


%% get reference roughness results

% file names
AvgSpec_fileName = [ wav_file '.Specific Roughness (Hearing Model).asc']; % channel 1 and 2
AvgSpecCombBinaural_fileName = [ wav_file '.Specific Roughness (Hearing Model)_combined_binaural.asc']; % combined binaural

TDep_fileName = [ wav_file '.Roughness (Hearing Model) vs. Time.asc']; % channel 1 and 2
TDepCombBinaural_fileName = [ wav_file '.Roughness (Hearing Model) vs. Time_combined_binaural.asc']; % combined binaural

Spec_TDep_fileName_channel_1 = [ wav_file '.Specific Roughness (Hearing Model) vs. Time_channel_1.asc']; % channel 1
Spec_TDep_fileName_channel_2 = [ wav_file '.Specific Roughness (Hearing Model) vs. Time_channel_2.asc']; % channel 2
Spec_TDep_fileName_combined_binaural = [ wav_file '.Specific Roughness (Hearing Model) vs. Time_combined_binaural.asc']; % combined binaural

% load files
OUT.roughness.AvgSpec = readmatrix( [R_dir AvgSpec_fileName] , 'FileType', 'text');
OUT.roughness.AvgSpecCombBinaural = readmatrix( [R_dir AvgSpecCombBinaural_fileName] , 'FileType', 'text');

OUT.roughness.TDep = readmatrix( [R_dir TDep_fileName] , 'FileType', 'text');
OUT.roughness.TDepCombBinaural = readmatrix( [R_dir TDepCombBinaural_fileName] , 'FileType', 'text');

OUT.roughness.Spec_TDep_channel_1 = readmatrix( [R_dir Spec_TDep_fileName_channel_1] , 'FileType', 'text');
OUT.roughness.Spec_TDep_channel_2 = readmatrix( [R_dir Spec_TDep_fileName_channel_2] , 'FileType', 'text');
OUT.roughness.Spec_TDep_combined_binaural = readmatrix( [R_dir Spec_TDep_fileName_combined_binaural] , 'FileType', 'text');

% single values
% Analysis Name    Channel Name    R/asper    AnalysisRange
% Specific Roughness (Hearing Model) vs. Time    Ch1    0,128    [0,3 - 29,8 s]
% Specific Roughness (Hearing Model) vs. Time    Ch2    0,103    [0,3 - 29,8 s]
% Specific Roughness (Hearing Model) vs. Time    CombinedBinaural Ch1&Ch2    0,115    [0,3 - 29,8 s]

% Analysis Name	Channel Name	Default	{Default}	[Default]	AnalysisRange
% Roughness (Hearing Model) vs. Time	Ch1	R	0.141	asper	[0.3 - 24.9 s]
% Roughness (Hearing Model) vs. Time	Ch2	R	0.107	asper	[0.3 - 24.9 s]
% Roughness (Hearing Model) vs. Time	CombinedBinaural Ch1&Ch2	R	0.125	asper	[0.3 - 24.9 s]

switch wav_file
    case 'TrainStation7-0100-0130'
        OUT.roughness.single_values = [0.128 0.103 0.115]; % [Ch1 Ch2 combBinaural]
    case 'Park3-0002-0027_UAS'
        OUT.roughness.single_values = [0.141 0.107 0.125]; % [Ch1 Ch2 combBinaural]
end

%% get reference tonality results

% file name
AvgSpec_fileName = [ wav_file '.Specific Tonality (Hearing Model).asc'];
TDep_fileName = [ wav_file '.Tonality (Hearing Model) vs. Time.asc'];
Spec_TDep_fileName_channel_1 = [ wav_file '.Specific Tonality (Hearing Model) vs. Time_channel_1.asc'];
Spec_TDep_fileName_channel_2 = [ wav_file '.Specific Tonality (Hearing Model) vs. Time_channel_2.asc'];

% load files
OUT.tonality.AvgSpec = readmatrix( [T_dir AvgSpec_fileName] , 'FileType', 'text');
OUT.tonality.TDep = readmatrix( [T_dir TDep_fileName] , 'FileType', 'text');
OUT.tonality.Spec_TDep_channel_1 = readmatrix( [T_dir Spec_TDep_fileName_channel_1] , 'FileType', 'text');
OUT.tonality.Spec_TDep_channel_2 = readmatrix( [T_dir Spec_TDep_fileName_channel_2] , 'FileType', 'text');

% single values
% Channel Name    T/tuHMS
% Ch1    0,660
% Ch2    0,314

% Analysis Name	Channel Name	Default	{Default}	[Default]	AnalysisRange
% Tonality (Hearing Model) vs. Time	Ch1	T	0.271	tuHMS	
% Tonality (Hearing Model) vs. Time	Ch2	T	0.208	tuHMS	

switch wav_file
    case 'TrainStation7-0100-0130'
        OUT.tonality.single_values = [0.660 0.314]; % [Ch1 Ch2]
    case 'Park3-0002-0027_UAS'
        OUT.tonality.single_values = [0.271 0.208]; % [Ch1 Ch2]
end


%% get reference loudness results

% file name
AvgSpec_fileName = [ wav_file '.Specific Loudness (Hearing Model).asc']; % channel 1 and 2
AvgSpecCombBinaural_fileName = [ wav_file '.Specific Loudness (Hearing Model)_combined_binaural.asc']; % combined binaural

TDep_fileName = [ wav_file '.Loudness (Hearing Model) vs. Time.asc']; % channel 1 and 2
TDepCombBinaural_fileName = [ wav_file '.Loudness (Hearing Model) vs. Time_combined_binaural.asc']; % combined binaural

Spec_TDep_fileName_channel_1 = [ wav_file '.Specific Loudness (Hearing Model) vs. Time_channel_1.asc']; % channel 1
Spec_TDep_fileName_channel_2 = [ wav_file '.Specific Loudness (Hearing Model) vs. Time_channel_2.asc']; % channel 2
Spec_TDep_fileName_combined_binaural = [ wav_file '.Specific Loudness (Hearing Model) vs. Time_combined_binaural.asc']; % combined binaural

% load files
OUT.loudness.AvgSpec = readmatrix( [L_dir AvgSpec_fileName] , 'FileType', 'text');
OUT.loudness.AvgSpecCombBinaural = readmatrix( [L_dir AvgSpecCombBinaural_fileName] , 'FileType', 'text');

OUT.loudness.TDep = readmatrix( [L_dir TDep_fileName] , 'FileType', 'text');
OUT.loudness.TDepCombBinaural = readmatrix( [L_dir TDepCombBinaural_fileName] , 'FileType', 'text');

OUT.loudness.Spec_TDep_channel_1 = readmatrix( [L_dir Spec_TDep_fileName_channel_1] , 'FileType', 'text');
OUT.loudness.Spec_TDep_channel_2 = readmatrix( [L_dir Spec_TDep_fileName_channel_2] , 'FileType', 'text');
OUT.loudness.Spec_TDep_combined_binaural = readmatrix( [L_dir Spec_TDep_fileName_combined_binaural] , 'FileType', 'text');

% single values
% Analysis Name    Channel Name    N/soneHMS    AnalysisRange
% Loudness (Hearing Model) vs. Time    Ch1    8,10    [0,304 - 29,9 s]
% Loudness (Hearing Model) vs. Time    Ch2    6,43    [0,304 - 29,9 s]
% Loudness (Hearing Model) vs. Time    CombinedBinaural Ch1&Ch2    7,32    [0,304 - 29,9 s]

% Analysis Name	Channel Name	Default	{Default}	[Default]	AnalysisRange
% Loudness (Hearing Model) vs. Time	Ch1	N	3.879	soneHMS	[0.304 - 25 s]
% Loudness (Hearing Model) vs. Time	Ch2	N	3.681	soneHMS	[0.304 - 25 s]
% Loudness (Hearing Model) vs. Time	CombinedBinaural Ch1&Ch2	N	3.787	soneHMS	[0.304 - 25 s]

switch wav_file
    case 'TrainStation7-0100-0130'
        OUT.loudness.single_values = [8.10 6.43 7.32]; % [Ch1 Ch2 combBinaural]
    case 'Park3-0002-0027_UAS'
        OUT.loudness.single_values = [3.879 3.681 3.787]; % [Ch1 Ch2 combBinaural]
end

end