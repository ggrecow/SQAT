% script pub_Lotinga2025_Forum_Acusticum_ECMA418_2
%
% Generates the figures of the following contribution:
%
% Lotinga, M., Torjussen, M., and Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS 
% OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED 
% SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS). 
% Forum Acusticum.
%
% The signals used are:
%
% Stereo signal #1: binaural audio recording of a 'train station' environment (30 seconds, 2-channel binaural)
% The signal 'TrainStation.7.wav' was extracted from the EigenScape database 
% (https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between
%  01m00s and 01m30s. The EigenScape database, which is described by 
% Green et al (https://doi.org/10.3390/app7111204), is licensed 
% under Creative Commons Attribution 4.0.
% - Signal label: <ExStereo_TrainStation7-0100-0130.wav>
%
% Stereo signal #2: ambisonic recording of a 'park' environment with
% unmanned aircraft system (UAS / drone) flight overhead (25 seconds,
% 2-channel binaural)
% The signal 'Park.3.wav' was extracted from the EigenScape database 
% and trimmed between 00m02s and 00m27s. The original 4th order ambisonic
% recording was reduced to 2nd order for playback over a 16-channel array.
% An auralised UAS was superimposed on the recording, as described by
% Lotinga et al (https://doi.org/10.1038/s44384-024-00001-6) using the
% software of Green
% (https://github.com/acoustics-code-salford/uas-sound-propagation).
% The overall sound was re-recorded in binaural using a
% head-and-torso-simulator, with pre-equalisation applied.
% - Signal label: <ExStereo_Park3-0002-0027_UAS.wav>
%
% - The signals are stored in the following folder: <sound_files\reference_signals\>. 
%
% Author: Gil Felix Greco, Braunschweig 27.02.2025
% Modified: 03.04.2025 Mike Lotinga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Load .wav file

dir_sound1 = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];
dir_sound2 = [basepath_SQAT 'publications' filesep 'pub_Lotinga2025_Forum_Acusticum_ECMA418_2' filesep 'data' filesep 'Audio' filesep];

fileTag = 'ExStereo_';

wav_fileTrain = 'TrainStation7-0100-0130';
[insigTrain, fs] = audioread([dir_sound1 fileTag wav_fileTrain '.wav']); 

wav_fileUAS = 'Park3-0002-0027_UAS';
[insigUAS, ~] = audioread([dir_sound2 fileTag wav_fileUAS '.wav']);

%% Compute sound quality metrics from .wav signals

% Common settings
fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')
% NOTE: fieldtype is set to match capability in reference software
% (recorded soundfields would not match free frontal incidence assumption)

% Train station
% Roughness
OUTTrain.roughness = Roughness_ECMA418_2( insigTrain, fs, fieldtype );

% Tonality 
OUTTrain.tonality = Tonality_ECMA418_2( insigTrain, fs, fieldtype);

% Loudness
OUTTrain.loudness = Loudness_ECMA418_2( insigTrain, fs, fieldtype );

% UAS over park
% Roughness
OUTUAS.roughness = Roughness_ECMA418_2( insigUAS, fs, fieldtype );

% Tonality 
OUTUAS.tonality = Tonality_ECMA418_2( insigUAS, fs, fieldtype);

% Loudness
OUTUAS.loudness = Loudness_ECMA418_2( insigUAS, fs, fieldtype );

%% Plot spectrograms - channel 1 only
label_fig = [wav_fileTrain ' (Channel 1)' '_Spectro'];
plt_Spectro(insigTrain(:, 1), fs, 4, 0.5, true, save_figs, label_fig)

label_fig = [wav_fileUAS ' (Channel 1)' '_Spectro'];
plt_Spectro(insigUAS(:, 1), fs, 4, 0.5, true, save_figs, label_fig)

%% Load reference results

ref_resultsTrain = get_reference_results('TrainStation7-0100-0130');
ref_resultsUAS = get_reference_results('Park3-0002-0027_UAS');

%% plot train station time-dependent quantities (only channel 1)

% Roughness 
metric = 'roughness';
xRef = ref_resultsTrain.roughness.TDep(:,1);
yRef = ref_resultsTrain.roughness.TDep(:,2);
xImplementation = OUTTrain.roughness.timeOut;
yImplementation = OUTTrain.roughness.roughnessTDep(:,1);
label_fig = [wav_fileTrain ' (Channel 1)' '_TDep_Roughness'];

plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)

% Tonality
metric = 'tonality';
xRef = ref_resultsTrain.tonality.TDep(:,1);
yRef = ref_resultsTrain.tonality.TDep(:,2);
xImplementation = OUTTrain.tonality.timeOut;
yImplementation = OUTTrain.tonality.tonalityTDep(:,1);
label_fig = [wav_fileTrain ' (Channel 1)' '_TDep_Tonality'];

plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)

% Loudness
metric = 'loudness';
xRef = ref_resultsTrain.loudness.TDep(:,1);
yRef = ref_resultsTrain.loudness.TDep(:,2);
xImplementation = OUTTrain.loudness.timeOut;
yImplementation = OUTTrain.loudness.loudnessTDep(:,1);
label_fig = [wav_fileTrain ' (Channel 1)' '_TDep_Loudness'];

% call plot function
plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)

%% plot train station time-dependent specific tonality (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_resultsTrain.tonality.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_resultsTrain.tonality.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_resultsTrain.tonality.Spec_TDep_channel_1(2:end, 2:end) ; % specific tonality

% implementation (Channel 1)
x2Axis =  OUTTrain.tonality.timeOut; % time vector
y2Axis =  OUTTrain.tonality.bandCentreFreqs; % freq vector
z2Axis =  OUTTrain.tonality.specTonality(:,:,1) ; % specific tonality

label_fig = [wav_fileTrain ' (Channel 1)' '_tDep_Specific_Tonality'];

% call plot function
plt_tDepSpecific(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', 'tonality', save_figs, label_fig)

%% plot train station time-dependent specific loudness (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_resultsTrain.loudness.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_resultsTrain.loudness.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_resultsTrain.loudness.Spec_TDep_channel_1(2:end, 2:end) ; % specific loudness

% implementation (Channel 1)
x2Axis =  OUTTrain.loudness.timeOut; % time vector
y2Axis =  OUTTrain.loudness.bandCentreFreqs; % freq vector
z2Axis =  OUTTrain.loudness.specLoudness(:,:,1) ; % specific loudness

label_fig = [wav_fileTrain ' (Channel 1)' '_tDep_Specific_Loudness'];

% call plot function
plt_tDepSpecific(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', 'loudness', save_figs, label_fig)

%% plot train station time-dependent specific roughness (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_resultsTrain.roughness.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_resultsTrain.roughness.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_resultsTrain.roughness.Spec_TDep_channel_1(2:end, 2:end) ; % specific roughness

% implementation (Channel 1)
x2Axis =  OUTTrain.roughness.timeOut; % time vector
y2Axis =  OUTTrain.roughness.bandCentreFreqs; % freq vector
z2Axis =  OUTTrain.roughness.specRoughness(:,:,1) ; % specific roughness

label_fig = [wav_fileTrain ' (Channel 1)' '_tDep_Specific_Roughness'];

% call plot function
plt_tDepSpecific(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', 'roughness', save_figs, label_fig)

%% plot train station time-averaged specific tonality (only channel 1)

yRef = ref_resultsTrain.tonality.AvgSpec(:,2);
yImplementation = OUTTrain.tonality.specTonalityAvg(:,1);
label_fig = [wav_fileTrain ' (Channel 1)' '_avgSpecific_Tonality'];

% call plot function
plt_AggSpecific( yRef, yImplementation, 'tonality', save_figs, label_fig)

%% plot train station time-averaged specific loudness (only channel 1)

yRef = ref_resultsTrain.loudness.AvgSpec(:,2);
yImplementation = OUTTrain.loudness.specLoudnessPowAvg(:,1);
label_fig = [wav_fileTrain ' (Channel 1)' '_avgSpecific_Loudness'];

% call plot function
plt_AggSpecific( yRef, yImplementation, 'loudness', save_figs, label_fig)

%% plot train station time-averaged specific roughness (only channel 1)

yRef = ref_resultsTrain.roughness.AvgSpec(:,2);
yImplementation = OUTTrain.roughness.specRoughnessAvg(:,1);
label_fig = [wav_fileTrain ' (Channel 1)' '_avgSpecific_Roughness'];

% call plot function
plt_AggSpecific( yRef, yImplementation, 'roughness', save_figs, label_fig)

%% plot - train station overall values

% tonality
ch1_ref_tonality = ref_resultsTrain.tonality.single_values(1);
ch1_implementation_tonality = OUTTrain.tonality.tonalityAvg(1);

single_values1 = [ch1_ref_tonality, ...
                             ch1_implementation_tonality];
% loudness
ch1_ref_loudness = ref_resultsTrain.loudness.single_values(1);
ch1_implementation_loudness = OUTTrain.loudness.loudnessPowAvg(1);

single_values2 = [ch1_ref_loudness, ...
                             ch1_implementation_loudness];

% roughness
ch1_ref_roughness = ref_resultsTrain.roughness.single_values(1);
ch1_implementation_roughness = OUTTrain.roughness.roughness90Pc(1);

single_values3 = [ch1_ref_roughness, ...
                             ch1_implementation_roughness];

label_fig = [wav_fileTrain '_singleValues'];

% call plot function
plt_singleValues(single_values1, single_values2, single_values3, save_figs, label_fig)


%% plot UAS over park time-dependent quantities (only channel 1)

% Roughness 
metric = 'roughness';
xRef = ref_resultsUAS.roughness.TDep(:,1);
yRef = ref_resultsUAS.roughness.TDep(:,2);
xImplementation = OUTUAS.roughness.timeOut;
yImplementation = OUTUAS.roughness.roughnessTDep(:,1);
label_fig = [wav_fileUAS ' (Channel 1)' '_TDep_Roughness'];

plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)

% Tonality
metric = 'tonality';
xRef = ref_resultsUAS.tonality.TDep(:,1);
yRef = ref_resultsUAS.tonality.TDep(:,2);
xImplementation = OUTUAS.tonality.timeOut;
yImplementation = OUTUAS.tonality.tonalityTDep(:,1);
label_fig = [wav_fileUAS ' (Channel 1)' '_TDep_Tonality'];

plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)

% Loudness
metric = 'loudness';
xRef = ref_resultsUAS.loudness.TDep(:,1);
yRef = ref_resultsUAS.loudness.TDep(:,2);
xImplementation = OUTUAS.loudness.timeOut;
yImplementation = OUTUAS.loudness.loudnessTDep(:,1);
label_fig = [wav_fileUAS ' (Channel 1)' '_TDep_Loudness'];

% call plot function
plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs, label_fig)

%% plot UAS over park time-dependent specific tonality (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_resultsUAS.tonality.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_resultsUAS.tonality.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_resultsUAS.tonality.Spec_TDep_channel_1(2:end, 2:end) ; % specific tonality

% implementation (Channel 1)
x2Axis =  OUTUAS.tonality.timeOut; % time vector
y2Axis =  OUTUAS.tonality.bandCentreFreqs; % freq vector
z2Axis =  OUTUAS.tonality.specTonality(:,:,1) ; % specific tonality

label_fig = [wav_fileUAS ' (Channel 1)' '_tDep_Specific_Tonality'];

% call plot function
plt_tDepSpecific(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', 'tonality', save_figs, label_fig)

%% plot UAS over park time-dependent specific loudness (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_resultsUAS.loudness.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_resultsUAS.loudness.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_resultsUAS.loudness.Spec_TDep_channel_1(2:end, 2:end) ; % specific loudness

% implementation (Channel 1)
x2Axis =  OUTUAS.loudness.timeOut; % time vector
y2Axis =  OUTUAS.loudness.bandCentreFreqs; % freq vector
z2Axis =  OUTUAS.loudness.specLoudness(:,:,1) ; % specific loudness

label_fig = [wav_fileUAS ' (Channel 1)' '_tDep_Specific_Loudness'];

% call plot function
plt_tDepSpecific(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', 'loudness', save_figs, label_fig)

%% plot UAS over park time-dependent specific roughness (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_resultsUAS.roughness.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_resultsUAS.roughness.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_resultsUAS.roughness.Spec_TDep_channel_1(2:end, 2:end) ; % specific roughness

% implementation (Channel 1)
x2Axis =  OUTUAS.roughness.timeOut; % time vector
y2Axis =  OUTUAS.roughness.bandCentreFreqs; % freq vector
z2Axis =  OUTUAS.roughness.specRoughness(:,:,1) ; % specific roughness

label_fig = [wav_fileUAS ' (Channel 1)' '_tDep_Specific_Roughness'];

% call plot function
plt_tDepSpecific(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', 'roughness', save_figs, label_fig)

%% plot UAS over park time-averaged specific tonality (only channel 1)

yRef = ref_resultsUAS.tonality.AvgSpec(:,2);
yImplementation = OUTUAS.tonality.specTonalityAvg(:,1);
label_fig = [wav_fileUAS ' (Channel 1)' '_avgSpecific_Tonality'];

% call plot function
plt_AggSpecific( yRef, yImplementation, 'tonality', save_figs, label_fig)

%% plot UAS over park time-averaged specific loudness (only channel 1)

yRef = ref_resultsUAS.loudness.AvgSpec(:,2);
yImplementation = OUTUAS.loudness.specLoudnessPowAvg(:,1);
label_fig = [wav_fileUAS ' (Channel 1)' '_avgSpecific_Loudness'];

% call plot function
plt_AggSpecific( yRef, yImplementation, 'loudness', save_figs, label_fig)

%% plot UAS over park time-averaged specific roughness (only channel 1)

yRef = ref_resultsUAS.roughness.AvgSpec(:,2);
yImplementation = OUTUAS.roughness.specRoughnessAvg(:,1);
label_fig = [wav_fileUAS ' (Channel 1)' '_avgSpecific_Roughness'];

% call plot function
plt_AggSpecific( yRef, yImplementation, 'roughness', save_figs, label_fig)

%% plot - UAS over park overall values

% tonality
ch1_ref_tonality = ref_resultsUAS.tonality.single_values(1);
ch1_implementation_tonality = OUTUAS.tonality.tonalityAvg(1);

single_values1 = [ch1_ref_tonality, ...
                             ch1_implementation_tonality];
% loudness
ch1_ref_loudness = ref_resultsUAS.loudness.single_values(1);
ch1_implementation_loudness = OUTUAS.loudness.loudnessPowAvg(1);

single_values2 = [ch1_ref_loudness, ...
                             ch1_implementation_loudness];

% roughness
ch1_ref_roughness = ref_resultsUAS.roughness.single_values(1);
ch1_implementation_roughness = OUTUAS.roughness.roughness90Pc(1);

single_values3 = [ch1_ref_roughness, ...
                             ch1_implementation_roughness];

label_fig = [wav_fileUAS '_singleValues'];

% call plot function
plt_singleValues(single_values1, single_values2, single_values3, save_figs, label_fig)
