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
% Stereo signal: binaural audio recording of a 'train station' environment (30 seconds, 2-channel binaural)
% The signal 'TrainStation.7.wav' was extracted from the EigenScape database 
% (https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between
%  01m00s and 01m30s. The EigenScape database, which is described by 
% Green et al (https://doi.org/10.3390/app7111204), is licenced 
% under Creative Commons Attribution 4.0. 
% - Signal label: <ExStereo_TrainStation7-0100-0130.wav>
% - The signal is stored in the following folder: <sound_files\reference_signals\>. 
%
% Author: Gil Felix Greco, Braunschweig 13.02.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % save figure flag

%% Load .wav file

dir_sound = [basepath_SQAT 'sound_files' filesep 'reference_signals' filesep];

fileTag = 'ExStereo_';
wav_file = 'TrainStation7-0100-0130';

[insig, fs] = audioread([dir_sound fileTag wav_file '.wav']); 

%% Compute sound quality metrics from .wav signal

% Common settings
fieldtype = 'free-frontal'; % string (default: 'free-frontal'; or 'diffuse')

% Roughness  
OUT.roughness = Roughness_ECMA418_2( insig, fs, fieldtype );

% Tonality 
OUT.tonality = Tonality_ECMA418_2( insig, fs, fieldtype);

% Loudness
OUT.loudness = Loudness_ECMA418_2( insig, fs, fieldtype );

%% Load reference results

ref_results = get_reference_results;

%% plot time-dependent quantities (only channel 1)

% Roughness 
metric = 'roughness';
xRef = ref_results.roughness.TDep(:,1);
yRef = ref_results.roughness.TDep(:,2);
xImplementation = OUT.roughness.timeOut;
yImplementation = OUT.roughness.roughnessTDep(:,1);
label_fig = [wav_file ' (Channel 1)' '_TDep_Roughness'];

plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs,  label_fig)

% Tonality
metric = 'tonality';
xRef = ref_results.tonality.TDep(:,1);
yRef = ref_results.tonality.TDep(:,2);
xImplementation = OUT.tonality.timeOut;
yImplementation = OUT.tonality.tonalityTDep(:,1);
label_fig = [wav_file ' (Channel 1)' '_TDep_Tonality'];

plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs,  label_fig)

% Loudness
metric = 'loudness';
xRef = ref_results.loudness.TDep(:,1);
yRef = ref_results.loudness.TDep(:,2);
xImplementation = OUT.loudness.timeOut;
yImplementation = OUT.loudness.loudnessTDep(:,1);
label_fig = [wav_file ' (Channel 1)' '_TDep_Loudness'];

% call plot function
plt_tDep(xRef, yRef, xImplementation, yImplementation, metric, save_figs,  label_fig)

%% plot time-dependent specific tonality (only channel 1)

% commercial software (Channel 1)
x1Axis =  ref_results.tonality.Spec_TDep_channel_1(2:end,1) ; % time vector
y1Axis =  ref_results.tonality.Spec_TDep_channel_1(1, 2:end) ; % freq vector
z1Axis =  ref_results.tonality.Spec_TDep_channel_1(2:end, 2:end) ; % specific loudness

% implementation (Channel 1)
x2Axis =  OUT.tonality.timeOut; % time vector
y2Axis =  OUT.tonality.bandCentreFreqs; % freq vector
z2Axis =  OUT.tonality.specTonality(:,:,1) ; % specific loudness

label_fig = [wav_file ' (Channel 1)' '_tDep_Specific_Tonality'];

% call plot function
plt_TonalitySpectrogram(x1Axis, y1Axis, z1Axis', x2Axis, y2Axis, z2Axis', save_figs,  label_fig)

%% plot time-averaged specific tonality (only channel 1)

yRef = ref_results.tonality.AvgSpec(:,2);
yImplementation = OUT.tonality.specTonalityAvg(:,1);
label_fig = [wav_file ' (Channel 1)' '_avgSpecific_Tonality'];

% call plot function
plt_TonalityAvgSpecific( yRef, yImplementation, save_figs, label_fig)

%% plot - overall values

% tonality
ch1_ref_tonality = ref_results.tonality.single_values(1);
ch1_implementation_tonality = OUT.tonality.tonalityAvg(1);

single_values1 = [ch1_ref_tonality, ...
                             ch1_implementation_tonality];
% loudness
ch1_ref_loudness = ref_results.loudness.single_values(1);
ch1_implementation_loudness = OUT.loudness.loudnessPowAvg(1);

single_values2 = [ch1_ref_loudness, ...
                             ch1_implementation_loudness];

% roughness
ch1_ref_roughness = ref_results.roughness.single_values(1);
ch1_implementation_roughness = OUT.roughness.roughness90Pc(1);

single_values3 = [ch1_ref_roughness, ...
                             ch1_implementation_roughness];

label_fig = [wav_file '_singleValues'];

% call plot function
plt_singleValues(single_values1, single_values2, single_values3, save_figs,  label_fig)
