# About this code

This code generates the tables and figures of the following publication: **Felix Greco, G., Merino-Martinez, R., Osses, A., Langer, S. C.** SQAT: a MATLAB-based toolbox for quantitative sound quality. In Proceedings of the INTER–NOISE and NOISE–CON Congress, the 52nd International Congress and Exhibition on Noise Control Engineering, Chiba, Japan, 20–23 August 2023. 

# How to use this code
This code is a wraper to the codes used to verify the psychoacoustic metrics implemented in SQAT v1.0. The following codes are called by this fuction:
- `Loudness_ISO532_1`: verification using synthetic stationary signals (see code [here](../../validation/Loudness_ISO532_1/1_synthetic_signals_stationary_loudness) and time-varying technical signals (see code [here](../../validation/Loudness_ISO532_1/3_technical_signals_time_varying_loudness)).
- `Sharpness_DIN45692`: verification using narrowband and broadband signals (see code [here](../../validation/Sharpness_DIN45692)).
- `Roughness_Daniel1997`: verification using AM tones as a function of the modulation depth (see code [here](../../validation/Roughness_Daniel1997/2_AM_modulation_depth) and modulation frequency (see code [here](../../validation/Roughness_Daniel1997/1_AM_modulation_freq)).  
- `FluctuationStrength_Osses2016`: verification using AM tones as a function of the modulation frequency (see code [here](../../validation/FluctuationStrength_Osses2016/1_AM_tones_fmod) and AM broadband noises as a function of the modulation frequency (see code [here](../../validation/FluctuationStrength_Osses2016/2_AM_BBN_fmod)).   
- `Tonality_Aures1985`: verification using pure tones with different SNR from narrowband noises (see code [here](../../validation/Tonality_Aures1985).

In order to run the validation codes, the user needs to download the dataset of sound files from zenodo <a href="https://doi.org/10.5281/zenodo.7933206" target="_blank">here</a>. The obtained folder called `validation_SQAT_v1_0` has to be included in the `sound_files` folder of the toolbox. 

This code is part of SQAT v1.0, released 14.05.2023

