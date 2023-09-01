# About this code

The `pub_Greco2023_Internoise.m` code generates the figures of the following publication:

- Felix Greco, G., Merino-Martínez, R., Osses, A., & Langer, S. C. (2023). SQAT: a MATLAB-based toolbox for quantitative sound quality analysis. In Proceedings of the INTER–NOISE and NOISE–CON Congress, the 52nd International Congress and Exhibition on Noise Control Engineering, Chiba, Greater Tokyo, Japan, 20–23 August 2023. Institute of Noise Control Engineering: Reston, VA, USA, 2023. [[Link]](https://www.researchgate.net/publication/373334884_SQAT_a_MATLAB-based_toolbox_for_quantitative_sound_quality_analysis). 

The `pub_Greco2023_Internoise.m` is a wraper to the codes used to verify the psychoacoustic metrics implemented in SQAT v1.0. The following codes are called by this wraper:
- `Loudness_ISO532_1`: verification of loudness implementation according to ISO 532-1:2016 using synthetic stationary signals (see code [here](../../validation/Loudness_ISO532_1/1_synthetic_signals_stationary_loudness)) and time-varying technical signals (see code [here](../../validation/Loudness_ISO532_1/3_technical_signals_time_varying_loudness)).
- `Sharpness_DIN45692`: verification of sharpness implementation according to DIN 45692:2009 using narrowband and broadband signals (see code [here](../../validation/Sharpness_DIN45692)).
- `Roughness_Daniel1997`: verification of Daniel & Weber roughness model implementation using AM tones as a function of the modulation depth (see code [here](../../validation/Roughness_Daniel1997/2_AM_modulation_depth)) and modulation frequency (see code [here](../../validation/Roughness_Daniel1997/1_AM_modulation_freq)).  
- `FluctuationStrength_Osses2016`: verification of Osses et al. fluctuation strength model implementation using AM tones as a function of the modulation frequency (see code [here](../../validation/FluctuationStrength_Osses2016/1_AM_tones_fmod)) and AM broadband noises as a function of the modulation frequency (see code [here](../../validation/FluctuationStrength_Osses2016/2_AM_BBN_fmod)).   
- `Tonality_Aures1985`: verification of Aures' tonality model implementation using pure tones with different SNR from narrowband noises (see code [here](../../validation/Tonality_Aures1985)).

# How to use this code
To run the `pub_Greco2023_Internoise.m` code, the user needs to download the dataset of sound files used by the verification codes from zenodo <a href="https://doi.org/10.5281/zenodo.7933206" target="_blank">here</a>. The obtained folder called `validation_SQAT_v1_0` has to be included in the `sound_files` folder of the toolbox. 

# Log
Code released in SQAT v1.0, 14.05.2023

