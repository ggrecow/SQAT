# SQAT: a sound quality analysis toolbox for MATLAB
This is the repository of SQAT, an open-source **S**ound **Q**uality **A**nalysis **T**oolbox for MATLAB.

# Toolbox structure
The toolbox has the following directories:
- `psychoacoustic_metrics`: this directory contains a number of algorithms implementing a specific psychoacoustic metric. 
- `sound_level_meter`: contains scripts to obtain sound pressure levels in dB(A), dB(C), or dB(Z), using fast, slow or impulse time weightings.
- `utilities`: contains some scripts that are complementary to any of the toolbox functions.
- `examples`: an example script is provided for each metric.
- `sound_files`: this directory hosts reference sounds in .wav format that are used mainly by the `examples` codes. 
- `validation`: this directory contains scripts used to validate each algorithm. In order to reproduce the validation codes, the dataset of test sounds needs to be downloaded from Zenodo ([link](https://doi.org/10.5281/zenodo.7933206)). The obtained folder called `validation_SQAT_v1_0` has to be included in the `sound_files` folder of the toolbox. 
- `publications`: scripts reproducing figures and/or tables for publications using SQAT are provided in this directory.

# How to use the toolbox
After downloading this repository, you just need to add the toolbox into the path of your MATLAB. The `startup_SQAT` code provided can be used to automatically include all folders to the MATLAB path, until the MATLAB session ends. However, in order to avoid conflicts, the `startup_SQAT` needs to be used every time MATLAB is (re)started. If you just want to use the metrics and no `sound_files`, you can add manually only the relevant folders to the path (e.g., `psychoacoustic_metrics`, `sound_level_meter` and `utilities`). 

# Sound quality metrics available in SQAT v1.0
The folowing psychoacoustic-based metrics are available in the `psychoacoustic_metrics` folder :
- `Loudness_ISO532_1`: Zwicker loudness model according to ISO 532-1:2016 (see validation [here](validation/Loudness_ISO532_1) and example [here](examples/Loudness_ISO532_1/ex_Loudness_ISO532_1.m)).
- `Sharpness_DIN45692`: Sharpness according to DIN 45692:2009 (see validation [here](validation/Sharpness_DIN45692) and example [here](examples/Sharpness_DIN45692/ex_Sharpness_DIN45692.m)). 
- `Roughness_Daniel1997`: Roughness model from Daniel & Weber (see validation [here](validation/Roughness_Daniel1997) and example [here](examples/Roughness_Daniel1997/ex_Roughness_Daniel1997.m)).  
- `FluctuationStrength_Osses2016`: Fluctuation strength model from Osses et al. (see validation [here](validation/FluctuationStrength_Osses2016) and example [here](examples/FluctuationStrength_Osses2016/ex_FluctuationStrength_Osses2016.m)).   
- `Tonality_Aures1985`: Tonality model from Aures (see validation [here](validation/Tonality_Aures1985) and example [here](examples/Tonality_Aures1985/ex_Tonality_Aures1985.m)).
- `PsychoacousticAnnoyance_Zwicker1999`: psychoacoustic annoyance model from Zwicker et al. (see example [here](examples/PsychoacousticAnnoyance_Zwicker1999/ex_PsychoacousticAnnoyance_Zwicker1999.m)).
- `PsychoacousticAnnoyance_More2010`: psychoacoustic annoyance model from More (see example [here](examples/PsychoacousticAnnoyance_More2010/ex_PsychoacousticAnnoyance_More2010.m)).
- `PsychoacousticAnnoyance_Di2016`: psychoacoustic annoyance model from Di et al. (see example [here](examples/PsychoacousticAnnoyance_Di2016/ex_PsychoacousticAnnoyance_Di2016.m)).

The following SPL-based metrics using different frequency weightings (A, C, or Z) and time weightings (fast, slow, or impulse) can be calculated using the codes available in `sound_level_meter` folder (see example [here](examples/sound_level_meter/ex_sound_level_meter.m)):

- Sound pressure level over time.
- Equivalent sound pressure level.
- Maximum sound pressure level.
- Sound exposure level.

# How to cite this repository
This repository can be cited as follows: 

G. Felix Greco, R. Merino-Martínez, & A. Osses (2023). "SQAT: A sound quality analysis toolbox for MATLAB." url: [https://github.com/ggrecow/sqat](https://github.com/ggrecow/sqat)

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.


