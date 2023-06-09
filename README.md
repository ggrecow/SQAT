# SQAT: a sound quality analysis toolbox for MATLAB
This is the repository of **SQAT**, an open-source **S**ound **Q**uality **A**nalysis **T**oolbox for MATLAB. It contains a collection of codes implementing key metrics for **quantitative** sound quality analysis. With **SQAT** you can conduct quick quantitative sound quality analysis on any calibrated input sound file, in Pascal units. To give a transparent indication of how close the implementations are from the original models, we provide a detailed set of verification routines. Moreover, a number of example codes and exemplary sound files is provided in order to facilitate the initial use of the algorithms.  

# Toolbox structure
The toolbox has the following directories:
- `psychoacoustic_metrics`: this directory contains a number of algorithms implementing a specific psychoacoustic metric (see [folder](psychoacoustic_metrics)). 
- `sound_level_meter`: contains scripts to obtain sound pressure levels in dB(A), dB(C), or dB(Z), using fast, slow or impulse time weightings (see [folder](sound_level_meter)). 
- `utilities`: contains some scripts that are complementary to any of the toolbox functions (see [folder](utilities)).
- `examples`: an example script is provided for each metric (see [folder](examples)).
- `sound_files`: this directory hosts reference sounds in .wav format that are used mainly by the `examples` codes (see [folder](sound_files)). 
- `validation`: this directory contains scripts used to validate each algorithm. Instructions on how to run these codes are provided in each respective folder and directly on the header of the codes (see [folder](validation)). 
- `publications`: contains scripts to reproduce figures and/or tables of publications from the toolbox authors (see [folder](publications)). 

# How to use the toolbox
After downloading this repository, you just need to add the toolbox into the path of your MATLAB. The `startup_SQAT` code provided can be used to automatically include all folders to the MATLAB path, until the MATLAB session ends. In order to avoid conflicts, the `startup_SQAT` needs to be used every time MATLAB is (re)started. If you just want to use the metrics, you can add manually only the relevant folders to the MATLAB path (e.g., `psychoacoustic_metrics`, `sound_level_meter` and `utilities`). 

# Sound quality metrics available in SQAT v1.0

<!---The folowing psychoacoustic-based metrics are available in the `psychoacoustic_metrics` folder :
- `Loudness_ISO532_1`: Zwicker loudness model according to ISO 532-1:2017 (see [validation](validation/Loudness_ISO532_1) and [example](examples/Loudness_ISO532_1/ex_Loudness_ISO532_1.m)).
- `Sharpness_DIN45692`: Sharpness according to DIN 45692:2009 (see [validation](validation/Sharpness_DIN45692) and [example](examples/Sharpness_DIN45692/ex_Sharpness_DIN45692.m)). 
- `Roughness_Daniel1997`: Roughness model from Daniel & Weber (see [validation](validation/Roughness_Daniel1997) and [example](examples/Roughness_Daniel1997/ex_Roughness_Daniel1997.m)).  
- `FluctuationStrength_Osses2016`: Fluctuation strength model from Osses *et al.* (see [validation](validation/FluctuationStrength_Osses2016) and [example](examples/FluctuationStrength_Osses2016/ex_FluctuationStrength_Osses2016.m)).   
- `Tonality_Aures1985`: Tonality model from Aures (see [validation](validation/Tonality_Aures1985) and [example](examples/Tonality_Aures1985/ex_Tonality_Aures1985.m)).
- `PsychoacousticAnnoyance_Zwicker1999`: psychoacoustic annoyance model from Zwicker *et al.* (see [example](examples/PsychoacousticAnnoyance_Zwicker1999/ex_PsychoacousticAnnoyance_Zwicker1999.m)).
- `PsychoacousticAnnoyance_More2010`: psychoacoustic annoyance model from More (see [example](examples/PsychoacousticAnnoyance_More2010/ex_PsychoacousticAnnoyance_More2010.m)).
- `PsychoacousticAnnoyance_Di2016`: psychoacoustic annoyance model from Di *et al.* (see [example](examples/PsychoacousticAnnoyance_Di2016/ex_PsychoacousticAnnoyance_Di2016.m)). --->

The implemented metrics available in the `psychoacoustic_metrics` folder are listed in the table below:

| Metric  | Model | Implementation | Validation | Example |
| :---: | :---: |:---: |:---: |:---: |
| Loudness  | ISO 532-1:2017 [1]  |  [link](psychoacoustic_metrics/Loudness_ISO532_1) | [link](validation/Loudness_ISO532_1) | [link](examples/Loudness_ISO532_1/ex_Loudness_ISO532_1.m) |
| Sharpness  |  DIN 45692:2009 [2]  |  [link](psychoacoustic_metrics/Sharpness_DIN45692) | [link](validation/Sharpness_DIN45692) | [link](examples/Sharpness_DIN45692/ex_Sharpness_DIN45692.m) |
| Roughness  | Daniel & Weber [3]  |  [link](psychoacoustic_metrics/Roughness_Daniel1997) |  [link](validation/Roughness_Daniel1997) | [link](examples/Roughness_Daniel1997/ex_Roughness_Daniel1997.m) |
| Fluctuation Strength  | Osses *et al.* [4]  |  [link](psychoacoustic_metrics/FluctuationStrength_Osses2016) |  [link](validation/FluctuationStrength_Osses2016) | [link](examples/FluctuationStrength_Osses2016/ex_FluctuationStrength_Osses2016.m) |
| Tonality  | Aures [5]  |   [link](psychoacoustic_metrics/Tonality_Aures1985) |  [link](validation/Tonality_Aures1985) | [link](examples/Tonality_Aures1985/ex_Tonality_Aures1985.m) |
| Psychoacoustic Annoyance  | Zwicker & Fastl [6]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Zwicker1999) | - | [link](examples/PsychoacousticAnnoyance_Zwicker1999/ex_PsychoacousticAnnoyance_Zwicker1999.m) |
| Psychoacoustic Annoyance  | More [7]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_More2010) | - | [link](examples/PsychoacousticAnnoyance_More2010/ex_PsychoacousticAnnoyance_More2010.m)|
| Psychoacoustic Annoyance  | Di *et al.* [8]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Di2016) | - |  [link](examples/PsychoacousticAnnoyance_Di2016/ex_PsychoacousticAnnoyance_Di2016.m)|

The following SPL-based metrics using different frequency weightings (A, C, or Z) and time weightings (fast, slow, or impulse) can be calculated using the codes available in `sound_level_meter` folder (see [example](examples/sound_level_meter/ex_sound_level_meter.m)):

- Sound pressure level over time.
- Equivalent sound pressure level.
- Maximum sound pressure level.
- Sound exposure level.

**References**

[1] International Organization for Standardization. (2017). Acoustics - Methods for calculating loudness - Part 1: Zwicker method (ISO Standard No. 532-1).

[2] Deutsches Institut für Normung. (2009). Measurement technique for the simulation of the auditory sensation of sharpness (DIN Standard No. 45692).

[3] Daniel, P., & Weber, R. (1997). Psychoacoustical Roughness: Implementation of an Optimized Model. [Acta Acustica united with Acustica](https://www.ingentaconnect.com/content/dav/aaua/1997/00000083/00000001/art00020), 83(1), 113-123.

[4] Osses Vechi, A., García León, R., & Kohlrausch, A. (2016). Modelling the sensation of fluctuation strength. [Proceedings of Meetings on Acoustics](https://doi.org/10.1121/2.0000410), 28(1), 050005.  

[5] Aures, W. (1985). Berechnungsverfahren für den sensorischen Wohlklang beliebiger Schallsignale (A model for calculating the sensory euphony of various sounds). [Acta Acustica united with Acustica](https://www.ingentaconnect.com/content/dav/aaua/1985/00000059/00000002/art00008), 59(2), 130-141.

[6] Zwicker, E., & Fastl, H. (1999). Psychoacoustics: facts and models, Second edition. [Springer-Verlag](https://doi.org/10.1007/978-3-662-09562-1).

[7] More, S. R. (2010). Aircraft noise characteristics and metrics. PhD thesis, [Purdue University](https://docs.lib.purdue.edu/dissertations/AAI3453255/).

[8] Di, G. Q., Chen, X. W., Song, K., Zhou, B., & Pei, C. M. (2016). Improvement of Zwicker’s psychoacoustic annoyance model aiming at tonal noises. [Applied Acoustics](https://doi.org/10.1016/j.apacoust.2015.12.006), 105, 164-170.


# How to cite this repository
If you use this toolbox in your research, please help us to gain visibility by citing SQAT:

Gil Felix Greco, Roberto Merino-Martínez, & Alejandro Osses. (2023). SQAT: a sound quality analysis toolbox for MATLAB. Zenodo. [https://doi.org/10.5281/zenodo.7934709]( https://doi.org/10.5281/zenodo.7934709)

If you need to cite  the current release of SQAT, please refer to the "Cite this repository" feature in the "About" section of this GitHub repository.


<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.


