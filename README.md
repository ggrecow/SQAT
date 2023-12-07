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

# Sound quality metrics available in SQAT 

<!---The folowing psychoacoustic-based metrics are available in the `psychoacoustic_metrics` folder :
- `Loudness_ISO532_1`: Zwicker loudness model according to ISO 532-1:2017 (see [validation](validation/Loudness_ISO532_1) and [example](examples/Loudness_ISO532_1/ex_Loudness_ISO532_1.m)).
- `Sharpness_DIN45692`: Sharpness according to DIN 45692:2009 (see [validation](validation/Sharpness_DIN45692) and [example](examples/Sharpness_DIN45692/ex_Sharpness_DIN45692.m)). 
- `Roughness_Daniel1997`: Roughness model from Daniel & Weber (see [validation](validation/Roughness_Daniel1997) and [example](examples/Roughness_Daniel1997/ex_Roughness_Daniel1997.m)).  
- `FluctuationStrength_Osses2016`: Fluctuation strength model from Osses *et al.* (see [validation](validation/FluctuationStrength_Osses2016) and [example](examples/FluctuationStrength_Osses2016/ex_FluctuationStrength_Osses2016.m)).   
- `Tonality_Aures1985`: Tonality model from Aures (see [validation](validation/Tonality_Aures1985) and [example](examples/Tonality_Aures1985/ex_Tonality_Aures1985.m)).
- `PsychoacousticAnnoyance_Zwicker1999`: psychoacoustic annoyance model from Zwicker *et al.* (see [example](examples/PsychoacousticAnnoyance_Zwicker1999/ex_PsychoacousticAnnoyance_Zwicker1999.m)).
- `PsychoacousticAnnoyance_More2010`: psychoacoustic annoyance model from More (see [example](examples/PsychoacousticAnnoyance_More2010/ex_PsychoacousticAnnoyance_More2010.m)).
- `PsychoacousticAnnoyance_Di2016`: psychoacoustic annoyance model from Di *et al.* (see [example](examples/PsychoacousticAnnoyance_Di2016/ex_PsychoacousticAnnoyance_Di2016.m)). --->

The implemented metrics available in the `psychoacoustic_metrics` folder are listed in the table below, including the release on which each metric was first introduced:

| Metric  | Model | Implementation | Validation | Example | Release |
| :---: | :---: |:---: |:---: |:---: |:---: |
| Loudness  | ISO 532-1:2017 [1]  |  [link](psychoacoustic_metrics/Loudness_ISO532_1) | [link](validation/Loudness_ISO532_1) | [link](examples/Loudness_ISO532_1/ex_Loudness_ISO532_1.m) | v1.0 |
| Sharpness  |  DIN 45692:2009 [2]  |  [link](psychoacoustic_metrics/Sharpness_DIN45692) | [link](validation/Sharpness_DIN45692) | [link](examples/Sharpness_DIN45692/ex_Sharpness_DIN45692.m) | v1.0 |
| Roughness  | Daniel & Weber [3]  |  [link](psychoacoustic_metrics/Roughness_Daniel1997) |  [link](validation/Roughness_Daniel1997) | [link](examples/Roughness_Daniel1997/ex_Roughness_Daniel1997.m) | v1.0 |
| Fluctuation Strength  | Osses *et al.* [4]  |  [link](psychoacoustic_metrics/FluctuationStrength_Osses2016) |  [link](validation/FluctuationStrength_Osses2016) | [link](examples/FluctuationStrength_Osses2016/ex_FluctuationStrength_Osses2016.m) | v1.0 |
| Tonality  | Aures [5]  |   [link](psychoacoustic_metrics/Tonality_Aures1985) |  [link](validation/Tonality_Aures1985) | [link](examples/Tonality_Aures1985/ex_Tonality_Aures1985.m) | v1.0 |
| Psychoacoustic Annoyance  | Zwicker & Fastl [6]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Zwicker1999) | - | [link](examples/PsychoacousticAnnoyance_Zwicker1999/ex_PsychoacousticAnnoyance_Zwicker1999.m) | v1.0 |
| Psychoacoustic Annoyance  | More [7]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_More2010) | - | [link](examples/PsychoacousticAnnoyance_More2010/ex_PsychoacousticAnnoyance_More2010.m) | v1.0 |
| Psychoacoustic Annoyance  | Di *et al.* [8]  |  [link](psychoacoustic_metrics/PsychoacousticAnnoyance_Di2016) | - |  [link](examples/PsychoacousticAnnoyance_Di2016/ex_PsychoacousticAnnoyance_Di2016.m) | v1.0 |
| EPNL  | FAR Part 36 [9]  |  [link](psychoacoustic_metrics/EPNL_FAR_Part36) | [link](validation/EPNL_FAR_Part36) |  [link](examples/EPNL_FAR_Part36/ex_EPNL_FAR_Part36.m) | v1.1 |

The following SPL-based metrics using different frequency weightings (A, B, C, D or Z) and time weightings (fast, slow, or impulse) can be calculated using the codes available in `sound_level_meter` folder (see [examples](examples/sound_level_meter)):

- Sound pressure level over time.
- Equivalent sound pressure level.
- Maximum sound pressure level.
- Sound exposure level.
- Sound spectrum in 1/3 octave bands.

**References**

[1] International Organization for Standardization. (2017). Acoustics - Methods for calculating loudness - Part 1: Zwicker method (ISO Standard No. 532-1).

[2] Deutsches Institut für Normung. (2009). Measurement technique for the simulation of the auditory sensation of sharpness (DIN Standard No. 45692).

[3] Daniel, P., & Weber, R. (1997). Psychoacoustical Roughness: Implementation of an Optimized Model. [Acta Acustica united with Acustica](https://www.ingentaconnect.com/content/dav/aaua/1997/00000083/00000001/art00020), 83(1), 113-123.

[4] Osses, A., García, R., & Kohlrausch, A. (2016). Modelling the sensation of fluctuation strength. [Proceedings of Meetings on Acoustics](https://doi.org/10.1121/2.0000410), 28(1), 050005.  

[5] Aures, W. (1985). Berechnungsverfahren für den sensorischen Wohlklang beliebiger Schallsignale (A model for calculating the sensory euphony of various sounds). [Acta Acustica united with Acustica](https://www.ingentaconnect.com/content/dav/aaua/1985/00000059/00000002/art00008), 59(2), 130-141.

[6] Zwicker, E., & Fastl, H. (1999). Psychoacoustics: facts and models, Second edition. [Springer-Verlag](https://doi.org/10.1007/978-3-662-09562-1).

[7] More, S. R. (2010). Aircraft noise characteristics and metrics. PhD thesis, [Purdue University](https://docs.lib.purdue.edu/dissertations/AAI3453255/).

[8] Di, G. Q., Chen, X. W., Song, K., Zhou, B., & Pei, C. M. (2016). Improvement of Zwicker’s psychoacoustic annoyance model aiming at tonal noises. [Applied Acoustics](https://doi.org/10.1016/j.apacoust.2015.12.006), 105, 164-170.

[9] Federal Aviation Regulations, 14 CFR Parts 36 and 91, Docket No. FAA-2003-16526; Amendment No. 36-26, 91-288, (2005). [https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036](https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036) (Last viewed 30 Oct 2023)

# How to cite this repository
If you use this toolbox in your research, we would be grateful if you help us to gain visibility by citing SQAT. This is the main citation if you need to cite the toolbox repository itself:

> Felix Greco, G., Merino-Martínez, R., &  Osses, A. (2023). SQAT: a sound quality analysis toolbox for MATLAB. Zenodo. doi: [10.5281/zenodo.7934709](https://doi.org/10.5281/zenodo.7934709)

> [!TIP]
> **Apart from the toolbox repository itself (see above), each version released has its own doi in Zenodo. As differences between releases may occur, it is a good practice to cite the specific SQAT version being used. If you need to cite the current release of SQAT, please refer to the "Cite this repository" feature in the "About" section of this GitHub repository.**

The following paper is the main work describing SQAT and the metrics available in the first release:

> Felix Greco, G., Merino-Martínez, R., Osses, A., & Langer, S. C. (2023). SQAT: a MATLAB-based toolbox for quantitative sound quality analysis. INTER-NOISE and NOISE-CON Congress and Conference Proceedings, InterNoise23, Chiba, Japan. doi: [10.3397/IN_2023_1075](https://doi.org/10.3397/IN_2023_1075)

Additionaly, here's a paper by the members of the SQAT team showing three case studies where the SQAT toolbox was used to perform all analyses:

> Osses, A., Felix Greco, G., & Merino-Martínez, R. (2023). Considerations for the perceptual evaluation of steady-state and time-varying sounds using psychoacoustic metrics. Forum Acusticum, Turin, Italy, 11-15 September 2023. [[Link]](https://hal.science/hal-04186316). 
Raw data and extra scripts to reproduce all the paper figures can be found [here](https://doi.org/10.5281/zenodo.7933489).

# Works using SQAT

We would be very happy to know that you find SQAT useful and have used it in your own work. In this case, please reach out so that we can feature your work here. 

> Brandetti, L., Mulders, S. P., Merino-Martinez, R., Watson, S., & van Wingerden, J.-W. (2023). Multi-objective calibration of vertical-axis wind turbine controllers: balancing aero-servo-elastic performance and noise, [Wind Energy Science Discussions](https://doi.org/10.5194/wes-2023-154) [preprint].

> Pockelé, J. (2023). Auralisation of modelled wind turbine noise for psychoacoustic listening experiments: development and validation of the wind turbine auralisation tool WinTAur. Master's thesis, [Delft University of Technology](http://resolver.tudelft.nl/uuid:cc9e67b4-6bde-4114-97c0-43b11b4a48ef).

> Knuth, D., Ring, T. P., & Langer, S. C. (2023). Utilizing auralization to investigate psychoacoustic perception of vibrating structures. [Proceedings of 49. Jahrestagung für Akustik (DAGA)](https://pub.dega-akustik.de/DAGA_2023/data/articles/000414.pdf).

# Toolbox history

Prompted by the lack of well-documented and validated open-source implementations of sound quality metrics, the development of SQAT was initiated in 2019 by [Gil Felix Greco](https://www.linkedin.com/in/gil-felix-greco-363985101/) with the main objective of applying it to the sound quality analysis of environmental noise generated by virtual aircraft designs. In 2021, a preliminary version of SQAT, which works only with specific inputs provided by an aircraft noise prediction research software (PANAM-DLR), was first applied to analyze the sound quality of a novel aircraft concept (see publication below). 

> Felix Greco, G., Bertsch, L., Ring, T. P., & Langer, S. C. (2021). Sound quality assessment of a medium-range aircraft with enhanced fan-noise shielding design. [CEAS Aeronautical Journal](https://doi.org/10.1007/s13272-021-00515-9) 12, 481–493.

In October 2022, after inspiring discussions during an acoustics conference in South Korea, [Alejandro Osses](https://www.linkedin.com/in/alejandro-osses-10039883/) and [Roberto Merino-Martínez](https://www.linkedin.com/in/roberto-merino-martinez/) joined the SQAT team in order to foster the development and release of an open-source version of the toolbox. The first release of the open-source version of SQAT (v1.0), which works with any calibrated sound file as input, is freely available for use since May 2023 in this GitHub repository. 

The toolbox is expected to be carefully maintained by the SQAT team in the long-term, thus supporting its use as a reliable and readily accessible tool. Moreover, additional metrics are expected to be included in future releases as long as enough verification and documentation related to their implementations are available. With this strong requirement, we hope to only make publicly available implementations that are already mature enough for general use and therefore won't need to be updated after released. Although this shall allow for the backward compatibility between different releases, this may not always be the case as unforeseen issues may be reported by the users after a release. Thus, please always make sure to double check the description of the most up-to-date release in case you are using different versions of the toolbox.       

# Contact

If you would like to get in touch to report a bug, make suggestions or ask a question, please contact us on Github by opening an issue.

# Licensing
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.


