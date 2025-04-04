# About this script

The `pub_Lotinga2025_Forum_Acusticum_ECMA418_2.m` code generates the figures of the following publication:

- Lotinga, M., Torjussen, M., & Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS). **To be presented at Forum Acusticum in June 2025**.

The study presented in the paper compares results obtained using a commercial software (ref. results) and using the implementations in SQAT. To perform the study, two signals are used: 

1) a binaural file recorded in a 'train station' environment. The signal 'TrainStation.7.wav' was extracted from the EigenScape database [(link)](https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between 01m00s and 01m30s. The EigenScape database, which is described by Green et al. [2], is licensed under Creative Commons Attribution 4.0.  

2)  ambisonic recording of a 'park' environment with unmanned aircraft system (UAS / drone) flight overhead (25 seconds, 2-channel binaural). The signal 'Park.3.wav' was extracted from the EigenScape database and trimmed between 00m02s and 00m27s. The original 4th order ambisonic recording was reduced to 2nd order for playback over a 16-channel array. An auralised UAS was superimposed on the recording, as described by Lotinga et al (https://doi.org/10.1038/s44384-024-00001-6) using the software of Green (https://github.com/acoustics-code-salford/uas-sound-propagation). The overall sound was re-recorded in binaural using a head-and-torso-simulator, with pre-equalisation applied.

In this publication, selected results from the verification studies of the sound quality metrics from ECMA-418-2 [1] are presented. Additional verification studies are available in the verification folder of SQAT [(link)](../../validation), and in the refmap repository [(link)](https://github.com/acoustics-code-salford/refmap-psychoacoustics/tree/main/validation/ECMA-418-2).

# How to use this code

To run the `pub_Lotinga2025_Forum_Acusticum_ECMA418_2.m` script, in addition to have initialised the SQAT toolbox, the two signals are required. 

Signal 1, called `ExStereo_TrainStation7-0100-0130.wav`, is stored in the `sound_files/reference_signals` folder [(here)](../../sound_files/reference_signals), while the reference results are stored in the respective validation folders of each ECMA-418-1 metric implementation in the `validation` folder (for example, see reference results from loudness [here](../../validation/Loudness_ECMA418_2/2_Loudness_ECMA418_2_software_comparison)). 

Signal 2, called `ExStereo_Park3-0002-0027_UAS.wav` and associated reference results are stored in a dedicated zenodo repository (DOI: [10.5281/zenodo.15132459](https://doi.org/10.5281/zenodo.15132460)). After downloading the `pub_Lotinga2025_Forum_Acusticum_ECMA418_2` file from zenodo, unzip it and run the provided `copy_data_for_pub_Lotinga2025_to_SQAT.m` script to automatically place the files necessary data in the correct folders. Alternatively, the `data` folder inside the downloaded `pub_Lotinga2025_Forum_Acusticum_ECMA418_2` file can be placed manually in this folder to run the `pub_Lotinga2025_Forum_Acusticum_ECMA418_2.m` script.


# References
[1] Ecma International. (2024). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 3rd Edition/December 2024). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf) (Last viewed 22 Jan 2025)

[2] Green, M. C., & Murphy, D. (2017). EigenScape: A Database of Spatial Acoustic Scene Recordings. [Applied sciences](https://doi.org/10.3390/app7111204), 7(11), 1004.  

# Log

Last checked: Gil Felix Greco (13.02.2025), Mike Lotinga (19.03.2025), Gil Felix Greco (04.04.2025) 

