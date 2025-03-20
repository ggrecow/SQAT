# About this script

The `pub_Lotinga2025_Forum_Acusticum_ECMA418_2.m` code generates the figures of the following publication:

- Lotinga, M., Torjussen, M., & Felix Greco, G. (2025) VERIFIED IMPLEMENTATIONS OF THE SOTTEK PSYCHOACOUSTIC HEARING MODEL STANDARDISED SOUND QUALITY METRICS (ECMA-418-2 LOUDNESS, TONALITY AND ROUGHNESS). **To be presented at Forum Acusticum in June 2025**.

In this publication, selected results from the verification studies of the sound quality metrics from ECMA-418-2 [1] are presented. The complete study is available in the verification folder of the toolbox [(link)](../../validation).
 
The study presented in the paper compares results obtained using a commercial software (ref. results) and using the implementations in SQAT. To perform the study, a binaural file recorded in a 'train station' environment is used. The signal 'TrainStation.7.wav' was extracted from the EigenScape database [(link)](https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between 01m00s and 01m30s. The EigenScape database, which is described by 
Green et al. [2], is licensed under Creative Commons Attribution 4.0.  

# How to use this code
To run the `pub_Lotinga2025_Forum_Acusticum_ECMA418_2.m` code, the following sound file is required: `ExStereo_TrainStation7-0100-0130.wav`, which is stored in the `sound_files\reference_signals` folder. 

# References
[1] Ecma International. (2024). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 3rd Edition/December 2024). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf) (Last viewed 22 Jan 2025)

[2] Green, M. C., & Murphy, D. (2017). EigenScape: A Database of Spatial Acoustic Scene Recordings. [Applied sciences](https://doi.org/10.3390/app7111204), 7(11), 1004.  

# Log

Last checked: Gil Felix Greco (13.02.2025) Mike Lotinga (19.03.2025

