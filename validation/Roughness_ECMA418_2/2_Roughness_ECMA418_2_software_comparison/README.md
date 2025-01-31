# About this code 

The `Roughness_ECMA418_2_software_comparison.m` code compares roughness results (ECMA-418-2 model [1]) obtained in a commercial software and using the implementation in SQAT (see `Roughness_ECMA418_2.m` code [here](../../../psychoacoustic_metrics/Roughness_ECMA418_2/Roughness_ECMA418_2.m)). To perform this study, a binaural file recorded in a 'train station' environment is used. The signal 'TrainStation.7.wav' was extracted from the EigenScape database [(link)](https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between 01m00s and 01m30s. The EigenScape database, which is described by 
Green et al. [2], is licenced under Creative Commons Attribution 4.0.   

<!-- Artemis Suite 15.7-->

# How to use this code
This is a standalone code. Therefore, no extra steps are required to run it.

# Results

## Time-dependent roughness

<img src='figs/TrainStation7-0100-0130 (Channel 1)_TDep_Roughness.png' width=500> <img src='figs/TrainStation7-0100-0130 (Channel 2)_TDep_Roughness.png' width=500> <img src='figs/TrainStation7-0100-0130 (Combined binaural)_TDep_Roughness.png' width=500>

## Time-averaged specific roughness

<img src='figs/TrainStation7-0100-0130 (Channel 1)_avgSpecific_Roughness.png' width=1000> 

<img src='figs/TrainStation7-0100-0130 (Channel 2)_avgSpecific_Roughness.png' width=1000> 

<img src='figs/TrainStation7-0100-0130 (Combined binaural)_avgSpecific_Roughness.png' width=1000>

## Overall roughness

<img src='figs/TrainStation7-0100-0130_singleValues_Roughness.png' width=500>

# References
[1] Ecma International. (2024). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 3rd Edition/December 2024). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf) (Last viewed 22 Jan 2025)

[2] Green, M. C., & Murphy, D. (2017). EigenScape: A Database of Spatial Acoustic Scene Recordings. [Applied sciences](https://doi.org/10.3390/app7111204), 7(11), 1004.  


# Log
Created by Gil Felix Greco (31.01.2025)



