# About this code 

The `Roughness_ECMA418_2_software_comparison.m` code compares roughness results (ECMA-418-2 model [1]) obtained using a commercial software (ref. results) and using the implementation in SQAT (see `Roughness_ECMA418_2.m` code [here](../../../psychoacoustic_metrics/Roughness_ECMA418_2/Roughness_ECMA418_2.m)). To perform this study, a binaural file recorded in a 'train station' environment is used. The signal 'TrainStation.7.wav' was extracted from the EigenScape database [(link)](https://zenodo.org/doi/10.5281/zenodo.1012808), and trimmed between 01m00s and 01m30s. The EigenScape database, which is described by 
Green et al. [2], is licenced under Creative Commons Attribution 4.0.   

<!-- Artemis Suite 17.1-->

# How to use this code
The following sound file is required: `ExStereo_TrainStation7-0100-0130.wav`, which is stored in the `sound_files\reference_signals` folder. 

# Results

## Time-dependent roughness

| Channel 1       | Channel 2          |
| -------------- | -------------- |
| <img src='figs/TrainStation7-0100-0130 (Channel 1)_TDep_Roughness.png' width=500>  | <img src='figs/TrainStation7-0100-0130 (Channel 2)_TDep_Roughness.png' width=500>  |

Combined binaural |  
:-------------------------:| 
| <img src='figs/TrainStation7-0100-0130 (Combined binaural)_TDep_Roughness.png' width=500>|  
    

## Time-averaged specific roughness

Channel 1  |  
:-------------------------:| 
| <img src='figs/TrainStation7-0100-0130 (Channel 1)_avgSpecific_Roughness.png' width=1000>  | 

Channel 2 |  
:-------------------------:| 
| <img src='figs/TrainStation7-0100-0130 (Channel 2)_avgSpecific_Roughness.png' width=1000>  | 

Combined binaural |  
:-------------------------:| 
| <img src='figs/TrainStation7-0100-0130 (Combined binaural)_avgSpecific_Roughness.png' width=1000> | 
 
## Time-dependent specific roughness

| Reference (channel 1)       | Implementation (channel 1)          |
| -------------- | -------------- |
| <img src='figs/TrainStation7-0100-0130 (Channel 1)_tDep_Specific_Roughness_ref.png' width=500>  | <img src='figs/TrainStation7-0100-0130 (Channel 1)_tDep_Specific_Roughness_implementation.png' width=500>   |

|   Reference (channel 2)       | Implementation (channel 2)          |
| -------------- | -------------- |
| <img src='figs/TrainStation7-0100-0130 (Channel 2)_tDep_Specific_Roughness_ref.png' width=500>  | <img src='figs/TrainStation7-0100-0130 (Channel 2)_tDep_Specific_Roughness_implementation.png' width=500>  |

## Overall roughness

<img src='figs/TrainStation7-0100-0130_singleValues_Roughness.png' width=500>

# References
[1] Ecma International. (2025). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 4th Edition/June 2025). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf) (Last viewed 16 Nov 2025)

[2] Green, M. C., & Murphy, D. (2017). EigenScape: A Database of Spatial Acoustic Scene Recordings. [Applied sciences](https://doi.org/10.3390/app7111204), 7(11), 1004.  


# Log
Last checked: Mike Lotinga (16.11.2025)



