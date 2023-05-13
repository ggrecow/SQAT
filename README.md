# SQAT: A sound quality analysis toolbox for MATLAB
This is the repository of SQAT, a **S**ound **Q**uality **A**nalysis **T**oolbox for MATLAB.

# Toolbox structure
The toolbox has the following directories:
- **psychoacoustic_metrics**: This is the main directory of the toolbox. It contains a number of algorithms implementing a specific psychoacoustic metric. 
- **examples**: An example script is provided for each metric.
- **validation**: This directory contains scripts used to validate each algorithm. In order to reproduce the validation codes, the dataset of test sounds needs to be downloaded from Zenodo ([link](https://doi.org/10.5281/zenodo.7933206)) and the paste <validation_SQAT_v1_0> has to be included in the <sound_files> folder of the toolbox.
- **sound_files**: This directory hosts specific wav files that are used in the toolbox. Due to the size of the subfolders, some of the sounds need to be downloaded from Zenodo and placed in the corresponding folders (instructions will be uploaded soon)
- **sound_level_meter**: Contains scripts to obtain sound pressure levels in dB(A), dB(C), or dB(Z), using fast or slow time weightings.
- **utilities**: Contains some scripts that are complementary to any of the toolbox functions.


# How to cite this repository
This repository can be cited as follows: 

G. Felix Greco, R. Merino-Martínez, & A. Osses (2023). "SQAT: A sound quality analysis toolbox for MATLAB." url: [https://github.com/ggrecow/sqat](https://github.com/ggrecow/sqat)

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.


