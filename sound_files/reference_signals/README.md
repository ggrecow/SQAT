# About this folder

The `reference_signals` folder contains the sound files used by the examples (see [folder](../../examples)) and by some validation scripts (see [folder](../../validation)) of SQAT. Most of them are the reference signals of the psychoacoustic models implemented in SQAT, i.e. the signal that, by definition of each model, yields a value of 1 in its unit.

# Calibration of the sound files

The sound files are stored in .wav format, and therefore their amplitudes are limited to the range [-1, +1]. To obtain sound pressure values in pascals, the full-scale level (dBFS) of each file needs to be known, that is, the sound pressure level (in dB SPL) corresponding to a sample with an amplitude of 1. The SQAT convention is a full-scale level of 94 dB SPL, meaning that an amplitude of 1 corresponds to 1 Pa. A sound file stored according to this convention can be used directly as input to the SQAT functions, without any calibration.

The table below gives the full-scale level of each file. If you use these files outside the provided example scripts (for instance in the SQAT graphical user interface), set the calibration of each file to the dBFS value given below.

| Sound file | Content | Level | dBFS | fs (Hz) | Duration (s) |
| :--- | :--- | :---: | :---: | :---: | :---: |
| `RefSignal_Loudness_ISO532_1.wav` | 1 kHz pure tone | 40 dB SPL | 94 | 48000 | 5 |
| `RefSignal_Loudness_ECMA418_2.wav` | 1 kHz pure tone | 40 dB SPL | 94 | 48000 | 5 |
| `RefSignal_Tonality_ECMA418_2.wav` | 1 kHz pure tone | 40 dB SPL | 94 | 48000 | 5 |
| `RefSignal_Tonality_Aures1985.wav` | 1 kHz pure tone | 60 dB SPL | 94 | 44100 | 5 |
| `RefSignal_Sharpness_DIN45692.wav` | 1 kHz pure tone | 60 dB SPL | **90** | 48000 | 5 |
| `RefSignal_Roughness_Daniel1997.wav` | 1 kHz tone, 100% amplitude-modulated at 70 Hz | 60 dB SPL | 94 | 48000 | 5 |
| `RefSignal_Roughness_ECMA418_2.wav` | 1 kHz tone, 100% amplitude-modulated at 70 Hz | 60 dB SPL | 94 | 48000 | 5 |
| `RefSignal_FluctuationStrength_Osses2016.wav` | 1 kHz tone, 100% amplitude-modulated at 4 Hz | 60 dB SPL | 94 | 44100 | 5 |
| `ExSignal_A320_auralized_departure_104dBFS.wav` | Auralised departure of an A320 aircraft | - | 104 | 44100 | 40 |
| `ExStereo_TrainStation7-0100-0130.wav` | Binaural recording of a train station environment (2 channels), from [1] | - | see notes | 48000 | 30 |

Notes:

- `ExSignal_A320_auralized_departure_104dBFS.wav` was auralised with a full-scale level of 104 dB SPL to avoid clipping (see `ex_EPNL_FAR_Part36.m`).
- `ExStereo_TrainStation7-0100-0130.wav` has **no known calibration**. It was obtained from the recording `TrainStation.7.wav` of the EigenScape database [(link)](https://zenodo.org/doi/10.5281/zenodo.1012808), trimmed between 01m00s and 01m30s. The EigenScape recordings are 4th-order ambisonic recordings that were not calibrated to absolute sound pressure levels [1], and the file provided here is a binaural (2-channel) version of it. In the examples, validation and publication scripts of SQAT, the file is used as it is, i.e. with the SQAT convention of 94 dBFS. The resulting levels (about 76 and 68 dB SPL for channels 1 and 2) are therefore a consequence of this convention, and not the level of the original sound scene. The EigenScape database is licensed under Creative Commons Attribution 4.0.

# References

[1] Green, M. C., & Murphy, D. (2017). EigenScape: A Database of Spatial Acoustic Scene Recordings. Applied Sciences, 7(11), 1204. DOI: [10.3390/app7111204](https://doi.org/10.3390/app7111204)

[2] Deutsches Institut für Normung. (2009). Measurement technique for the simulation of the auditory sensation of sharpness (DIN Standard No. 45692).

# Log

- This README was written by Gil Felix Greco with the aid of an AI assistant (Claude Opus 5.5, Anthropic), 15.10.2026. The levels, durations and sampling frequencies in the table were measured from the sound files themselves.

- Reminder: before the release of SQAT v2.0, renaming the files to include their full-scale level (e.g. `RefSignal_Loudness_ISO532_1_94dBFS.wav`) was considered. Documenting the full-scale levels in this README was chosen instead, to avoid renaming the files and breaking the scripts of users that load them by name.
