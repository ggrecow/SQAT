# About this code 
The `validation_synthetic_signals_time_varying.m` code is used to verify the loudness implementation according to ISO 532-1 [1] (see `Loudness_ISO532_1`code [here](../../../psychoacoustic_metrics/Loudness_ISO532_1/Loudness_ISO532_1.m)). The verification of the time-varying loudness method is performed considering the synthetic test sounds provided in the Annex B.4 of the ISO standard:

- Test signal 6 (pure tone, $f_{\mathrm{c}}=250~\mathrm{Hz}$ and $L_{\mathrm{p}}=30-80~\mathrm{dB~SPL}$)
- Test signal 7 (pure tone, $f_{\mathrm{c}}=1~\mathrm{kHz}$ and $L_{\mathrm{p}}=30-80~\mathrm{dB~SPL}$)
- Test signal 8 (pure tone, $f_{\mathrm{c}}=4~\mathrm{kHz}$ and $L_{\mathrm{p}}=30-80~\mathrm{dB~SPL}$)
- Test signal 9 (pink noise, $L_{\mathrm{p}}=0-50~\mathrm{dB~SPL}$)
- Test signal 10 (10 ms tone pulse, $f_{\mathrm{c}}=1~\mathrm{kHz}$ and $L_{\mathrm{p}}=70~\mathrm{dB~SPL}$)
- Test signal 11 (50 ms tone pulse, $f_{\mathrm{c}}=1~\mathrm{kHz}$ and $L_{\mathrm{p}}=70~\mathrm{dB~SPL}$)
- Test signal 12 (500 ms tone pulse, $f_{\mathrm{c}}=1~\mathrm{kHz}$, $L_{\mathrm{p}}=70~\mathrm{dB~SPL}$)
- Test signal 13 (Combined tone pulses, $f_{\mathrm{c}}=1~\mathrm{kHz}$)

Here, the center frequency $f_{\mathrm{c}}$ and sound pressure level $L_{\mathrm{p}}$ are used to describe the signals.

# How to use this code
In order to run this code and reproduce the figures available in the `figs` folder, the user needs to download the dataset of sound files from zenodo <a href="https://doi.org/10.5281/zenodo.7933206" target="_blank">here</a>. The obtained folder called `validation_SQAT_v1_0` has to be included in the `sound_files` folder of the toolbox. 

# Results
The figures below compare the results obtained using the `Loudness_ISO532_1` implementation in SQAT with the tolerance reference values from the ISO standard. 

| Test signal 6   | Test signal 6 (specific loudness @ 2.5 Bark)         |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_6.png)   | ![](figs/validation_time_varying_loudness_signal_6_specific_loudness.png)  |

| Test signal 7   | Test signal 7 (specific loudness @ 8.5 Bark)        |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_7.png)   | ![](figs/validation_time_varying_loudness_signal_7_specific_loudness.png)  |

| Test signal 8    | Test signal 8 (specific loudness @ 17.5 Bark)         |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_8.png)   | ![](figs/validation_time_varying_loudness_signal_8_specific_loudness.png)  |

| Test signal 9    | Test signal 9 (specific loudness @ 17.5 Bark)         |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_9.png)   | ![](figs/validation_time_varying_loudness_signal_9_specific_loudness.png)  |

| Test signal 10   | Test signal 10 (specific loudness @ 8.5 Bark)         |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_10.png)   | ![](figs/validation_time_varying_loudness_signal_10_specific_loudness.png)  |

| Test signal 11   | Test signal 11 (specific loudness @ 8.5 Bark)       |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_11.png)   | ![](figs/validation_time_varying_loudness_signal_11_specific_loudness.png)  |

| Test signal 12    | Test signal 12 (specific loudness @ 8.5 Bark)         |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_12.png)   | ![](figs/validation_time_varying_loudness_signal_12_specific_loudness.png)  |

| Test signal 13    | Test signal 13 (specific loudness @ 8.5 Bark)        |
| -------------- | -------------- |
| ![](figs/validation_time_varying_loudness_signal_13.png)   | ![](figs/validation_time_varying_loudness_signal_13_specific_loudness.png)  |

A summary presenting the differences (max. total loudness and loudness value exceed 5\% of the time) between calculated (SQAT) and reference values provided by ISO 532-1 is presented below for signals 06-13. With the reformulation of the code (released from v2.x), results computed using prior version (released in v1.x) are provided for comparison. Results of max. total loudness are close to zero with the new version of the code. Differences in exceeded value are no bigger than 0.5 sone. Despite providing reference values, the ISO standard do not stipulates any tolerance values for technical time-varying signals in terms of those single-value indicators.

| Summary of loudness differences for signals 06-13 (Max. total loudness)   | Summary of loudness differences for signals 06-13 (loudness value exceed 5\% of the time)         |
| -------------- | -------------- |
| ![](figs/validation_time_varying_synthetic_signals_loudness_difference_Nmax.png)   | ![](figs/validation_time_varying_synthetic_signals_loudness_difference_N5.png) |


# References
[1] International Organization for Standardization. (2017). Acoustics - Methods for calculating loudness - Part 1: Zwicker method (ISO Standard No. 532-1).

# Log
- This code was released in SQAT v1.0 (Gil Felix Greco, 14.05.2023)

- Summary of loudness differences included (Gil Felix Greco, 08.12.2024)

- Summary of loudness differences include now comparison between code versions (Gil Felix Greco, 21.08.2026)

