# Loudness according to ISO 532-1: verification of the implementation in SQAT
The Zwicker method for calculating the loudness of stationary and time-varying sounds is standardized in ISO 532-1 [1]. In SQAT, the implementation of this loudness model is named `Loudness_ISO532_1` (see code [here](../../psychoacoustic_metrics/Loudness_ISO532_1/Loudness_ISO532_1.m)). For the purpose of validating a specific implementation, the ISO standard provides a set of 25 test signals and the following requirements:

- Criteria for stationary loudness: for all stationary test signals given in the Annex B.2 and B.3, the specific loudness values shall differ by no more than $\pm5~\\%$ or $\pm0.1~\mathrm{sone/Bark}$ from the reference values, and the total loudness shall not deviate from the reference values by more than $\pm5~\\%$ or $\pm1~\mathrm{sone}$.

- Criteria for time-varying loudness: for all test signals given in the Annex B.4, the specific loudness vs. time shall not differ by more than $\pm5~\\%$ or $\pm0.1~\mathrm{sone/Bark}$ from the reference values within a temporal tolerance of $\pm2~\mathrm{ms}$. The tolerance can be extended to $\pm10~\\%$ or $\pm0.2~\mathrm{sone/Bark}$ within a temporal tolerance of $\pm2~\mathrm{ms}$, but only for maximum $1~\\%$ of the sampled specific loudness vs. time function using a time resolution of $\pm2~\mathrm{ms}$.

- Criteria for time-varying loudness: for all test signals given in the Annex B.4 and B.5, the deviation of the total loudness vs. time shall differ by no more than $\pm5~\\%$ or $\pm0.1~\mathrm{sone/Bark}$ from the reference values within a temporal tolerance of $\pm2~\mathrm{ms}$. The tolerance can be extended to $\pm10~\\%$ or $\pm0.2~\mathrm{sone/Bark}$ within a temporal tolerance of $\pm2~\mathrm{ms}$, but only for maximum $1~\\%$ of the sampled total loudness vs. time function using a time resolution of $\pm2~\mathrm{ms}$.

The dataset of test sounds provided by the ISO standard is freely available and can be downloaded from the following link: [http://standards.iso.org/iso/532/-1/ed-1/en](http://standards.iso.org/iso/532/-1/ed-1/en). The verification of the `Loudness_ISO532_1` implementation in SQAT is presented within this folder:    

- [Stationary loudness: synthetic signals provided in Annex B.2 and B.3](1_synthetic_signals_stationary_loudness)
- [Time-varying loudness: synthetic signals provided in Annex B.4](2_synthetic_signals_time_varying_loudness)
- [Time-varying loudness: technical signals provided in Annex B.5](3_technical_signals_time_varying_loudness)

# References
[1] International Organization for Standardization. (2017). Acoustics - Methods for calculating loudness - Part 1: Zwicker method (ISO Standard No. 532-1).
