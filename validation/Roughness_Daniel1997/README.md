# Roughness model from Daniel & Weber: verification of the implementation in SQAT
In SQAT, an implementation of the roughness model from Daniel & Weber [1] is provided. The implementation is named `Roughness_Daniel1997` (see code [here](../../psychoacoustic_metrics/Roughness_Daniel1997/Roughness_Daniel1997.m)). In this folder, a thorough verification of the model implementation in SQAT is provided. The verification concerns the roughness dependence of amplitude-modulated (AM) and frequency-modulated (FM) sounds on different parameters, such as:

- [Roughness dependence of AM tones on the modulation frequency](1_AM_modulation_freq)
- [Roughness dependence of AM tones on the modulation depth](2_AM_modulation_depth)
- [Roughness dependence of FM tones on the modulation depth](3_FM_modulation_depth)
- [Roughness dependence of FM tones on the level](4_FM_level)
- [Roughness of unmodulated and AM white noise](5_white_noise)

The implementation was corrected in September 2026 against the revised code that Dik Hermes provided to the SQAT team (see issue 47 and the log of `Roughness_Daniel1997.m`), and the results in this folder were recomputed. Over a grid of 105 AM tones the two implementations agree within 0.041 asper.

# References
[1] Daniel, P., & Weber, R. (1997). Psychoacoustical Roughness: Implementation of an Optimized Model. [Acta Acustica united with Acustica](https://www.ingentaconnect.com/content/dav/aaua/1997/00000083/00000001/art00020), 83(1), 113-123.

