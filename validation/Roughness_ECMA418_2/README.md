# Roughness model according to ECMA-418-2:2024: verification of the implementation in SQAT

In SQAT, an implementation of the roughness model according to ECMA-418-2:2024 [1] is provided.  The implementation was developed within the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license as detailed in the original code repository [(link)](https://github.com/acoustics-code-salford/refmap-psychoacoustics). The original implementation, called `acousticSHMRoughness.m` was ported into SQAT, where it is named  `Roughness_ECMA418_2` (see code [here](../../psychoacoustic_metrics/Roughness_ECMA418_2/Roughness_ECMA418_2.m)). The following verification studies are available for this implementation:

- [Roughness dependence of AM tones on the modulation frequency](1_AM_modulation_freq)

# References
[1] Ecma International. (2024). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 3rd Edition/December 2024). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf) (Last viewed 22 Jan 2025)

