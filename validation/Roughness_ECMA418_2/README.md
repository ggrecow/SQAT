# Roughness model according to ECMA-418-2:2025: verification of the implementation in SQAT

In SQAT, an implementation of the roughness model according to ECMA-418-2:2025 [1] is provided. The implementation was developed within the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license as detailed in the original code repository [(link)](https://github.com/acoustics-code-salford/refmap-psychoacoustics). The original implementation, called `acousticSHMRoughness.m` was ported into SQAT, where it is renamed according to the SQAT convention as `Roughness_ECMA418_2.m` (see code [here](../../psychoacoustic_metrics/Roughness_ECMA418_2/Roughness_ECMA418_2.m)). The following verification studies are available for this implementation:

- [Roughness dependence of AM tones on the modulation frequency](1_AM_modulation_freq)
- [Comparison with commercial software](2_Roughness_ECMA418_2_software_comparison)

# References
[1] Ecma International. (2025). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 4th Edition/June 2025). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf) (Last viewed 16 Nov 2025)

