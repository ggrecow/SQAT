# Loudness model according to ECMA-418-2:2024: verification of the implementation in SQAT

In SQAT, an implementation of the loudness model according to ECMA-418-2:2024 [1] is provided. The implementation was developed within the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license as detailed in the original code repository [(link)](https://github.com/acoustics-code-salford/refmap-psychoacoustics). The original implementation, called `acousticSHMLoudness.m` was ported into SQAT, where it is renamed according to the SQAT convention as `Loudness_ECMA418_2.m` (see code [here](../../psychoacoustic_metrics/Loudness_ECMA418_2/Loudness_ECMA418_2.m)). The following verification studies are available for this implementation:

- [Equal-loudness-contours](1_equal_loudness_level_contours)
- [Comparison with commercial software](2_Loudness_ECMA418_2_software_comparison)

# References
[1] Ecma International. (2024). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 3rd Edition/December 2024). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf) (Last viewed 22 Jan 2025)

