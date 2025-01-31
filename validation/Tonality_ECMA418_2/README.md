# Tonality model according to ECMA-418-2:2024: verification of the implementation in SQAT

In SQAT, an implementation of the tonality model according to ECMA-418-2:2024 [1] is provided. The implementation was developed within the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license as detailed in the original code repository [(link)](https://github.com/acoustics-code-salford/refmap-psychoacoustics). The original implementation, called `acousticSHMTonality.m` was ported into SQAT, where it is renamed according to the SQAT convention as `Tonality_ECMA418_2.m` (see code [here](../../psychoacoustic_metrics/Tonality_ECMA418_2/Tonality_ECMA418_2.m)). The following verification studies are available for this implementation:

- [Comparison with commercial software](Tonality_ECMA418_2_software_comparison)

# References
[1] Ecma International. (2024). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 3rd Edition/December 2024). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_3rd_edition_december_2024.pdf) (Last viewed 22 Jan 2025)

