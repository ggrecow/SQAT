# Tonality model according to ECMA-418-2:2025: verification of the implementation in SQAT

In SQAT, an implementation of the tonality model according to ECMA-418-2:2025 [1] is provided. The implementation was developed within the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license as detailed in the original code repository [(link)](https://github.com/acoustics-code-salford/refmap-psychoacoustics). The original implementation, called `acousticSHMTonality.m` was ported into SQAT, where it is renamed according to the SQAT convention as `Tonality_ECMA418_2.m` (see code [here](../../psychoacoustic_metrics/Tonality_ECMA418_2/Tonality_ECMA418_2.m)). The following verification studies are available for this implementation:

- [Comparison with commercial software](Tonality_ECMA418_2_software_comparison)

# References
[1] Ecma International. (2025). Psychoacoustic metrics for ITT equipment - Part 2 (methods for describing human perception based on the Sottek Hearing Model) (Standard No. 418-2, 4th Edition/June 2025). [https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf](https://ecma-international.org/wp-content/uploads/ECMA-418-2_4th_edition_june_2025.pdf) (Last viewed 16 Nov 2025)

