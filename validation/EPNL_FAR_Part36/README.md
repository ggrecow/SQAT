# Effective Perceived Noise Level (EPNL): verification of the implementation in SQAT

In SQAT, an implementation of the EPNL is provided according to Ref. [1]. The implementation is named `EPNL_FAR_Part36` (see code [here](../../psychoacoustic_metrics/EPNL_FAR_Part36/EPNL_FAR_Part36.m)). The EPNL is a metric mainly used in the context of environmental aircraft noise and is the most important acoustic indicator employed for aircraft noise certification. Despite its use for legal purposes, to the best knowledge of the Author of this text, a reference implementation or test signals tthat can be used to validate/verify a particular implementation of the EPNL calculation does not exist. Nevertheless, it is possible to verify if particular elements (or parts) of the EPNL calculation were correctly implemented. In this folder, we provide a (partial) verification of the EPNL code provided in SQAT, which concerns the following aspects:   

- [Implementation of the Perceived Noisiness](1_Perceived_Noisiness)
- [Implementation of the tone correction factor](2_Tone_Correction_Factor)

# References

[1] Federal Aviation Regulations, 14 CFR Parts 36 and 91, Docket No. FAA-2003-16526; Amendment No. 36-26, 91-288, (2005). [https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036](https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036) (Last viewed 30 Oct 2023)

# Log
README.md created on 30.10.2023 by Gil Felix Greco
