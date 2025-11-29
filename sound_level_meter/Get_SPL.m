function output = Get_SPL(insig, fs, weight_freq, weight_time)
% function output = Get_SPL(insig, fs, weight_freq, weight_time)
%
% function : calculate SPL using the sound_level_meter in SQAT
%
% Get SPL using the  <Do_SLM> function of SQAT
%
%   Inputs
%   insig : input signal (Pa)
%   fs : sampling frequency (Hz)
%   weight_freq : frquency-weighting
%   weight_time : frquency-weighting
%
% type <Help Do_SLM> for more info about the input format
%
%   Output
%   output : Struct containing the following results
%   t : time vector
%   lvl_dB : sound pressure level over time (vector)
%   Lmax : maximum sound pressure level (scalar)
%   SEL : sound exposure level (scalar)
%
% Author: Gil Felix Greco
% Date: 15/11/2024
% Date: 29/11/2025 (Alejandro Osses: Stored as a separate function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dBFS = 94;

[output.lvl_dB, ~] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS);
output.t = (1:length(output.lvl_dB))/fs;

output.Leq = Get_Leq(output.lvl_dB,fs); % Make sure you enter only mono signals
T = length(output.lvl_dB)/fs;

output.Lmax = max(output.lvl_dB);
output.SEL  = output.Leq + 10*log10(T);

