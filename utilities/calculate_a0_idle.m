function [B, freqs, a0] = calculate_a0_idle(fs,N)
% function [B, freqs, a0] = calculate_a0_idle(fs,N)
%
% The B coefficients contain an approximate transfer function of the outer 
% and middle ear. The 'idle' suffix indicates that the resonance of the
% ear canal was removed from the approximation. Actually, for mid and low
% frequencies, the B coefficient introduce no weighting at all (unit weigthing
% or 0 dB).
%
% % Stand-alone example:
% N = 4096; % defines the frequency resolution: delta_f = fs/N;
% fs = 44100; % Hz, sampling frequency in Hz
% calculate_a0_idle(fs,N);
% ylim([-103 3]);
%
% Author: Alejandro Osses
% Date: 13.11.2024. Created as a separate function (calculate_a0_idle.m)
% Date: 2014-2017. Original implementation as a subfunction or 'inline
%   function' (il_calculate_a0_idle, within the fluctuation strength code).

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = Get_Bark(N,qb,freqs);

a0tab = [
    0       0
    10      0
    19      0
    20      -1.43
    21		-2.59
    21.5	-3.57
    22		-5.19
    22.5	-7.41
    23		-11.3
    23.5	-20
    24		-40
    25		-130
    26		-999
];

a0            = zeros(1,N);
a0(qb)        = from_dB(interp1(a0tab(:,1),a0tab(:,2),Bark(qb)));
a0(isnan(a0)) = 0;

B = create_a0_FIR(freqs,a0(qb),N,fs);

if nargout == 0
    create_a0_FIR(freqs,a0(qb),N,fs);
end

a0 = a0(qb); % to have the same dimensions as freqs