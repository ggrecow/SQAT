function [B, freqs, a0] = calculate_a0(fs,N)
% function [B, freqs, a0] = calculate_a0(fs,N)
%
% Transmission factor for free-field, according to Fig 8.18 (page 226) in 
% Fastl & Zwicker Book, Psychoacoustics: facts and models 3rd edition
%
% % Stand-alone example:
% N = 4096; % defines the frequency resolution: delta_f = fs/N;
% fs = 44100; % Hz, sampling frequency in Hz
% [B_more_accurate, freqs, a0_more_accurate] = calculate_a0(fs,N);
% [B_to_compare, freqs, a0_to_compare] = calculate_a0_idle(fs,N); % as used in the FluctuationStrength code
%
% figure; 
% plot(freqs,20*log10(abs(a0_more_accurate)),'b-'); hold on; grid on
% plot(freqs,20*log10(abs(a0_to_compare))   ,'r--'); 
% ylim([-103 3]);
% xlabel('Linear frequency (Hz)');
% ylabel('Gain factor due to the a0 transmission curve (dB)');
%
% Author: Alejandro Osses
% Date: 2014-2017

if nargin == 0
    clc
    help calculate_a0;
    return;
end

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = Get_Bark(N,qb,freqs);

a0tab = [ % lower slope from middle ear (fig_a0.c, see Figure_Psychoacoustics_tex)
    0       -999
    0.5     -34.7
    1       -23
    1.5     -17
    2       -12.8
    2.5     -10.1
    3       -8
    3.5     -6.4
    4       -5.1
    4.5     -4.2
    5       -3.5
    5.5     -2.9
    6       -2.4
    6.5     -1.9
    7       -1.5
    7.5     -1.1 % 850 Hz
    8       -0.8
    8.5     0
    10      0     % 1.2 kHz
    12      1.15
    13      2.31
    14      3.85
    15      5.62
    16      6.92
    16.5    7.38
    17      6.92  % 3.5 kHz
    18      4.23
    18.5    2.31
    19      0     % 5.4 kHz
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
