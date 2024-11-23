function [B, freqs, a0] = calculate_a0(fs,N,a0_type)
% function [B, freqs, a0] = calculate_a0(fs,N,a0_type)
%
% Compensation of the outer and middle ear transmission effects (Free-field), called 
% a0 compensation factor. The default method is a0_type = 'fastl2007' is to use the 
% compensation as defined in Fastl2007, Fig. 8.18, page 226 (doi: 10.1007/978-3-540-68888-4). 
%
% A simplified a0 compensation can be adopted if a0_type is set to
% 'fluctuationstrength_osses2016', where the ear canal resonance of Fastl's a0 curve is 
% removed. In other words, the a0 curve is roughly approximated as a low-pass filter.
% Although not explicitly stated by Osses et al. 2016 (doi: 10.1121/2.0000410), the simplified
% a0 transmission curve leads to very similar results during the validation
% of their fluctuation strength algorithm.
%
% Run the stand-alone example, below, to compare both types of curves.
%
% % Stand-alone example:
% N = 4096; % defines the frequency resolution: delta_f = fs/N;
% fs = 44100; % Hz, sampling frequency in Hz
% [B_more_accurate, freqs, a0_more_accurate] = calculate_a0(fs,N);
% a0_type = 'FluctuationStrength_Osses2016';
% [B_to_compare, freqs, a0_to_compare] = calculate_a0(fs,N,a0_type); % as used in the FluctuationStrength code
%
% figure; 
% plot(freqs,20*log10(abs(a0_more_accurate)),'b-'); hold on; grid on
% plot(freqs,20*log10(abs(a0_to_compare))   ,'r--'); 
% ylim([-93 13]);
% xlabel('Linear frequency (Hz)');
% ylabel('Gain factor due to the a0 transmission curve (dB)');
%
% Author: Alejandro Osses
% Date: 22 November 2024 (extension to contain the default Fastl2007 curve
%       and its simplified version as used by Osses et al. 2016).
% Date: 2014-2017 (Implementation)

if nargin == 0
    clc
    help calculate_a0;
    return;
end

if nargin < 3
    a0_type = 'fastl2007';
end

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = Get_Bark(N,qb,freqs);

switch lower(a0_type)
    case 'fastl2007'

% <il_calculate_a0> in the original <Fluctuation_strength_Osses2016> code
% a0tab = [ 
%     0       -999
%     0.5     -34.7
%     1       -23
%     1.5     -17
%     2       -12.8
%     2.5     -10.1
%     3       -8
%     3.5     -6.4
%     4       -5.1
%     4.5     -4.2
%     5       -3.5
%     5.5     -2.9
%     6       -2.4
%     6.5     -1.9
%     7       -1.5
%     7.5     -1.1 % 850 Hz
%     8       -0.8
%     8.5     0
%     10      0     % 1.2 kHz
%     12      1.15
%     13      2.31
%     14      3.85
%     15      5.62
%     16      6.92
%     16.5    7.38
%     17      6.92  % 3.5 kHz
%     18      4.23
%     18.5    2.31
%     19      0     % 5.4 kHz
%     20      -1.43
%     21		-2.59
%     21.5	-3.57
%     22		-5.19
%     22.5	-7.41
%     23		-11.3
%     23.5	-20
%     24		-40
%     25		-130
%     26		-999];

% Compensation of the transmission factor from Free-field, as defined 
% in Fastl2007, Fig. 8.18, page 226 (doi: 10.1007/978-3-540-68888-4). 
% This is the same curve used in the <Roughness_Daniel1997> function
a0tab =	[ 0	     0
    10	 0
    12	 1.15
    13	 2.31
    14	 3.85
    15	 5.62
    16	 6.92
    16.5	 7.38
    17	 6.92
    18	 4.23
    18.5	 2.31
    19	 0
    20	-1.43
    21	-2.59
    21.5	-3.57
    22	-5.19
    22.5	-7.41
    23	-11.3
    23.5	-20
    24	-40
    25	-130
    26	-999];

    case 'fluctuationstrength_osses2016'
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
    26		-999]; % minus infinity
end

a0            = zeros(1,N);
a0(qb)        = from_dB(interp1(a0tab(:,1),a0tab(:,2),Bark(qb)));
a0(isnan(a0)) = 0;

B = create_a0_FIR(freqs,a0(qb),N,fs);

if nargout == 0
    create_a0_FIR(freqs,a0(qb),N,fs);
end

a0 = a0(qb); % to have the same dimensions as freqs
