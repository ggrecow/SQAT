function B = create_a0_FIR(f,a0,N,fs)
% function B = create_a0_FIR(f,a0,N,fs)
%
% It creates a FIR filter according to the frequency f and curve defined in a0.
%
% This function is used by calculate_a0.m
%
% Author: Alejandro Osses
% Date: 2014-2017

f = [0 f fs/2];
a0 = [a0(1) a0 a0(end)];

B = fir2(N,f/(fs/2),a0);

if nargout == 0
    [H1,Fn]=freqz(B,1,N/2);
    
    figure;
    plot(fs/2*Fn/pi, 20*log10(abs(H1)));
    xlabel('Frequency [Hz]')
    legend([num2str(N) ' taps']);
    title('FIR filter resulting from the input (curve) parameter a0')
    xlim([0 fs/2])
end
