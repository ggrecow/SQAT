function test_Gen_weighting
% function test_Gen_weighting
%
% 1. Description:
%   In this test function, the A, B, C, and D curves are plotted as a
%   function of frequency.
%
% 2. Stand-alone example:
%    test_Gen_weighting
%
% 3. Additional info:
%    MATLAB / GNU Octave compatible: Yes / Yes
%
% Programmed by Alejandro Osses
% Date: 16/08/2023
% Date: 29/11/2025, GNU Octave compatible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 44100;
K = fs/2; % to get df=1 Hz
[bA,aA] = Gen_weighting_filters(fs,'A');
[bB,aB] = Gen_weighting_filters(fs,'B');
[bC,aC] = Gen_weighting_filters(fs,'C');
[bD,aD] = Gen_weighting_filters(fs,'D');

[hA, f] = freqz(bA,aA,K,fs);
hB = freqz(bB,aB,K,fs);
hC = freqz(bC,aC,K,fs);
hD = freqz(bD,aD,K,fs);

figure;
semilogx(f,20*log10(abs(hA)),'b'); hold on;
plot(    f,20*log10(abs(hB)),'g'); hold on;
semilogx(f,20*log10(abs(hC)),'r'); hold on;
semilogx(f,20*log10(abs(hD)),'k'); hold on;
xlim([20 20000]); grid on;
ylim([-55 15])

legend({'A','B','C','D'},'Location','SouthEast');
title(sprintf('Frequency-weighting curves\n(note that B and D are no longer standardised)'));
XT = [32 63 125 250 500 1000 2000 4000 8000 16000]; % Hz
set(gca,'XTick',XT);
set(gca,'XTickLabel',XT);

YT = -54:3:12; % dB, gain
set(gca,'YTick',YT);

xlabel('Frequency (Hz)');
ylabel('Relative magnitude (dB)');
