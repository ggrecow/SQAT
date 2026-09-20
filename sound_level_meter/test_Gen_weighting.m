function test_Gen_weighting
% function test_Gen_weighting
%
% 1. Description:
% 2. Stand-alone example:
% 3. Additional info:
%
% Programmed by Alejandro Osses
% Created on    : 16/08/2023
% Last update on: 16/08/2023 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
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

legend({'A','B','C','D'},'Location','SouthEast');
title(sprintf('Frequency-weighting curves\n(note that B and D are no longer standardised)'));
XT = [32 63 125 250 500 1000 2000 4000 8000 16000];
set(gca,'XTick',XT);
set(gca,'XTickLabels',XT);

xlabel('Frequency (Hz)');
ylabel('Relative magnitude (dB)');
