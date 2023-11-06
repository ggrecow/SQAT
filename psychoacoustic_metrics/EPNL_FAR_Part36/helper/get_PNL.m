function [PN, PNL, PNLM, PNLM_idx] = get_PNL(input)
% function [PN, PNL, PNLM, PNLM_idx] = get_PNL(input)
%
% This function calculates the Perceived Noise Levels in two steps:
%
%       1) For each one-third octave band from 50 through 10,000 Hz, convert SPL(i, k) to perceived noisiness n(i, k),
%           by using the mathematical formulation of the noy table given in Section A36.4.7. (see Ref. [1] in the main function
%           or type <help EPNL_FAR_Part36> for more info)
%
%       2) The noy values are combined and then converted to instantaneous
%           perceived noise levels, PNL(k), using the formulation provided in Section A36.4.2.1 (see Ref. [1] in the main function
%           or type <help EPNL_FAR_Part36> for more info)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
%       input :
%       SPL[nTime,nFreq] matrix with nFreq=24 columns containing unweighted
%       SPL values for each third-octave band from 50 Hz to 10 kHz, and
%       nTime rows corresponding to time-steps
%
% OUTPUTS:
%       PN : [nTimes,1] vector
%       Perceived Noisiness versus time vector, in Noys
%
%       PNL : [nTimes,1] vector
%       Perceived Noise Level versus time vector, in PNdB
%
%       PNLM : scalar
%       max. Perceived Noise Level, in PNdB
%
%       PNLM_idx : scalar
%       index in the PNL(t) vector where max(PNL) occurs
%
% Function author: Gil Felix Greco, Braunschweig 27.10.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_times = size(input,1); % number of time steps
num_freqs = size(input,2); % number of freq bands

%% Read Table A36-3 from Ref. [1], which provides the constants for mathematically formulated NOY values

noy_tab = load('NOY_FORMULATING_TABLE.m');
% Band = noy_tab(:,1);
% f = noy_tab(:,2);
SPLa = noy_tab(:,3);
SPLb = noy_tab(:,4);
SPLc = noy_tab(:,5);
SPLd = noy_tab(:,6);
SPLe = noy_tab(:,7);
Mb = noy_tab(:,8);
Mc = noy_tab(:,9);
Md = noy_tab(:,10);
Me = noy_tab(:,11);

%% Convert SPL to Perceived Noisiness, nn

% DEFINITIONS
% i is index for the octave bands
% k is time-index vector

nn = zeros (num_times,num_freqs);

for i = 1:num_freqs
    for k = 1:num_times
        SPL = input(k,i);
        if (SPL >= SPLa(i))
            nn(k,i) = 10^(Mc(i) * (SPL - SPLc(i)));
        elseif (SPL >= SPLb(i) && SPL < SPLa(i))
            nn(k,i) = 10^(Mb(i) * (SPL - SPLb(i)));
        elseif (SPL >= SPLe(i) && SPL < SPLb(i))
            nn(k,i) = 0.3 * 10^(Me(i) * (SPL - SPLe(i)));
        elseif (SPL >= SPLd(i) && SPL < SPLe(i))
            nn(k,i) = 0.1 * 10^(Md(i) * (SPL - SPLd(i)));
        else
            nn(k,i) = 0;
        end
        clear SPL
    end
end

%% Combine the Perceived Noisiness, nn, to get PN and PNL

PNL = zeros(num_times,1);
PN = zeros(num_times,1);

for k = 1:num_times
    
    nmax = max(nn(k,:));
    PN(k) = 0.85 * nmax + 0.15 * (sum(nn(k,:))); % Perceived Noisiness, unit is Noys
    
    % Convert the total Perceived Noisiness, N(k) into Perceived Noise Level, PNL(k):
    PNL(k) = 40 + (10/log10(2)) * log10(PN(k)); % Perceived Noise Level, unit is PNdB
    
end

[PNLM,PNLM_idx] = max(PNL); % Maximum Perceived Noise Level (PNLM), unit is PNdB

end

