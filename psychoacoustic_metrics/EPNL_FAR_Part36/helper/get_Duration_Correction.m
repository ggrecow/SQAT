function [D, idx_t1, idx_t2] = get_Duration_Correction( PNLT, PNLTM, PNLTM_idx, dt, threshold )
% function [D, idx_t1, idx_t2] = get_Duration_Correction( PNLT, PNLTM, PNLTM_idx, dt, threshold )
%
% This function calculates the duration correction factor, D, by
% integration under the curve of the tone corrected perceived noise level
% versus time, during a time period of T=t(2)-t(1), where
%
%  - t(1) is the first point of time after which PNLT becomes greater than PNLTM–<threshold>
%  - t(2) is the point of time after which PNLT remains constantly less than PNLTM–<threshold>
%
% For aircraft certification process, the duration correction is calculated by integrating
% the PNLT based on a PNLTM-10(TPNdB) drop-down, here a different down point
% value can be considered by using the <threshold> input value.
%
% Moreover, for the aircraft certification process, dt should be 0.5 seconds and
% that is why the original formula has a -13 factor = 10 log10 (0.5/10) due
% to the normalization for the T0=10 seconds. If dt is not 0.5, a
% different factor applies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
%       PNLT : [nTimes,1] vector
%       Tone-corrected Perceived Noise Level versus time vector, in TPNdB
%
%       PNLTM : scalar
%       max. Tone-corrected Perceived Noise Level versus time vector, in TPNdB
%
%       PNLTM_idx : scalar
%       index in the PNLT(t) vector where max(PNLT) occurs
%
%       dt : integer
%       time-step, in seconds, in which the PNLT(t) is computed
%
%       threshold : integer
%       threshold value used to calculate the PNLT decay from PNLTM, in TPNdB
%
% OUTPUTS:
%
%       D : scalar
%       Duration correction factor, in TPNdB
%
%       idx_t1 : scalar
%       index in the PNLT(t) vector where t(1) occurs
%
%       idx_t2 : scalar
%       index in the PNLT(t) vector where t(2) occurs
%
% Function author: Gil Felix Greco, Braunschweig 27.10.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the PNLTM-threshold down points (t1 and t2)
Decay = PNLTM - threshold;
K=1; % find only first idx

% t(1) is the first point of time after which PNLT becomes greater than PNLTM minus threshold
idx_t1 = find( PNLT(1:PNLTM_idx)>Decay, K ); % idx_t1 is the first point where PNLT becomes >(PNLTM - threshold)

% t(2) is the point of time after which PNLT remains constantly less than PNLTM minus threshold
idx_t2 = find( PNLT(PNLTM_idx:end)<Decay, K); % idx_t2 is the first point where PNLT becomes <(PNLTM - threshold)
idx_t2 = idx_t2 + (PNLTM_idx-1); % correct for PNLTM_idx number because full vector is trimmed in the previous line

% if case for idx_t2 not found (PNLT never becomes lower than Decay)
if isempty(idx_t2)
    idx_t2 = size(PNLT, 1);  % take idx_t2 as last index in PNLT
	warning("The signal does not decay by more than the threshold within the available duration. An indicative EPNL value is calculated from the available duration, but should not be used for aircraft noise certification.");
end

% Calculate duration correction factor
D = 10*log10(sum(10.^(PNLT(idx_t1:idx_t2)/10))) - PNLTM + 10*log10(dt/10);

end
