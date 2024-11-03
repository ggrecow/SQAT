function [PNLT, PNLTM, PNLTM_idx, OUT] = get_PNLT( input, freq_bands, PNL )
% function [PNLT, PNLTM, PNLTM_idx, OUT] = get_PNLT( input, freq_bands, PNL )
%
% This function calculates the Tone-corrected Perceived Noise Level, PNLT(k), by:
%
%       1) Calculation the 10-steps tone-correction factor described in
%          Refs. [1,2,3] (type <help EPNL_FAR_Part36> for more info about the references)
%
%       2) The largest of the tone correction factors, Cmax(k), is added to
%       corresponding PNL(k) values in order to determine the Tone-corrected perceived noise levels, PNLT(k) 
%
%       3) Bandsharing adjustment to PNLTM added in 20.06.2024
%
% A detailed definition and discussion about the tone-correction factor is
% provided in:
%
% Correction Procedures for Aircraft Noise Data. Volume IV. Tone
% Perception. May,D. N.; Watson,E. E.; WYLE LABS EL SEGUNDO CALIFORNIA;
% 1980-02-01 - Available at: https://apps.dtic.mil/sti/citations/ADA083075
% (Last viewed: 26 Oct 2023)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
%       input : matrix
%       SPL[nTime,nFreq] matrix with nFreq=24 columns containing unweighted
%       SPL values for each third-octave band from 50 Hz to 10 kHz, and
%       nTime rows corresponding to time-steps
%
%       freq_bands : vector
%       [1, nFreq=24] vector corresponding to the 1/3 octave bands center frequency        
%
%       PNL : [nTimes,1] vector
%       Perceived Noise Level versus time vector, in PNdB
%
% OUTPUTS:
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
%       OUT : struct
%       Struct containing the S, diff, delS, SPLP, SP, SB, SPLPP, F, and C variables (only used for verifying the tone-correction factor implementation)
%
% Function author: Gil Felix Greco, Braunschweig 27.10.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINITIONS
% i is index for the octave bands
% k is time-index vector

num_times = size(input,1); % number of time steps
num_freqs = size(input,2); % number of freq bands

%% Step 1: start with the SPL in the 80 Hz TOB (band number 3 in this case),
% and calculate the changes in SPL (or "SLOPES) in the remainder TOB

index_80 = find(freq_bands>=80,1);

S = zeros(num_freqs,num_times);
for k = 1:num_times
    for i = index_80+1:num_freqs
        S(i,k) = input(k,i) - input(k,i-1);
    end
end

%% Step 2: Encircle the value of the SLOPE, S(i,k),
% where the absolute value of the change in slope is greater than 5 dB
% STEP 3 is also included here

delS = zeros(size(S));
SPLs = delS;
diff = zeros(1,num_freqs);

for k = 1:num_times
    for i = index_80:num_freqs
        diff(i) = abs(S(i,k) - S(i-1,k));
        if (diff(i) > 5) % <--Step 2
            delS(i,k) = S(i,k);
            if (S(i,k) > 0 && (S(i,k) > S(i-1,k))) % Step 3.1: If the encircled value of S(i,k) is positive and algebraically greater than S(i-1,k), then encircle SPL(i,k)
                SPLs(i,k) = input(k,i);
            elseif (S(i,k) <= 0 && (S(i-1,k) > 0)) % Step 3.2: if the encircled value of S(i,k) is zero or negative and S(i-1,k) is positive, then encircle SPL(i-1,k)
                SPLs(i-1,k) = input(k,i-1);
            end
        end % Step 3.3: for all other cases, no SPL value is encircled
    end
end

%% Step 4: Omit all SPL(i,k) encircled in Step 3 and
% compute the new SPL levels, SPLP(i,k):

SPLP = zeros(num_freqs,num_times);

for k = 1:num_times
    for i = 1:num_freqs
        if (SPLs(i,k) == 0) % Step 4.1: for non-encircled SPL, set SPLP equal to the original SPL, SPLP(i,k) = SPL(i,k)
            SPLP(i,k) = input(k,i);
        elseif (SPLs(i,k) > 0 && i < num_freqs) % Step 4.2: for encircled SPL in bands 1 (50 Hz) till 23 (8 kHz) inclusive, set SPLP equal to the arithmetic average of the preceding and following SPL
            SPLP(i,k) = 0.5 * (input(k,i-1) + input(k,i+1));
        elseif (SPLs(i,k) > 0 && i == num_freqs) % Step 4.3: if the SPL in the highest freq band (10 kHz) is encircled, set the SPLP in that band equal to SPLP (k,24) = SPL (k,23) + S(k,23).
            SPLP(i,k) = (input(k,i-1) + S(i-1,k));
        elseif (SPLs(i,k) <= 0 && i == num_freqs)
            SPLP(i,k) = input(k,i);
        end
    end
end

%% STEP 5: Recompute the new SLOPE, SP(i,k),
% including one for an imaginary 25th freq band (i.e. 12.5 kHz)

SP = zeros(num_freqs+1,num_times);

for k = 1:num_times
    for i = index_80+1:num_freqs-1
        SP(i,k) = SPLP(i,k) - SPLP(i-1,k);
    end
    SP(index_80,k) = SP(index_80+1,k);
    SP(num_freqs,k) = SPLP(num_freqs,k) - SPLP(num_freqs-1,k);
    SP(num_freqs+1,k) = SP(num_freqs,k);
end

%% STEP 6: for i, from 3 (80 Hz) till 23 (8 kHz),
% compute the arithmetic average of the three adjacent slopes (i.e. next 2)

SB = zeros(num_freqs+1,num_times);
for k = 1:num_times
    for i = index_80:num_freqs-1
        SB(i,k) = (1/3) * (SP(i,k) + SP(i+1,k) + SP(i+2,k));
    end
end

%%  STEP 7: compute the final adjusted TOB SPL, SPLPP(i,k),
% by beginning with the band number 3 (80 Hz) and proceeding to band number 24 (10 kHz)

SPLPP = zeros(num_freqs,num_times);

for k = 1:num_times
    SPLPP(index_80,k) = input(k,index_80);
    for i = index_80+1:num_freqs-1
        SPLPP(i,k) = (SPLPP(i-1,k) + SB(i-1,k));
    end
    SPLPP(num_freqs,k) = SPLPP(num_freqs-1,k) + SB(num_freqs-1,k);
end

%% Step 8: calculate the differences, F(i,k),
% between the original SPL and the final background SPL, SPLPP(i,k)

F = zeros(num_freqs,num_times);

for k = 1:num_times
    for i = index_80:num_freqs
        F(i,k) = (input(k,i) - SPLPP(i,k));
        if F(i,k)<=1.5 % note only values equal to or greater than 1.5.
            F(i,k)=0;
        end
    end
end

%% Steps 9 & 10: for each TOB from 80 Hz till 10 kHz (i.e. band 3 through 24),
% determine the tone correction factor C(i,k) from the SPL differences F(i,k), and Table A36-2
% (Step 10 is included here also)

C = zeros(num_freqs,num_times);
Cmax = zeros(num_times,1);

for k = 1:num_times
    for i = index_80:num_freqs
        if (freq_bands(i) >= 50 && freq_bands(i) < 500 && F(i,k) >= 1.5 && F(i,k) < 3)
            C(i,k) = (F(i,k)/3 - 0.5);
        elseif (freq_bands(i) >= 50 && freq_bands(i) < 500 && F(i,k) >= 3 && F(i,k) < 20)
            C(i,k) = (F(i,k)/6);
        elseif (freq_bands(i) >= 50 && freq_bands(i) < 500 && F(i,k) >= 20)
            C(i,k) = 10/3;
        elseif (freq_bands(i) >= 500 && freq_bands(i) <= 5000 && F(i,k) >= 1.5 && F(i,k) < 3)
            C(i,k) = ((2 * F(i,k)/3) - 1);
        elseif (freq_bands(i) >= 500 && freq_bands(i) <= 5000 && F(i,k) >= 3 && F(i,k) < 20)
            C(i,k) = (F(i,k)/3);
        elseif (freq_bands(i) >= 500 && freq_bands(i) <= 5000 && F(i,k) >= 20)
            C(i,k) = (20/3);
        elseif (freq_bands(i) > 5000 && freq_bands(i) <= 10000 && F(i,k) >= 1.5 && F(i,k)<3)
            C(i,k) = (F(i,k)/3 - 0.5);
        elseif (freq_bands(i) > 5000 && freq_bands(i) <= 10000 && F(i,k) >= 3 && F(i,k) <20)
            C(i,k) = (F(i,k)/6);
        elseif (freq_bands(i) > 5000 && freq_bands(i) <= 10000 && F(i,k) >= 20)
            C(i,k) = (10/3);
        end
    end
    Cmax(k) = max(C(:,k)); % <--- STEP 10: designate the largest of the tone correction factors determined in STEP 9 as Cmax(k).
end

%% PNLT calculation
% The Tone-corrected perceived noise levels, PNLT(k), must be determined by
% adding Cmax(k) values to corresponding PNL(k) values

PNLT = zeros(num_times,1);

for k = 1:num_times
    PNLT(k) = PNL(k) + Cmax(k); % TONE-CORRECTED PERCEIVED NOISE LEVEL, unit is TPNdB
end

[PNLTM,PNLTM_idx] = max(PNLT); % MAXIMUM TONE-CORRECTED PERCEIVED NOISE LEVEL (PNLTM)

%% Bandsharing adjustment to PNLTM

if size(Cmax,1)~=1 % workaround to run the <run_validation_tone_correction.m> code, where only one time-step is considered

    Cavg = sum( [Cmax(PNLTM_idx - 2), Cmax(PNLTM_idx - 1), Cmax(PNLTM_idx),...
                           Cmax(PNLTM_idx + 1), Cmax(PNLTM_idx + 2)] )/5;

    if Cavg > Cmax(PNLTM_idx)
        DeltaB = Cavg*Cmax(PNLTM_idx);
    else
        DeltaB = 0;
    end

    % apply adjustment
    PNLTM = PNLTM + DeltaB;

else
end

%% Output variables for verification of the tone correction implementation

OUT.S = S;
OUT.delS = delS;
OUT.SPLP = SPLP;
OUT.SP = SP;
OUT.SB = SB;
OUT.SPLPP = SPLPP;
OUT.F = F;
OUT.C =C; 
OUT.diff = diff;

end
