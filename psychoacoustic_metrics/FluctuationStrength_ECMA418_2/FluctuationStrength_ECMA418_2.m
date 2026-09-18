function OUT = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip, show)
% OUT = FluctuationStrength_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%
% Returns fluctuation strength values according to ECMA-418-2:2025
% (using the Sottek Hearing Model) for an input calibrated
% single mono or single stereo audio (sound pressure) time-series
% signal, insig. For stereo signals, Fluctuation Strength is calculated
% for each channel [left ear, right ear], and also for the combination
% of both, denominated as "combined binaural" (see ECMA-418-2:2025,
% Section 9.1.15).
%
% According to ECMA-418-2:2025 (Section 9.1.14), the 90th percentile
% (a.k.a. value exceeded 10% of the time) of the time-dependent
% fluctuation strength must be used as the representative single value
% (overall fluctuation strength). This value is provided here by the
% <fluctStrength90Pc> output variable.
%
% NOTE: according to ECMA-418-2:2025 (Section 9.1.12), the fluctuation
% strength values corresponding to the first 36 output samples
% (l50 = 0 to 35, approximately the first 683 ms of the input signal)
% must be discarded due to the transient responses of the digital
% filters. Therefore, <time_skip> must be greater or equal to 700 ms
% to compute any time-aggregated quantity.
%
% Reference signal: 60 dBSPL 1 kHz tone 100% modulated at 4 Hz yields
% 1 vacilHMS.
%
% Inputs
% ------
% insig : column vector [Nx1] mono or [Nx2] binaural
%   the input signal as single mono or stereo audio (sound
%   pressure) signals
%
% fs : integer
%   the sample rate (frequency) of the input signal(s)
%
% fieldtype : keyword string (default: 'free-frontal')
%   determines whether the 'free-frontal' or 'diffuse' field stages
%   are applied in the outer-middle ear filter
%
% time_skip : double (default: 700 ms - see Section 9.1.12 ECMA-418-2:2025)
%   skip start of the signal in <time_skip> seconds so that
%   the transient response of the digital filters is avoided.
%   Best-practice: <time_skip> must be equal or higher than default value
%
% show : Boolean true/false (default: false)
%   flag indicating whether to generate a figure from the output
%
% Returns
% -------
%
% OUT : structure
%   contains the following fields:
%
% specFluctStrength : matrix
%   time-dependent specific fluctuation strength for each critical band
%   arranged as [time, bands(, channels)]
%
% specFluctStrengthAvg : matrix
%   time-averaged specific fluctuation strength for each critical band
%   arranged as [bands(, channels)]
%   OBS: takes <time_skip> into consideration
%
% fluctStrengthTDep : vector or matrix
%   time-dependent overall fluctuation strength arranged as
%   [time(, channels)]
%
% fluctStrength90Pc : number or vector
%   90th percentile (a.k.a. value exceeded 10% of the time) of the
%   time-dependent fluctuation strength. According to ECMA-418-2:2025
%   (Section 9.1.14), this quantity must be used as the representative
%   single value (overall fluctuation strength).
%   OBS: takes <time_skip> into consideration
%
% bandCentreFreqs : vector
%   centre frequencies corresponding with each critical band rate
%
% timeOut : vector
%   time (seconds) corresponding with time-dependent outputs
%
% timeInsig : vector
%   time (seconds) of insig
%
% soundField : string
%   identifies the soundfield type applied (the input argument fieldtype)
%
% Several statistics based on fluctStrengthTDep
%         ** FSmean : mean value of instantaneous fluctuation strength
%               (vacilHMS)
%         ** FSstd : standard deviation of instantaneous fluctuation
%               strength (vacilHMS)
%         ** FSmax : maximum of instantaneous fluctuation strength
%               (vacilHMS)
%         ** FSmin : minimum of instantaneous fluctuation strength
%               (vacilHMS)
%         ** FSx : fluctuation strength value exceeded during x percent
%               of the time (vacilHMS)
%             in case of binaural input, FSx(1,3), being 1st, 2nd and 3rd
%             column corresponding to [left ear, right ear, comb. binaural]
%             OBS: all quantities here take <time_skip> into consideration
%
% In case of stereo inputs, the following additional fields are provided
% separately for the "comb. binaural" case (i.e. combination of left and
% right ears)
%
% specFluctStrengthBin : matrix
%   time-dependent specific fluctuation strength for each critical band
%   arranged as [time, bands]
%
% specFluctStrengthAvgBin : matrix
%   time-averaged specific fluctuation strength for each critical band
%   arranged as [bands]
%   OBS: takes <time_skip> into consideration
%
% fluctStrengthTDepBin : vector
%   time-dependent overall fluctuation strength arranged as [time]
%
% fluctStrength90PcBin : number
%   90th percentile of the time-dependent fluctuation strength
%   (combined binaural). OBS: takes <time_skip> into consideration.
%
% If show==true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall fluctuation strength, with the latter also indicating the
% time-aggregated value. In case of stereo signals, a set of plots is
% returned for each input channel, with another set for the combined
% binaural fluctuation strength. For the latter, the indicated
% time-averaged A-weighted sound level corresponds with the channel with
% the highest sound level.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in
% Pascals (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Sergio Aguirre &
%          Gil Felix Greco
%
% Date created: 17.09.2026
% Date last modified: 17.09.2026
% MATLAB version: 2026a
%
% Copyright statement: This file is part of the SQAT toolbox and is
% subject to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt>
% in the SQAT repository root. Some files in SQAT carry a different
% license, always stated in their own header; where this file depends on
% them, the combined work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is",
% WITHOUT WARRANTY OF ANY KIND, express or implied, including but not
% limited to the warranties of MERCHANTABILITY and FITNESS FOR A
% PARTICULAR PURPOSE.
%
% This code calls sub-component file 'cmap_inferno.txt'. The contents of
% the file includes a copy of data obtained from the repository
% https://github.com/BIDS/colormap, and is CC0 1.0 licensed for modified
% use, see https://creativecommons.org/publicdomain/zero/1.0 for
% information. See also the <licenses/mpl-colormaps LICENCE.txt> file in
% the SQAT repository root.
%
% Checked by:
% Date last checked:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Arguments validation
    arguments (Input) % Matlab R2018b or newer
        insig (:, :) double {mustBeReal}
        fs (1, 1) double {mustBePositive, mustBeInteger}
        fieldtype (1, :) string {mustBeMember(fieldtype,...
                                                       {'free-frontal',...
                                                        'diffuse'})} = 'free-frontal'
        time_skip (1, 1) double {mustBeReal} = 700e-3
        show {mustBeNumericOrLogical} = false
    end

%% Input checks

% define time threshold value from which all values before must be
% dropped (Section 9.1.12: first 36 values at 50 Hz)
t_threshold = 700e-3;

% check insig dimension (only [Nx1] or [Nx2] are valid)
if size(insig,1) > 2 & size(insig,2) > 2 % insig has more than 2 channels
    error('Error: Input signal has more than 2 channels. ')
elseif size(insig, 2) > 2  % insig is [1xN] or [2xN]
    insig = insig';
    fprintf('\nWarning: Input signal is not [Nx1] or [Nx2] and was transposed.\n');
end

% Check the length of the input data (block size is sb = 65536 at 48 kHz,
% i.e. 1.3653 s; one full block beyond the start padding is required)
if size(insig, 1) < ceil(1.5*fs)
    error("Error: Input signal is too short along the specified axis to calculate fluctuation strength (must be at least 1.5 s)")
end

% Check the channel number of the input data
if size(insig, 2) > 2
    error("Error: Input signal comprises more than two channels")
else
    chansIn = size(insig, 2);
    if chansIn == 2
        chans = ["Stereo left"; "Stereo right"];
        binaural = true;
    else
        chans = "Mono";
        binaural = false;
    end
end

if time_skip < t_threshold
    warning("Time_skip must be at least 700 ms to avoid transient responses of the digital filters (see ECMA-418-2:2025, Section 9.1.12). Setting time_skip to 700 ms!!!")
    time_skip = t_threshold;
end

%% Define constants

signalT = size(insig, 1)/fs;  % duration of input signal
sampleRate48k = 48e3;  % processing sample rate (Section 5.1.1) [r_s]

% Critical band structure (Section 5.1.4.1 ECMA-418-2:2025)
deltaFreq0 = 81.9289;  % [deltaf(f=0)]
c = 0.1618;  % centre-frequency constant
dz = 0.5;  % critical band overlap [deltaz]
halfBark = 0.5:dz:26.5;  % half-overlapping critical band rate scale [z]
nBands = length(halfBark);  % number of bands [CBF = 53]
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Equation 9 [F(z)]

% Block and hop sizes (Section 9.1.1 ECMA-418-2:2025)
overlap = 0.75;  % block overlap proportion
blockSize = 65536;  % block size [s_b]
hopSize = (1 - overlap)*blockSize;  % hop size [s_h] = 16384

% Downsampled rates (Section 9.1.2 ECMA-418-2:2025)
downSample = 32;  % downsampling factor
sampleRate1500 = sampleRate48k/downSample;  % [r_s_tilde] = 1500 Hz
blockSize1500 = blockSize/downSample;  % [s_b_tilde] = 2048
hopSize1500 = hopSize/downSample;  % [s_h_tilde] = 512
resDFT1500 = sampleRate1500/blockSize1500;  % DFT resolution [deltaf]

% Envelope analysis window defaults (Section 9.1.3 ECMA-418-2:2025)
winMargin = 64;  % default zeros at beginning/end [nzb = nze = sb~/32]
medianWinLen = 33;  % moving median length (Section 9.1.3.2)
roundDigits = 8;  % decimal rounding of envelopes (Section 9.1.3.2)
roundScale = 10^roundDigits;  % rounding scale
pEmaxMin = 5e-6;  % Pa, block gate (Section 9.1.3.2)
qpThrFactor = 0.01;  % quieter period threshold factor (Section 9.1.3.3)
minQuietLen = 320;  % minimum quieter period length (Section 9.1.3.4)
n2Min = 511;  % validity bound (Section 9.1.3.6)
regTrim = 10;  % regression trim (Section 9.1.3.6)
relStdMin = 1e-3;  % relative std bound (Section 9.1.3.6)

% HSA parameters (Sections 9.1.4 and 9.1.5 ECMA-418-2:2025)
kMax = 48;  % spectral lines used for positive modulation rates
phiEmin = 0.15;  % absolute maximum criterion, Equation 143
phi0Factor = 0.001;  % relative maximum criterion, Equation 143
fiGrid = 0.25*2.^(((1:16) - 2)/3);  % error minima candidate rates [f_i], Section 9.1.5
dupTol = 1.25;  % duplicate tolerance in deltaf units, Equation 145
dupTolHz = dupTol*resDFT1500;  % duplicate tolerance (Hz)
preSelFactor = 0.05;  % preselection factor, Equation 146

% Newton fine tuning parameters (Section 9.1.7 ECMA-418-2:2025)
newtonDx = 1e-5;  % differential quotient step
newtonDamp = 0.25;  % damping factor, Equation 152
newtonStepCap = 2e-4;  % step cap, Equation 152
newtonTol = 1e-7;  % stop tolerance, Equation 152
newtonMaxIter = 40;  % iteration cap, Equation 152
minModRate = 0.125;  % Hz, block discard rate, Section 9.1.7

% Harmonic analysis parameters (Section 9.1.8 ECMA-418-2:2025)
harmOrders = 1:3;  % assumed orders of f_c,1,opt
harmMaxOrder = 5;  % highest order considered
harmTol = 0.04;  % harmonic tolerance, Equation 154

% Sum weighting parameters (Section 9.1.9 ECMA-418-2:2025)
wbwA = 0.79577;  % Equation 158
wbwExp = 0.43461;  % Equation 158

% Aural nonlinearity constants (Section 5.1.8 Equation 23 and Table 2,
% Section 5.1.9 Table 3 ECMA-418-2:2025; same values as in
% <shmBasisLoudness.m>)
cal_N = 0.0211668;  % calibration factor [c_N]
cal_Nx = 1.00132;  % calibration multiplier
aNL = 1.5;  % nonlinearity exponent [alpha]
p_threshold = 2e-5*10.^((15:10:85)/20).';  % Table 2
vNL = [1, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434];  % Table 2
LTQz = [0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,...
        0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,...
        0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,...
        0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,...
        0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,...
        0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,...
        0.0296, 0.0339, 0.0398, 0.0485, 0.0622];  % Table 3 [LTQ(z)]

% Scaling threshold (Section 9.1.10 ECMA-418-2:2025)
aMinThreshold = 5.2519;

% Output stage constants (Section 9.1.11 ECMA-418-2:2025)
sampleRate50 = 50;  % output sampling rate [r_s50]
cal_F = 0.003840572;  % calibration factor [c_F], Equation 163
eCoef = 0.37106;  % exponent coefficient, Equation 164
eSlope = 1.6407;  % exponent slope, Equation 164
eCenter = 2.5804;  % exponent centre, Equation 164
eBase = 0.58449;  % exponent offset, Equation 164
medianWinB = 71;  % moving median length for B, Equation 165
tauFS = 0.75;  % lowpass time constant [tau], Equation 168

% Standardised epsilon (footnote 14) and machine epsilon
epsilon = 1e-12;
eps0 = eps;

%% Signal processing

% Input pre-processing
% --------------------
if fs ~= sampleRate48k  % Resample signal
    [p_re, ~] = shmResample(insig, fs);
else  % don't resample
    p_re = insig;
end

% get time vector of input signal (at the processing sample rate)
timeInsig = (0 : length(p_re(:,1))-1) ./ sampleRate48k;

% Input signal samples
n_samples = size(p_re, 1);

% Section 5.1.2 ECMA-418-2:2025 Fade in weighting and zero-padding
% (only the start is zero-padded, with n_zeros,start = s_b = 65536 as
% required for fluctuation strength in Section 5.1.2.2)
pn = shmPreProc(p_re, max(blockSize), max(hopSize), true, false);

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2025 Outer and middle/inner ear signal filtering
pn_om = shmOutMidEarFilter(pn, fieldtype);

% Loop through channels in file
% -----------------------------
for chan = size(pn_om, 2):-1:1

    % Apply auditory filter bank
    % --------------------------
    %
    % Section 5.1.4.2 ECMA-418-2:2025
    pn_omz = shmAuditoryFiltBank(pn_om(:, chan), false);

    % Note: at this stage, memory imposes a loop over the critical bands
    % until the downsampling of the envelopes is applied
    for zBand = nBands:-1:1
        % Segmentation into blocks
        % ------------------------
        %
        % Section 5.1.5 ECMA-418-2:2025 (block size and hop size for
        % fluctuation strength as defined in Section 9.1.1)
        i_start = 1;
        [pn_lz, iBlocksOut] = shmSignalSegment(pn_omz(:, zBand), 1,...
                                               blockSize,...
                                               overlap, i_start, true);

        % Envelope calculation and downsampling
        % -------------------------------------
        %
        % Section 9.1.2 ECMA-418-2:2025: magnitude of Hilbert transform
        % analytic signal with downsampling, Equation 119
        % [p_E,l,z(n_tilde)]
        envelopes(:, :, zBand) = downsample(abs(hilbert(pn_lz)),...
                                            downSample, 0);

    end  % end of for loop over critical bands

    % Analysis of the envelopes of each block and band
    % -------------------------------------------------
    % Sections 9.1.3 to 9.1.10 ECMA-418-2:2025

    nBlocks = size(envelopes, 2);
    specLoudHSA = zeros(nBlocks, nBands);  % [N'_HSA(l,z)]
    powerHSA = zeros(nBlocks, nBands);  % p_0^2 + 2*sum(A_i) + eps
    modAmpHat = zeros(nBlocks, nBands);  % [A_hat(l,z)]

    for lBlock = 1:nBlocks
        for zBand = 1:nBands

            pE = envelopes(:, lBlock, zBand);  % [p_E,l,z(n_tilde)]

            % Section 9.1.3.2 ECMA-418-2:2025: smoothing by moving median,
            % rounding to 8 digits and weighting with the default window
            pEs = round(movmedian(pE, medianWinLen)*roundScale)/roundScale;
            nzb = winMargin;  % [n_zb]
            nze = winMargin;  % [n_ze]
            wInit = zeros(blockSize1500, 1);
            wInit(nzb + 1:blockSize1500 - nze) = 1;
            pEw = pEs.*wInit;

            % block gate on the maximum of the weighted envelope
            pEmax = max(pEw);
            if pEmax <= pEmaxMin
                continue  % the entire block is a quieter period
            end

            % Section 9.1.3.3 ECMA-418-2:2025: threshold (rounded to 8
            % digits) and update of the border parameters with the first
            % and the last sample above the threshold
            pEthr = round(qpThrFactor*pEmax*roundScale)/roundScale;
            aboveThr = pEw >= pEthr;
            firstIdx = find(aboveThr, 1, 'first');
            lastIdx = find(aboveThr, 1, 'last');
            nzb = firstIdx - 1;
            nze = blockSize1500 - lastIdx;

            % Section 9.1.3.4 ECMA-418-2:2025: longest quieter period inside
            % the updated interval
            belowThr = ~aboveThr(firstIdx:lastIdx);
            runEdges = diff([false; belowThr; false]);
            runStart = find(runEdges == 1);
            runEnd = find(runEdges == -1) - 1;
            runLength = runEnd - runStart + 1;
            isLong = runLength > minQuietLen;
            if any(isLong)
                runLengthLong = runLength(isLong);
                runStartLong = runStart(isLong);
                runEndLong = runEnd(isLong);
                [~, iLong] = max(runLengthLong);
                nqpmb = firstIdx - 1 + runStartLong(iLong) - 1;  % 0-based index
                nqpme = firstIdx - 1 + runEndLong(iLong) - 1;  % 0-based index
                % Section 9.1.3.5 ECMA-418-2:2025: keep the longer part with
                % ones, with the margin of s_b_tilde/32 zeros
                if (nqpmb - (nzb + winMargin)) > ((blockSize1500 - 1 - nze - winMargin) - nqpme)
                    nze = blockSize1500 - 1 - nqpmb + winMargin;
                else
                    nzb = nqpme + winMargin;
                end
            end

            % Section 9.1.3.6 ECMA-418-2:2025: validity checks of the final
            % analysis interval
            n1 = nzb;
            n2 = blockSize1500 - 1 - nze;
            if (n2 - n1 + 1) < minQuietLen || n2 < n2Min
                continue  % the entire block is a quieter period
            end
            % relative standard deviation of a linear regression of the
            % smoothed envelope over [n1 + 10, n2 - 10]
            i1 = n1 + regTrim;
            i2 = n2 - regTrim;
            segY = pEs(i1 + 1:i2 + 1);
            segX = (i1:i2).';
            regCoef = polyfit(segX, segY, 1);
            regResid = segY - polyval(regCoef, segX);
            relStd = std(regResid)/(mean(segY) + eps0);
            if relStd < relStdMin
                continue  % the entire block is a quieter period
            end
            wE = zeros(blockSize1500, 1);  % Equation 120 [w_E,l,z(n_tilde)]
            wE(n1 + 1:n2 + 1) = 1;

            % Section 9.1.4 ECMA-418-2:2025: spectrum and power spectrum,
            % Equations 121 and 122, for k = 0..48
            PE = fft(pE.*wE);
            PE = PE(1:kMax + 1);
            PhiE = abs(PE).^2;

            % Section 9.1.5 ECMA-418-2:2025: local maxima fulfilling
            % Equation 143; the end bins k = 0 and k = 48 cannot be maxima
            [PhiPks, locPks] = findpeaks(PhiE);
            kPks = locPks - 1;  % 0-based bin index
            kPks = kPks(PhiPks >= max(phi0Factor*PhiE(1), phiEmin));

            % Equation 144: modulation rates of the maxima. The printed
            % "- 1" corresponds to one-based bin indices; with the zero-based
            % k of the DFT definition (footnote 38) the rate is the weighted
            % centroid of the three bins times deltaf
            fp = zeros(numel(kPks), 1);
            for iPk = 1:numel(kPks)
                kCentroid = kPks(iPk) + (-1:1);
                PhiCentroid = PhiE(kCentroid + 1);
                PhiCentroid = PhiCentroid(:).';
                fp(iPk) = sum(kCentroid.*PhiCentroid)/sum(PhiCentroid)*resDFT1500;
            end

            % local minima of the error function E((0, f_i)) on the grid f_i
            errGrid = zeros(numel(fiGrid), 1);
            for iGrid = 1:numel(fiGrid)
                errGrid(iGrid) = shmHighResSpecAnalysis(PE, nzb, nze, fiGrid(iGrid));
            end
            isLocMin = false(size(errGrid));
            isLocMin(2:end - 1) = errGrid(2:end - 1) < errGrid(1:end - 2)...
                                  & errGrid(2:end - 1) < errGrid(3:end);

            if ~any(isLocMin)
                if isempty(kPks)
                    continue  % no local minimum and no local maximum: no modulation
                end
                candRates = fp(:).';  % all local maxima are the candidates
            else
                idxLocMin = find(isLocMin);
                [~, iLocMin] = min(errGrid(idxLocMin));
                fmin = fiGrid(idxLocMin(iLocMin));  % [f_min(l,z)]
                % Equation 145: duplicates of spectral line pairs
                isDup = abs(fmin - fp(:)) < dupTolHz;
                if any(isDup)
                    candRatesI = [fmin, fp(~isDup).'];  % case I
                    candRatesII = fp(:).';  % case II
                    errI = shmHighResSpecAnalysis(PE, nzb, nze, candRatesI);
                    errII = shmHighResSpecAnalysis(PE, nzb, nze, candRatesII);
                    if errI <= errII
                        candRates = candRatesI;
                    else
                        candRates = candRatesII;
                    end
                else
                    candRates = [fmin, fp(:).'];
                end
            end
            candRates = sort(candRates);  % [f_c]

            % HSA of the candidates and preselection, Equation 146
            [~, ~, lineAmp] = shmHighResSpecAnalysis(PE, nzb, nze, candRates);
            modAmp = abs(lineAmp).^2;  % [A_i(l,z)]
            keepLines = modAmp > preSelFactor*max(modAmp);
            modRateTilde = candRates(keepLines).';  % [f_tilde_c,i(l,z)]
            if isempty(modRateTilde)
                continue
            end

            % Section 9.1.6 ECMA-418-2:2025: weighted power spectrum,
            % Equations 147 and 148 [A_tilde_i(l,z)]
            modAmpTilde = modAmp(keepLines).*shmFluctWeight(modRateTilde, bandCentreFreqs(zBand));
            [~, iMax] = max(modAmpTilde);

            % Section 9.1.7 ECMA-418-2:2025: fine tuning by the damped Newton
            % method, Equations 149 to 152. With the rate in Hz, the step
            % is limited to 5e-5 Hz, so 40 iterations move the rate by at
            % most 2e-3 Hz
            modRateStart = modRateTilde(iMax);
            modRateOpt = modRateStart;  % [f_c,1,opt(l,z)]
            for kIter = 1:newtonMaxIter
                errMinus = shmHighResSpecAnalysis(PE, nzb, nze, modRateOpt - newtonDx);
                errCentre = shmHighResSpecAnalysis(PE, nzb, nze, modRateOpt);
                errPlus = shmHighResSpecAnalysis(PE, nzb, nze, modRateOpt + newtonDx);
                errDeriv1 = (errPlus - errMinus)/(2*newtonDx);  % Equation 149
                errDeriv2 = (errPlus - 2*errCentre + errMinus)/(newtonDx^2);  % Equation 150
                newtonStep = newtonDamp*sign(errDeriv1)*min(abs(errDeriv1)/(abs(errDeriv2) + eps0),...
                                                            newtonStepCap);  % Equation 152
                modRateOpt = modRateOpt - newtonStep;  % Equation 151
                if abs(newtonStep) <= newtonTol
                    break
                end
            end
            tuneValid = abs(modRateOpt - modRateStart) <= dupTolHz;
            if ~tuneValid
                modRateOpt = modRateStart;  % optimisation failed: tuning cancelled
            end
            if modRateOpt < minModRate
                continue  % the modulation in this block is discarded
            end
            if tuneValid
                modRateTilde(iMax) = modRateOpt;
                [~, ~, lineAmpOpt] = shmHighResSpecAnalysis(PE, nzb, nze, modRateOpt);
                modAmpTilde(iMax) = abs(lineAmpOpt).^2 ...
                                    *shmFluctWeight(modRateOpt, bandCentreFreqs(zBand));
            end

            % Section 9.1.8 ECMA-418-2:2025: harmonic analysis, Equations
            % 153 to 156, for the assumed orders 1 to 3 of f_c,1,opt
            energyOrder = zeros(numel(harmOrders), 1);
            inSetOrder = false(numel(modRateTilde), numel(harmOrders));
            ratioOrder = zeros(numel(modRateTilde), numel(harmOrders));
            for iOrder = 1:numel(harmOrders)
                modRateFund = modRateOpt/harmOrders(iOrder);  % [f_c,1,o(l,z)]
                intRatio = round(modRateTilde/modRateFund);  % Equation 153
                intRatio(intRatio > harmMaxOrder) = 0;
                ratioOrder(:, iOrder) = intRatio;
                inSet = false(size(modRateTilde));
                hasRatio = intRatio > 0;
                inSet(hasRatio) = abs(modRateTilde(hasRatio)./(intRatio(hasRatio).*modRateFund) - 1) < harmTol;  % Equation 154
                inSetOrder(:, iOrder) = inSet;
                energyOrder(iOrder) = sum(modAmpTilde(inSet));  % Equation 155
            end
            [~, iOrderMax] = max(energyOrder);
            modRateFundMax = modRateOpt/harmOrders(iOrderMax);  % Equation 156 [f_1(l,z)]

            % rates of the harmonic complex corrected to integer multiples of
            % f_1, re-fitted with the constant part and one spectral line pair
            % per order; the constant part is the mean of the predictions
            inSetMax = inSetOrder(:, iOrderMax);  % [I_max(l,z)]
            modRateCorr = ratioOrder(inSetMax, iOrderMax)*modRateFundMax;
            modRateOrders = unique(modRateCorr);
            constPartOrders = zeros(numel(modRateOrders), 1);
            lineAmpCorr = zeros(numel(modRateCorr), 1);
            for iOrder = 1:numel(modRateOrders)
                [~, constPartOrders(iOrder), lineAmpOrder] = shmHighResSpecAnalysis(PE, nzb, nze, modRateOrders(iOrder));
                lineAmpCorr(modRateCorr == modRateOrders(iOrder)) = lineAmpOrder;
            end
            constPart = mean(constPartOrders);  % [p_0(l,z)]
            modAmpCorr = abs(lineAmpCorr).^2;
            modAmpTildeCorr = modAmpCorr.*shmFluctWeight(modRateCorr, bandCentreFreqs(zBand));

            % Section 9.1.9 ECMA-418-2:2025: weighting of the sum of the
            % harmonic complex, Equations 157 and 158; the exponent
            % applies to the distance between the centre of gravity of the
            % components and f_c,1,opt
            sumModAmpTilde = sum(modAmpTildeCorr);
            centreGravity = sum(modRateCorr.*modAmpTildeCorr)/(sumModAmpTilde + eps0);
            weightBW = 1 + wbwA*abs(centreGravity - modRateOpt)^wbwExp;  % [w_bw]
            modAmpHat(lBlock, zBand) = weightBW*sumModAmpTilde;

            % Section 9.1.10 ECMA-418-2:2025: HSA-based loudness, Equations
            % 160 and 161, with the nonlinearity of Section 5.1.8
            sumModAmp = sum(modAmpCorr);
            powerHSA(lBlock, zBand) = constPart^2 + 2*sumModAmp + eps0;
            rmsHSA = sqrt(0.5*(constPart^2 + 2*sumModAmp));
            specLoudHSATilde = cal_N*cal_Nx*(rmsHSA/20e-6).*prod((1 + (rmsHSA./p_threshold).^aNL).^((diff(vNL)/aNL).'), 1);
            specLoudHSA(lBlock, zBand) = max(specLoudHSATilde - LTQz(zBand), 0);

        end  % end of for loop over critical bands
    end  % end of for loop over blocks

    % Equation 159: scaling with the HSA-based loudness [A(l,z)]
    modAmpScaled = zeros(nBlocks, nBands);
    for lBlock = 1:nBlocks
        specLoudMax = max(specLoudHSA(lBlock, :));
        if specLoudMax > 0
            scaleNum = specLoudHSA(lBlock, :).^2.*blockSize1500.*modAmpHat(lBlock, :);
            scaleDen = powerHSA(lBlock, :).*(specLoudMax + eps0);
            modAmpScaled(lBlock, :) = scaleNum./scaleDen;
            modAmpScaled(lBlock, scaleDen == 0) = 0;  % bands without analysis
        end
    end
    modAmpScaled(modAmpScaled < aMinThreshold) = 0;  % threshold of Section 9.1.10

    % Time-dependent specific fluctuation strength
    % --------------------------------------------
    % Section 9.1.11 ECMA-418-2:2025

    % interpolation to 50 Hz sampling rate using a piecewise cubic
    % Hermite function; Equation 162
    l_50Last = floor(n_samples/sampleRate48k*sampleRate50) + 1;
    x = (iBlocksOut - 1)/sampleRate48k;
    xq = (0:l_50Last - 1)/sampleRate50;  % equidistant grid l_50/50
    F_est = zeros(l_50Last, nBands);
    for zBand = 1:nBands
        F_est(:, zBand) = pchip(x, modAmpScaled(:, zBand), xq);
    end  % end of for loop for interpolation
    F_est(F_est < 0) = 0;  % [F'_est(l_50,z)]

    % Equations 166 and 167: squared and linear mean over bands
    F_estRMS = rms(F_est, 2);  % [F_tilde_est]
    F_estAvg = mean(F_est, 2);  % [F_bar_est]

    % Equation 165: band distribution ratio [B_hat(l_50)]
    Bl50hat = F_estRMS./(F_estAvg + epsilon);

    % Equation 165: moving median smoothing of B_hat [B(l_50)]
    Bl50 = movmedian(Bl50hat, medianWinB);

    % Equation 164: exponent [E(l_50)]
    El50 = eCoef*(tanh(eSlope*(Bl50 - eCenter)) + 1)*0.5 + eBase;

    % Equation 163: calibration and nonlinear transformation
    % [F_hat'(l_50,z)]
    F_hat = cal_F*(F_est.^El50);

    % Equation 168: first-order lowpass filtering with tau = 0.75 s
    % along the time axis [F'(l_50,z)], with F'(0,z) = F_hat'(0,z) as
    % the initial condition
    aLP = 1 - exp(-1/(sampleRate50*tauFS));
    zLP = (1 - aLP)*F_hat(1, :);
    specFluctStrength(:, :, chan) = filter(aLP, [1, -(1 - aLP)],...
                                           F_hat, zLP, 1);

end  % end of for loop over channels

% Binaural fluctuation strength
% Section 9.1.15 ECMA-418-2:2025 [F'_B(l_50,z)]
if chansIn == 2 && binaural
    specFluctStrength(:, :, 3) = sqrt(sum(specFluctStrength(:, :, 1:2).^2, 3)/2);  % Equation 170
    chansOut = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Combined binaural"];
else
    chansOut = chansIn;  % assign number of output channels
end

% time (s) corresponding with results output [t]
timeOut = (0:(size(specFluctStrength, 1) - 1))/sampleRate50;

[~, time_skip_idx] = min( abs(timeOut-time_skip) ); % find idx of time_skip on timeOut
time_skip_idx = max(time_skip_idx, 37); % Section 9.1.12: discard l_50 = 0..35 (rows 1 to 36) at least
[~, idx_insig] = min( abs(timeInsig - time_skip) ); % find idx of time_skip on timeInsig

% Section 9.1.12 ECMA-418-2:2025
% Time-averaged specific fluctuation strength [F'(z)]
specFluctStrengthAvg = mean(specFluctStrength(time_skip_idx:end, :, :), 1); %<--- time index takes <time_skip> into consideration

% Section 9.1.13 ECMA-418-2:2025
% Time-dependent fluctuation strength Equation 169 [F(l_50)]
% Discard singleton dimensions
if chansOut == 1
    fluctStrengthTDep = sum(specFluctStrength.*dz, 2);
    specFluctStrengthAvg = transpose(specFluctStrengthAvg);
else
    fluctStrengthTDep = squeeze(sum(specFluctStrength.*dz, 2));
    specFluctStrengthAvg = squeeze(specFluctStrengthAvg);
end

% Section 9.1.14 ECMA-418-2:2025
% Overall fluctuation strength [F]
fluctStrength90Pc = prctile(fluctStrengthTDep(time_skip_idx:end, :), 90, 1); %<--- time index takes <time_skip> into consideration

%% Output assignment

% Assign outputs to structure
if chansOut == 3 % stereo case ["Stereo left"; "Stereo right"; "Combined binaural"];

    % outputs only with ["Stereo left"; "Stereo right"]
    OUT.specFluctStrength = specFluctStrength(:, :, 1:2);
    OUT.specFluctStrengthAvg = specFluctStrengthAvg(:, 1:2);
    OUT.fluctStrengthTDep = fluctStrengthTDep(:, 1:2);
    OUT.fluctStrength90Pc = fluctStrength90Pc(:, 1:2);

    % outputs only with  "combined binaural"
    OUT.specFluctStrengthBin = specFluctStrength(:, :, 3);
    OUT.specFluctStrengthAvgBin = specFluctStrengthAvg(:, 3);
    OUT.fluctStrengthTDepBin = fluctStrengthTDep(:, 3);
    OUT.fluctStrength90PcBin = fluctStrength90Pc(:, 3);

    % general outputs
    OUT.bandCentreFreqs = bandCentreFreqs;
    OUT.timeOut = timeOut;
    OUT.timeInsig = timeInsig;
    OUT.soundField = fieldtype;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fluctuation strength statistics based on fluctStrengthTDep ["Stereo left"; "Stereo right"; "Combined binaural"];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    metric_statistics = 'FluctuationStrength_ECMA418_2';
    OUT_statistics = get_statistics( fluctStrengthTDep(time_skip_idx:end,1:chansOut), metric_statistics ); % get statistics

    % copy fields of <OUT_statistics> struct into the <OUT> struct
    fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

    for i = 1:numel(fields_OUT_statistics)
        fieldName = fields_OUT_statistics{i};
        if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
            OUT.(fieldName) = OUT_statistics.(fieldName);
        end
    end

    clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else % mono case

    OUT.specFluctStrength = specFluctStrength;
    OUT.specFluctStrengthAvg = specFluctStrengthAvg;
    OUT.fluctStrengthTDep = fluctStrengthTDep;
    OUT.fluctStrength90Pc = fluctStrength90Pc;

    OUT.bandCentreFreqs = bandCentreFreqs;
    OUT.timeOut = timeOut;
    OUT.timeInsig = timeInsig;
    OUT.soundField = fieldtype;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fluctuation strength statistics based on fluctStrengthTDep (mono case)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    metric_statistics = 'FluctuationStrength_ECMA418_2';
    OUT_statistics = get_statistics( fluctStrengthTDep(time_skip_idx:end), metric_statistics ); % get statistics

    % copy fields of <OUT_statistics> struct into the <OUT> struct
    fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

    for i = 1:numel(fields_OUT_statistics)
        fieldName = fields_OUT_statistics{i};
        if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
            OUT.(fieldName) = OUT_statistics.(fieldName);
        end
    end

    clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% Output plotting

if show

    % colormap
    cmap_inferno = load('cmap_inferno.txt');

    % generate A-weighting filter for LAeq calculation
    [b, a] = Gen_weighting_filters(fs, 'A');
    insig_A = filter(b, a, insig);  % filter signal
    LAeq_all = 20*log10(rms(insig_A(idx_insig:end, :))./2e-5);  % calculate LAeq

    for chan = chansOut:-1:1
        % Plot results
        fig = figure('name', sprintf( 'Fluctuation strength analysis - ECMA-418-2 (%s signal)', chans(chan) ) );
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specFluctStrength(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(specFluctStrength(:, :, chan), [], 'all')*500)/500];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3];
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"];
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency (Hz)';
        ax1.XLabel.String = 'Time (s)';
        ax1.FontName =  'Times';
        ax1.FontSize = 11;
        colormap(cmap_inferno);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific fluctuation strength,'; '(vacil_{HMS}/Bark_{HMS})'});

        if chan == 3 % the binaural channel

             % take the higher channel level as representative (PD ISO/TS 12913-3:2019 Annex D)
            [LAeq, LR] = max(LAeq_all);

            % if branch to identify which channel is higher
            if LR == 1
                whichEar = 'left ear';
            else
                whichEar = 'right ear';
            end  % end of if branch

        elseif chan == 2 % Stereo right

            LAeq = LAeq_all(chan);
            whichEar = 'right ear';

        elseif chan == 1 % Stereo left or mono

            LAeq = LAeq_all(chan);
            if chansOut~=1
                whichEar = 'left ear';
            else
                whichEar = 'mono';
            end
        end

        titleString = sprintf('%s signal, $L_{\\textrm{Aeq,%s}} =$ %.3g (dB SPL)', chans(chan), whichEar, LAeq);

        title(titleString, 'Interpreter','Latex' );

        ax2 = nexttile(2);
        plot(ax2, timeOut, fluctStrengthTDep(:, chan), 'color', cmap_inferno(166, :),...
            'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold on;
        plot(ax2, timeOut, fluctStrength90Pc(1, chan)*ones(size(timeOut)),'--', 'color',...
            cmap_inferno(34, :), 'LineWidth', 1, 'DisplayName', "90th" + string(newline) + "percentile");
        hold off;
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];

        if max(fluctStrengthTDep(:, chan)) > 0
            ax2.YLim = [0, 1.1*ceil(max(fluctStrengthTDep(:, chan))*10)/10];
        end

        ax2.XLabel.String = 'Time (s)';
        ax2.YLabel.String = 'Fluctuation strength (vacil_{HMS})';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Times';
        ax2.FontSize = 11;
        legend('Location', 'eastoutside', 'FontSize', 8);
        set(gcf,'color','w');

    end  % end of for loop for plotting over channels
end  % end of if branch for plotting if outplot true

end %of function
