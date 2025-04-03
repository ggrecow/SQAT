function OUT = Tonality_ECMA418_2(insig, fs, fieldtype, time_skip, show)
% function OUT = Tonality_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%
% Returns tonality values and frequencies according to ECMA-418-2:2024
% (using the Sottek Hearing Model) for an input calibrated single mono
% or single stereo audio (sound pressure) time-series signal, insig. For stereo
% signals, Tonality is calculated for each channel [left ear, right ear].  
% 
% According to  ECMA-418-2:2024 (Section 6.2.11), the representative 
% single value to express the
% overall tonality is the time-averaged overall tonality, provided here by
% the <tonalityAvg> output variable.
%
% NOTE: according to ECMA-418-2:2024 (Section 6.2.9), the tonality
% values corresponding to the first 300 ms of the input signal must be discarded 
% due to the transient responses of the digital filters. Considering the temporal resolution
% of the tonality model, this means any values below 304 ms must be discarded.
% Therefore, <time_skip> must be greater or equal to 304 ms to compute any
% time-aggregated quantity. 
%
%  Reference signal: pure tone with center frequency of 1 kHz and RMS value
%  of 40 dBSPL equals 1 tu_HMS
%
% Inputs
% ------
% insig : column vector [Nx1] mono or [Nx2] binaural
%            the input signal as single mono or stereo audio (sound
%            pressure) signals
%
% fs : integer
%                the sample rate (frequency) of the input signal(s)
%
% fieldtype : keyword string (default: 'free-frontal')
%                 determines whether the 'free-frontal' or 'diffuse' field stages
%                 are applied in the outer-middle ear filter
%
% time_skip : integer (default: 304 ms - see Section 6.2.9 ECMA-418-2:2024)
%                   skip start of the signal in <time_skip> seconds so that
%                   the transient response of the digital filters is avoided.
%                   Best-practice: <time_skip> must be equal or higher than
%                   default value
%
% show : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
% 
% Returns
% -------
%
% OUT : structure
%            contains the following fields:
%
% specTonality : matrix
%                time-dependent specific tonality for each (half) critical
%                band
%                arranged as [time, bands(, channels)]
%
% specTonalityFreqs : matrix
%                     time-dependent frequencies of the dominant tonal
%                     components corresponding with each of the
%                     time-dependent specific tonality values in each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% specTonalityAvg : matrix
%                   time-averaged specific tonality for each (half)
%                   critical band
%                   arranged as [bands(, channels)]
%                   OBS: takes <time_skip> into consideration
%
% specTonalityAvgFreqs : matrix
%                        frequencies of the dominant tonal components
%                        corresponding with each of the
%                        time-averaged specific tonality values in each
%                        (half) critical band
%                        arranged as [bands(, channels)]
%                        OBS: takes <time_skip> into consideration
%
% specTonalLoudness : matrix
%                     time-dependent specific tonal loudness for each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% specNoiseLoudness : matrix
%                     time-dependent specific noise loudness for each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% tonalityTDep : vector or matrix
%                time-dependent overall tonality
%                arranged as [time(, channels)]
%
% tonalityTDepFreqs : vector or matrix
%                     time-dependent frequencies of the dominant tonal
%                     components corresponding with the
%                     time-dependent overall tonality values
%                     arranged as [time(, channels)]
%
% tonalityAvg : number or vector
%               time-averaged overall tonality
%               arranged as [tonality(, channels)]
%               OBS: takes <time_skip> into consideration
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each (half)
%                   critical band rate scale width
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% timeInsig : vector
%           time (seconds) of insig
%
% soundField : string
%              identifies the soundfield type applied (the input argument
%              fieldtype)
%
% Several statistics based on tonalityTDep
%         ** Tmean : mean value of instantaneous tonality (tu_HMS)
%         ** Tstd : standard deviation of instantaneous tonality (tu_HMS)
%         ** Tmax : maximum of instantaneous tonality (tu_HMS)
%         ** Tmin : minimum of instantaneous tonality (tu_HMS)
%         ** Tx : tonality value exceeded during x percent of the time (tu_HMS)
%             in case of binaural input, Tx(1,3), being 1st and 2nd column 
%             corresponding to [left ear, right ear] 
%             OBS: all quantities here take <time_skip> into consideration
%
% If show==true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall tonality, with the latter also indicating the time-aggregated
% value. A set of plots is returned for each input channel.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 07/08/2023
% Date last modified: 02/04/2025
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license,
% as detailed in the original code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics). 
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% This code calls sub-component file 'cmap_plasma.txt'. The contents of
% the file includes a copy of data obtained from the repository 
% https://github.com/BIDS/colormap, and is CC0 1.0 licensed for modified
% use, see https://creativecommons.org/publicdomain/zero/1.0 for
% information.
%
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by: Gil Felix Greco
% Date last checked: 16.02.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Arguments validation
    arguments (Input) % Matlab R2018b or newer
        insig (:, :) double {mustBeReal}
        fs (1, 1) double {mustBePositive, mustBeInteger}
        fieldtype (1,:) string {mustBeMember(fieldtype,...
                                                       {'free-frontal',...
                                                        'diffuse'})} = 'free-frontal'
        time_skip (1, 1) double {mustBeReal} = 304e-3
        show {mustBeNumericOrLogical} = false
    end

%% Input checks

% define time threshold value from which all values before must be dropped.
t_threshold =  304e-3;

% check insig dimension (only [Nx1] or [Nx2] are valid)
if  size(insig,1) > 2 & size(insig,2) > 2 % insig has more than 2 channels
    error('Error: Input signal has more than 2 channels. ')
elseif  size(insig, 2) > 2  % insig is [1xN] or [2xN]
    insig = insig';
    fprintf('\nWarning: Input signal is not [Nx1] or [Nx2] and was transposed.\n');
end

% Check the length of the input data (must be at least 304 ms)
if size(insig, 1) <  t_threshold*fs
    error('Error: Input signal is too short along the specified axis to calculate tonality (must be at least 304 ms)')
end

% Check the channel number of the input data
if size(insig, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    inchans = size(insig, 2);
    if inchans > 1
        chans = ["Stereo left"; "Stereo right"];
    else
        chans = "Mono";
    end
end

if time_skip<t_threshold
        warning("Time_skip must be at least 304 ms to avoid transient responses of the digital filters (see ECMA-418-2:2024 (Section 6.2.9)). Setting time_skip to 304 ms!!!")
        time_skip = t_threshold;
end

%% Define constants

sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2024 [r_s]
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2024 [deltaf(f=0)]
c = 0.1618;  % Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2024

halfBark = 0.5:0.5:26.5;  % half-critical band rate scale [z]
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2024 [F(z)]
dfz = sqrt(deltaFreq0^2 + (c*bandCentreFreqs).^2);  % Section 5.1.4.1 Equation 10 ECMA-418-2:2024 [deltaf(z)]

% Block and hop sizes Section 6.2.2 Table 4 ECMA-418-2:2024
overlap = 0.75;  % block overlap proportion
% block sizes [s_b(z)]
blockSize = [8192*ones(1, 3), 4096*ones(1, 13), 2048*ones(1, 9), 1024*ones(1, 28)];
% hop sizes (section 5.1.2 footnote 3 ECMA 418-2:2022) [s_h(z)]
hopSize = (1 - overlap)*blockSize;

% Output sample rate based on hop sizes - Resampling to common time basis
% Section 6.2.6 ECMA-418-2:2024 [r_sd]
sampleRate1875 = sampleRate48k/min(hopSize);

% Number of bands that need averaging. Section 6.2.3 Table 5
% ECMA-418-2:2024 [NB]
NBandsAvg = [0, 1, 2*ones(1,14), ones(1,9), zeros(1,28);...
             1, 1, 2*ones(1,14), ones(1,9), zeros(1,28)];

% Critical band interpolation factors from Section 6.2.6 Table 6
% ECMA-418-2:2024 [i]
i_interp = blockSize/min(blockSize);

% Noise reduction constants from Section 6.2.7 Table 7 ECMA-418-2:2024
alpha = 20;
beta = 0.07;

% Sigmoid function factor parameters Section 6.2.7 Table 8 ECMA-418-2:2024
% [c(s_b(z))]
csz_b = [18.21*ones(1, 3), 12.14*ones(1, 13), 417.54*ones(1, 9),...
         962.68*ones(1, 28)]; 
% [d(s_b(z))]
dsz_b = [0.36*ones(1, 3), 0.36*ones(1, 13), 0.71*ones(1, 9),...
         0.69*ones(1, 28)]; 

% Scaling factor constants from Section 6.2.8 Table 9 ECMA-418-2:2024
A = 35;
B = 0.003;

cal_T = 2.8758615;  % calibration factor in Section 6.2.8 Equation 51 ECMA-418-2:2024 [c_T]
cal_Tx = 1/0.9999043734252;  % Adjustment to calibration factor (Footnote 22 ECMA-418-2:2024)

%% Signal processing

% Input pre-processing
% --------------------
if fs ~= sampleRate48k  % Resample signal
    [p_re, ~] = shmResample(insig, fs);
else  % don't resample
    p_re = insig;
end

% get time vector of input signal
timeInsig = (0 : length(p_re(:,1))-1) ./ fs;

% Section 5.1.2 ECMA-418-2:2024 Fade in weighting and zero-padding
pn = shmPreProc(p_re, max(blockSize), max(hopSize));

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2024 Outer and middle/inner ear signal filtering
pn_om = shmOutMidEarFilter(pn, fieldtype);

% Loop through channels in file
% -----------------------------
for chan = size(pn_om, 2):-1:1

    % Apply auditory filter bank
    % --------------------------
    
    % Filter equalised signal using 53 1/2Bark ERB filters according to 
    % Section 5.1.4.2 ECMA-418-2:2024
    pn_omz = shmAuditoryFiltBank(pn_om(:, chan), false);

    % Autocorrelation function analysis
    % ---------------------------------
    % Duplicate Banded Data for ACF
    % Averaging occurs over neighbouring bands, to do this the segmentation
    % needs to be duplicated for neigbouring bands. 'Dupe' has been added
    % to variables to indicate that the vectors/matrices have been modified
    % for duplicated neigbouring bands.
    
    pn_omzDupe = [pn_omz(:, 1:5), pn_omz(:, 2:18), pn_omz(:, 16:26),...
                   pn_omz(:, 26:53)];
    blockSizeDupe = [8192*ones(1, 5), 4096*ones(1, 17), 2048*ones(1, 11),...
                     1024*ones(1, 28)];
    bandCentreFreqsDupe = [bandCentreFreqs(1:5),...
                           bandCentreFreqs(2:18),...
                           bandCentreFreqs(16:26),...
                           bandCentreFreqs(26:53)];
    % (duplicated) indices corresponding with the NB bands around each z band
    i_NBandsAvgDupe = [1, 1, 1, 6:18, 23:31, 34:61;
                       2, 3, 5, 10:22, 25:33, 34:61];

    for zBand = 61:-1:1

        % Segmentation into blocks
        % ------------------------
        % Section 5.1.5 ECMA-418-2:2024
        i_start = blockSizeDupe(1) - blockSizeDupe(zBand) + 1;
        [pn_lz, ~] = shmSignalSegment(pn_omzDupe(:, zBand), 1,...
                                      blockSizeDupe(zBand), overlap, i_start);
 
        % Transformation into Loudness
        % ----------------------------
        % Sections 5.1.6 to 5.1.9 ECMA-418-2:2024 [N'_basis(z)]
        [pn_rlz, bandBasisLoudness, ~]...
            = shmBasisLoudness(pn_lz, bandCentreFreqsDupe(zBand));
        basisLoudnessArray{zBand} = bandBasisLoudness;

        % Apply ACF
        % ACF implementation using DFT
        % Section 6.2.2 Equations 27 & 28 ECMA-418-2:2024
        % [phi_unscaled,l,z(m)]
        unscaledACF = ifft(abs(fft(pn_rlz, 2*blockSizeDupe(zBand), 1)).^2,...
                           2*blockSizeDupe(zBand), 1);
        % Section 6.2.2 Equation 29 ECMA-418-2:2024 [phi_l,z(m)]
        denom = sqrt(cumsum(pn_rlz.^2, 1, 'reverse').*flipud(cumsum(pn_rlz.^2)))...
                + 1e-12;

        % note that the block length is used here, rather than the 2*s_b,
        % for compatability with the remaining code - beyond 0.75*s_b is
        % assigned (unused) zeros in the next line
        unbiasedNormACF = unscaledACF(1:blockSizeDupe(zBand), :)./denom;
        unbiasedNormACF((0.75*blockSizeDupe(zBand) + 1):blockSizeDupe(zBand), :) = 0;

        % Section 6.2.2 Equation 30 ECMA-418-2:2024 [phi_z'(m)
        unbiasedNormACFDupe{zBand} = basisLoudnessArray{zBand}.*unbiasedNormACF;

    end
    
    % Average the ACF over nB bands - Section 6.2.3 ECMA-418-2:2024        
    for zBand = 53:-1:1  % Loop through 53 critical band filtered signals
        
        NBZ = NBandsAvg(1, zBand) + NBandsAvg(2, zBand) + 1; % Total number of bands to average over

        % Averaging of frequency bands
        meanScaledACF = mean(reshape(cell2mat(unbiasedNormACFDupe(i_NBandsAvgDupe(1, zBand):i_NBandsAvgDupe(2, zBand))),...
                                     blockSize(zBand), [], NBZ), 3);

        % Average the ACF over adjacent time blocks [phibar_z'(m)]
        if zBand <= 16 
            meanScaledACF(:, 2:end - 1) = movmean(meanScaledACF, 3, 2, 'omitnan',...
                                                  'EndPoints', 'discard');
        end
        
        % Application of ACF lag window Section 6.2.4 ECMA-418-2:2024
        tauz_start = max(0.5/dfz(zBand), 2e-3);  % Equation 31 ECMA-418-2:2024 [tau_start(z)]
        tauz_end = max(4/dfz(zBand), tauz_start + 1e-3);  % Equation 32 ECMA-418-2:2024 [tau_end(z)]
        % Equations 33 & 34 ECMA-418-2:2024
        mz_start = ceil(tauz_start*sampleRate48k);  % Starting lag window index [m_start(z)]
        mz_end = floor(tauz_end*sampleRate48k);  % Ending lag window index [m_end(z)]
        M = mz_end - mz_start + 1;
        % Equation 35 ECMA-418-2:2024
        % lag-windowed, detrended ACF [phi'_z,tau(m)]
        lagWindowACF = zeros(size(meanScaledACF));
        lagWindowACF(mz_start:mz_end, :) = meanScaledACF(mz_start:mz_end, :)...
                                           - mean(meanScaledACF(mz_start:mz_end, :));
        
        % Estimation of tonal loudness
        % ----------------------------
        % Section 6.2.5 Equation 36 ECMA-418-2:2024
        % ACF spectrum in the lag window [Phi'_z,tau(k)]
        magFFTlagWindowACF = abs(fft(lagWindowACF, 2*max(blockSize), 1));
        magFFTlagWindowACF(isnan(magFFTlagWindowACF)) = 0;
    
        % Section 6.2.5 Equation 37 ECMA-418-2:2024 [Nhat'_tonal(z)]
        % first estimation of specific loudness of tonal component in critical band
        bandTonalLoudness = meanScaledACF(1, :);
        mask = 2*max(magFFTlagWindowACF, [], 1)/(M/2) <= meanScaledACF(1, :);
        bandTonalLoudness(mask) = 2*max(magFFTlagWindowACF(:, mask), [], 1)/(M/2);
    
        % Section 6.2.5 Equation 38 & 39 ECMA-418-2:2024
        % [k_max(z)]
        [~, kz_max] = max(magFFTlagWindowACF, [], 1);
        % frequency of maximum tonal component in critical band [f_ton(z)]
        bandTonalFreqs = (kz_max - 1)*(sampleRate48k/(2*max(blockSize)));

        % Section 6.2.7 Equation 41 ECMA-418-2:2024 [N'_signal(l,z)]
        % specific loudness of complete band-pass signal in critical band
        bandLoudness = meanScaledACF(1, :);
        
        % Resampling to common time basis Section 6.2.6 ECMA-418-2:2024
        if i_interp(zBand) > 1
            % Note: use of interpolation function avoids rippling caused by
            % resample function, which otherwise affects specific loudness 
            % calculation for tonal and noise components
            l_n = size(meanScaledACF, 2);
            x = linspace(1, l_n, l_n);
            xq = linspace(1, l_n, i_interp(zBand)*(l_n - 1) + 1);
            bandTonalLoudness = interp1(x, bandTonalLoudness, xq);
            bandLoudness = interp1(x, bandLoudness, xq);
            bandTonalFreqs = interp1(x, bandTonalFreqs, xq);

        end

        % Remove end zero-padded samples Section 6.2.6 ECMA-418-2:2024
        l_end = ceil(size(p_re, 1)/sampleRate48k*sampleRate1875) + 1;  % Equation 40 ECMA-418-2:2024

        bandTonalLoudness = bandTonalLoudness(1:l_end);
        bandLoudness = bandLoudness(1:l_end);
        bandTonalFreqs = bandTonalFreqs(1:l_end);

        % Noise reduction Section 6.2.7 ECMA-418-2:2020
        % ---------------------------------------------
        % Equation 42 ECMA-418-2:2024 signal-noise-ratio first approximation
        % (ratio of tonal component loudness to non-tonal component loudness in critical band)
        % [SNRhat(l,z)]
        SNRlz1 = bandTonalLoudness./((bandLoudness - bandTonalLoudness) + 1e-12);

        % Equation 43 ECMA-418-2:2024 low pass filtered specific loudness
        % of non-tonal component in critical band [Ntilde'_tonal(l,z)]
        bandTonalLoudness = shmNoiseRedLowPass(bandTonalLoudness, sampleRate1875);

        % Equation 44 ECMA-418-2:2024 lowpass filtered SNR (improved estimation)
        % [SNRtilde(l,z)]
        SNRlz = shmNoiseRedLowPass(SNRlz1, sampleRate1875);

        % Equation 46 ECMA-418-2:2024 [g(z)]
        gz = csz_b(zBand)/(bandCentreFreqs(zBand)^dsz_b(zBand));

        % Equation 45 ECMA-418-2:2024 [nr(l,z)]
        crit = exp(-alpha*((SNRlz/gz) - beta));
        nrlz = 1 - crit;  % sigmoidal weighting function
        nrlz(crit >= 1) = 0;

        % Equation 47 ECMA-418-2:2024 [N'_tonal(l,z)]
        bandTonalLoudness = nrlz.*bandTonalLoudness;

        % Section 6.2.8 Equation 48 ECMA-418-2:2024 [N'_noise(l,z)]
        bandNoiseLoudness = shmNoiseRedLowPass(bandLoudness, sampleRate1875) - bandTonalLoudness;  % specific loudness of non-tonal component in critical band
    
        % Store critical band results
        % ---------------------------
        % specific time-dependent signal-noise-ratio in each critical band
%         specSNR(:, zBand, chan) = SNRlz1;

        % specific time-dependent loudness of signal in each critical band
%         specLoudness(:, zBand, chan) = bandLoudness;

        % specific time-dependent loudness of tonal component in each critical band  [N'_tonal(l,z)]
        specTonalLoudness(:, zBand, chan) = bandTonalLoudness;

        % specific time-dependent loudness of non-tonal component in each critical band [N'_noise(l,z)]
        specNoiseLoudness(:, zBand, chan) = bandNoiseLoudness;

        % time-dependent frequency of tonal component in each critical band [f_ton(z)]
        specTonalityFreqs(:, zBand, chan) = bandTonalFreqs;

    end

    % Calculation of specific tonality
    % --------------------------------
    % Section 6.2.8 Equation 49 ECMA-418-2:2024 [SNR(l)]
    overallSNR = max(specTonalLoudness, [], 2)./(1e-12 + sum(specNoiseLoudness, 2));  % loudness signal-noise-ratio
    
    % Section 6.2.8 Equation 50 ECMA-418-2:2024 [q(l)]
    crit = exp(-A*(overallSNR - B));
    ql = 1 - crit;  % sigmoidal scaling factor
    ql(crit >= 1) = 0;
    
    % Section 6.2.8 Equation 51 ECMA-418-2:2024 [T'(l,z)]
    specTonality = cal_T*cal_Tx*ql.*specTonalLoudness;  % time-dependent specific tonality
    
    % Section 6.2.8 Equation 52 ECMA-418-2:2024
    % time (s) corresponding with results output [t]
    timeOut = (0:(size(specTonality, 1) - 1))/sampleRate1875;

    [~, time_skip_idx] = min( abs(timeOut-time_skip) ); % find idx of time_skip on timeOut
    [~, idx_insig] = min( abs(timeInsig - time_skip) ); % find idx of time_skip on timeInsig

    % Calculation of time-averaged specific tonality Section 6.2.9
    % ECMA-418-2:2024 [T'(z)]
    for zBand = 53:-1:1
        mask = specTonality(:, zBand, chan) > 0.02;  % criterion Section 6.2.9 point 2
        mask(1:(time_skip_idx - 1)) = 0;  % criterion Section 6.2.9 point 1 %<--- time index takes <time_skip> into consideration

        % Section 6.2.9 Equation 53 ECMA-418-2:2024  
        specTonalityAvg(1, zBand, chan)...
            = sum(specTonality(mask, zBand, chan), 1)./(nnz(mask) + 1e-12); %<--- time index takes <time_skip> into consideration
        specTonalityAvgFreqs(1, zBand, chan)...
            = sum(specTonalityFreqs(mask, zBand, chan), 1)./(nnz(mask) + 1e-12); %<--- time index takes <time_skip> into consideration
    end

    % Calculation of total (non-specific) tonality Section 6.2.10
    % -----------------------------------------------------------
    % Further update can add the user input frequency range to determine
    % total tonality - not yet incorporated

    % Section 6.2.10 Equation 61 ECMA-418-2:2024
    % Time-dependent total tonality [T(l)]
    [tonalityTDep(:, chan), zmax] = max(specTonality(:, :, chan), [], 2);

    for ll = size(specTonalityFreqs, 1):-1:1
        tonalityTDepFreqs(ll, chan) = specTonalityFreqs(ll, zmax(ll), chan);
    end
    
    % Calculation of representative values Section 6.2.11 ECMA-418-2:2024
    % Time-averaged total tonality
    mask = tonalityTDep(:, chan) > 0.02;  % criterion Section 6.2.9 point 2
    mask(1:(time_skip_idx - 1)) = 0;    % criterion Section 6.2.9 point 1 %<--- time index takes <time_skip> into consideration

    % Section 6.2.11 Equation 63 ECMA-418-2:2024
    % Time-averaged total tonality [T]
    tonalityAvg(chan) = sum(tonalityTDep(mask, chan))/(nnz(mask) + eps); %<--- time index takes <time_skip> into consideration

%% Output plotting

    if show

        % colormap
        cmap_plasma = load('cmap_plasma.txt');

        % generate A-weighting filter for LAeq calculation
        [b, a] = Gen_weighting_filters(fs, 'A');
        insig_A = filter(b, a, insig);  % filter signal
        LAeq_all = 20*log10(rms(insig_A(idx_insig:end, :))./2e-5);  % calculate LAeq

        % Plot results
        fig = figure('name', sprintf( 'Tonality analysis - ECMA-418-2 (%s signal)', chans(chan) ) );
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specTonality(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(tonalityTDep(:, chan))*10)/10];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"]; 
        ax1.YScale = 'log';
        ax1.YLabel.String = "Frequency (Hz)";
        ax1.XLabel.String = "Time (s)";
        ax1.FontName =  'Times';
        ax1.FontSize = 11;
        colormap(cmap_plasma);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific Tonality,'; '(tu_{HMS}/Bark_{HMS})'});

        LAeq = LAeq_all(chan);

        titleString = sprintf('%s signal, $L_{\\textrm{Aeq}} =$ %.3g (dB SPL)', chans(chan), LAeq);

        title(titleString, 'Interpreter','Latex' );
        
        ax2 = nexttile(2); 

        plot(ax2, timeOut, tonalityTDep(:, chan), 'color',  cmap_plasma(166, :),...
            'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold on
        plot(ax2, timeOut, tonalityAvg(1, chan)*ones(size(timeOut)), '--', 'color', cmap_plasma(34, :),...
            'LineWidth', 1, 'DisplayName', "Time-" + string(newline) + "average");
        hold off

        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax2.YLim = [0, 1.1*ceil(max(tonalityTDep(:, chan))*10)/10];
        ax2.XLabel.String = "Time (s)";
        ax2.YLabel.String = "Tonality (tu_{HMS})";
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Times';
        ax2.FontSize = 11;
        legend('Location', 'eastoutside', 'FontSize', 8);
        set(gcf,'color','w');
    end

end

%% Output assignment

% Discard singleton dimensions
if inchans > 1
    specTonalityAvg = squeeze(specTonalityAvg);
    specTonalityAvgFreqs = squeeze(specTonalityAvgFreqs);
else
    specTonalityAvg = transpose(specTonalityAvg);
    specTonalityAvgFreqs = transpose(specTonalityAvgFreqs);
end

OUT.specTonality = specTonality;
OUT.specTonalityAvg = specTonalityAvg;
OUT.specTonalityFreqs = specTonalityFreqs;
OUT.specTonalityAvgFreqs = specTonalityAvgFreqs;

OUT.specTonalLoudness = specTonalLoudness;
OUT.specNoiseLoudness = specNoiseLoudness;

OUT.tonalityTDep = tonalityTDep;
OUT.tonalityAvg = tonalityAvg;
OUT.tonalityTDepFreqs = tonalityTDepFreqs;
OUT.bandCentreFreqs = bandCentreFreqs;

OUT.timeOut = timeOut;
OUT.timeInsig = timeInsig;
OUT.soundField = fieldtype;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tonality statistics based on tonalityTDep ["Stereo left"; "Stereo right"]; for stereo case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metric_statistics = 'Tonality_ECMA418_2';
OUT_statistics = get_statistics( tonalityTDep(time_skip_idx:end,1:size(tonalityTDep, 2)), metric_statistics ); % get statistics

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

end %of function
