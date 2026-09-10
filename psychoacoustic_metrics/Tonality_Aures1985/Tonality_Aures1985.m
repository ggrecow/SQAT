function OUT = Tonality_Aures1985(insig,fs,LoudnessField,time_skip,show)
% function OUT = Tonality_Aures1985(insig,fs,LoudnessField,time_skip,show)
%
%   This function calculates tonality metric by:
%
%   [1] Aures, Wilhelm (1985). "Berechnungsverfahren fuer den sensorischen Wohlklang 
%       beliebiger Schallsignale." Acta Acustica united with Acustica 59: p. 130-141.
%
%   The Aures' tonality is based on Terhard's virtual pitch theory, given by:
%
%   [2] Terhardt, E., Stoll, G. and Seewann, M. (1982). Algorithm for 
%       extraction of pitch and pitch salience from complex tonal signals. 
%       J. Acoust. Soc. Am., 71, 679-688. doi:10.1121/1.387544
%
%  Loudness calculation is conducted according to ISO 532:1-2017
%  (type <help Loudness_ISO532_1> for more info)
%
%   Reference: a pure tone with 1000 Hz and 60 dBSPL has a tonality of 1 t.u.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT ARGUMENTS
%   insig : array
%   acoustic signal, monophonic (Pa)
%
%   fs : integer
%   sampling frequency (Hz).
%
%   LoudnessField : integer
%   chose field for loudness calculation; free field = 0; diffuse field = 1;
%   type <help Loudness_ISO532_1> for more info

%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistic calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUT:
%   OUT : struct containing the following fields
%
%       * InstantaneousTonality: instantaneous tonality (t.u.)  vs time
%       * TonalWeighting: tonal weighting as a function of time
%       * LoudnessWeighting: loudness weighting as a function of time
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousTonality
%         ** Kmean : mean value of InstantaneousTonality (t.u.)
%         ** Kstd : standard deviation of InstantaneousTonality (t.u.)
%         ** Kmax : maximum of InstantaneousTonality (t.u.)
%         ** Kmin : minimum of InstantaneousTonality (t.u.)
%         ** Kx : Tonality value exceeded during x percent of the time (t.u.)
%
% Author: Gil Felix Greco, Braunschweig 13/07/2020 (updated 14.04.2023)
% Author: Gil Felix Greco, Braunschweig 16.02.2025 - introduced get_statistics function
% Author: Sergio Aguirre, September 2026 - the tones are removed from the
%   windowed spectrum and the notch is at least one main lobe wide, so the
%   result no longer depends on where the tone falls between two FFT bins
% Author: Sergio Aguirre, September 2026 - the sound pressure excess is now
%   stored per tonal component, and the bins that replace a tone keep the
%   phase they already had, so the function is deterministic
% Author: Sergio Aguirre, September 2026 - the regions of the spectrum that
%   are narrower than a critical band are extracted as tonal components as
%   well, as Aures asks in section 2.3.2, and the noise term of the level
%   excess is summed over the spectrum those extractions leave behind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

%% resampling
% resampling audio to 44.1 kHz or 48kHz
if ~(fs == 44100 || fs == 48000)
    gcd_fs = gcd(44100,fs); % greatest common denominator
    insig = resample(insig,44100/gcd_fs,fs/gcd_fs);
    fs = 44100;
end

%% window parameters

% time_resolution=80e-3;    % window length fixed in 80 ms (Terhard), gives a df=12.5 Hz
time_resolution=250e-3; % window length, chosen from the scale over which w1
                        % varies; gives df = 4 Hz at both sampling rates

N=round(fs*time_resolution); % define window length, N bins
window = hann(N);

fftgain = 2^0.5/(N*mean(hann(N))); % gain to be applied based on the FFT length
ENBW = N*sum(window.^2)/sum(window)^2; % equivalent noise bandwidth of the window, in bins (1.5 for Hann)

%% freq vectors based on window input signals

% from Terhardt [3]: Aurally relevant tonal information of any signal is
%  confined in the frequency region of about 20 Hz to 5 kHz.

MinFrequency=20;  
MinFrequencyindex = ceil( 1 + ( MinFrequency*(N/fs) ) ); % index corresponding to min frequency (20 Hz) for tone extraction

MaxFrequency=5000;   
MaxFrequencyIndex = ceil( 1 + ( MaxFrequency*(N/fs) ) ); % index corresponding to max frequency (5 kHz) for tone extraction

Freq = fs*((1:round(N))'-1)/N;  % freq vector
FreqCrop = Freq(MinFrequencyindex:MaxFrequencyIndex); % croped freq vector from MinFrequencyindex till MaxFrequencyIndex
df=FreqCrop(2)-FreqCrop(1); % freq discretization

%% initialize windowed vectors

t_b = ( 1:length(insig) )/fs; % time vector
 
overlap = round(0.5*N);       % overlap 

insig = buffer(insig,N,overlap,'nodelay');
t_b = buffer(t_b,N,overlap,'nodelay');

nFrames = size(insig,2)-1;

tone=cell(nFrames,1);            % Memory allocation: tone cell per time frame
tonality=zeros(nFrames,1);       % Memory allocation for tonality computation
t=zeros(nFrames,1);              % Memory allocation: time vector for iFrames
w_gr=zeros(nFrames,1);           % Memory allocation: loudness weighting function per time frame
w_tonal=zeros(nFrames,1);        % Memory allocation: tonal weighting function per time frame
TINY_VALUE = 1e-99;

%% Here we go ...

for iFrame = 1:nFrames
    
    %% windowed time-frame
    
    Winsig = insig(:,iFrame);     % cut insig for each iFrames
    
    t(iFrame,1) = t_b(1,iFrame);  % output time vector for iFrames
    
    Winsig = ( window.*Winsig );  % Apply window to frame
    
    %% compute SPL for each time-frame
    
    SpectralEnergy = abs( fft(Winsig.*fftgain) ).^2;
    SPL = 10.*log10( (SpectralEnergy+TINY_VALUE)./4e-10 ); % dBSPL    
        
    %%%% check plot (only for debugging) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      figure; plot(Freq,SPL) % check plot   
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Find peaks according to Terhard's criteria for each time-frame
    
    SPLcrop = SPL(MinFrequencyindex:MaxFrequencyIndex); % crop SPL vector from MinFrequencyindex to MaxFrequencyIndex
    
    threshold = 7;  % condition for tonal component, in dBSPL

    % The criterion of Terhardt asks for the level two and three bins away. Those
    % two bins were 25 and 37.5 Hz at the resolution of the original method, so
    % the distances are kept in hertz and the criterion no longer moves with the
    % length of the analysis window.
    k2 = max(round(25/df),1);
    k3 = max(round(37.5/df),2);
    
    ToneIdx = zeros(length(SPLcrop),1); % initialize vector, tonal components idx
    k = 1; % initialize counter
    
    % find tones... The bins the criterion compares with are read from the
    % full spectrum, so that a tone near the lower edge of the range still has
    % neighbours below it: with the distances in hertz the criterion reaches
    % 37.5 Hz down, and read from the cropped vector alone it could not fire
    % below 57.5 Hz.
    for i = 1:length(SPLcrop)
        
        j = i + MinFrequencyindex - 1; % index of the same bin on the full spectrum
        
        if j-k3 < 1 || j+k3 > length(SPL)
            continue
        end
        
        if SPL(j) > SPL(j-1) && ... % first condition
           SPL(j) >= SPL(j+1) && ...
           SPL(j) - SPL(j-k3) >= threshold && ... % second condition
           SPL(j) - SPL(j-k2) >= threshold && ...
           SPL(j) - SPL(j+k2) >= threshold && ...
           SPL(j) - SPL(j+k3) >= threshold
            
           ToneIdx(k) = i; % get the idx of the tones on Lcrop
           k = k+1;
        end
    end
           
    % save tone information
    ToneIdx(ToneIdx==0) = [];   % if no tones were found, ToneIdx shall remain empty
    ToneL = SPLcrop(ToneIdx);   % SPL of the tones
    NTones = find(ToneIdx);     % number of tones
    ToneF = FreqCrop(ToneIdx);  % central freq of the tone
    
    % estimate bandwidth of the i-th tone using half-power (-3 dB decay) criteria (this analysis is made on the full SPL and freq vectors)
    flow=zeros(1,length(NTones));  % declare variable for memory allocation
    fhigh=zeros(1,length(NTones));
    BW=zeros(length(NTones),1);
       
    for i=1:length(NTones) %Source: https://de.mathworks.com/matlabcentral/answers/1441689-i-am-trying-to-find-the-full-width-at-half-max-value-and-plot-the-waveform-with-markers?s_tid=srchtitle
                
        ymx = ToneL(i); % SPL of the i-th tone
        [~,idx] = min( abs(Freq-ToneF(i)) ); % index of the i-th tone 
        hafmax = ymx-3; % half power, three decibels below the peak
        % hafmax = ymx-3; % target value (-3 dB decay)
        
        idxrng1 = find(SPL(1:idx)<hafmax, 1, 'last');
        
        if isempty(idxrng1) || idxrng1<4 % if idxrng1 is empty, it means hafmax is below the 1st bin of the signal (probably due to a low freq tone with large bandwidth)
            idxrng1 = 4; % in this case, truncate idxrng1 to 4
        end

        idxrng2 = find(SPL(idx+1:numel(Freq))<hafmax,1,'first')+idx;
               
        flow(i) = interp1(SPL(idxrng1:idxrng1+1), Freq(idxrng1:idxrng1+1), hafmax);  % low freq of the band
        fhigh(i) = interp1(SPL(idxrng2-1:idxrng2), Freq(idxrng2-1:idxrng2), hafmax); % high freq of the band
               
        BW(i,1) = fhigh(i) - flow(i); % tone's bandwidth
        
        if BW(i,1)==0 % if BW is zero, truncate BW to 1
           BW(i,1)=1;
        end
               
        clear idxrng1 idxrng2 idx
    end
    
    BW( isinf(BW) | isnan(BW) ) = 1;  % replace inf and NaN 
    BWnotch = max(BW, 4*df);  % the notch has to cover the main lobe of the window,
                              % four bins wide for a Hann window. BW itself stays the
                              % measured width, which is what the weighting w1 needs
        
    % Aures (1985), section 2.3.2: besides the sinusoidal components, the
    % spectrum holds regions narrower than a critical band whose critical band
    % stands at least 7 dB above each of the two neighbouring critical bands.
    % Those regions count as tonal components as well, and what is left once
    % they are taken out is the noise spectrum of the model. The geometry
    % follows FindBand_V, appendix A.3.8 of the Purdue thesis of Hastings
    % (2004), the only published operational form of this step.
    %
    % The paper states the criterion and leaves open how it ranks against the
    % sinusoidal extraction, and that has to be settled here: the criterion of
    % Terhardt compares a peak with the bins 25 and 37.5 Hz away, so inside a
    % band narrower than that reach it measures the edges of the band and every
    % local maximum of the noise passes it. A region therefore stands as one
    % component and the sinusoids inside it are dropped, since the region
    % already carries their power; a sinusoid stands where no region covers it.
    % For a pure tone the two paths agree: the region around it has the width
    % of the analysis window, which the weighting below takes out again, and
    % its level is the SPL of the tone.

    nHalf = ceil( (N+1)/2 );  % bins of the single sided spectrum
    [nbF,nbBW,nbL,nbIdx,nbBridge] = il_find_narrowband( SPL(1:nHalf),...
                              Freq(1:nHalf), threshold, MinFrequency, MaxFrequency );

    % the level of a region is a sum of bin powers of a windowed spectrum, so
    % it is divided by the equivalent noise bandwidth of the window to read as
    % the level of the component; a pure tone then reads its own SPL
    nbL = nbL - 10*log10(ENBW);

    for i = 1:numel(nbF)   % a region carries the power of the sinusoids in it
        ToneIdx( ToneF >= Freq(nbIdx{i}(1)) & ToneF <= Freq(nbIdx{i}(end)) ) = 0;
    end
    absorbed = (ToneIdx==0);
    ToneIdx(absorbed) = [];  ToneL(absorbed) = [];  ToneF(absorbed) = [];
    BW(absorbed) = [];  BWnotch(absorbed) = [];  NTones(absorbed) = [];

    if isempty(ToneIdx)==1 && isempty(nbF)  % if ToneRef is empty, then there are no tones for this time-frame
        
        %% OUTPUTS for this case
        
        w_tonal(iFrame,1) = 0;  % Tonal weighting
        w_gr(iFrame,1) = 0;     % loudness weighting  
        tonality(iFrame,1) = 0; % tonality
        
    else     % if tones were found ...
        
        idx = find(ToneL>0);    % find idx of only positive levels (i.e.,
                                %   tones with SPL above 0 dB) - necessary 
                                %   because resampling may introduce several 
                                %   tones with very low amplitude
        ToneIdx = ToneIdx(idx); % idx of the tone
        ToneL = ToneL(idx);     % SPL of the tones
        NTones = NTones(idx);   % number of tones
        ToneF = ToneF(idx);     % central freq of the tone
        BW = BW(idx);           % bandwidth
        BWnotch = BWnotch(idx);
        
        if isempty(ToneIdx)==1 && isempty(nbF)  % if ToneRef is empty (there are no tonal
                                % components with SPL>0 dB), then there are 
                                % no tones for this time-frame
            %% OUTPUTS for this case
            w_tonal(iFrame,1) = 0;  % Tonal weighting
            w_gr(iFrame,1) = 0;     % loudness weighting
            tonality(iFrame,1) = 0; % tonality
            
        else     % if tones were found and their SPL is above 0 dB ...
                                    
            %% filtering out the tones from the signal
            
            y=window.*insig(:,iFrame);  % windowed frame: the notch below removes a
                                        % tone from the same spectrum that found it
            
            insigSpectrum=fft(y);  % spectrum of insig for each iFrames
            
            %%%% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure; semilogy(Freq,abs(insigSpectrum).^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            SingleSidedinsigSpectrum = insigSpectrum(1:ceil((length(insigSpectrum)+1)/2)); % single-sided spectrum of insig for each iFrames
            
            FreqSingleSidedinsigSpectrum=0:fs/length(y):fs/2;  % freq vector of single-sided spectrum of insig for each iFrames
            
            %%%% check plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure; semilogy(FreqSingleSidedinsigSpectrum,abs(SingleSidedinsigSpectrum).^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=1:length(NTones) % loop across tones
                
                index_low = find (FreqSingleSidedinsigSpectrum>=(ToneF(i)-(BWnotch(i)./2)),1,'first'); % find idx of i-th tone's lower freq
                index_up = find (FreqSingleSidedinsigSpectrum>=(ToneF(i)+(BWnotch(i)./2)),1,'first');  % find idx of i-th tone's upper freq
                
                if isempty(index_low) 
                    index_low = 1;
                end
                
                if isempty(index_up)
                    index_up = numel(FreqSingleSidedinsigSpectrum);
                end
                
                if index_low==1 % may happen with low-freq tones with large bandwidth
                    magn=0.5.*(abs(SingleSidedinsigSpectrum(index_low))+abs(SingleSidedinsigSpectrum(index_up+1))); % create a magnitude vector
                else
                    magn=0.5.*(abs(SingleSidedinsigSpectrum(index_low-1))+abs(SingleSidedinsigSpectrum(index_up+1))); % create a magnitude vector
                end
                
                phase = angle( SingleSidedinsigSpectrum(index_low:index_up) ).'; % keep the phase the replaced bins already had
                SingleSidedinsigSpectrum(index_low:index_up) = magn.*exp(1j.*phase); % replace tones
                
            end
            
            %%%% check plot (only for debugging) %%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure; semilogy(FreqSingleSidedinsigSpectrum,abs(SingleSidedinsigSpectrum).^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Footprint of the estimator that measured each width, taken out of
            % it before the width is turned into Bark. The half power estimator
            % of the sinusoids reads 0.997 to 1.43 bins for a pure tone, the
            % quantile estimator of the regions reads 2.082 to 2.667 bins over
            % 300 Hz to 3 kHz and every position between two bins. The larger of
            % the two is taken for each, so that a pure tone always comes out
            % with no width left and w1 stays at one for it.
            
            ToneF = ToneF(:); ToneL = ToneL(:); BW = BW(:);
            Wfoot = repmat( 1.43*df, numel(ToneF), 1 );
            
            % the regions are replaced by the floor the detector measured, so
            % that what is left of the spectrum is the noise of the model
            for i = 1:numel(nbF)
                rr = nbIdx{i};
                magn = sqrt( 10.^(nbBridge{i}(:)./10).*4e-10 )./fftgain;
                phase = angle( SingleSidedinsigSpectrum(rr) );
                SingleSidedinsigSpectrum(rr) = magn.*exp(1j.*phase);
            end
            
            if ~isempty(nbF)
                ToneF = [ToneF; nbF];
                ToneL = [ToneL; nbL];
                BW = [BW; nbBW];
                Wfoot = [Wfoot; repmat( 2.667*df, numel(nbF), 1 )];
                [ToneF,iSort] = sort(ToneF);
                ToneL = ToneL(iSort); BW = BW(iSort); Wfoot = Wfoot(iSort);
            end
            
            % the noise spectrum of the model, which the level excess reads
            FreqNoise = FreqSingleSidedinsigSpectrum(:);
            SPLnoise = 10.*log10( (abs(SingleSidedinsigSpectrum(:)).*fftgain).^2./4e-10 + TINY_VALUE );
            
            doubleSideFilteredSpectrum = [SingleSidedinsigSpectrum; conj(flipud(SingleSidedinsigSpectrum(2:end-1)))]; % double-side the filtered spectrum
            
            %%%% check plot (only for debugging) %%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure; semilogy(Freq,abs(doubleSideFilteredSpectrum).^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            filtered_signal=ifft(doubleSideFilteredSpectrum,'symmetric');  % get filtered signal in time-domain
            
            %% Compute w_gr (loudness weighting)
            
            % compute loudness from input signal 
            % assume a stationary loudness within iFrame
            
            L_total = Loudness_ISO532_1(y, fs,...   % input signal and sampling freq.
                                LoudnessField,...   % field; free field = 0; diffuse field = 1;
                                            1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                         time_resolution*0.05,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                            0);     % show results; 0=no, 1=yes
            
            % compute loudness of the filtered signal (i.e. input signal with tones removed) 
            % assume a stationary loudness within the iFrame
            
            L_filtered = Loudness_ISO532_1(filtered_signal,fs,...   % input signal and sampling freq.
                                                LoudnessField,...   % field; free field = 0; diffuse field = 1;
                                                            1,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                                         time_resolution*0.05,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                                                            0);     % show results; 0=no, 1=yes
            
            % loudness weighting per time frame
            w_gr(iFrame,1)= 1 - ( L_filtered.Loudness/L_total.Loudness );
            
            %	Note: On rare occasions, it is possible for the Loudness of Noise to be greater
            %   than the total Loudness.  This occurs because filtering the tones may slightly
            %	elevate the noise.  If the signal is almost all noise, then this may push it
            %	higher.  If this happens, then the signal should not be considered tonal,
            %	therefore, for this case set Wgr == 0.
            
            if w_gr(iFrame,1)<0
               w_gr(iFrame,1)=0;
            end
            
            clear y insigSpectrum SingleSidedinsigSpectrum
            clear FreqSingleSidedinsigSpectrum doubleSideSpectrum filtered_signal
            
            %% Compute tonal weighting
            
            tone{iFrame,1}.Lcrop = SPLcrop;  %  SPL of the spectrum - SPLcrop = SPL(MinFrequencyindex:MaxFrequencyIndex);
            tone{iFrame,1}.freq = FreqCrop;  %  frequency vector - freq = freq_all(MinFrequencyindex:MaxFrequencyIndex);
            tone{iFrame,1}.ToneF = ToneF;    %  ToneF: central frequency of the tones
            tone{iFrame,1}.ToneL = ToneL;    %  ToneL: SPL of the tones
            tone{iFrame,1}.BW = BW;          %  bandwidth of the tones
            tone{iFrame,1}.df = df;          %  freq discretization
            tone{iFrame,1}.Wfoot = Wfoot;    %  footprint of the width estimator of each component, Hz
            tone{iFrame,1}.Lnoise = interp1( FreqNoise, SPLnoise,...
                                             FreqCrop, 'linear', 'extrap' ); % noise spectrum on the cropped axis
                          
            tone{iFrame,1}.LX=il_SPL_excess(tone{iFrame,1}); %  Sound pressure excess calculation (define aurally relevance of the tones)
                               
            w_tonal(iFrame,1)=il_tonal_weighting(tone{iFrame,1});  % Tonal weighting
            
            %% TONALITY
            
            C=1.1055;  % is a constant such that 1 kHz pure tone with a level of 60 dB would have a tonalness of 1, which for an ideal implementaiton should be =1.09
            
            tonality(iFrame,1) = abs( C.*w_tonal(iFrame,1).^(0.29).*w_gr(iFrame,1).^(0.79) );
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Data

% main output results
OUT.InstantaneousTonality = tonality;  % instantaneous tonality
OUT.TonalWeighting = w_tonal;          % instantaneous tonal weighting
OUT.LoudnessWeighting = w_gr;          % instantaneous loudness weighting
OUT.time = t;                          % time vector

% get statistics from Time-varying tonality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector

metric_statistics = 'Tonality_Aures1985';
OUT_statistics = get_statistics( tonality(idx:end), metric_statistics ); % get statistics

% copy fields of <OUT_statistics> struct into the <OUT> struct
fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

for i = 1:numel(fields_OUT_statistics)
    fieldName = fields_OUT_statistics{i};
    if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
        OUT.(fieldName) = OUT_statistics.(fieldName);
    end
end

clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots

if show == true
    
    figure('NAME','Aures tonality analysis',...
        'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
    %%%
    subplot(3,1,1)
    plot(t,tonality);
    title('Instantaneous tonality','Interpreter','Latex');
    ylabel('Aures tonality, $K$ (t.u.)','Interpreter','Latex');
    xlabel('Time, $t$ (s)','Interpreter','Latex');
    ylim([0 1.1]);
    
    %%%
    subplot(3,1,2)
    plot(t,w_gr,'k');
    title('Loudness weighting','Interpreter','Latex');
    ylabel('Loudness weighting, $W_{\mathrm{Loudness}}$','Interpreter','Latex');
    xlabel('Time, $t$ (s)','Interpreter','Latex');
    ylim([0 1.1]);
    
    %%%
    subplot(3,1,3)
    
    plot(t,w_tonal,'k');
    title('Tonal weighting','Interpreter','Latex');
    ylabel('Tonal weighting, $W_{\mathrm{Tonal}}$','Interpreter','Latex');
    xlabel('Time, $t$ (s)','Interpreter','Latex');
    ylim([0 1.1]);
    
    set(gcf,'color','w')
    
end
end

% End-of-file main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of inline functions:
function LX=il_SPL_excess(input)
% function LX=il_SPL_excess(input)
%
%   INPUT: tone struct containing
%          * tone.freq - freq vector; FreqCrop = Freq(MinFrequencyindex:MaxFrequencyIndex); 
%          * tone.Lcrop - SPL vector; SPLcrop = SPL(MinFrequencyindex:MaxFrequencyIndex);
%          * tone.ToneF - vector containing the central frequency of each tone
%          * tone.ToneL - vector containing SPl of each tone%
%
%   OUTPUT
%          * LX (sound pressure level excess of each tonal component)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main source: https://github.com/densilcabrera/aarae/blob/master/Analysers/Pitch%20and%20Frequency/Terhardt_VirtualPitch.m
% original source: See reference [2], Terhardt et al. (1982)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gil Felix Greco - Braunschweig 10.06.2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Terhardt et al. (1982), eq. (4), asks for a dimensionless intensity summed
% over the spectrum that holds the noise. The sum below runs over the spectrum
% the extraction stages left behind, and it carries no reference pressure: the
% terms it is added to further down are dimensionless as well.
Intensity = 10.^(input.Lnoise./10);

freq_Lx = input.freq;   % freq vector of the tone
ToneF = input.ToneF;    % tone(s) central frequency
ToneL = input.ToneL;    % tone(s) level
NTones = size(ToneF,1); % number of tones

toneBark = il_Fq2Bark(ToneF);       % convert central freq of tones to Bark scale
spectrumBark = il_Fq2Bark(freq_Lx); % convert freq vector to Bark scale

LX = deal(zeros(NTones,1)); % initialize sound pressure level excess vector

for i = 1:NTones

    % Intensity of noise for each tone paragraph after eq 7b in Ref. [3] (Terhard's papers)

    idx_cb = spectrumBark >= round( toneBark(i)-0.5 )...
           & spectrumBark <= round( toneBark(i)+0.5 ); % idx of the critical band around the tonal component

    idx_toneBark = find( round( spectrumBark==toneBark(i) )); % find idx of the tone on the Bark vector

    idx_cb(idx_toneBark-2:idx_toneBark+2) = 0; % skip the five central samples around the tonal component

    EGR = sum( Intensity(idx_cb) ); % Masking intensity of broadband noise

    % Secondary excitation level
    sumlo = 1e-99;
    sumhi = 1e-99;

    for j = 1:NTones

        if (j < i)

            s = -24 - (230./(ToneF(j))) + (0.2.*ToneL(j)); % eq 7b from Ref. [3]
            Lji = ToneL(j) - s .* (toneBark(j) - toneBark(i));
            sumlo = sumlo + 10.^(Lji./20);

        elseif (j > i)

            s=27;
            Lji = ToneL(j) - s .* (toneBark(j) - toneBark(i));
            sumhi = sumhi + 10.^(Lji./20);

        end

    end

    AEK = sumlo + sumhi;

    % Intensity at threshold of hearing
    EHS = il_Threshold(ToneF(i));
    EHS = 10.^(EHS/10);

    % Sound pressure level excess - NOTE: in the original paper from Terhard [3]
    % -10log10 is used while in the paper of Aures [1] simply -log10 is used

    if NTones==1 % if there is only one tone
        LXi = ToneL(i) - 10.*log10( EGR  + EHS ); %eq 4 from Ref. [3]
    else
        LXi = ToneL(i) - 10.*log10( AEK.^2 + EGR  + EHS ); %eq 4 from Ref. [3]
    end

    if LXi > 0
        LX(i) = LXi;
    end

end
end % end il_SPL_excess

function [w_tonal]=il_tonal_weighting(input)

bw=input.BW;      % bandwidth of the tones [Hz]
fc=input.ToneF;   % central frequency of the tonal components
delta_L=input.LX; % SPL excess for each tonal component

%% w1 accounts for each tonal component bandwidth

% The measured width carries the width of the analysis window. For two smooth
% kernels the widths add in quadrature, so the window is removed the same way.
% The width to remove is the one the estimator reads on a pure tone, which
% depends on where the tone falls between two bins and reaches 1.43 bins for a
% Hann window. Taking that worst case makes any pure tone give dz = 0.
Wfoot = input.Wfoot(:);
Wtrue = sqrt( max(bw(:).^2 - Wfoot.^2, 0) );
zup   = il_Fq2Bark(fc+(Wtrue./2));
zlow  = il_Fq2Bark(fc-(Wtrue./2));
dz    = zup - zlow;

w1 = ( 0.13./(dz+0.13) );

%% w2 accounts for each tonal component's center frequency

w2 =  ( 1./( sqrt (1+0.2.*(fc./700 + 700./fc).^2) ) ).^(0.29);

%% w3 accounts for each tonal component SPL excess

w3 =( 1-exp(-delta_L/15) ).^(0.29);

%% prime weightings

ww1 = w1.^(1./0.29);
ww2 = w2.^(1./0.29);
ww3 = w3.^(1./0.29);

%% total tonal weighting

w_tonal= sqrt(sum( (ww1 .* ww2 .* ww3).^2 ) );

end % End il_tonal_weightin

%% function: find the narrow band components of a spectrum

function [nbF,nbBW,nbL,nbIdx,nbBridge] = il_find_narrowband(SPL,f,threshold,fmin,fmax)
% function [nbF,nbBW,nbL,nbIdx,nbBridge] = il_find_narrowband(SPL,f,threshold,fmin,fmax)
%
%   Narrow band components of the residual spectrum. Aures (1985b), section
%   2.3.2, asks for the residue to be searched for regions narrower than a
%   critical band whose neighbouring critical bands are at least 7 dB lower.
%   The paper states the criterion and leaves the geometry open. The geometry
%   used here follows FindBand_V, appendix A.3.8 of the Purdue thesis of
%   Hastings (2004), which is the only published operational form: a band sum
%   over three bins, a half power walk that sets the width and the centre, and
%   an extension out to a noise floor taken as the level exceeded by 90 % of
%   the half critical band beyond the half power point.
%
%   INPUT
%     SPL       : [Nx1] sound pressure level of the residual spectrum, dB
%     f         : [Nx1] frequency vector, Hz
%     threshold : level a region must have over each neighbouring critical
%                 band to count as a component, dB
%     fmin,fmax : frequency range of the search, Hz
%
%   OUTPUT
%     nbF      : centre frequency of each component, Hz, geometric mean of
%                the two half power points
%     nbBW     : half power width of each component, Hz
%     nbL      : level of each component over the floor, dB
%     nbIdx    : index range each component occupies
%     nbBridge : level of the floor over that range, dB

SPL = SPL(:); f = f(:);
nBins = numel(SPL);
df    = f(2)-f(1);

MinPeakLevel = 0;     % dB, the search stops below this peak level
MaxToneCount = 20;    % the default of 1 of the reference is not usable here
Fraction     = 0.5;   % of a critical band, the span of the noise floor estimate
PercentBelow = 0.1;   % the noise floor is the level exceeded by 90 % of it
BWt          = 3;     % bins of the smoothing band sum, must be odd

side = (BWt-1)/2;

nbF = zeros(MaxToneCount,1); nbBW = nbF; nbL = nbF;
nbIdx = cell(MaxToneCount,1); nbBridge = cell(MaxToneCount,1);
nFound = 0;

% band sum over BWt bins, expressed as a mean per bin (step 1 of the reference)
Yav = SPL;
for k = (1+side):(nBins-side)
    Yav(k) = 10*log10( sum(10.^(SPL(k-side:k+side)./10)) ) - 10*log10(BWt);
end

% -100 marks a bin the search must not enter, either out of range or already taken
Ysearch = Yav;
Ysearch( f < fmin | f > fmax ) = -100;

Ywork = SPL;   % the residue, which the floors below are written into
zAll  = il_Fq2Bark(f);

for iComp = 1:MaxToneCount

    [Ymax,PeakIndex] = max(Ysearch);
    if Ymax < MinPeakLevel
        break
    end

    % half power point on each side, or the edge of an identified region
    LToneFlag = 0; RToneFlag = 0;
    ink = 0;
    while true
        ink = ink+1;
        if PeakIndex-ink <= 0
            LeftPowerIndex = 1; break
        elseif Yav(PeakIndex-ink)+3 < Ymax
            LeftPowerIndex = PeakIndex-ink; break
        elseif Ysearch(PeakIndex-ink) == -100
            LeftPowerIndex = PeakIndex-ink; LToneFlag = 1; break
        end
    end
    ink = 0;
    while true
        ink = ink+1;
        if PeakIndex+ink > nBins
            RightPowerIndex = nBins; break
        elseif Yav(PeakIndex+ink)+3 < Ymax
            RightPowerIndex = PeakIndex+ink; break
        elseif Ysearch(PeakIndex+ink) == -100
            RightPowerIndex = PeakIndex+ink; RToneFlag = 1; break
        end
    end

    % noise floor on the left, over a fraction of a critical band
    if f(LeftPowerIndex) < 500
        StartLeft = LeftPowerIndex - round(Fraction*100/df);
    else
        StartLeft = LeftPowerIndex - round(Fraction*0.2*f(LeftPowerIndex)/df);
    end
    StartLeft = max(StartLeft,1);

    if LToneFlag == 1
        LeftRegion = LeftPowerIndex;
    else
        tmp = sort( Yav(StartLeft:LeftPowerIndex) );
        ti  = max( floor(PercentBelow*numel(tmp))-1, 1 );
        LeftNoiseFloor = tmp(ti);
        LeftRegion = StartLeft;
        for k = LeftPowerIndex:-1:StartLeft
            if Ysearch(k) == -100, LeftRegion = k+1; break; end
            if Yav(k) <= LeftNoiseFloor
                LeftRegion = k; break
            end
        end
    end

    % noise floor on the right, same construction
    if f(RightPowerIndex) < 500
        StopRight = RightPowerIndex + round(Fraction*100/df);
    else
        StopRight = RightPowerIndex + round(Fraction*0.2*f(RightPowerIndex)/df);
    end
    StopRight = min(StopRight,nBins);

    if RToneFlag == 1
        RightRegion = RightPowerIndex;
    else
        tmp = sort( Yav(RightPowerIndex:StopRight) );
        ti  = max( floor(PercentBelow*numel(tmp))-1, 1 );
        RightNoiseFloor = tmp(ti);
        RightRegion = StopRight;
        for k = RightPowerIndex:StopRight
            if Ysearch(k) == -100, RightRegion = k-1; break; end
            if Yav(k) <= RightNoiseFloor
                RightRegion = k; break
            end
        end
    end

    % the region is taken out of the search whether or not it turns out tonal
    Ysearch(LeftRegion:RightRegion) = -100;

    idx    = (LeftRegion:RightRegion).';
    span   = max(RightRegion-LeftRegion,1);
    bridge = Ywork(LeftRegion) + ...
             ( Ywork(RightRegion)-Ywork(LeftRegion) ).*(idx-LeftRegion)./span;

    fc  = sqrt( f(LeftPowerIndex)*f(RightPowerIndex) );
    bwp = f(RightPowerIndex) - f(LeftPowerIndex);
    CBW = 25 + 75*(1+1.4*(fc/1000)^2)^0.69;

    if fc <= 0 || bwp <= 0 || bwp >= CBW
        continue   % a region as wide as a critical band is noise
    end

    % the two conditions of Aures, the region against its neighbouring bands
    zc   = il_Fq2Bark(fc);
    own  = zAll >= zc-0.5  & zAll <= zc+0.5;
    low  = zAll >= zc-1.5  & zAll <  zc-0.5;
    upp  = zAll >  zc+0.5  & zAll <= zc+1.5;
    if ~any(low) || ~any(upp)
        continue
    end
    Lown = 10*log10( sum(10.^(Ywork(own)./10)) );
    Llow = 10*log10( sum(10.^(Ywork(low)./10)) );
    Lupp = 10*log10( sum(10.^(Ywork(upp)./10)) );
    if ~( Lown-Llow >= threshold && Lown-Lupp >= threshold )
        continue
    end

    % level of the component over the floor, taken over the half power band
    inHalf = idx >= LeftPowerIndex & idx <= RightPowerIndex;
    Pold   = sum( 10.^(Ywork(LeftPowerIndex:RightPowerIndex)./10) );
    Pnew   = sum( 10.^(bridge(inHalf)./10) );
    if Pold <= Pnew
        continue
    end

    % Width of the component. The half power walk of the reference measures the
    % width of a fluctuation peak once the component is noise, and saturates:
    % ideal bands of 33, 66 and 133 Hz all read 11 Hz at a resolution of 1 Hz.
    % The width used here is the span of the middle 90 % of the power of the
    % region, divided by 0.9 so that a rectangular band reads its own width. It
    % reads those same three bands within 12 %, and it does not move when the
    % component sits in broadband noise, where the moments of the distribution
    % are carried by the tails. Aures states the criterion for these components
    % and leaves the measurement of their width open, so this is a choice made
    % here and not something the paper prescribes.
    Preg = 10.^(Ywork(idx)./10) - 10.^(bridge./10);
    Preg(Preg<0) = 0;
    if sum(Preg) <= 0
        continue
    end
    freg = f(idx);
    cw   = cumsum(Preg)./sum(Preg) + (1:numel(Preg)).'.*1e-12;  % strictly rising
    flo90 = interp1(cw,freg,0.05,'linear','extrap');
    fhi90 = interp1(cw,freg,0.95,'linear','extrap');
    wrms  = (fhi90-flo90)./0.9;

    Ywork(LeftRegion:RightRegion) = bridge;   % the residue keeps only the floor

    nFound = nFound+1;
    nbF(nFound)      = fc;
    nbBW(nFound)     = wrms;
    nbL(nFound)      = 10*log10( Pold-Pnew );
    nbIdx{nFound}    = idx;
    nbBridge{nFound} = bridge;

end

nbF = nbF(1:nFound); nbBW = nbBW(1:nFound); nbL = nbL(1:nFound);
nbIdx = nbIdx(1:nFound); nbBridge = nbBridge(1:nFound);

end

%% function: convert frequency to bark

function B = il_Fq2Bark(f)

% critical band rate corresponding to a given frequency
% input f is frequency in Hz
% output B is critical band rate in Barks

f=f./1000;
B = 13 .* atan(0.76 .* f) + 3.5 .* atan ((f./7.5).^2);

end % end il_Fq2Bark

%% function: hearing threshold

function L = il_Threshold(f)

    % hearing threshold
    % input f is frequency in Hz
    % output L is threshold in dB

    f=f/1000;
    L = 3.64 * f.^-0.8 ...
        - 6.5 * exp(-0.6 * (f - 3.3).^2) ...
        + 1e-3 * f.^4;
end % end il_Threshold

%**************************************************************************
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
%    may be used to endorse or promote products derived from this software 
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%**************************************************************************
