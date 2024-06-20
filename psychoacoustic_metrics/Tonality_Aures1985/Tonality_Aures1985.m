function OUT = Tonality_Aures1985(insig, fs, LoudnessField, start_skip, end_skip, show)
% function OUT = Tonality_Aures1985(insig, fs, LoudnessField, start_skip, end_skip, show)
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

%   start_skip : number
%   end_skip : number
%   skip start/end of the signal in seconds for statistic calculations
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
%         ** Kx : percentile InstantaneousTonality exceeded during x percent of the signal (t.u.)
%
% Author: Gil Felix Greco, Braunschweig 13/07/2020 (updated 14.04.2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
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
time_resolution=160e-3; % window length fixed in 160 ms, gives a df=6.25 Hz

N=round(fs*time_resolution); % define window length, N bins
window = hann(N);

fftgain = 2^0.5/(N*mean(hann(N))); % gain to be applied based on the FFT length

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
    
    ToneIdx = zeros(length(SPLcrop),1); % initialize vector, tonal components idx
    k = 1; % initialize counter
    
    % find tones...
    for i = 4:(length(SPLcrop)-3)
        
        if SPLcrop(i) > SPLcrop(i-1) && ... % first condition
           SPLcrop(i) >= SPLcrop(i+1) && ...
           SPLcrop(i) - SPLcrop(i-3) >= threshold && ... % second condition
           SPLcrop(i) - SPLcrop(i-2) >= threshold && ...
           SPLcrop(i) - SPLcrop(i+2) >= threshold && ...
           SPLcrop(i) - SPLcrop(i+3) >= threshold
            
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
        hafmax = ymx.*0.707; % target value 
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
        
    if isempty(ToneIdx)==1  % if ToneRef is empty, then there are no tones for this time-frame
        
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
        
        if isempty(ToneIdx)==1  % if ToneRef is empty (there are no tonal
                                % components with SPL>0 dB), then there are 
                                % no tones for this time-frame
            %% OUTPUTS for this case
            w_tonal(iFrame,1) = 0;  % Tonal weighting
            w_gr(iFrame,1) = 0;     % loudness weighting
            tonality(iFrame,1) = 0; % tonality
            
        else     % if tones were found and their SPL is above 0 dB ...
                                    
            %% filtering out the tones from the signal
            
            y=insig(:,iFrame);     % get insig for each iFrames
            
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
                
                index_low = find (FreqSingleSidedinsigSpectrum>=(ToneF(i)-(BW(i)./2)),1,'first'); % find idx of i-th tone's lower freq
                index_up = find (FreqSingleSidedinsigSpectrum>=(ToneF(i)+(BW(i)./2)),1,'first');  % find idx of i-th tone's upper freq
                
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
                
                phase = (rand(1,index_up-index_low+1)-0.5).*pi.*2; % create random phase vector
                SingleSidedinsigSpectrum(index_low:index_up) = magn.*exp(1j.*phase); % replace tones
                
            end
            
            %%%% check plot (only for debugging) %%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure; semilogy(FreqSingleSidedinsigSpectrum,abs(SingleSidedinsigSpectrum).^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
                          
            tone{iFrame,1}.LX=il_SPL_excess(tone{iFrame,1}); %  Sound pressure excess calculation (define aurally relevance of the tones)
                               
            w_tonal(iFrame,1)=il_tonal_weighting(tone{iFrame,1});  % Tonal weighting
            
            %% TONALITY
            
            C=1.125;  % is a constant such that 1 kHz pure tone with a level of 60 dB would have a tonalness of 1, which for an ideal implementaiton should be =1.09
            
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tonality statistics

[~, idx] = min( abs(OUT.time - start_skip) ); % find idx of start_skip on time vector
[~, idxEnd] = min( abs(OUT.time - (max(OUT.time) - end_skip)) ); % find idx of end_skip on time vector

OUT.Kmean = mean(tonality(idx:idxEnd));
OUT.Kstd = std(tonality(idx:idxEnd));
OUT.Kmax = max(tonality(idx:idxEnd));
OUT.Kmin = min(tonality(idx:idxEnd));
OUT.K1 = get_percentile(tonality(idx:idxEnd),1);
OUT.K2 = get_percentile(tonality(idx:idxEnd),2);
OUT.K3 = get_percentile(tonality(idx:idxEnd),3);
OUT.K4 = get_percentile(tonality(idx:idxEnd),4);
OUT.K5 = get_percentile(tonality(idx:idxEnd),5);
OUT.K10 = get_percentile(tonality(idx:idxEnd),10);
OUT.K20 = get_percentile(tonality(idx:idxEnd),20);
OUT.K30 = get_percentile(tonality(idx:idxEnd),30);
OUT.K40 = get_percentile(tonality(idx:idxEnd),40);
OUT.K50 = median(tonality(idx:idxEnd));
OUT.K60 = get_percentile(tonality(idx:idxEnd),60);
OUT.K70 = get_percentile(tonality(idx:idxEnd),70);
OUT.K80 = get_percentile(tonality(idx:idxEnd),80);
OUT.K90 = get_percentile(tonality(idx:idxEnd),90);
OUT.K95 = get_percentile(tonality(idx:idxEnd),95);

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

pref = 2e-5; % reference pressure, Pa
Intensity = pref.*10.^(input.Lcrop./10);

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

    NTonesM = 0;
    if LXi > 0
        NTonesM = NTonesM + 1;
        LX(NTonesM) = LXi;
    end

end
end % end il_SPL_excess

function [w_tonal]=il_tonal_weighting(input)

bw=input.BW;      % bandwidth of the tones [Hz]
fc=input.ToneF;   % central frequency of the tonal components
delta_L=input.LX; % SPL excess for each tonal component
df=input.df;      % freq discretization

%% w1 accounts for each tonal component bandwidth

zup = il_Fq2Bark(fc+(bw./2));
zlow = il_Fq2Bark(fc-(bw./2));
dz = (zup-zlow)/df^2;

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
