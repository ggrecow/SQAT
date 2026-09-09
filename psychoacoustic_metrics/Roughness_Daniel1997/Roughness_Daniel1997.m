function OUT = Roughness_Daniel1997(insig,fs,time_skip,show)
% function OUT = Roughness_Daniel1997(insig,fs,time_skip,show)
%
%   This function calculates time-varying roughness and time-averaged specific
%     roughness using the roughness model by Daniel & Weber:
%     Daniel, P., & Weber, R. (1997). Psychoacoustical roughness: implementation
%     of an optimized model. Acustica(83), 113-123.
%
%   Reference signal: 1 kHz tone, 100% amplitude modulated at 70 Hz, with a
%   sound pressure level of 60 dB (rms of the modulated signal), yields
%   1 asper (see the calibration below).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   insig : array [Nx1]
%   acoustic signal, monophonic (Pa)
%
%   fs : integer
%   sampling frequency (Hz)
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUT:
%   OUT : struct containing the following fields
%
%       * InstantaneousRoughness: instantaneous roughness (asper) as a 
%         function of time
%       * InstantaneousSpecificRoughness: specific roughness(asper/Bark) as
%         a function of time and frequency (Bark scale)
%       * TimeAveragedSpecificRoughness: time-averaged specific roughness 
%         (asper/Bark) as a function of frequency (Bark scale)
%       * barkAxis : vector of Bark band numbers used for the computation
%         of specific roughness computation
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousRoughness
%         ** Rmean : mean value of instantaneous roughness (asper)
%         ** Rstd : standard deviation of instantaneous roughness (asper)
%         ** Rmax : maximum of instantaneous roughness (asper)
%         ** Rmin : minimum of instantaneous roughness (asper)
%         ** Rx : roughness value exceeded during x percent of the time (asper)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure of the implementation (see also FluctuationStrength_Osses2016,
% which Alejandro Osses derived from this model):
%   1. Analysis windows of 200 ms (Blackman) with 50 % overlap, at 48 kHz
%      (the input is resampled unless it is sampled at 40960, 44100 or
%      48000 Hz).
%   2. Spectrum of each window, weighted by the a0 transmission factor of
%      the outer and middle ear (utilities/calculate_a0.m, curve of Fastl &
%      Zwicker 2007), with the level conversion described below.
%   3. Excitation patterns of the 47 half-Bark channels following Terhardt
%      (1979): utilities/Terhardt_filterbank.m, shared with the fluctuation
%      strength, with the parameters of utilities/Terhardt_filterbank_params.m.
%   4. Modulation depth of the temporal envelope of each channel, weighted
%      by the modulation transfer functions H(fmod) of Daniel & Weber
%      (private/Get_Hweight_roughness.m).
%   5. Cross-correlation of the envelopes of channels 1 Bark apart, and
%      specific roughness r_i = (g(z_i)^0.5 * m_i * k_i)^2 with the g(z)
%      weighting of private/Get_gzi_roughness.m.
%   6. Total roughness R = cal * sum(r_i), Daniel & Weber Eq. 9, and the
%      statistics of get_statistics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code history
%
% The first MATLAB implementation of the Daniel & Weber model was written by
% J. E. Schrader at TU Eindhoven, supervised by Dik Hermes:
%   Schrader, J. E. (2002). A MATLAB implementation of a model of auditory
%   roughness. Technische Universiteit Eindhoven. Available at 
%   (last viewed September 2026):
%   https://research.tue.nl/en/publications/a-matlab-implementation-of-a-model-of-auditory-roughness/
%   Hermes distributed and maintained that code, from which the versions
%   available in PsySound3 and AARAE derive.
%
% - SQAT v1.x: <roughnessDW.m> obtained from
%   https://github.com/densilcabrera/aarae/ (accessed 11.02.2020), where it
%   is attributed to Dik Hermes. Adapted and verified for SQAT by Gil Felix
%   Greco (2023); scaling of the specific roughness corrected by Alejandro
%   Osses (10.05.2023); <get_statistics> introduced by Gil Felix Greco
%   (16.02.2025).
%
% - SQAT v2.0: implementation restructured by Sergio Aguirre and 
%   Gil Felix Greco (September 2026, see issue #47) to follow the structure
%   of <FluctuationStrength_Osses2016.m>, so that shared components are 
%   reused from <utilities>. The corrections deployed were taken from a 
%   revised implementation that Dik Hermes wrote
%   for his book "The Perceptual Structure of Sound" (Springer, 2023, Chap. 6,
%   Sec. 2) and sent to the SQAT team in 2025. Results change with respect 
%   to SQAT v1.x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file and code is subject to the MIT license in
% its entirety, as detailed in the license text reproduced at the end of
% this file.
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help Roughness_Daniel1997;
    return;
end

if nargin < 4
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

%% window settings

time_resolution = 0.2;    % time-step for the windowing (s)

audio = insig;
if size(audio,2)~=1 % if the insig is not a [Nx1] array
    audio = audio';   % correct the dimension of the insig
end

%% resampling input signal

if ~(fs == 44100 || fs == 40960 || fs == 48000)
    gcd_fs = gcd(48000,fs); % greatest common denominator
    audio = resample(audio,48000/gcd_fs,fs/gcd_fs);
    fs = 48000;
end

N = round(fs*time_resolution); % window length, in bins
hopsize = N/2;                 % number of samples hop between successive windows
window = blackman(N,'periodic');
samples = size(audio,1);
n = floor((samples-N)/hopsize); % number of analysis windows

%% model parameters

% frequency axis, Bark scale, hearing threshold: shared with the
% <FluctuationStrength_Osses2016.m> (utilities)
params = Terhardt_filterbank_params(N,fs);
Chno   = params.Chno;  % number of critical-band channels, half-Bark spacing

% a0 transmission factor of the outer and middle ear (Fastl & Zwicker
% 2007), applied on the spectrum of each window as in the reference code
% of Hermes; utilities/calculate_a0.m returns the curve over the audible
% bins params.qb
[~, ~, a0_qb] = calculate_a0(fs,N,'fastl2007');
a0 = ones(1,N);
a0(params.qb) = a0_qb;

Hweight = Get_Hweight_roughness(N,fs);  % modulation weighting functions (private)
gzi     = Get_gzi_roughness(Chno);      % g(z) weighting, Hermes 2025 table (private)

dz = 0.5;            % Bark, integration step of the specific roughness
zi = (1:Chno)'/2;    % Bark axis of the specific roughness

% Level conversion: the spectrum of each window must carry the sound
% pressure level of every component in dB SPL, since the hearing threshold
% and the upper excitation slope of Terhardt depend on it. With the
% magnitude of the windowed spectrum scaled by 2/(N*mean(window)), a
% sinusoid of amplitude A lands at A, so the constant that turns a
% component of L dB SPL (rms) into a magnitude of 10^(L/20) is
% 20*log10(1/20e-6) - 20*log10(sqrt(2)) = 90.97 dB (reference pressure and
% peak to rms of a sinusoid). This is the conversion used by Hermes; the
% previous value of 91.2 dB carried an empirical adjustment of 0.23 dB.
L_cal  = 20*log10(1/20e-6) - 20*log10(sqrt(2));
AmpCal = 10^(L_cal/20)*2/(N*mean(window));

% Calibration of the asper scale: Daniel & Weber (1997), Eq. 9, R =
% cal*sum(r_i), with cal = 0.25 chosen so that the reference signal (1 kHz
% tone, 100 % amplitude modulated at 70 Hz, 60 dB SPL as the rms of the
% modulated signal) gives 1 asper. With the corrected excitation slope,
% the revised g(z) table and the level conversion above, cal = 0.25 gives
% R_ref asper on the reference signal of the toolbox
% (sound_files/reference_signals/RefSignal_Roughness_Daniel1997.wav), so
% cal is re-derived as 0.25/R_ref, keeping the definition of the paper.
% The specific roughness is scaled by Cal = cal/dz, so that its integral
% over the Bark axis gives R.
R_ref = 1.006197; % asper, measured on the reference signal with cal = 0.25 (MATLAB R2026a, 48 kHz)
cal   = 0.25/R_ref;
Cal   = cal/dz;

%% process window

startIndex = 1;
endIndex = N;
[TimePoints,R_mat] = deal(zeros(n,1));
ri_mat = zeros(Chno,n);
clampFrames = 0; clampLdB = -Inf; clampFreq = NaN; % windows with a clamped Terhardt slope

for windowNum = 1:n  % for each window

    dataIn = audio(startIndex:endIndex,1).*window;
    currentTimePoint = startIndex/fs;

    % 1. Spectrum of the window with the level conversion and the a0
    %    transmission factor
    FreqIn = a0.*fft( transpose(dataIn*AmpCal) );

    % 2. Excitation patterns of the critical-band channels (Terhardt),
    %    shared with the fluctuation strength. <info> is requested, so the
    %    filterbank does not warn per window; one warning per call is
    %    raised below
    [ei, info] = Terhardt_filterbank(FreqIn, params);
    if info.clamp.n > 0
        clampFrames = clampFrames + 1;
        if info.clamp.LdB > clampLdB
            clampLdB  = info.clamp.LdB;
            clampFreq = info.clamp.freq;
        end
    end

    % 3. Modulation depth of the temporal envelope of each channel
    [mdept,hBPi] = il_modulation_depths(ei,Hweight);

    % 4. Cross-correlation coefficients between channels 1 Bark apart
    ki = il_cross_correlation(hBPi);

    % 5. Specific roughness and total roughness
    ri = il_specific_roughness(mdept,ki,gzi,Cal,Chno);
    R  = dz*sum(ri);  % total R = integration of the specific R pattern

    % matrices to return
    R_mat(windowNum) = R;
    ri_mat(1:Chno,windowNum) = ri;
    TimePoints(windowNum,1) = currentTimePoint;

    startIndex = startIndex+hopsize;
    endIndex = endIndex+hopsize;

end

% One warning per call when the Terhardt upper slope was clamped to zero
% in any window (component level above 120 + 1150/f dB, about 121 dB at
% 1 kHz), as in FluctuationStrength_Osses2016
if clampFrames > 0
    warning('SQAT:Roughness:TerhardtSlopeClamped', ...
        ['Terhardt upper slope clamped to zero in %d of %d window(s): at least one ' ...
         'component exceeds 120 + 1150/f dB (highest: %.1f dB at %.0f Hz). The ' ...
         'roughness of those windows is an extrapolation outside the range over ' ...
         'which the model was validated; check the dBFS calibration.'], ...
        clampFrames, n, clampLdB, clampFreq);
end

%% ************************************************************************
% output struct
% *************************************************************************

% main output results
OUT.InstantaneousRoughness = R_mat;                       % instantaneous roughness
OUT.InstantaneousSpecificRoughness = ri_mat;              % time-varying specific roughness
OUT.TimeAveragedSpecificRoughness = mean(ri_mat,2);       % mean specific roughness
OUT.time = TimePoints;                                    % time
OUT.barkAxis = zi;                                        % critical band rate (for specific roughness)
OUT.dz = dz;

% Roughness statistics based on InstantaneousRoughness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector

metric_statistics = 'Roughness_Daniel1997';
OUT_statistics = get_statistics( R_mat(idx:end), metric_statistics ); % get statistics

% copy fields of <OUT_statistics> struct into the <OUT> struct
fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

for i = 1:numel(fields_OUT_statistics)
    fieldName = fields_OUT_statistics{i};
    if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
        OUT.(fieldName) = OUT_statistics.(fieldName);
    end
end

clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plots

if show == true
    
    figure('name','Roughness analysis',...
        'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
    
    % Time-varying roughness
    subplot(2,2,1:2)
    
    plot(TimePoints,R_mat,'r-');
    
    title ('Instantaneous roughness','Interpreter','Latex');
    xlabel('Time (s)','Interpreter','Latex');
    ylabel('Roughness, $R$ (asper)','Interpreter','Latex');
    
    % Time-averaged roughness as a function of critical band
    subplot(2,2,3)
    
    plot(zi, mean(ri_mat,2),'r-');
    
    title('Time-averaged specific roughness','Interpreter','Latex');
    xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');
    ylabel('Specific roughness, $R^{\prime}$ (asper/Bark)','Interpreter','Latex');
    
    % Specific roughness spectrogram
    subplot(2,2,4)
    
    [xx,yy]=meshgrid(TimePoints,OUT.barkAxis);
    pcolor(xx,yy,ri_mat);
    shading interp; colorbar; axis tight;
    
    set(gca,'YDir','normal');
    title('Instantaneous specific roughness','Interpreter','Latex');
    xlabel('Time (s)','Interpreter','Latex');
    ylabel('Critical band, $z$ (Bark)','Interpreter','Latex');
    ylabel(colorbar, 'Specific roughness, $R^{\prime}$ ($\mathrm{asper}/\mathrm{Bark}$)','Interpreter','Latex');
    
    set(gcf,'color','w')
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdept,hBPi] = il_modulation_depths(ei,Hweight)
% function [mdept,hBPi] = il_modulation_depths(ei,Hweight)
%
% Modulation depth of the temporal envelope of each channel: the envelope
% is the full-wave rectified excitation, its fluctuation is band-passed by
% the modulation weighting Hweight in the frequency domain, and the depth
% is the rms of the band-passed fluctuation over the DC component h0.

[Chno,N] = size(ei);

% channel by channel, in the operation order of the original code: the
% correlation coefficient of channels whose envelope sits at the rounding
% floor of the FFT depends on that rounding pattern, so a batched FFT
% along the channel dimension moves the roughness of a few signals by up
% to 5e-4 asper. The loop keeps the results identical to the previous
% implementation.
hBPi   = zeros(Chno,N);
h0     = zeros(1,Chno);
hBPrms = zeros(1,Chno);
mdept  = zeros(1,Chno);

for k = 1:Chno
    etmp      = abs(ei(k,:));
    h0(k)     = mean(etmp);
    Fei       = fft(etmp-h0(k));
    hBPi(k,:) = 2*real(ifft(Fei.*Hweight(k,:)));
    hBPrms(k) = rms(hBPi(k,:));
    if h0(k) > 0
        mdept(k) = min(hBPrms(k)/h0(k), 1);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ki = il_cross_correlation(hBPi)
% function ki = il_cross_correlation(hBPi)
%
% Pearson correlation coefficient between the band-passed envelopes of
% channels 1 Bark (two channels) apart.

Chno = size(hBPi,1);
ki = zeros(1,Chno-2);

for k=1:1:Chno-2
    cfac = cov(hBPi(k,:),hBPi(k+2,:));
    den  = diag(cfac);
    den  = sqrt(den*den');
    if den(2,1)>0
        ki(k) = cfac(2,1)/den(2,1);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ri = il_specific_roughness(mdept,ki,gzi,Cal,Chno)
% function ri = il_specific_roughness(mdept,ki,gzi,Cal,Chno)
%
% Specific roughness of each channel, Daniel & Weber (1997), Eq. 12. The
% g(z) weighting enters linearly: gzi carries the square root of the
% tabulated g(z), so squaring the product applies the table once.

ri = zeros(1,Chno);

ri(1) = (gzi(1)*mdept(1)*ki(1))^2;
ri(2) = (gzi(2)*mdept(2)*ki(2))^2;

for k = 3:1:Chno-2
    ri(k) = (gzi(k)*mdept(k)*ki(k-2)*ki(k))^2;
end

ri(Chno-1) = (gzi(Chno-1)*mdept(Chno-1)*ki(Chno-3))^2;
ri(Chno)   = (gzi(Chno)*mdept(Chno)*ki(Chno-2))^2;

ri = Cal*ri; % appropriately scaled specific roughness

end

%**************************************************************************
%
% copyright © 1999-2023 Dik Hermes
% 
% This software is licenced under the MIT license:
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the “Software”), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, a
% nd/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%
%**************************************************************************