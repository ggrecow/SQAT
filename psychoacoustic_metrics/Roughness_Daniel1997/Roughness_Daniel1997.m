function OUT = Roughness_Daniel1997(insig, fs, start_skip, end_skip, show)
% function OUT = Roughness_Daniel1997(insig, fs, start_skip, end_skip, show)
%
%   This function calculates time-varying roughness and time-averaged specific
%     roughness using the roughness model by Daniel & Weber:
%     Daniel, P., & Weber, R. (1997). Psychoacoustical roughness: implementation
%     of an optimized model. Acustica(83), 113-123.
%
%   Reference signal: 60 dB 1 kHz tone 100% modulated at 70 Hz should yield 1 asper.
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
%   start_skip : number
%   end_skip : number
%   skip start/end of the signal in seconds for statistics calculations
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
%         ** Rx : percentile roughness exceeded during x percent of the signal (asper)
%
% Original file name: roughnessDW.m obtained from 
%   https://github.com/densilcabrera/aarae/ (accessed 11/02/2020)
% Author: Dik Hermes (2000-2005)
% Author: Matt Flax (2006) and Farhan Rizwi (2007), adapted for the PsySound3 toolbox
% Author: Gil Felix Greco (2023). Adapted (and verified) for SQAT. 
% Author: Alejandro Osses, 10/05/2023. Appropriate scaling for the specific roughness.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help Roughness_Daniel1997;
    return;
end

if nargin < 5
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end
%% window settings

time_resolution=0.2;  % time-step for the windowing
N=fs*time_resolution; % window length

audio=insig;

if size(audio,2)~=1 % if the insig is not a [Nx1] array
    audio=audio';   % correct the dimension of the insig
end

hopsize=N/2; % hopsize is the number of samples hop between successive windows (window length is 8192).
window = blackman(N);

%% resampling input signal

if ~(fs == 44100 || fs == 40960 || fs == 48000)
    gcd_fs = gcd(48000,fs); % greatest common denominator
    audio = resample(audio,48000/gcd_fs,fs/gcd_fs);
    fs = 48000;
end

samples = size(audio,1);
n = floor((samples-N)/hopsize);

%% %%%%%%%%%%%%%%%
% BEGIN InitAll %
%%% %%%%%%%%%%%%%%

Bark = [0     0	   50	 0.5
    1   100	  150	 1.5
    2   200	  250	 2.5
    3   300	  350	 3.5
    4   400	  450	 4.5
    5   510	  570	 5.5
    6   630	  700	 6.5
    7   770	  840	 7.5
    8   920	 1000	 8.5
    9  1080	 1170	 9.5
    10  1270 1370	10.5
    11  1480 1600	11.5
    12  1720 1850	12.5
    13  2000 2150	13.5
    14  2320 2500	14.5
    15  2700 2900	15.5
    16  3150 3400	16.5
    17  3700 4000	17.5
    18  4400 4800	18.5
    19  5300 5800	19.5
    20  6400 7000	20.5
    21  7700 8500	21.5
    22  9500 10500	22.5
    23 12000 13500	23.5
    24 15500 20000	24.5];

N2	= N/2+1;
dFs	= fs/N;
Bark2	= [
    sort([Bark(:,2);Bark(:,3)]),...
    sort([Bark(:,1);Bark(:,4)])
    ];
N0	= round(20*N/fs)+1;         % low frequency index @ 20 Hz
N01	= N0-1;
Ntop	= round(20000*N/fs)+1; % high frequency index @ 20 kHz?


%% Make list with Barknumber of each frequency bin

Barkno	  = zeros(1,N2);
f	      = N0:1:Ntop;
Barkno(f) = interp1(Bark2(:,1),Bark2(:,2),(f-1)*dFs);

%% Make list of frequency bins closest to Cf's

Cf = ones(2,24);
for a=1:1:24
    Cf(1,a)=round(Bark((a+1),2)*N/fs)+1-N0;
    Cf(2,a)=Bark(a+1,2);
end

%% Make list of frequency bins closest to Critical Band Border frequencies

Bf = ones(2,24);
Bf(1,1)=round(Bark(1,3)*N/fs);

for a=1:1:24
    Bf(1,a+1)=round(Bark((a+1),3)*N/fs)+1-N0;
    Bf(2,a)=Bf(1,a)-1;
end

Bf(2,25)=round(Bark((25),3)*N/fs)+1-N0;

%% Make list of minimum excitation (Hearing Treshold)

HTres= [0		  130
    0.01      70
    0.17	  60
    0.8	      30
    1	      25
    1.5	      20
    2		  15
    3.3	      10
    4		  8.1
    5		  6.3
    6		  5
    8		  3.5
    10		  2.5
    12		  1.7
    13.3	  0
    15		 -2.5
    16		 -4
    17		 -3.7
    18		 -1.5
    19		  1.4
    20		  3.8
    21		  5
    22		  7.5
    23 	      15
    24 	      48
    24.5 	  60
    25		  130];

k = (N0:1:Ntop);
MinExcdB = interp1(HTres(:,1),HTres(:,2),Barkno(k));

%% Initialize constants and variables

dz   = 0.5; % Barks
z    = (0.5:dz:23.5)'; % frequency in Barks
zb    = sort([Bf(1,:),Cf(1,:)]);
MinBf = MinExcdB(zb);
ei    = zeros(47,N);
Fei   = zeros(47,N);

gr  = [
    0 1 2.5 4.9  6.5 8 9 10 11 11.5 13 17.5 21 24
    0 0.35 0.7 0.7 1.1 1.25 1.26 1.18 1.08 1 0.66 0.46 0.38 0.3
    ];

gzi    = zeros(1,47);
h0     = zeros(1,47);
k      = 1:1:47;
gzi(k) = sqrt(interp1(gr(1,:)',gr(2,:)',k/2,'spline'));


% calculate a0
a0tab =	[ 0	     0
    10	 0
    12	 1.15
    13	 2.31
    14	 3.85
    15	 5.62
    16	 6.92
    16.5	 7.38
    17	 6.92
    18	 4.23
    18.5	 2.31
    19	 0
    20	-1.43
    21	-2.59
    21.5	-3.57
    22	-5.19
    22.5	-7.41
    23	-11.3
    23.5	-20
    24	-40
    25	-130
    26	-999];

a0    = ones(1,N);
k     = (N0:1:Ntop);
a0(k) = db2mag(interp1(a0tab(:,1),a0tab(:,2),Barkno(k)));

%%%%%%%%%%%%%%%
% END InitAll %
%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%
% BEGIN Hweights %
%%% %%%%%%%%%%%%%%%
% weights for freq. bins < N/2

DCbins	= 2;

H2 = [
    0       0
    17      0.8
    23		0.95
    25		0.975
    32		1
    37		0.975
    48		0.9
    67      0.8
    90		0.7
    114     0.6
    171     0.4
    206     0.3
    247     0.2
    294     0.1
    358     0
    ];

H5 = [
    0       0
    32      0.8
    43      0.95
    56      1
    69      0.975
    92      0.9
    120     0.8
    142     0.7
    165     0.6
    231     0.4
    277     0.3
    331     0.2
    397     0.1
    502     0
    ];

H16 = [
    0		0
    23.5	0.4
    34		0.6
    47		0.8
    56		0.9
    63		0.95
    79		1
    100     0.975
    115     0.95
    135     0.9
    159     0.85
    172     0.8
    194     0.7
    215     0.6
    244     0.5
    290     0.4
    348     0.3
    415     0.2
    500     0.1
    645     0
    ];

H21 = [
    0		0
    19		0.4
    44		0.8
    52.5	0.9
    58		0.95
    75		1
    101.5	0.95
    114.5	0.9
    132.5	0.85
    143.5	0.8
    165.5	0.7
    197.5   0.6
    241     0.5
    290     0.4
    348     0.3
    415     0.2
    500     0.1
    645     0
    ];


H42 = [
    0		0
    15		0.4
    41		0.8
    49		0.9
    53		0.965
    64		0.99
    71		1
    88		0.95
    94		0.9
    106     0.85
    115     0.8
    137     0.7
    180     0.6
    238     0.5
    290     0.4
    348     0.3
    415     0.2
    500     0.1
    645     0
    ];

Hweight	= zeros(47,N);

% weighting function H2
last	= floor((358/fs)*N) ;
k	= DCbins+1:1:last;
f	= (k-1)*fs/N;
Hweight(2,k) = interp1(H2(:,1),H2(:,2),f(k-DCbins));

% weighting function H5
last	=	floor((502/fs)*N);
k	=	DCbins+1:1:last;
f	=	(k-1)*fs/N;
Hweight(5,k)	= interp1(H5(:,1),H5(:,2),f(k-DCbins));

% weighting function H16
last	=	floor((645/fs)*N);
k	=	DCbins+1:1:last;
f	=	(k-1)*fs/N;
Hweight(16,k)	= interp1(H16(:,1),H16(:,2),f(k-DCbins));

% weighting function H21
Hweight(21,k)	= interp1(H21(:,1),H21(:,2),f(k-DCbins));

% weighting function H42
Hweight(42,k)	= interp1(H42(:,1),H42(:,2),f(k-DCbins));

% H1-H4
Hweight(1,:) = Hweight(2,:);
Hweight(3,:) = Hweight(2,:);
Hweight(4,:) = Hweight(2,:);

% H5-H15
for l =	6:1:15
    Hweight(l,:) = Hweight(5,:);
end

% H17-H20
for l =	17:1:20
    Hweight(l,:) = Hweight(16,:);
end

% H22-H41
for l =	22:1:41
    Hweight(l,:) = Hweight(21,:);
end

% H43-H47
for l =	43:1:47
    Hweight(l,:) = Hweight(42,:);
end

%%%%%%%%%%%%%%%%
% END Hweights %
%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%
% BEGIN process window %
%%%%%%%%%%%%%%%%%%%%%%%%

AmpCal = db2mag(91.2)*2/(N*mean(blackman(N, 'periodic')));
%     AmpCal=length(window)/sum(window);

% Calibration between wav-level and loudness-level (assuming
% blackman window and FFT will follow)

Chno	=	47;     % number of channels
Cal	 	=	0.50;   % calibration factor, twice the old value (0.25)
qb		=	N0:1:Ntop;
freqs	=	(qb+1)*fs/N;
hBPi	=	zeros(Chno,N);
hBPrms	=	zeros(1,Chno);
mdept	=	zeros(1,Chno);
ki		=	zeros(1,Chno-2);
ri		=	zeros(1,Chno);

startIndex = 1;
endIndex = N;
[TimePoints,R_mat,SPL_mat] = deal(zeros(n,1));
ri_mat = zeros(Chno,n);

for windowNum = 1:n    %for each frame
    
    dataIn = audio(startIndex:endIndex,1).*window;
    currentTimePoint = startIndex/fs;
    
    % Calculate Excitation Patterns
    TempIn =  dataIn*AmpCal;
    [rt,~]=size(TempIn);
    [r,~]=size(a0);
    if rt~=r; TempIn=TempIn'; end
    
    TempIn	=	a0.*fft(TempIn);
    Lg		=	abs(TempIn(qb));    % get absolute value of fourier transform for  indices in range of human hearing
    LdB		=	mag2db(Lg);
    whichL	=	find(LdB>MinExcdB); % extract indices where FFT magnitudes exceed excitation threshold
    sizL	=	length(whichL);     % get number of frequencies where this holds
    
    % steepness of slopes (Terhardt)
    S1 = -27;
    S2 = zeros(1,sizL);             % preallocate
    
    for w = 1:1:sizL
        
        % Steepness of upper slope [dB/Bark] in accordance with Terhardt
        steep = -24-(230/freqs(w))+(0.2*LdB(whichL(w)));
        
        if steep < 0
            S2(w) = steep;      % set S2 with steepness value calculated earlier
        end
    end
    whichZ	= zeros(2,sizL);    % preallocate
    qd		= 1:1:sizL;         % indices of frequencies above excitation threshold
    whichZ(1,:)	= floor(2*Barkno(whichL(qd)+N01));  % get bark band numbers
    whichZ(2,:)	= ceil(2*Barkno(whichL(qd)+N01));
    
    ExcAmp = zeros(sizL,Chno);
    Slopes = zeros(sizL,Chno);
    
    for k=1:1:sizL    %loop over freq indices above threshold
        Ltmp = LdB(whichL(k)); % copy FFT magnitude (in dB) above threshold
        Btmp = Barkno(whichL(k)+N01); % and the bark number associat
        
        for l = 1:1:whichZ(1,k) % loop up to floored bark number of freq index k
            Stemp = (S1*(Btmp-(l*0.5)))+Ltmp;
            if Stemp>MinBf(l)
                Slopes(k,l)=db2mag(Stemp);
            end
        end
        
        for l = whichZ(2,k):1:Chno % loop up to ceil'd bark number
            Stemp =	(S2(k)*((l*0.5)-Btmp))+Ltmp;
            if Stemp>MinBf(l)
                Slopes(k,l)=db2mag(Stemp); % critical filterbank upper side
            end
        end
    end
    
    for k=1:Chno % loop over each channel
        etmp = zeros(1,N);
        for l=1:1:sizL   % for each l index of fft bin in human hearing freq range
            N1tmp = whichL(l); % get freq index of bin
            N2tmp = N1tmp + N01;
            if (whichZ(1,l) == k)
                ExcAmp(N1tmp, k) = 1;
            elseif (whichZ(2,l) == k)
                ExcAmp(N1tmp, k) = 1;
            elseif (whichZ(2,l) > k)
                ExcAmp(N1tmp,k) = Slopes(l,k+1)/Lg(N1tmp);
            else
                ExcAmp(N1tmp,k) = Slopes(l,k-1)/Lg(N1tmp);
            end
            etmp(N2tmp) = ExcAmp(N1tmp,k)*TempIn(N2tmp);
        end      % this is the specific excitation time function
        
        % ifft to get time domain blocks of signal
        ei(k,:)	= N*real(ifft(etmp));
        etmp	= abs(ei(k,:));
        h0(k)	= mean(etmp);
        Fei(k,:)	= fft(etmp-h0(k));
        hBPi(k,:)	= 2*real(ifft(Fei(k,:).*Hweight(k,:)));
        hBPrms(k)	= rms(hBPi(k,:));
        
        if h0(k)>0
            mdept(k) = hBPrms(k)/h0(k);
            
            if mdept(k)>1
                mdept(k)=1;
            end
        else
            mdept(k)=0;
        end
    end
    
    % find cross-correlation coefficients
    for k=1:1:45
        cfac	=	cov(hBPi(k,:),hBPi(k+2,:));
        den	=	diag(cfac);
        den	=	sqrt(den*den');
        if den(2,1)>0
            ki(k)	=	cfac(2,1)/den(2,1);
        else
            ki(k)	=	0;
        end
    end
    
    % Calculate specific roughness ri and total roughness R
    ri(1)	=	(gzi(1)*mdept(1)*ki(1))^2;
    ri(2)	=	(gzi(2)*mdept(2)*ki(2))^2;
    
    for k = 3:1:45
        ri(k)	=	(gzi(k)*mdept(k)*ki(k-2)*ki(k))^2;
    end
    
    ri(46)	=	(gzi(46)*mdept(46)*ki(44))^2;
    ri(47)	=	(gzi(47)*mdept(47)*ki(45))^2;
    
    ri      = Cal*ri; % appropriately scaled specific roughness
    R       = dz*sum(ri); % total R = integration of the specific R pattern
    
    SPL = mean(rms(dataIn));
    if SPL > 0
        SPL = mag2db(SPL)+83; % -20 dBFS <--> 60 dB SPL
    else
        SPL = -400;
    end
    
    % matrices to return
    R_mat(windowNum) = R;
    ri_mat(1:Chno,windowNum) = ri;
    SPL_mat(windowNum) = SPL;
    
    startIndex = startIndex+hopsize;
    endIndex = endIndex+hopsize;
    TimePoints(windowNum,1) = currentTimePoint;
end

%%%%%%%%%%%%%%%%%%%%%%
% END process window %
%%%%%%%%%%%%%%%%%%%%%%

%% *************************************************************************
% output struct
% *************************************************************************

% main output results
OUT.InstantaneousRoughness = R_mat;                       % instantaneous roughness
OUT.InstantaneousSpecificRoughness = ri_mat;              % time-varying specific roughness
OUT.TimeAveragedSpecificRoughness = mean(ri_mat,2);       % mean specific roughness
OUT.time = TimePoints;                                    % time
OUT.barkAxis = z;                                         % critical band rate (for specific roughness)
OUT.dz = dz;

% Roughness statistics based on InstantaneousRoughness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,idx] = min( abs(OUT.time - start_skip) ); % find idx of start_skip on time vector
[~,idxEnd] = min( abs(OUT.time - (max(OUT.time) - end_skip)) ); % find idx of end_skip on time vector

OUT.Rmax = max(R_mat(idx:idxEnd));
OUT.Rmin = min(R_mat(idx:idxEnd));
OUT.Rmean = mean(R_mat(idx:idxEnd));
OUT.Rstd = std(R_mat(idx:idxEnd));
OUT.R1 = get_percentile(R_mat(idx:idxEnd),1);
OUT.R2 = get_percentile(R_mat(idx:idxEnd),2);
OUT.R3 = get_percentile(R_mat(idx:idxEnd),3);
OUT.R4 = get_percentile(R_mat(idx:idxEnd),4);
OUT.R5 = get_percentile(R_mat(idx:idxEnd),5);
OUT.R10 = get_percentile(R_mat(idx:idxEnd),10);
OUT.R20 = get_percentile(R_mat(idx:idxEnd),20);
OUT.R30 = get_percentile(R_mat(idx:idxEnd),30);
OUT.R40 = get_percentile(R_mat(idx:idxEnd),40);
OUT.R50 = median(R_mat(idx:idxEnd));
OUT.R60 = get_percentile(R_mat(idx:idxEnd),60);
OUT.R70 = get_percentile(R_mat(idx:idxEnd),70);
OUT.R80 = get_percentile(R_mat(idx:idxEnd),80);
OUT.R90 = get_percentile(R_mat(idx:idxEnd),90);
OUT.R95 = get_percentile(R_mat(idx:idxEnd),95);

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
    
    plot((1:47)'/2, mean(ri_mat,2),'r-');
    
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