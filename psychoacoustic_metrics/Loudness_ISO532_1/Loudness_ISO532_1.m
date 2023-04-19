function OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
% function OUT = Loudness_ISO532_1(insig, fs, field, method, time_skip, show)
%
%  Zwicker Loudness model according to ISO 532-1 for stationary 
%  signals (Method A) and arbitrary signals (Method B)  
%
%  Reference signal: 40 dBSPL 1 kHz tone yields 1 sone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT ARGUMENTS
%   insig : array 
%   for method = 0 [1xN] array, insig is an array containing N=28 third octave unweighted SPL from 25 Hz to 12500 Hz
%   for method = 1 and method = 2 [Nx1] array, insig is a monophonic calibrated audio signal (Pa), 1 channel only as specified by the standard
%
%   fs : integer 
%   sampling frequency (Hz). For method = 0, provide a dummy scalar
%
%   field : integer 
%   free field = 0; diffuse field = 1;
%
%   method : integer 
%   0 = stationary (from input 1/3 octave unweighted SPL)
%   1 = stationary (from audio file) 
%   2 = time varying (from audio file)  
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUTS (method==0 and method==1; stationary method)
%   OUT : struct containing the following fields
%
%       * time_insig - time vector of the audio input, in seconds
%       * barkAxis - bark vector
%       * SpecificLoudness - time-averaged specific loudness (sone/Bark) 
%       * Loudness - loudness (sone) 
%       * LoudnessLevel - loudness level (phon)
%       * TimeAveragedSPL - time-averaged overall SPL (1/3 octave bands, DBSPL)
% 
% OUTPUTS (method==2; time-varying method)
%   OUT : struct containing the following fields
%
%       * barkAxis - vector of Bark band numbers used for specific loudness computation
%       * time - time vector of the final loudness calculation, in seconds
%       * time_insig - time vector of insig, in seconds
%       * InstantaneousLoudness - instantaneous loudness (sone) vs time
%       * InstantaneousSpecificLoudness - specific loudness (sone/Bark) vs time
%       * InstantaneousLoudnessLevel - instantaneous loudness level (phon) vs time
%       * SpecificLoudness - time-averaged specific loudness (sone/Bark) 
%       * InstantaneousSPL - overall SPL (1/3 octave bands) for each time step, in dBSPL
%       * Several statistics based on the InstantaneousLoudness 
%         ** Nmean : mean value of InstantaneousLoudness (sone)
%         ** Nstd : standard deviation of InstantaneousLoudness (sone)
%         ** Nmax : maximum of InstantaneousLoudness (sone)
%         ** Nmin : minimum of InstantaneousLoudness (sone)
%         ** Nx : percentile loudness exceeded during x percent of the signal (sone)
%         ** N_ratio : ratio between N5/N95 ( 1.1 (stationary)> N_ratio > 1.1 (time varying) )
%           *** HINT: loudness calculation takes some time to have a steady-response
%                     therefore, it is a good practice to consider a time_skip to compute the statistics 
%                     due to transient effects in the beginning of the loudness calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source: C code is provided in the ISO532 Annex A (2014).
%
% Author: Ella Manor - MATLAB implementation for AARAE (2015)
% Author: Gil Felix Greco, Braunschweig 22.02.2023 - adapted and validatec 
%                   for SQAT. The validation was based on the test signals
%                   provided from ISO 532-1:2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help Loudness_ISO532_1;
    return;
end

if nargin < 6
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

if method==0
    
    if size(insig,1)~=1 % if the insig is not a [Nx1] array
        insig=insig';   % correct the dimension of the insig
    end
    
elseif method==1 || method==2
    
    if size(insig,2)~=1 % if the insig is not a [Nx1] array
        insig=insig';   % correct the dimension of the insig
    end
    
end

% Time constants for non-linear temporal decay
Tshort = 0.005;
Tlong = 0.015;
Tvar = 0.075;

% Factors for virtual upsampling/inner iterations
NL_ITER = 24;
% Sampling rate to which third-octave-levels are downsampled
SR_LEVEL = 2000;
% Sampling rate to which output/total summed loudness is downsampled
SR_LOUDNESS = 500;
% Tiny value for adjusting intensity levels for stationary signals
TINY_VALUE = 1e-12;
% ref value for stationary signals
I_REF =4e-10;
% bark vector
barkAxis=(1:240)/10;

%% ***************************
% STEP 1 - Resample to 48 kHz
% ****************************

if method ~= 0 % if different from stationary (from input 1/3 octave unweighted SPL)
    
    audio=insig; % input signal
    
    if fs ~= 48000
        gcd_fs = gcd(48000,fs); % greatest common denominator
        audio = resample(audio,48000/gcd_fs,fs/gcd_fs);
        fs = 48000;
        len = size(audio,1);
    else
        len = size(audio,1);
    end

    %% ******************************************************************
    % Assign values to global variables according to the selected method
    % *******************************************************************
    
    if method == 1 % stationary from audio signal
        SampleRateLevel = 1;
        DecFactorLoudness = 1;
        NumSamplesLevel = 1;
    elseif method == 2 % time_varying from audio signal
        SampleRateLevel = SR_LEVEL;
        SampleRateLoudness = SR_LOUDNESS;
        DecFactorLevel = fs/SampleRateLevel;
        DecFactorLoudness = SampleRateLevel/SampleRateLoudness;
        NumSamplesLevel = ceil(len/DecFactorLevel);
        NumSamplesLoudness = ceil(NumSamplesLevel/DecFactorLoudness);
    end
    
    %% **************************************************
    % STEP 2 - Create filter bank and filter the signal
    % ***************************************************
    
    % reference
    br = [1,2,1;1,0,-1;1,-2,1];
    ar = [1,-2,1;1,-2,1;1,-2,1];
    
    CenterFrequency = zeros(28,1);
    for i = 1:28
        CenterFrequency(i) = 10^(((i-1)-16)/10.) * 1000; % calculate centre frequencies
    end
    
    % filter 'a' coefficient offsets TABLES A.1 A.2
    ad = cat(3,[0,-6.70260e-004,6.59453e-004;...
        0,-3.75071e-004,3.61926e-004;...
        0,-3.06523e-004,2.97634e-004],... % 25 Hz
        [0,-8.47258e-004,8.30131e-004;...
        0,-4.76448e-004,4.55616e-004;...
        0,-3.88773e-004,3.74685e-004],... % 31.5 Hz
        [0,-1.07210e-003,1.04496e-003;...
        0,-6.06567e-004,5.73553e-004;...
        0,-4.94004e-004,4.71677e-004],... % 40 Hz
        [0,-1.35836e-003,1.31535e-003;...
        0,-7.74327e-004,7.22007e-004;...
        0,-6.29154e-004,5.93771e-004],... % 50 Hz
        [0,-1.72380e-003,1.65564e-003;...
        0,-9.91780e-004,9.08866e-004;...
        0,-8.03529e-004,7.47455e-004],... % 63 Hz
        [0,-2.19188e-003,2.08388e-003;...
        0,-1.27545e-003,1.14406e-003;...
        0,-1.02976e-003,9.40900e-004],... % 80 Hz
        [0,-2.79386e-003,2.62274e-003;...
        0,-1.64828e-003,1.44006e-003;...
        0,-1.32520e-003,1.18438e-003],... % 100 Hz
        [0,-3.57182e-003,3.30071e-003;...
        0,-2.14252e-003,1.81258e-003;...
        0,-1.71397e-003,1.49082e-003],... % 125 Hz
        [0,-4.58305e-003,4.15355e-003;...
        0,-2.80413e-003,2.28135e-003;...
        0,-2.23006e-003,1.87646e-003],... % 160 Hz
        [0,-5.90655e-003,5.22622e-003;...
        0,-3.69947e-003,2.87118e-003;...
        0,-2.92205e-003,2.36178e-003],... % 200 Hz
        [0,-7.65243e-003,6.57493e-003;...
        0,-4.92540e-003,3.61318e-003;...
        0,-3.86007e-003,2.97240e-003],... % 250 Hz
        [0,-1.00023e-002,8.29610e-003;...
        0,-6.63788e-003,4.55999e-003;...
        0,-5.15982e-003,3.75306e-003],... % 315 Hz
        [0,-1.31230e-002,1.04220e-002;...
        0,-9.02274e-003,5.73132e-003;...
        0,-6.94543e-003,4.71734e-003],... % 400 Hz
        [0,-1.73693e-002,1.30947e-002;...
        0,-1.24176e-002,7.20526e-003;...
        0,-9.46002e-003,5.93145e-003],... % 500 Hz
        [0,-2.31934e-002,1.64308e-002;...
        0,-1.73009e-002,9.04761e-003;...
        0,-1.30358e-002,7.44926e-003],... % 630 Hz
        [0,-3.13292e-002,2.06370e-002;...
        0,-2.44342e-002,1.13731e-002;...
        0,-1.82108e-002,9.36778e-003],... % 800 Hz
        [0,-4.28261e-002,2.59325e-002;...
        0,-3.49619e-002,1.43046e-002;...
        0,-2.57855e-002,1.17912e-002],... % 1000 Hz
        [0,-5.91733e-002,3.25054e-002;...
        0,-5.06072e-002,1.79513e-002;...
        0,-3.69401e-002,1.48094e-002],... % 1250 Hz
        [0,-8.26348e-002,4.05894e-002;...
        0,-7.40348e-002,2.24476e-002;...
        0,-5.34977e-002,1.85371e-002],... % 1600 Hz
        [0,-1.17018e-001,5.08116e-002;...
        0,-1.09516e-001,2.81387e-002;...
        0,-7.85097e-002,2.32872e-002],... % 2000 Hz
        [0,-1.67714e-001,6.37872e-002;...
        0,-1.63378e-001,3.53729e-002;...
        0,-1.16419e-001,2.93723e-002],... % 2500 Hz
        [0,-2.42528e-001,7.98576e-002;...
        0,-2.45161e-001,4.43370e-002;...
        0,-1.73972e-001,3.70015e-002],... % 3150 Hz
        [0,-3.53142e-001,9.96330e-002;...
        0,-3.69163e-001,5.53535e-002;...
        0,-2.61399e-001,4.65428e-002],... % 4000 Hz
        [0,-5.16316e-001,1.24177e-001;...
        0,-5.55473e-001,6.89403e-002;...
        0,-3.93998e-001,5.86715e-002],... % 5000 Hz
        [0,-7.56635e-001,1.55023e-001;...
        0,-8.34281e-001,8.58123e-002;...
        0,-5.94547e-001,7.43960e-002],... % 6300 Hz
        [0,-1.10165e+000,1.91713e-001;...
        0,-1.23939e+000,1.05243e-001;...
        0,-8.91666e-001,9.40354e-002],... % 8000 Hz
        [0,-1.58477e+000,2.39049e-001;...
        0,-1.80505e+000,1.28794e-001;...
        0,-1.32500e+000,1.21333e-001],... % 10000 Hz
        [0,-2.50630e+000,1.42308e-001;...
        0,-2.19464e+000,2.76470e-001;...
        0,-1.90231e+000,1.47304e-001]); % 12500 Hz
    
    % filter gains
    filtgain = [4.30764e-011;... % 25 Hz
        8.59340e-011;... % 31.5 Hz
        1.71424e-010;... % 40 Hz
        3.41944e-010;... % 50 Hz
        6.82035e-010;... % 63 Hz
        1.36026e-009;... % 80 Hz
        2.71261e-009;... % 100 Hz
        5.40870e-009;... % 125 Hz
        1.07826e-008;... % 160 Hz
        2.14910e-008;... % 200 Hz
        4.28228e-008;... % 250 Hz
        8.54316e-008;... % 315 Hz
        1.70009e-007;... % 400 Hz
        3.38215e-007;... % 500 Hz
        6.71990e-007;... % 630 Hz
        1.33531e-006;... % 800 Hz
        2.65172e-006;... % 1000 Hz
        5.25477e-006;... % 1250 Hz
        1.03780e-005;... % 1600 Hz
        2.04870e-005;... % 2000 Hz
        4.05198e-005;... % 2500 Hz
        7.97914e-005;... % 3150 Hz
        1.56511e-004;... % 4000 Hz
        3.04954e-004;... % 5000 Hz
        5.99157e-004;... % 6300 Hz
        1.16544e-003;... % 8000 Hz
        2.27488e-003;... % 10000 Hz
        3.91006e-003];   % 12500 Hz
    
    filteredaudio = zeros(len,28);
    
    for n = 1:28
        filteredaudio(:,n) = filtgain(n) * filter(br(3,:),ar(3,:)-ad(3,:,n),...
            filter(br(2,:),ar(2,:)-ad(2,:,n),...
            filter(br(1,:),ar(1,:)-ad(1,:,n),audio)));
    end
    
    %% **************************************************************
    % STEP 3 - Squaring and smoothing by 3 1st order lowpass filters
    % ***************************************************************
    
    filteredaudio = filteredaudio.^2;
    
    ThirdOctaveLevel = zeros(NumSamplesLevel,28);
    
    for i = 1:28
        
        if method == 2 % time-varying from audio signal
            
            smoothedaudio = zeros(len,28);
            
            if CenterFrequency(i) <= 1000
                Tau = 2/(3*CenterFrequency(i));
            else
                Tau = 2/(3*1000.);
            end
            
            % 3x smoothing 1st order low-pass filters in series
            A1 = exp(-1 ./ (fs * Tau));
            B0 = 1 - A1;
            Y1 = 0;
            for k = 1:3
                for j = 1:length(filteredaudio)
                    %                 smoothedaudio(j,i) = A1*temp(j,i) + B0*Y1;
                    smoothedaudio(j,i)= (B0*filteredaudio(j,i))+(A1*Y1); % <----- modified from original by gfg
                    Y1 = smoothedaudio(j,i);
                end
            end
            
            c=1;
            for j = 1:NumSamplesLevel
                ThirdOctaveLevel(j,i) = 10*log10((smoothedaudio(c,i)+TINY_VALUE)/I_REF);
                c = c+DecFactorLevel;
            end
            
        elseif method == 1 % stationary from audio signal
            
            NumSkip = floor(time_skip * fs);
            smoothedaudio = zeros(len-NumSkip,28);
            
            if NumSkip > len/2
                warndlg('time signal too short');
            end
            
            if NumSkip == 0; NumSkip = 1;end
            smoothedaudio(1:len-NumSkip,i) = filteredaudio(NumSkip:len-1,i);
            %         ThirdOctaveLevel(NumSamplesLevel,i) = 10*log10((sum(smoothedaudio(:,i))/len+TINY_VALUE)/I_REF);
            ThirdOctaveLevel(NumSamplesLevel,i) = 10*log10((sum(smoothedaudio(:,i)/len)+TINY_VALUE)/I_REF); % <----- modified from original by gfg
        end
        
    end

else % if method == 0
    
    SampleRateLevel = 1;
    DecFactorLoudness = 1;
    NumSamplesLevel = 1;
    
    ThirdOctaveLevel=insig; % get 1/3 octave levels from insig if method = 0
end

%% ***********************************************************
% STEP 4 - Apply weighting factor to first 3 1/1 octave bands
% ************************************************************

% WEIGHTING BELLOW 315Hz TABLE A.3

% Ranges of 1/3 Oct bands for correction at low frequencies according to equal loudness contours
RAP = [45 55 65 71 80 90 100 120 ];

% Reduction of 1/3 Oct Band levels at low frequencies according to equal loudness contours
% within the eight ranges defined by RAP (DLL)
DLL = [-32 -24 -16 -10 -5 0 -7 -3 0 -2 0;
    -29 -22 -15 -10 -4 0 -7 -2 0 -2 0;
    -27 -19 -14 -9 -4 0 -6 -2 0 -2 0;
    -25 -17 -12 -9 -3 0 -5 -2 0 -2 0;
    -23 -16 -11 -7 -3 0 -4 -1 0 -1 0;
    -20 -14 -10 -6 -3 0 -4 -1 0 -1 0;
    -18 -12 -9 -6 -2 0 -3 -1 0 -1 0;
    -15 -10 -8 -4 -2 0 -3 -1 0 -1 0];

CorrLevel = zeros(NumSamplesLevel,11);
Intens = zeros(NumSamplesLevel,11);
CBI = zeros(NumSamplesLevel,3);
LCB = zeros(NumSamplesLevel,3);

for j = 1:NumSamplesLevel
    for i = 1:size(DLL,2)
        k=1;
        while ( ( ThirdOctaveLevel(j,i) > RAP(k)-DLL(k,i) ) ) && (k < 8)
            k=k+1;
        end
        CorrLevel(j,i) = ThirdOctaveLevel(j,i) + DLL(k,i); % attenuated levels
        Intens(j,i) = 10^(CorrLevel(j,i)/10); % attenuated 1/3 octave intensities
    end

    % *************************************************************
    % STEP 5 - Sumup intensity values of the first 3 critical bands
    % *************************************************************

    CBI(j,1) = sum(Intens(j,1:6)); % first critical band (sum of octaves (25Hz to 80Hz))
    CBI(j,2) = sum(Intens(j,7:9)); % second critical band (sum of octaves (100Hz to 160Hz))
    CBI(j,3) = sum(Intens(j,10:11)); % third critical band (sum of octaves (200Hz to 250Hz))

    FNGi = 10*log10(CBI);

    for i = 1:3
        if CBI(j,i)>0
            LCB(j,i) = FNGi(j,i);
        else
            LCB(j,i) = 0;
        end
    end
end

%% **********************************************************************
% STEP 6 - Calculate core loudness for each critical band
% ***********************************************************************

% LEVEL CORRECTIONS TABLE A.5 (LDF0) DDF
% Level correction to convert from a free field to a diffuse field (last critical band 12.5kHz is not included)
DDF = [0 0 0.5 0.9 1.2 1.6 2.3 2.8 3 2 0 -1.4 -2 -1.9 -1 0.5 3 4 4.3 4];

% LEVEL CORRECTIONS TABLE A.6 (LTQ)
% Critical band level at absolute threshold without taking into account the 
% transmission characteristics of the ear
LTQ =  [30 18 12 8 7 6 5 4 3 3 3 3 3 3 3 3 3 3 3 3]; % Threshold due to internal noise
% Hearing thresholds for the excitation levels (each number corresponds to a critical band 12.5kHz is not included)

% LEVEL CORRECTIONS TABLE A.7 DCB
% Correction factor because using third octave band levels (rather than critical bands)
DCB = [-0.25 -0.6 -0.8 -0.8 -0.5 0 0.5 1.1 1.5 1.7 1.8 1.8 1.7 1.6 1.4 1.2 0.8 0.5 0 -0.5]; 

% LEVEL CORRECTIONS TABLE A.4 (A0)
% % Attenuation due to transmission in the middle ear
A0  = [ 0 0 0 0 0 0 0 0 0 0 -0.5 -1.6 -3.2 -5.4 -5.6 -4 -1.5 2 5 12];
% Moore et al disagrees with this being flat for low frequencies

Le = zeros(NumSamplesLevel,20);
CoreL = zeros(NumSamplesLevel,21);

for j = 1:NumSamplesLevel
    
    for i = 1:20
        
        Le(j,i) = ThirdOctaveLevel(j,i+8);
        
        if i <= 3
            Le(j,i) = LCB(j,i);
        end
        
        Le(j,i) = Le(j,i) - A0(i);
        
        if field == 1
            Le(j,i) = Le(j,i) + DDF(i);
        end
        
        if Le(j,i) > LTQ(i)
            S = 0.25;
            Le(j,i) = Le(j,i) - DCB(i);
            MP1 = 0.0635 .* 10.^(0.025 .* LTQ(i));
            MP2 =  ( ( (1 - S) + S.*10^(0.1.*(Le(j,i)-LTQ(i)))).^0.25) - 1;
            CoreL(j,i) = MP1 .* MP2;
            
            if CoreL(j,i) <= 0
                CoreL(j,i) = 0;
            end
        end
    end
end

%% *************************************************************************
% STEP 7 - Correction of specific loudness within the lowest critical band
% **************************************************************************

for j = 1:NumSamplesLevel
    CorrCL = 0.4 + 0.32 .* CoreL(j,1).^(0.2);
    
    if CorrCL > 1
        CorrCL = 1;
    end
    
    CoreL(j,1) = CoreL(j,1)*CorrCL;
end

%% **********************************************************************
% STEP 8 - Implementation of NL Block
% ***********************************************************************

if method == 2 % time-varying from audio signal
    
    DeltaT = 1 / (SampleRateLevel*NL_ITER);
    P = (Tvar + Tlong) / (Tvar*Tshort);
    Q = 1/(Tshort*Tvar);
    Lambda1 =-P/2 + sqrt(P*P/4 - Q);
    Lambda2 =-P/2 - sqrt(P*P/4 - Q);
    Den = Tvar * (Lambda1 - Lambda2);
    E1 = exp(Lambda1 * DeltaT);
    E2 = exp(Lambda2 * DeltaT);
    
    NlLpB(1) = (E1 - E2) / Den;
    NlLpB(2) =((Tvar * Lambda2 + 1) * E1 - (Tvar * Lambda1 + 1) * E2) / Den;
    NlLpB(3) =((Tvar * Lambda1 + 1) * E1 - (Tvar * Lambda2 + 1) * E2) / Den;
    NlLpB(4) = (Tvar * Lambda1+1) * (Tvar * Lambda2 + 1) * (E1-E2) / Den;
    NlLpB(5) = exp(-DeltaT / Tlong);
    NlLpB(6) = exp(-DeltaT / Tvar);

    for i = 1:21
        
        NlLpUoLast = 0; % At beginning capacitors C1 and C2 are discharged
        NlLpU2Last = 0;
        
        for j = 1:NumSamplesLevel-1
            NextInput = CoreL(j+1,i);
            % interpolation steps between current and next sample
            Delta = (NextInput - CoreL(j,i)) / NL_ITER;
            Ui = CoreL(j,i);
            
            % f_nl_lp FUNCTION STARTS
            % case 1
            if Ui < NlLpUoLast
                if NlLpUoLast > NlLpU2Last
                    % case 1.1
                    U2 = NlLpUoLast*NlLpB(1) - NlLpU2Last*NlLpB(2);
                    Uo = NlLpUoLast*NlLpB(3) - NlLpU2Last*NlLpB(4);
                    if  Uo < Ui
                        Uo  = Ui;
                    end
                    if U2 > Uo
                        U2 = Uo;
                    end
                else
                    % case 1.2
                    Uo = NlLpUoLast*NlLpB(5);
                    if  Uo < Ui
                        Uo = Ui;
                    end
                    U2 = Uo;
                end
                % case 2
            elseif Ui == NlLpUoLast
                Uo = Ui;
                % case 2.1
                if Uo > NlLpUoLast
                    U2 = (NlLpUoLast - Ui)*NlLpB(6) + Ui;
                    % case 2.2
                else
                    U2 = Ui;
                end
                % case 3
            else
                Uo = Ui;
                U2 = (NlLpU2Last - Ui)*NlLpB(6) + Ui;
            end
            
            NlLpUoLast = Uo;
            NlLpU2Last = U2;
            
            CoreL(j,i) = Uo;
            % f_nl_lp FUNCTION ENDS
            
            Ui = Ui + Delta;
            
            % inner iteration
            for k = 1:NL_ITER
                % f_nl_lp FUNCTION STARTS
                % case 1
                if Ui < NlLpUoLast
                    if NlLpUoLast > NlLpU2Last
                        % case 1.1
                        U2 = NlLpUoLast*NlLpB(1) - NlLpU2Last*NlLpB(2);
                        Uo = NlLpUoLast*NlLpB(3) - NlLpU2Last*NlLpB(4);
                        if Ui > Uo
                            Uo  = Ui;
                        end
                        if U2 > Uo
                            U2 = Uo;
                        end
                    else
                        % case 1.2
                        Uo = NlLpUoLast*NlLpB(5);
                        if Ui > Uo
                            Uo = Ui;
                        end
                        U2 = Uo;
                    end
                    % case 2
                elseif Ui == NlLpUoLast
                    Uo = Ui;
                    % case 2.1
                    if Uo > NlLpUoLast
                        U2 = (NlLpUoLast - Ui)*NlLpB(6) + Ui;
                        % case 2.2
                    else
                        U2 = Ui;
                    end
                    % case 3
                else
                    Uo = Ui;
                    U2 = (NlLpU2Last - Ui)*NlLpB(6) + Ui;
                end
                
                NlLpUoLast = Uo;
                NlLpU2Last = U2;
                
                CoreL(j,i) = Uo;
                % f_nl_lp FUNCTION ENDS
                Ui = Ui + Delta;
            end
        end
    end
    
end

%% **********************************************************************
% STEP 9 - CALCULATE THE SLOPES
% ***********************************************************************

% Upper limits of the approximated critical bands in Bark
% TABLE A.8
ZUP  = [.9 1.8 2.8 3.5 4.4 5.4 6.6 7.9 9.2 10.6 12.3 13.8 15.2 16.7 18.1 19.3 20.6 21.8 22.7 23.6 24];

% TABLE A.9
% Range of specific loudness for the determination of the steepness of the upper slopes in the specific loudness 
% - critical band rate pattern (used to plot the correct USL curve)
RNS = [21.5 18 15.1 11.5 9 6.1 4.4 3.1 2.13 1.36 0.82 0.42 0.30 0.22 0.15 0.10 0.035 0];

% This is used to design the right hand slope of the loudness
USL = [13 8.2 6.3 5.5 5.5 5.5 5.5 5.5;
    9 7.5 6 5.1 4.5 4.5 4.5 4.5;
    7.8 6.7 5.6 4.9 4.4 3.9 3.9 3.9;
    6.2 5.4 4.6 4.0 3.5 3.2 3.2 3.2;
    4.5 3.8 3.6 3.2 2.9 2.7 2.7 2.7;
    3.7 3.0 2.8 2.35 2.2 2.2 2.2 2.2;
    2.9 2.3 2.1 1.9 1.8 1.7 1.7 1.7;
    2.4 1.7 1.5 1.35 1.3 1.3 1.3 1.3;
    1.95 1.45 1.3 1.15 1.1 1.1 1.1 1.1;
    1.5 1.2 0.94 0.86 0.82 0.82 0.82 0.82;
    0.72 0.67 0.64 0.63 0.62 0.62 0.62 0.62;
    0.59 0.53 0.51 0.50 0.42 0.42 0.42 0.42;
    0.40 0.33 0.26 0.24 0.24 0.22 0.22 0.22;
    0.27 0.21 0.20 0.18 0.17 0.17 0.17 0.17;
    0.16 0.15 0.14 0.12 0.11 0.11 0.11 0.11;
    0.12 0.11 0.10 0.08 0.08 0.08 0.08 0.08;
    0.09 0.08 0.07 0.06 0.06 0.06 0.06 0.05;
    0.06 0.05 0.03 0.02 0.02 0.02 0.02 0.02];

LN = zeros(NumSamplesLevel,1);
N_mat = zeros(NumSamplesLevel,1);
Spec_N = zeros(1,240);
ZUP = ZUP+0.0001; %<----- add constant factor to ZUP according to code provided by ISO 532-1 (see ISO 532-1 - Program etc\Annex A.4\ISO_532-1_LIB\src\ISO_532-1.c - line 862)
ns = zeros(NumSamplesLevel,240);   

for l = 1:NumSamplesLevel
    
    N = 0;
    z1 = 0; % critical band rate starts at 0
    n1 = 0; % loudness level starts at 0
    iz = 1;
    z = 0.1;
    j=18;
    
    for i = 1:21 % specific loudness
        
        % Determines where to start on the slope
        ig = i - 1;
        
        % steepness of upper slope (USL) for bands above 8th one are identical
        if ig > 8
            ig = 8;
        end
        
        while z1 < ZUP(i)
            
            if n1 <=  CoreL(l,i)      % Nm is the main loudness level
                % contribution of unmasked main loudness to total loudness
                % and calculation of values
                if n1 < CoreL(l,i)
                    j=1;
                    
                    while (RNS(j) > CoreL(l,i)) && (j < 18) % the value of j is used below to build a slope
                        j = j+1; % j becomes the index at which Nm(i)                        % to the range of specific loudness
                    end
                end
                
                z2 = ZUP(i);
                n2 = CoreL(l,i);
                N = N + n2*(z2-z1);
                k = z;                     % initialisation of k
                
                while (k <= z2)
                    ns(l,iz) = n2;
                    iz = iz + 1;
                    k = k+(1/10);
                end
                
                z = k;
                
            else %if N1 > NM(i)
                % decision wether the critical band in question is completely
                % or partly masked by accessory loudness
                
                n2 = RNS(j);
                
                if n2 < CoreL(l,i)
                    n2 = CoreL(l,i);
                end
                
                dz = (n1-n2) / USL(j,ig);
                z2 = z1 + dz;
                
                if z2 > ZUP(i)
                    z2 = ZUP(i);
                    dz = z2 - z1;
                    n2 = n1 - dz*USL(j,ig);
                end
                
                N = N + dz*(n1+n2)/2;
                k = z;                     % initialisation of k
                
                while (k <= z2)
                    ns(l,iz) = n1 - (k-z1)*USL(j,ig);
                    iz = iz + 1;
                    k = k+(1/10);
                end
                
                z = k;
                
            end
            
            if (n2 <= RNS(j)) && (j < 18)
                j = j + 1;
            end
            
            if (n2 <= RNS(j)) && (j >= 18)
                j = 18;
            end
            
            z1 = z2;     % N1 and Z1 for next loop
            n1 = n2;
            
        end
    end
    
    if N < 0
        N = 0;
    end
    
    if N <= 16
        N = (N*1000+.5)/1000;
    else
        N = (N*100+.5)/100;
    end
    
    LN(l) = 40*(N + .0005)^.35;
    
    if LN(l) < 3
        LN(l) = 3;
    end
    
    if N >= 1
        LN(l) = 10*log10(N)/log10(2) + 40;
    end
    
    if method==0 || method==1 % stationary method
        LN = 40 * N.^0.35;
        LN( N>=1 ) = 40 + 10*log2( N( N>=1 ) );
        LN( LN < 0 ) = 0;
        LN( LN < 3 ) = 3;
    end
    
    N_mat(l) = N; % total loudness at current timeframe l
end

% specific Loudness as a function of Bark number
for i = 1:240
    Spec_N(i) = mean(ns(:,i));
end

%% **********************************************************************
% STEP 10 - Apply Temporal Weighting to Arbitrary signals
% ***********************************************************************

if method == 2 % time-varying from audio signal
    
    Loudness_t1 = zeros(NumSamplesLevel,1);
    Loudness_t2 = zeros(NumSamplesLevel,1);
    Loudness = zeros(NumSamplesLevel,1);
    
    % 1st order low-pass A
    Tau = 3.5e-3;
    A1 = exp(-1 / (SampleRateLevel * DecFactorLevel * Tau));
    B0 = 1 - A1;
    Y1 = 0;
    
    for i = 1:NumSamplesLevel
        X0 = N_mat(i);
        Y1 = B0 * X0 + A1 * Y1;
        Loudness_t1(i) = Y1;
        
        if i < NumSamplesLevel - 1
            Xd = (N_mat(i) - X0) / DecFactorLevel;
            for j = 1:DecFactorLevel
                X0 = X0 + Xd;
                Y1 = B0 * X0 + A1 * Y1;
            end
        end
    end
    
    % 1st order low-pass B
    Tau = 70e-3;
    A1 = exp(-1 / (SampleRateLevel * DecFactorLevel * Tau));
    B0 = 1 - A1;
    Y1 = 0;
    
    for i = 1:NumSamplesLevel
        X0 = N_mat(i);
        Y1 = B0 * X0 + A1 * Y1;
        Loudness_t2(i) = Y1;
        if i < NumSamplesLevel - 1
            Xd = (N_mat(i) - X0) / DecFactorLevel;
            for j = 1:DecFactorLevel
                X0 = X0 + Xd;
                Y1 = B0 * X0 + A1 * Y1;
            end
        end
    end
    
    % combine the filters
    
    for i = 1:NumSamplesLevel
        Loudness(i) = (0.47 * Loudness_t1(i)) + (0.53 * Loudness_t2(i));
    end
    
    % Decimate signal for decreased computation time by factor of 24 (fs =
    % 2 Hz)
    
    Total_Loudness = zeros(NumSamplesLoudness,1);
    sC = 1;
    for i = 1:NumSamplesLoudness
        Total_Loudness(i) = Loudness(sC);
        sC = sC+DecFactorLoudness;
    end
    
    ns_dec = zeros(NumSamplesLoudness,240);
    sC = 1;
    for i = 1:NumSamplesLoudness
        ns_dec(i,:) = ns(sC,:);
        sC = sC+DecFactorLoudness;
    end
    
    %% **********************************************************************
    %  Compute loudness level - conversion from sone to phon
    % ***********************************************************************
    
    LN = 40 * Total_Loudness.^0.35;
    LN( Total_Loudness>=1 ) = 40 + 10*log2( Total_Loudness( Total_Loudness>=1 ) );
    LN( LN < 0 ) = 0;
    LN( LN < 3 ) = 3;

    %% **********************************************************************
    %  output struct for time-varying signals
    % ***********************************************************************
    
    OUT.barkAxis=barkAxis; % bark vector
    OUT.time=(0:length(Total_Loudness)-1)' * 2e-3; % time vector of the final loudness calculation, in seconds
    OUT.time_insig=(0 : length(audio)-1) ./ fs;  % time vector of the audio input, in seconds
    OUT.InstantaneousLoudness=Total_Loudness; % Time-varying Loudness, in sone
    OUT.SpecificLoudness=Spec_N; % time-averaged specific loudness (sone/Bark)
    OUT.InstantaneousSpecificLoudness=ns_dec; % specific loudness (sone/Bark) vs time
    OUT.InstantaneousLoudnessLevel=LN ; % Time-varying Loudness level, in phon
    OUT.InstantaneousSPL=10.*log10(sum(10.^(ThirdOctaveLevel(:,1:end)./10),2)); % total SPL (1/3 octave bands) for each time step, in dBSPL
    
    % statistics from Time-varying Loudness (sone)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector
    
    OUT.Nmax = max(Total_Loudness(idx:end));
    OUT.Nmin = min(Total_Loudness(idx:end));
    OUT.Nmean = mean(Total_Loudness(idx:end));
    OUT.Nstd = std(Total_Loudness(idx:end));
    OUT.N1 = prctile(Total_Loudness(idx:end),99);
    OUT.N2 = prctile(Total_Loudness(idx:end),98);
    OUT.N3 = prctile(Total_Loudness(idx:end),97);
    OUT.N4 = prctile(Total_Loudness(idx:end),96);
    OUT.N5 = prctile(Total_Loudness(idx:end),95);
    OUT.N10 = prctile(Total_Loudness(idx:end),90);
    OUT.N20 = prctile(Total_Loudness(idx:end),80);
    OUT.N30 = prctile(Total_Loudness(idx:end),70);
    OUT.N40 = prctile(Total_Loudness(idx:end),60);
    OUT.N50 = median(Total_Loudness(idx:end));
    OUT.N60 = prctile(Total_Loudness(idx:end),40);
    OUT.N70 = prctile(Total_Loudness(idx:end),30);
    OUT.N80 = prctile(Total_Loudness(idx:end),20);
    OUT.N90 = prctile(Total_Loudness(idx:end),10);
    OUT.N_ratio=prctile(Total_Loudness(idx:end),5)/prctile(Total_Loudness(idx:end),95); % ratio between N5/N95 (1.1 (stationary)> N_ratio>1.1 (time varying)
    
    %% **********************************************************************
    % show plots (time-varying)
    % ***********************************************************************
    
    if show == true
        
        figure('name','Loudness analysis (time-varying)',...
            'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
        
        xmax = OUT.time(end);
        
        % plot input signal
        subplot( 2, 6, [1,2])
        plot( OUT.time_insig, audio); hold on;
        %         a=yline(rms(audio),'k--'); % plot( OUT.time_insig,rms(audio).*ones(length(OUT.time_insig)),'k--');
        %         legend(a,sprintf('$p_{\\mathrm{rms}}=$%g (Pa)',rms(audio)),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
        ax = axis; axis([0 xmax max(audio)*-2 max(audio)*2]);
        title('Input signal','Interpreter','Latex');
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('Sound pressure, $p$ (Pa)','Interpreter','Latex'); %grid on;
        
        % plot instantaneous sound pressure level (dBSPL)
        subplot( 2, 6, [3,4])
        plot(linspace(0,OUT.time_insig(end),length(OUT.InstantaneousSPL)), OUT.InstantaneousSPL);
        ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
        title('Instantaneous overall SPL (1/3 octave)','Interpreter','Latex');
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex'); grid on;
        
        % plot instantaneous loudness level (phon)
        subplot( 2, 6, [5,6])
        plot( OUT.time, abs(OUT.InstantaneousLoudnessLevel));
        ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
        title('Instantaneous loudness level','Interpreter','Latex');
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('Loudness level, $L_{\mathrm{N}}$ (phon)','Interpreter','Latex'); grid on;
        
        % plot instantaneous loudness (sone)
        subplot( 2, 6, [7,8])
        plot( OUT.time, Total_Loudness);
        ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
        title('Instantaneous loudness','Interpreter','Latex');
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('Loudness, $N$ (sone)','Interpreter','Latex'); grid on;
        
        % plot specific loudness (sone/bark)
        subplot( 2, 6, [9,10])
        plot( OUT.barkAxis, OUT.SpecificLoudness);
        ax = axis; axis([0 24 ax(3) ax(4)*1.1]);
        title('Time-averaged specific loudness','Interpreter','Latex');
        xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');
        ylabel('Specific loudness, $N^{\prime}$ ($\mathrm{sone}/\mathrm{Bark}$)','Interpreter','Latex'); grid on;
        
        % plot instantaneous specific loudness (sone/bark)
        subplot( 2, 6, [11,12])
        [xx,yy]=meshgrid(OUT.time,OUT.barkAxis);
        pcolor(xx,yy,OUT.InstantaneousSpecificLoudness');
        shading interp; colorbar; axis tight;
        
        title('Instantaneous specific loudness','Interpreter','Latex');
        xlabel('Time, $t$','Interpreter','Latex');
        ylabel(colorbar, 'Specific Loudness, $N^{\prime}$ ($\mathrm{sone}/\mathrm{Bark}$)','Interpreter','Latex');
        
        %freq labels
        ax = gca;
        set(ax,'YTick',[0 4 8 12 16 20 24]);
        ylabel('Critical band, $z$ (Bark)','Interpreter','Latex');
        
        set(gcf,'color','w');
        
    end
    
elseif method==0 || method==1
    %% **********************************************************************
    %  output struct for stationary signals
    % ***********************************************************************
    
    if method==1 % stationary from audio signal
        OUT.time_insig=(0 : length(audio)-1) ./ fs;  % time vector of the audio input, in seconds
    else
    end
    
    OUT.barkAxis=(1:240)/10; % bark vector
    OUT.SpecificLoudness=Spec_N; % time-averaged specific loudness (sone/Bark)
    OUT.Loudness=N; % loudness (sone)
    OUT.LoudnessLevel=LN ; % loudness level (phon)
    OUT.TimeAveragedSPL=10.*log10(sum(10.^(ThirdOctaveLevel(:,1:end)./10),2)); % total SPL (1/3 octave bands) for each time step, in dBSPL
    
    %% **********************************************************************
    % show plots (stationary)
    % ***********************************************************************
    
    if show == true
        
        figure('name','Loudness analysis (stationary)',...
            'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
        
        xmax = OUT.time_insig(end);
        
        % plot input signal
        subplot( 2, 1, 1)
        plot( OUT.time_insig, audio); hold on;
        a=yline(rms(audio),'k--'); % plot( OUT.time_insig,rms(audio).*ones(length(OUT.time_insig)),'k--');
        legend(a,sprintf('$p_{\\mathrm{rms}}=$%g (Pa) \n SPL=%g (dBSPL)',rms(audio),20.*log10(rms(audio)/2e-5)),'Location','NorthEast','Interpreter','Latex'); %legend boxoff
        ax = axis; axis([0 xmax max(audio)*-2 max(audio)*2]);
        title({sprintf('Input signal')},'Interpreter','Latex');
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('Sound pressure, $p$ (Pa)','Interpreter','Latex'); %grid on;
        
        % plot specific loudnes (sone/bark)
        subplot( 2, 1, 2)
        plot( OUT.barkAxis, OUT.SpecificLoudness );
        
        annotation('textbox',...
            [0.774500000000006 0.382222222222222 0.113999999999998 0.0477777777777778],...
            'String',{sprintf('Loudness, $N$=%g (sone) \nLoudness level, $L_{\\mathrm{N}}$=%g (phon)',N,LN)},...
            'Interpreter','latex',...
            'BackgroundColor',[1 1 1]);
        
        ax = axis; axis([0 24 ax(3) ax(4)*1.1]);
        title('Specific loudness','Interpreter','Latex');
        xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');
        ylabel('Specific loudness, $N^{\prime}$ ($\mathrm{sone}/\mathrm{Bark}$)','Interpreter','Latex'); grid on;
        
        set(gcf,'color','w');
        
    end
    
end

end

%**************************************************************************
% Copyright (c) <2015>, <Ella Manor>
% All rights reserved.
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
%**************************************************************************