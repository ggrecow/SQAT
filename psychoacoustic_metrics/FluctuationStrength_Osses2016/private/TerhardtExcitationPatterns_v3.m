function [ei,ei_f,freq] = TerhardtExcitationPatterns_v3(insig,fs,dBFS)
% function [ei,ei_f,freq] = TerhardtExcitationPatterns_v3(insig,fs,dBFS)
%
% Author: Alejandro Osses, extracted from FluctuationStrength_Osses2016.m on 12/05/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dBFS = 100; % AMT toolbox convention
end
corr = dBFS + 3;

dB2calibrate = rmsdb(insig)+dBFS;

% General parameters
params = il_calculate_params(insig,fs);
N01 = params.N01;
freqs = params.freqs;

dfreq = fs/params.N;
freq = dfreq*(1:params.N); % freqs and freq are the same array, but freqs starts at bin N0

% Transforms input signal to frequency domain
corr_factor = 10^(corr/20); % il_From_dB(corr)
insig = corr_factor*fft(insig)/params.N; % 3 dB added to adjust the SPL values to be put into slope equations

% Use only samples that fall into the audible range
Lg  = abs(insig(params.qb));
LdB = 20*log10(Lg); % il_To_dB(Lg);

% Use only components that are above the hearing threshold
whichL = find(LdB > params.MinExcdB);
nL     = length(whichL);

% Steepness of slopes
S1 = -27;			
S2 = zeros(1,nL);
for w = 1:nL
    idx = whichL(w);
    steep = -24 - ( 230 / freqs(idx)) + (0.2 * LdB(idx) ); 
    if steep < 0
        S2(w) = steep;
    end
end

whichZ      = zeros(2,nL);
whichZ(1,:)	= floor(2 * params.Barkno(whichL+N01));
whichZ(2,:)	=  ceil(2 * params.Barkno(whichL+N01));

% Calculate slopes from steep values
Slopes = zeros(nL,params.Chno);
Slopes_dB = nan(nL,params.Chno);

for l = 1:nL
    Li = LdB(whichL(l));
    zi = params.Barkno(whichL(l)+N01);
    
    for k = 1:whichZ(1,l)
        zk = k * 0.5;
        delta_z = zi - zk;
        Stemp =	(S1 * delta_z) + Li;
        if Stemp > params.MinBf(k)
            Slopes(l,k) = 10^(Stemp/20);
            Slopes_dB(l,k) = Stemp;
        end
    end

    for k = whichZ(2,l):params.Chno
        zk = k * 0.5;
        delta_z = zk - zi;
        Stemp = S2(l) * delta_z + Li;
        if Stemp > params.MinBf(k)
            Slopes(l,k) = 10^(Stemp/20);
            Slopes_dB(l,k) = Stemp;
        end
    end 
end

% Excitation patterns:
%   Each frequency having a level above the absolute threshold is looked at.
%   The contribution of that level (and frequency) onto the other critical
%   band levels is computed and then assigned.
ExcAmp  = zeros(nL,params.Chno);
ei      = zeros(params.Chno,params.N);
for i = 1:params.Chno
    etmp = zeros(1,params.N);
    for l = 1:nL
        N1tmp = whichL(l);
        N2tmp = N1tmp + N01;

        if i == 24
            disp('')
        end
        if whichZ(1,l) == i
            ExcAmp(N1tmp,i)	= 1;
        elseif whichZ(2,l) == i
            ExcAmp(N1tmp,i)	= 1;
        elseif whichZ(2,l) > i
            ExcAmp(N1tmp,i) = Slopes(l,i+1)/Lg(N1tmp);
        else % whichZ(1,l) < k
            ExcAmp(N1tmp,i) = Slopes(l,i-1)/Lg(N1tmp);
        end

        etmp(N2tmp) = ExcAmp(N1tmp,i)*insig(N2tmp); % for each level, the level is projected to that in the respective critical band i
    end

    if nargout >= 2
        ei_f(i,:) = 20*log10(abs(etmp));
    end
    ei(i,:) = 2*params.N*real(ifft(etmp)); % figure; plot(To_dB(abs(etmp)),'x','LineWidth',4); xlim([1950 2050])
end

outsig = sum(ei,1);
gain = dB2calibrate - (rmsdb(outsig)+dBFS);
gain_factor = 10^(gain/20);
ei = gain_factor*ei;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = il_calculate_params(x,fs)

params      = struct;
params.N    = length(x);
params.Chno = 47;

% Defines audible range indexes and frequencies
df           = fs/params.N;
N0           = round(20/df)+1; % start at 20 Hz
Ntop         = round(20e3/df)+1; % start at 20 kHz
params.N01   = N0-1;
params.qb    = N0:Ntop;
params.freqs = (params.qb-1)*df;

[params.Barkno,Bark_raw] = Get_Bark(params.N,params.qb,params.freqs);
% Loudness threshold related parameters
params.MinExcdB = il_calculate_MinExcdB(params.N01,params.qb,params.Barkno);
params.MinBf    = il_calculate_MinBf(params.N01,df,Bark_raw,params.MinExcdB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinExcdB = il_calculate_MinExcdB(N01,qb,Barkno)

HTres = [
    0		130
    0.01    70
    0.17    60
    0.8     30
    1       25
    1.5     20
    2		15
    3.3     10
    4		8.1
    5		6.3
    6		5
    8		3.5
    10		2.5
    12		1.7
    13.3	0
    15		-2.5
    16		-4
    17		-3.7
    18		-1.5
    19		1.4
    20		3.8
    21		5
    22		7.5
    23      15
    24      48
    24.5 	60
    25		130
];

MinExcdB            = zeros(1,length(qb));
MinExcdB(qb-N01)    = interp1(HTres(:,1),HTres(:,2),Barkno(qb));
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinBf = il_calculate_MinBf(N01,df,Bark,MinExcdB)
    
Cf = round(Bark(2:25,2)'/df)-N01+1;
Bf = round(Bark(1:25,3)'/df)-N01+1;  

zb      = sort([Bf Cf]);
MinBf   = MinExcdB(zb);