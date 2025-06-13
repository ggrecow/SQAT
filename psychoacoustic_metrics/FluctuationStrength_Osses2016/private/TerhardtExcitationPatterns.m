function [ei,ei_f,freq] = TerhardtExcitationPatterns(insig,fs,dBFS)
% function [ei,ei_f,freq] = TerhardtExcitationPatterns(insig,fs,dBFS)
%
% Author: Alejandro Osses, extracted from FluctuationStrength_Osses2016.m on 12/05/2023
% Modified: Mike Lotinga, May 2025 (parallelised code to omit loop over
% whichL for improved performance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dBFS = 100; % AMT toolbox convention
end
corr = dBFS + 3;

dB2calibrate = rmsdb(insig)+dBFS;

% General parameters
params = il_calculate_params(insig,fs);
N01 = params.N01;
freqs = params.freqs.';

dfreq = fs/params.N;
freq = dfreq*(1:params.N); % freqs and freq are the same array, but freqs starts at bin N0

% Transforms input signal to frequency domain
corr_factor = 10^(corr/20); % il_From_dB(corr)
insig = corr_factor*fft(insig)/params.N; % 3 dB added to adjust the SPL values to be put into slope equations

% Use only samples that fall into the audible range
Lg  = abs(insig(params.qb)).';
LdB = 20*log10(Lg); % il_To_dB(Lg);

% Use only components that are above the hearing threshold
MinExcdB = params.MinExcdB(:);
whichL = find(LdB > MinExcdB);

if isempty(whichL)
    ei = zeros(params.Chno, params.N);  % Return silence
    ei_f = zeros(params.Chno, params.N);  % Return silence
    freq = zeros(1, params.N);  % Return silence
    return
end

nL = length(whichL);

% Steepness of slopes
S1 = -27;			
S2 = zeros(nL, 1);

steep = -24 - (230./freqs(whichL)) + (0.2*LdB(whichL));
S2(steep < 0) = steep;
S2 = repmat(S2, 1, params.Chno);

whichZ = zeros(nL, 2);
whichZ(:, 1)	= floor(2*params.Barkno(whichL + N01));
whichZ(:, 2)	=  ceil(2*params.Barkno(whichL + N01));

% Calculate slopes from steep values
Slopes = zeros(nL, params.Chno);
Stemp = Slopes;

Li = repmat(LdB(whichL), 1, params.Chno);

kk = zeros(nL, params.Chno);
kk1 = kk;
kk2 = kk;
delta_z = zeros(nL, params.Chno);
for l = nL:-1:1    
    for k = whichZ(l, 1):-1:1
        kk1(l, k) = k;
    end
    for k = params.Chno:-1:whichZ(l, 2)
        kk2(l, k) = k;
    end
end

kk1mask = kk1 > 0;
kk2mask = kk2 > 0;

kk(kk1mask) = kk1(kk1mask);
kk(kk2mask) = kk2(kk2mask);

zk = 0.5*(kk);
zi = repmat(params.Barkno(whichL + N01).', 1, size(zk, 2));
delta_z(kk1mask) = zi(kk1mask) - zk(kk1mask);
delta_z(kk2mask) = zk(kk2mask) - zi(kk2mask);
Stemp(kk1mask) = S1*delta_z(kk1mask) + Li(kk1mask);
Stemp(kk2mask) = S2(kk2mask).*delta_z(kk2mask) + Li(kk2mask);
maxk = max(kk, [], 'all');
MinBfRep = repmat(params.MinBf(1:maxk), nL, 1);
mask = Stemp > MinBfRep;
Slopes(mask) = 10.^(Stemp(mask)/20);

% Excitation patterns:
%   Each frequency having a level above the absolute threshold is looked at.
%   The contribution of that level (and frequency) onto the other critical
%   band levels is computed and then assigned.
ExcAmp  = zeros(max(whichL), params.Chno);
ei      = zeros(params.Chno, params.N);

for i = params.Chno:-1:1
    etmp = zeros(1,params.N);

    if i ~= 1
        ExcAmp(whichL, i) = Slopes(:, i - 1)./Lg(whichL);
    end

    if i ~= 47
        mask1 = whichZ(:, 2) > i;
        ExcAmp(whichL(mask1), i) = Slopes(mask1, i + 1)./Lg(whichL(mask1));
    end

    mask2 = whichZ(:, 2) == i;
    ExcAmp(whichL(mask2), i) = 1;

    mask3 = whichZ(:, 1) == i;
    ExcAmp(whichL(mask3), i) = 1;

    etmp(whichL + N01) = ExcAmp(whichL, i).*insig(whichL + N01).'; % for each level, the level is projected to that in the respective critical band i

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