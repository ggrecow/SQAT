function params = Terhardt_filterbank_params(N,fs)
% function params = Terhardt_filterbank_params(N,fs)
%
% Parameters of the Terhardt critical-band filterbank (see
% Terhardt_filterbank.m) for an analysis window of N samples at the sampling
% frequency fs. The fields are:
%
%   N        : window length in samples
%   Chno     : number of half-Bark channels (47)
%   N01      : index offset of the audible range, so that bin q of the
%              spectrum corresponds to element q-N01 of the range vectors
%   qb       : bins of the audible range, 20 Hz to 20 kHz
%   freqs    : frequency of each bin of qb, (qb-1)*fs/N
%   Barkno   : critical-band rate of each bin (see Get_Bark.m)
%   MinExcdB : hearing threshold, in dB, at each bin of qb
%   MinBf    : hearing threshold at the critical-band edges and centres
%
% This function is shared by FluctuationStrength_Osses2016 (through
% TerhardtExcitationPatterns.m) and, from the rewrite discussed in issue 47
% on, by Roughness_Daniel1997.
%
% Author: Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016, as the local
%   functions il_calculate_params, il_calculate_MinExcdB and il_calculate_MinBf
%   of TerhardtExcitationPatterns.m (FluctuationStrength_Osses2016)
% Author: Sergio Aguirre, September 2026 - moved to the utilities folder
%   unchanged, so that both modulation metrics build the same parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params      = struct;
params.N    = N;
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
