function [b, a] = Gen_weighting_filters(fs, weightingType)
% function [b, a] = Gen_weighting_filters(fs, weightingType)
%
% 1. Description:
%		Generates the weighting filters with the following poles:
%	f1 =    20.5990 Hz -> w1 =   129.4273... rad/s % low frequency pole, see e.g. IEC 61672-1:2013, Annex E
%	f2 =   107.6526 Hz -> w2 =   676.4015... rad/s
% f3 =   737.8622 Hz -> w3 =  4636.1251... rad/s
%	f4 = 12194.2171 Hz -> w4 = 76618.5260... rad/s % high frequency pole
% f5 =   158.5    Hz -> w5 =    995.88 rad/s (source Osses2010)
%   A-curve uses: f1(x2),2(x1),3(x1),4(x2)
%   B-curve uses: f1(x2),5(x1),4(x2)
%   C-curve uses: f1(x2),4(x2)
%
% 2. Stand-alone example:
%       fs = 44100; % Hz
%       Gen_weighting_filters(fs,'A');
%
%       fs = 44100; % Hz
%       Gen_weighting_filters(fs,'C');
%
% 3. Additional info:
%       Tested cross-platform: No
%       References: Osses2010 (thesis), section 2.2.6; IEC 61672-1:2002;
%       IEC 61672-1:2013
%
% Author: Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Date: 14/07/2016
% Date: 03/04/2025 (Mike Lotinga, University of Salford: Use MATLAB bilinear
%       function and precise frequencies)
% Date: 29/11/2025 (Alejandro Osses: Adding GNU Octave compatibility)
% MATLAB / GNU Octave compatible: Yes / Yes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 16384; % arbitrary value
% The frequency of each Fourier bin
f = (0:fs/2/(N-1):fs/2);

% Determine the weighting filter frequency responses convienient to
% accurately set the desired filter orders (n,m)

w1 =   129.42731565506293; % rad/s
w2 =   676.4015402329549; % rad/s
w3 =  4636.125126885012; % rad/s
w4 = 76618.52601685845; % rad/s
w5 =   995.88; % rad/s
switch weightingType
    case 'A' % A-weighting filter
        K = 7.3901e9;
        n = 4; % at most we need a 4th order filter
        m = n;

        zrs =  [0; 0; 0; 0];
        pls = -[w1; w1; w2; w3; w4; w4];

    case 'B' % B-weighting filter
        K = 5.9862e9;
        n = 4; % at most we need a 4'th order filter
        m = n;

        zrs =  [0; 0; 0];
        pls = -[w1; w1; w5; w4; w4];

    case 'C' % C-weighting filter
        K = 5.9124e9;
        n = 2; % at most we need a 4'th order filter
        m = 4;

        zrs =  [0; 0];
        pls = -[w1; w1; w4; w4];

    case 'D' % D-weighting filter
        K = 91103.49;
        n = 3; % at most we need a 4'th order filter
        m = 4;

        zrs = [0; roots([1 6532 4.0975e7])];
        pls = [-1776.3; -7288.5; roots([1 21514 3.8836e8])];

    case 'Z' % un-weighted
        a = 1;
        b = 1;
        return

    otherwise % unknown request
        error(['weightingType=''' weightingType ''' is unknown. Options are ''A'', ''B'', ''C'' or ''D'''])
end

% Generate the filter

% Use the bilinear transformation to discretize the above transfer function.

warnState = warning('off', 'MATLAB:nearlySingularMatrix');

bOctave = exist('OCTAVE_VERSION','builtin');
if bOctave
    % Loading signal package to get bilinear.m
    pkg load signal
end

if ~bOctave
    % In MATLAB:
    [Zd, Pd, Kd] = bilinear(zrs, pls, K, fs);
else
    % In Octave, the fourth input argument is T=1/fs, instead of fs
    [Zd, Pd, Kd] = bilinear(zrs, pls, K, 1/fs);
end
[b, a] = zp2tf(Zd, Pd, Kd);

warning(warnState);

if nargout == 0
    if (2*f ~= fs) % correct small frequency error on the last fourier sample.
        f(end) = fs/2;
    end

    [H, w] = freqz(b,a,N/2);
    f = w/pi*(fs/2);
    figure;
    semilogx(f,20*log10(abs(H)));
    grid on;
    xlim([20 fs/2])
    ylim([-60 12])
    fticks = [32 63 125 250 500 1000 2000 4000 8000 16000];
    set(gca,'XTick',fticks);
    set(gca,'XTickLabel',fticks);
    xlabel('Frequency [Hz]');
    ylabel('Attenuation [dB]');
    title(sprintf('Frequency response %s-weighting curve',weightingType));
end

