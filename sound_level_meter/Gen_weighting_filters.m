function [b,a] = Gen_weighting_filters(fs,weightingType)
% function [b,a] = Gen_weighting_filters(fs,weightingType)
% function [b,a] = il_gen_weighting_Filters(fs,weightingType)
%
% 1. Description:
%		Generates the weighting filters with the follwing poles:
%	f1 =    20.6 Hz -> w1 =   129.43 rad/s % low frequency pole, see e.g. IEC 61672-1:2002, 5.4.6, 5.4.11
%	f2 =   107.7 Hz -> w2 =   676.70 rad/s
% 	f3 =   737.9 Hz -> w3 =  4636.36 rad/s
%	f4 = 12194   Hz -> w4 = 76617.16 rad/s % high frequency pole
%  	f5 =   158.5 Hz -> w5 =   995.88 rad/s (source Osses2010)
% 	A-curve uses: f1(x2),2(x1),3(x1),4(x2)
% 	B-curve uses: f1(x2),5(x1),4(x2)
%	C-curve uses: f1(x2),4(x2)
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
%       References: Osses2010 (thesis), section 2.2.6; IEC 61672-1:2002
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 14/07/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 16384; % arbitrary value
% the frequency of each Fourier bin
f = (0:fs/2/(N-1):fs/2);

% take a perceptual sampling of the frequency domain - for design purposes
eMin   = il_freqtoaud(min(f)); % ERB min
eMax   = il_freqtoaud(max(f)); % ERB max
eScale = eMin:(eMax-eMin)/(length(f)-1):eMax; % ERB scale
fScale = audtofreq(eScale); % frequencies sample according to a linear ERB scale

% fLinear = f; % save the linear frequency scale
f = fScale;  % switch the reference frequencies to be f

s = 1i*2*pi*f; % set up the s-plane variable

% determine the weighting filter frequency responses convienient to
% accurately set the desired filter orders (n,m)

w1 =   129.43; % rad/s
w2 =   676.70; % rad/s
w3 =  4636.36; % rad/s
w4 = 76617.16; % rad/s
w5 =   995.88; % rad/s
switch weightingType
    case 'A' % A-weighting filter
        K = 7.3901e9;
        freqResp = K*s.^4./((s+w1).^2 .* (s+w2).*(s+w3).*(s+w4).^2);

        n = 4; % at most we need a 4th order filter
        m = n;

        zrs =  [0; 0; 0; 0];
        pls = -[w1; w1; w2; w3; w4; w4];

    case 'B' % B-weighting filter
        K = 5.9862e9;
        freqResp = K*s.^3./((s+w1).^2 .* (s+w5).*(s+w4).^2);

        n = 4; % at most we need a 4'th order filter
        m = n;

        zrs =  [0; 0; 0];
        pls = -[w1; w1; w5; w4; w4];

    case 'C' % C-weighting filter
        K = 5.9124e9;
        freqResp = K*s.^2./((s+w1).^2 .*(s+w4).^2);

        n = 2; % at most we need a 4'th order filter
        m = 4;

        zrs =  [0; 0];
        pls = -[w1; w1; w4; w4];

    case 'D' % D-weighting filter
        K = 91103.49;
        freqResp = K*s.*(s.^2+6532*s+4.0975e7)./((s+1773.6) .*(s+7288.5).*(s.^2+21514*s+3.8836e8));

        n = 3; % at most we need a 4'th order filter
        m = 4;

        zrs = [0; roots([1 6532 4.0975e7])];
        pls = [-1776.3; -7288.5; roots([1 21514 3.8836e8])];

	 case 'R'
		% Filter weightings from [1], pg 12
		% These are defined for 48k
		b = [1 -2 1];
		a = [1 -1.99004745483398 0.99007225036621];

		% Use direct substituition of the definition of the z-transform
		% (z=exp(s*T)) to recalculate coeffecients for a different sampling
		% rate
		% Note: This could be another option for pre-filtering

        if fs ~= 48000;
            poles = roots(a);

            % Make polynomial after fixing up the roots
            %
            % z = exp(s*T) --> s = ln(z)/T
            % s = ln(z1)/T1 = ln(z2)/T2  -->  z2 = exp(ln(z1)*T2/T1)

            a = poly(exp(log(poles)*48000/fs));

            % Note that the two zeros at 1 remain there.
            % Note also, that the negligible high frequency gain adjustment
            % is ignored.
        end
        return

    case 'Z' % un-weighted
        a = 1;
        b = 1;
        return

    otherwise % unknown request
        error(['weightingType=''' weightingType ''' is unknown. Options are ''A'', ''B'', ''C'' or ''D'''])
end

m = m+1;
n = n*2;
m = m*2;

% % the total frequency response
% totalResp = freqResp;

% generate the filter
if (2*f ~= fs) % correct small frequency error on the last fourier sample.
    f(end) = fs/2;
end

% Use the bilinear transformation to discretize the above
% transfer function.

warnState = warning('off', 'MATLAB:nearlySingularMatrix');

[Zd, Pd, Kd] = bilinear_local(zrs, pls, K, fs);
[b, a] = zp2tf(Zd, Pd, Kd);
warning(warnState);

if nargout == 0
    [H, w] = freqz(b,a,N/2);
    f = w/pi*(fs/2);
    figure;
    semilogx(f,20*log10(abs(H)));
    grid on;
    xlim([20 fs/2])
    ylim([-60 12])
    fticks = [63 125 250 500 1000 2000 4000 8000 16000];
    set(gca,'XTick',fticks);
    xlabel('Frequency [Hz] (audible range)');
    ylabel('Attenuation [dB]');
    title(sprintf('Frequency response %s-weighting curve',weightingType));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f_ERB = il_freqtoaud(f)
% function f_ERB = il_freqtoaud(f)
%
% Taken and adapted from LTFAT 2.4 (from flags.do_erb)

%%% Comment from the original code:
% There is a round-off error in the Glasberg & Moore paper, as
% 1000/(24.7*4.37)*log(10) = 21.332 and not 21.4 as they state.
% The error is tiny, but may be confusing.
% On page 37 of the paper, there is Fortran code with yet another set
% of constants:
%     2302.6/(24.673*4.368)*log10(1+freq*0.004368);
f_ERB = 9.2645*sign(f).*log(1+abs(f)*0.00437);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = audtofreq(f_ERB)

f = (1/0.00437)*sign(f_ERB).*(exp(abs(f_ERB)/9.2645)-1);
