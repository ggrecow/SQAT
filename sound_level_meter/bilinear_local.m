function [zd, pd, kd, dd] = bilinear_local(z, p, k, fs)
%BILINEAR Bilinear transformation with optional frequency prewarping.
%   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
%   function specified by Z, P, and K to a z-transform discrete
%   equivalent obtained from the bilinear transformation:
%
%      H(z) = H(s) |
%                  | s = 2*Fs*(z-1)/(z+1)
%
%   where column vectors Z and P specify the zeros and poles, scalar
%   K specifies the gain, and Fs is the sample frequency in Hz.
%
%   % Example: see test_Gen_weighting
%
%   Author: J.N. Little, 28/04/1987, The MathWorks, Inc.
%   Author (modification): Alejandro Osses, 21/10/2023, simplification and
%     compatibility MATLAB/Octave

zs = z(:);
ps = p(:);
ks = k(1);

sampleFreq = 2*double(fs(1));
zs = zs(isfinite(zs));  % Strip infinities from zeros
                        % Do bilinear transformation
prodzs = prod(sampleFreq-zs,1);
zd1 = (1 + zs/sampleFreq)./(1 - zs/sampleFreq);

prodps = prod(sampleFreq - ps,1);
pd = (1 + ps/sampleFreq)./(1 - ps/sampleFreq);

kd = (ks*prodzs./prodps);
zd = [zd1;-ones(length(pd)-length(zd1),1)];
% end isZeroPoleGain

% LocalWords:  prewarping Zd Kd Fs Fp th numd dend
