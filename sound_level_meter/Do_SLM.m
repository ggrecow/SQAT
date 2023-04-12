function [outsig_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS)
% function [outsig_dB, dBFS] = Do_SLM(insig,fs,weight_freq,weight_time,dBFS)
%
% 1. Description:
%       outsig is the weighted pressure in [Pa]
% 
% 2. Stand-alone example:
%   %%%%%%
%   % Example 2.1:
%   %%% Creating a sinuoid with an RMS of 60 dB SPL:
%   fs = 44100;
%   freq = 1000; % Hz
%   dur = 1; % s
%   t = 0:1/fs:dur-1/fs;
%   dBFS = 94; % full scale convention
%   target_lvl = 60; % 60 dB SPL
%   cal_factor = 10^((target_lvl-dBFS)/20);
%   A = 1*cal_factor*sqrt(2); % 1=full scale; sqrt(2)=relationship between RMS and peak value for a sinusoid
%   insig = A*sin(2*pi*freq*t); % zero phase sinusoid centred at f=1000;
%
%   %%% Obtaining its A-weighted amplitude, with a fast weighting:
%   outsig_dB = Do_SLM(insig,fs,'A','f',dBFS); % A-weighting
%   figure;
%   plot(outsig_dB); grid on;
%   xlabel('Time [samples]');
%   ylabel('Amplitude [dB(A)]');
% 
%   %%%%%%
%   % Example 2.2:
%   %%% This example assumes that the file is found on disk and that the dB 
%   % full scale convention (dBFS) is as the default:
%
%   file = 'D:\Databases\dir03-Speech\dutch\LISTman\jwz551.wav';
%   [insig, fs] = audioread(file);
%   Do_SLM(insig,fs,'Z','f'); % Z-weighting, fast response, with an automatic output figure
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also: Get_Leq
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/07/2016
% Last update on: 12/07/2016 
%               : 22/03/2023, AO: Stylised output figure and independency 
%                             of codes with respect to the LTFAT toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    weight_freq = 'A';
end

if nargin < 4
    weight_time = 'f';
end

if nargin < 5
    dBFS = 100;
end
if size(insig,2)~=1
    if size(insig,1) == 1
        fprintf('%s: Row array detected as an input to this function. The default is a column-wise array.\n', mfilename);
        insig = insig(:);
    end
end
[b,a] = Gen_weighting_filters(fs,weight_freq);

% same processing, but without AMT, as done in PsySound:
dBoffset = 0.93; % determined empirically on 13/07/2016 to obtain the same values
                 % with this implementation and the one in PsySound.
calCoeff = 10.^((dBFS+dBoffset-94)/20);
insig = calCoeff*insig; % same as: insig = setdbspl(insig,rmsdb(insig)+dBFS,'dboffset',94);

outsig = filter(b,a,insig);
outsig = il_integrator(outsig,fs,weight_time);

outsig_dB = 20*log10(abs(outsig)/2e-5);

try
    idx = find(outsig_dB<0);
    if ~isempty(idx)
        fprintf('\t%s There were %.0f samples (out of %.0f) with levels below 0 dB SPL. They were set to 0 dB SPL\n',mfilename,length(idx),length(outsig_dB));
        outsig_dB(idx) = 0;
    end
end

if nargout == 0
    Leq = Get_Leq(outsig_dB);
    
    t = (1:length(outsig_dB))/fs;
    figure;
    plot(t, outsig_dB ); grid on;
    switch weight_freq
        case 'Z'
            suff = '';
        otherwise
            suff = ['(' weight_freq ')'];
    end
    xlabel('Time [s]');
    xlim([0 max(t)]);
    
    unit = sprintf('[dB%s]',suff);
    ylabel(sprintf('Amplitude %s',unit));
    title(sprintf('Level Leq=%.1f %s',Leq,unit));
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline functions:
function outsig = il_integrator(insig,fs, weightingType)
% function outsig = il_integrator(insig,fs, weightingType)
%
% 1. Description: This function generates and applies an integration filter
%      to the input data. This code is based on an implementation taken from
%      the PsySound toolbox by Matt Flax (January 2007) and Farhan Rizwi 
%      (July 2007).
%
%      Inputs:
%         insig         - Incoming input signal, as a data vector
%         fs            - Sampling rate of the data
%         weightingType - RC time constant is 'f' fast (125 ms), 's' slow
%                       (1 s), or 'i' impulsive. The time constant for the 
%                       leaky integrator is tau. This is basically a low-pass 
%                       filter whose transfer function is:
%                      1
%        H(s) =  ---------------
%                tau s  +   1
%
% Code adapted by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter coeffecients

switch weightingType
    case 'f' % fast leak - time constant = 125 ms
        tau = 125e-3;
    case 's' % slow leak - time constant = 1 s
        tau = 1;
    case 'i'
        tau = 35e-3; % impulse
    % case 'p'
    %     tau = 50e-6;	
    % case 'l'
    %     tau = str2double(fastOrSlow(2:end));	
    % case 'r'
    %     tau = str2double(fastOrSlow(2:end));	
    otherwise
        error(['integrator: unknown leak case ' char(fastOrSlow)]);
end

% State vector
Z = [];

% Exponential term
E = exp(-1/(tau*fs));

% Filter numerator - with gain adjustment
b = 1 - E;

% Filter denominator
a = [1 -E];

% Create run function handle
% Use filter to perform the integration
outsig = filter(b, a, abs(insig), Z, 1);
