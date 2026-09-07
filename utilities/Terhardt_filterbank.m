function [ei,info,ei_f] = Terhardt_filterbank(insig_f,params)
% function [ei,info,ei_f] = Terhardt_filterbank(insig_f,params)
%
% Critical-band filterbank with the excitation patterns of Terhardt (1979),
% as used by the roughness model of Daniel & Weber (1997) and by the
% fluctuation strength model of Osses et al. (2016). Each spectral component
% above the hearing threshold is spread over the 47 half-Bark channels with
% a lower slope of -27 dB/Bark and an upper slope that depends on its level
% and frequency (Terhardt 1979, Eqs. 4a and 4b); the excitation of each
% channel is returned in the time domain.
%
% This function is shared by FluctuationStrength_Osses2016 (through its
% wrapper TerhardtExcitationPatterns.m, which applies the a0 transmission
% factor in the time domain, calibrates and transforms the frame) and, from
% the rewrite discussed in issue 47 on, by Roughness_Daniel1997 (which
% applies a0 on the spectrum of the frame, as in the reference code of
% Hermes). Windowing, calibration, a0 and any renormalisation of the
% output stay with the caller.
%
% INPUT:
%   insig_f : [1 x N] complex spectrum of one analysis window, with the a0
%       transmission factor already applied and calibrated so that a
%       sinusoidal component of L dB SPL has a magnitude of 10^(L/20)
%   params : struct from Terhardt_filterbank_params.m (N, Chno, N01, qb,
%       freqs, Barkno, MinExcdB, MinBf)
%
% OUTPUT:
%   ei : [Chno x N] excitation pattern of each channel in the time domain,
%       scaled as 2*N*real(ifft(.)); the scale cancels in the modulation
%       depth of both metrics
%   info : struct with
%       n_components : number of spectral components above the hearing
%           threshold (zero means the window is silent and ei is all zeros)
%       clamp : struct (n, LdB, freq) reporting how many components had
%           the upper slope clamped to zero, and the level and frequency of
%           the one farthest into the clamped regime (see below). A warning
%           is raised here only when <info> is not requested, so that a
%           caller running the filterbank once per frame can aggregate the
%           frames and warn once per call
%   ei_f : [Chno x N] level, in dB, of the excitation pattern of each
%       channel in the frequency domain (for inspection only)
%
% Author: Dik Hermes, TU/e (2000-2005), original implementation of the
%   Daniel & Weber model
% Author: Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016; extracted
%   from FluctuationStrength_Osses2016.m on 12/05/2023
% Modified: Mike Lotinga, May 2025 (parallelised code to omit loop over
%   whichL for improved performance)
% Modified: Sergio Aguirre, September 2026 (masked both sides of the S2
%   assignment, which crashed above a component level of about 121 dB;
%   report the clamping in <info>); moved to the utilities folder as the
%   filterbank shared by the two modulation metrics, with the spectrum as
%   input. The arithmetic is unchanged: FluctuationStrength_Osses2016
%   returns bitwise identical results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N01   = params.N01;
freqs = params.freqs.';

% Use only samples that fall into the audible range
Lg  = abs(insig_f(params.qb)).';
LdB = 20*log10(Lg); % il_To_dB(Lg);

% Use only components that are above the hearing threshold
MinExcdB = params.MinExcdB(:);
whichL = find(LdB > MinExcdB);
info = struct('n_components', numel(whichL), ...
              'clamp', struct('n', 0, 'LdB', [], 'freq', [])); % Terhardt upper slope clamped to zero: none so far

if isempty(whichL)
    ei = zeros(params.Chno, params.N);  % Return silence
    ei_f = zeros(params.Chno, params.N);  % Return silence
    return
end

nL = length(whichL);

% Steepness of slopes
S1 = -27;			
S2 = zeros(nL, 1);

steep = -24 - (230./freqs(whichL)) + (0.2*LdB(whichL));
% Terhardt's upper slope is defined for steep < 0 only. A component whose
% level exceeds 120 + 1150/f dB (about 121 dB at 1 kHz) gives steep >= 0,
% which would mean an excitation that does not decay towards higher
% frequencies. S2 stays at zero for such a component (flat spread), so the
% result is an extrapolation outside the range over which the metric was
% validated (60 to 70 dB SPL); a wrong dBFS is the most likely cause.
% The clamping is reported in <clamp>: number of components, and level and
% frequency of the component farthest into the clamped regime. The warning
% is raised here only when the caller does not request <info>, so that a
% caller running the filterbank once per frame (the wrapper of the
% fluctuation strength, the frame loop of the roughness) can aggregate the
% frames and warn once per call.
info.clamp.n = nnz(steep >= 0);
if info.clamp.n > 0
    [~, iw] = max(steep);
    info.clamp.LdB  = LdB(whichL(iw));
    info.clamp.freq = freqs(whichL(iw));
    if nargout < 2
        warning('SQAT:Terhardt_filterbank:SlopeClamped', ...
            ['Terhardt upper slope clamped to zero for %d component(s) whose level ' ...
             'exceeds 120 + 1150/f dB (highest: %.1f dB at %.0f Hz). The excitation ' ...
             'of this window is an extrapolation outside the range over which the ' ...
             'model was validated; check the dBFS calibration.'], ...
            info.clamp.n, info.clamp.LdB, info.clamp.freq);
    end
end
% Both sides are masked on purpose: this reproduces the element-wise guard
% of the scalar loop in TerhardtExcitationPatterns_v3.m (S2 stays 0 whenever
% steep >= 0) and keeps the two sides the same size once any component has
% steep >= 0. With the right-hand side unmasked the assignment raised a size
% mismatch error as soon as one component crossed that level.
S2(steep < 0) = steep(steep < 0);
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

    etmp(whichL + N01) = ExcAmp(whichL, i).*insig_f(whichL + N01).'; % for each level, the level is projected to that in the respective critical band i

    if nargout >= 3
        ei_f(i,:) = 20*log10(abs(etmp));
    end
    ei(i,:) = 2*params.N*real(ifft(etmp)); % figure; plot(To_dB(abs(etmp)),'x','LineWidth',4); xlim([1950 2050])
end
