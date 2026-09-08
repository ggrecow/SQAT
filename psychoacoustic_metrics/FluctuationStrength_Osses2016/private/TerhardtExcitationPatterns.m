function [ei,ei_f,freq,clamp] = TerhardtExcitationPatterns(insig,fs,dBFS)
% function [ei,ei_f,freq,clamp] = TerhardtExcitationPatterns(insig,fs,dBFS)
%
% Excitation patterns of one analysis window of the fluctuation strength
% model: level calibration, transform to the frequency domain and the
% call to the shared critical-band filterbank of Terhardt, followed by the
% renormalisation of the summed excitation to the energy of the input. The
% filterbank itself lives in utilities/Terhardt_filterbank.m, shared with
% Roughness_Daniel1997; its parameters come from
% utilities/Terhardt_filterbank_params.m. The a0 transmission factor is
% applied by the caller, in the time domain (il_PeripheralHearingSystem_t).
%
% The optional fourth output <clamp> reports whether the Terhardt upper
% slope was clamped to zero for any component (see Terhardt_filterbank.m);
% when it is requested, no warning is raised here.
%
% Author: Alejandro Osses, extracted from FluctuationStrength_Osses2016.m on 12/05/2023
% Modified: Mike Lotinga, May 2025 (parallelised code to omit loop over
% whichL for improved performance)
% Modified: Sergio Aguirre, September 2026 (masked both sides of the S2
% assignment, which crashed above a component level of about 121 dB, and
% report in <clamp> when the upper slope is clamped to zero, warning only
% when that output is not requested); the filterbank core and its
% parameter builder moved to the utilities folder, this file keeps the
% calibration and the renormalisation. Results are bitwise identical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dBFS = 94; % SQAT toolbox convention
end
corr = dBFS + 3;

dB2calibrate = rmsdb(insig)+dBFS;

% General parameters
params = Terhardt_filterbank_params(length(insig),fs);

dfreq = fs/params.N;
freq = dfreq*(1:params.N); % params.freqs is the same array, but starts at bin N0

% Transforms input signal to frequency domain
corr_factor = 10^(corr/20); % il_From_dB(corr)
insig = corr_factor*fft(insig)/params.N; % 3 dB added to adjust the SPL values to be put into slope equations

% Shared filterbank (utilities). <info> is requested, so the filterbank does
% not warn; the warning of this wrapper is below
if nargout >= 2
    [ei, info, ei_f] = Terhardt_filterbank(insig, params);
else
    [ei, info] = Terhardt_filterbank(insig, params);
end
clamp = info.clamp;

if info.n_components == 0
    ei_f = zeros(params.Chno, params.N);  % Return silence (ei already is)
    freq = zeros(1, params.N);  % Return silence
    return
end

if clamp.n > 0 && nargout < 4
    warning('SQAT:FluctuationStrength:TerhardtSlopeClamped', ...
        ['Terhardt upper slope clamped to zero for %d component(s) whose level ' ...
         'exceeds 120 + 1150/f dB (highest: %.1f dB at %.0f Hz). The fluctuation ' ...
         'strength of this frame is an extrapolation outside the range over which ' ...
         'the metric was validated; check the dBFS calibration.'], ...
        clamp.n, clamp.LdB, clamp.freq);
end

outsig = sum(ei,1);
gain = dB2calibrate - (rmsdb(outsig)+dBFS);
gain_factor = 10^(gain/20);
ei = gain_factor*ei;
