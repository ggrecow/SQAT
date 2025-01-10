function OUT = FluctuationStrength_Osses2016(insig,fs,method,time_skip,show,struct_opt)
% function OUT = FluctuationStrength_Osses2016(insig,fs,method,time_skip,show,struct_opt)
%
%  This function calculates the fluctuation strength using the model 
%    developed by: [1] Osses, A., Garcia A., and Kohlrausch, A.. 
%    "Modelling the sensation of fluctuation strength." Proceedings of 
%    Meetings on Acoustics 22 ICA. Vol. 28, 050005. doi:10.1121/2.0000410
%
%  Reference signal: 60 dBSPL 1 kHz tone 100% modulated at 4 Hz should yield 1 vacil.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   insig : [Nx1] array
%   insig is a monophonic calibrated audio signal (Pa)
%
%   fs : sampling frequency (Hz) - Defaults of 48 kHz or 44.1 kHz (pre-computed
%                                    filters.
%   method : integer
%   method=0, stationary analysis - window size=length(insig) (s) kind of 
%             an rms value
%   method=1, time_varying analysis - window size=2 (s)
%             NOTE: if the signal's length is smaller than 2s, the analysis
%             is automatically changed to method=0
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
%   struct_opt: struct where some specific model parameters can
%   be set to a different value. If not specified, the default values are used. 
%   Currently, the only parameter that can be changed is the a0
%   (outer and middle ear transmission factor). For that, the struct <struct_opt> needs to 
%   contain the parameter <a0_type>, meaning <struct_opt.a0_type = 'string_input'> 
%   should be defined with one of the following <string_input>:
%
%  string_input = 'fluctuationstrength_osses2016' :  A simplified a0 factor can be adopted if 
%  <a0_type> is set to this option, where the ear canal resonance of Fastl's a0 curve is removed. 
%  In other words, the a0 curve is roughly approximated as a low-pass filter. Although not 
%  explicitly stated by Osses et al. 2016 (doi: 10.1121/2.0000410), the simplified a0 transmission 
%  curve leads to very similar results during the validation of their fluctuation strength algorithm.
%  This is the default which is used if no input is given at all, as defined by the model's author
%  
%  string_input = 'fastl2007' : uses the transmission factor for free-field,
%  according to Fig 8.18 (page 226) in Fastl & Zwicker Book, Psychoacoustics: facts and
%  models 3rd edition (doi: 10.1007/978-3-540-68888-4)
%
% OUTPUT:
%   OUT : struct containing the following fields
%
%       * InstantaneousFluctuationStrength: instantaneous fluctuation 
%           strength (vacil) vs time
%       * InstantaneousSpecificFluctuationStrength: specific fluctuation 
%           strength (vacil/Bark) vs time and Bark scale
%       * TimeAveragedSpecificFluctuationStrength: time-averaged specific 
%           fluctuation strength (vacil/Bark) vs Bark scale
%       * barkAxis : vector of Bark band numbers used for specific 
%           fluctuation strength computation
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousFS
%         ** FSmean : mean value of InstantaneousFS (vacil)
%         ** FSstd : standard deviation of InstantaneousFluctuationStrength
%              (vacil)
%         ** FSmax : maximum of InstantaneousFluctuationStrength (vacil)
%         ** FSmin : minimum of InstantaneousFluctuationStrength (vacil)
%         ** FSx : fluctuation strength value exceeded during x percent of the time (vacil)
%
% Original file name: FluctuationStrength_TUe.m from 
%   https://github.com/aosses-tue/mb/tree/master/FluctuationStrength_TUe (accessed 04/03/2020)
%
% Author: Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Author: Rodrigo Garcia, HTI, TU/e, the Netherlands, 2014-2016
% Author: Gil Felix Greco, Braunschweig 04.03.2020 - Modifications
%     1) includes resampling to 44100 Hz, which is preferible because it 
%        takes less time to compute than 48 kHz because of the filtering 
%        process of IIR filters for modeling the Hweigth parameter
%     2) include possibility to choose method (stationary or time-varying) 
%        which affects the window size
% Author: Alejandro Osses, 10/05/2023. Appropriate scaling for the specific 
%            fluctuation strength.
% Author: Alejandro Osses, 11/05/2023. Moving TerhardtExcitationPatterns_v3, 
%            Get_Bark to the private folder (old il_* functions)
% Author: Alejandro Osses, 13/11/2024. Included <struct_opt> input to allow
%            for changing the a0 transmission factor. the a0 transmission
%            factor were moved to the <utilities> folder of the toolbox as
%            standalone functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help FluctuationStrength_Osses2016;
    return;
end

if nargin < 5
    % Default for show
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end

if nargin < 6
    struct_opt = [];
end

if size(insig,2)~=1 % if the insig is not a [Nx1] array
    insig=insig';   % correct the dimension of the insig
end

%% Resampling audio to 44.1 kHz or 48 kHz
if ~(fs == 44100 || fs == 48000)
    gcd_fs = gcd(44100,fs); % greatest common denominator
    insig = resample(insig,44100/gcd_fs,fs/gcd_fs);
    fs = 44100;
end

if ~isfield(struct_opt,'a0_type')
    struct_opt.a0_type = 'fluctuationstrength_osses2016'; % this is the default of this model
end

%% Checking which method
if method==1 % 'time_varying'
    
    % This is the default from the original authors.
    time_resolution = 2;  % window length fixed in 2s (Osses et al., 2016)
    N=round(fs*time_resolution);
    
    if N>=length(insig) % if the signal's length is smaller than the window size, force method==0
        warning('The signal is shorter than 2 seconds. The analysis will be automatically changed to ''stationary'', i.e. method=0 and window size=length(insig). This analysis window may lead to inaccurate fluctuation-strength estimates, especially if the modulation components are low (below 10 Hz).');
        method=0;
    end
end

if method==0 %'stationary'
    N=length(insig); %  window size (N)=length(signal) kind of an rms value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
model_par = il_Get_fluctuation_strength_params(N,fs);
model_par.debug = 'none';

% model_par = ef(model_par,'window_type','cosine');

t_b = ( 1:length(insig) )/fs;

overlap = round(0.9*N);
insig = buffer(insig,N,overlap,'nodelay');
t_b     = buffer(t_b,N,overlap,'nodelay');
nFrames = size(insig,2);
fluct   = zeros(1,nFrames); % Memory allocation

%% ei = peripheral_stage(insig,fs,N);
% 1. Cosine window:

window = ones(N,1);
attackrelease = 50;

window = il_Do_cos_ramp(window,fs,attackrelease,attackrelease);

for iFrame = 1:nFrames
    
    signal = insig(:,iFrame);
    t(iFrame,1) = t_b(1,iFrame);
    
    % Apply window to frame
    signal = transpose(window .* signal);

    %% 2. Peripheral stages
    % 2.1 Peripheral hearing system (transmission factor a0)
    %     (see 'model_par.a0_in_time' == 1, in _debug version):
    %
    % 4096th order FIR filter:
    signal = il_PeripheralHearingSystem_t(signal,fs,struct_opt); 
    
    % 2.2 Excitation patterns
    %     (see model_par.filterbank == 'terhardt', in _debug version):
    
    dBFS = 94; % corresponds to 1 Pa (new default in SQAT)
    % dBFS = 100; % unit amplitude corresponds to 100 dB (AMT Toolbox 
                  % convention, default by the original authors)
    ei   = TerhardtExcitationPatterns_v3(signal,fs,dBFS);
    dz   = 0.5; % Barks, frequency step
    z    = 0.5:dz:23.5; % Bark
    fc   = bark2hz(z);
    flow = bark2hz(z-.5); flow(1) = 0.01;
    fup  = bark2hz(z+.5);
    BWHz = fup - flow;
        
    %% 3. Modulation depth (estimation)
    [mdept,hBPi] = il_modulation_depths(ei,model_par.Hweight);
    
    %% 4. Cross-correlation coefficient:
    %     (see model_par.dataset == 0, in _debug version)
                    
    % % here cross-correlation is computed before band-pass filtering:
    % Ki = il_cross_correlation(inoutsig); % with hBPi Ki goes down but not as much as 'it should'
    Ki = il_cross_correlation(hBPi);
    [fi_,mdept,kp,gzi] = il_specific_fluctuation(mdept,Ki,model_par);

    kp_fr(iFrame,:)= kp;
    gzi_fr(iFrame,:) = gzi;
    md_fr(iFrame,:) = mdept;
    fi(iFrame,:)  = model_par.cal * fi_;
    fluct(iFrame) = dz*sum(fi(iFrame,:)); % total fluct = integration of the specific fluct. strength pattern
    
end


%% ************************************************************************
% output struct
% *************************************************************************

% main output results
OUT.InstantaneousFluctuationStrength = fluct;             % instantaneous fluctuation strength
OUT.InstantaneousSpecificFluctuationStrength = fi;        % time-varying specific fluctuation strength
OUT.TimeAveragedSpecificFluctuationStrength = mean(fi,1); % mean specific fluctuation strength
OUT.time = t;                                             % time
OUT.barkAxis = transpose(z) ;                             % critical band rate (for specific fluctuation strength)
OUT.dz = dz;

%% Fluctuation Strength statistics based on InstantaneousFS:

[~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector

OUT.FSmax = max(fluct(idx:end));
OUT.FSmin = min(fluct(idx:end));
OUT.FSmean = mean(fluct(idx:end));
OUT.FSstd = std(fluct(idx:end));
OUT.FS1 = get_percentile(fluct(idx:end),1);
OUT.FS2 = get_percentile(fluct(idx:end),2);
OUT.FS3 = get_percentile(fluct(idx:end),3);
OUT.FS4 = get_percentile(fluct(idx:end),4);
OUT.FS5 = get_percentile(fluct(idx:end),5);
OUT.FS10 = get_percentile(fluct(idx:end),10);
OUT.FS20 = get_percentile(fluct(idx:end),20);
OUT.FS30 = get_percentile(fluct(idx:end),30);
OUT.FS40 = get_percentile(fluct(idx:end),40);
OUT.FS50 = median(fluct(idx:end));
OUT.FS60 = get_percentile(fluct(idx:end),60);
OUT.FS70 = get_percentile(fluct(idx:end),70);
OUT.FS80 = get_percentile(fluct(idx:end),80);
OUT.FS90 = get_percentile(fluct(idx:end),90);
OUT.FS95 = get_percentile(fluct(idx:end),95);

%% plots

if show == true && method==1
    
    figure('name','Fluctuation strength analysis',...
           'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
       
    % Time-varying Fluctuation Strength
    subplot(2,2,1:2)
    
    plot(t,fluct,'r-');
    
    title ('Instantaneous fluctuation strength','Interpreter','Latex');
    xlabel('Time (s)','Interpreter','Latex');
    ylabel('Fluctuation strength, $\mathrm{FS}$ (vacil)','Interpreter','Latex');
    
    % Time-averaged Fluctuation strength as a function of critical band
    subplot(2,2,3)
    
    plot(transpose(z), mean(fi,1),'r-');
    
    title('Time-averaged specific fluctuation strength','Interpreter','Latex');
    xlabel('Critical band, $z$ (Bark)','Interpreter','Latex');
    ylabel('Specific fluctuation strength, $\mathrm{FS}^{\prime}$ (vacil/Bark)','Interpreter','Latex');
    
    % Specific fluctuation strength spectrogram
    subplot(2,2,4)
    
    [xx,yy]=meshgrid(t,z); 
    pcolor(xx,yy,transpose(fi)); 
    shading interp; colorbar; axis tight; 
    
    set(gca,'YDir','normal');
    title('Instantaneous specific fluctuation strength','Interpreter','Latex');
    xlabel('Time (s)','Interpreter','Latex');
    ylabel('Critical band, $z$ (Bark)','Interpreter','Latex');
    ylabel(colorbar, 'Specific fluctuation strength, $\mathrm{FS}^{\prime}$ ($\mathrm{vacil}/\mathrm{Bark}$)','Interpreter','Latex');
           
    set(gcf,'color','w')
           
end


disp('')
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdept,hBPi,ei] = il_modulation_depths(ei,Hweight)
    
[Chno,Nc] = size(ei);
mdept   = zeros(1,Chno);

ei      = transpose(abs(ei));
h0      = mean(ei);
ei      = ei - repmat(h0,Nc,1);

if ~isnumeric( Hweight )
    % older versions of MATLAB
    hBPi = filter(Hweight,ei); % getting the envelopes
else
    hBPi = sosfilt(Hweight,ei);
end

try 
    % In case LTFAT toolbox is installed (overloads rms from signal processing toolbox)
    hBPrms = rms(hBPi,'dim',1); 
catch
    % uses the default rms calculation from the signal processing toolbox
    hBPrms = rms(hBPi,1);
end
hBPi = transpose(hBPi);

idx = find(h0>0);
mdept(idx) = hBPrms(idx)./h0(idx);

idx = find(h0==0);
mdept(idx) = 0;

idx = find(h0<0);
if length(idx) ~= 0
    error('There is an error in the algorithm')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ki = il_cross_correlation(hBPi)
[Chno,~] = size(hBPi);

ki = zeros(2,Chno);
for k=1:Chno-2
    try
        cfac = cov(hBPi(k,:),hBPi(k+2,:));
    catch
        error('You do not have the function cov (stats toolbox). Contact me at ale.a.osses@gmail.com to solve this problem');
    end
    den  = diag(cfac);
    den  = sqrt(den*den');

    if den(2,1) > 0 % Pearson correlation
        ki(1,k) = cfac(2,1)/den(2,1);
    elseif den(2,1) == 0
        ki(1,k) = 0;
    else
        warning('Cross correlation factor less than 1')
        ki(1,k) = 0;
    end
end

try
    ki(1,Chno-1) = interp1([0 0.5], ki(1,Chno-3:Chno-2),1,'spline');
    ki(1,Chno  ) = interp1([0 0.5], ki(1,Chno-2:Chno-1),1,'spline');
    ki(2,2     ) = interp1([0.5 1], ki(1,3:4),0,'spline');
    ki(2,1     ) = interp1([0.5 1], ki(1,2:3),0,'spline');
catch
    ki(1,Chno-1) = ki(1,Chno-2);
    ki(1,Chno  ) = ki(1,Chno-2);
    ki(2,1) = ki(1,3);
    ki(2,2) = ki(1,3);
end

ki(2,3:Chno) = ki(1,1:Chno-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fi,mdept,kp,gzi] = il_specific_fluctuation(mdept,Ki,model_par,dataset)

if nargin <4
    dataset = 0;
end

gzi = model_par.gzi; 
p_g = model_par.p_g;
p_m = model_par.p_m;
p_k = model_par.p_k;

Chno = length(gzi);

fi = zeros(1,Chno);

switch dataset
    case {0,90,99}

        % Version 3: % Improves approximation for FM tones
        thres = 0.7;
        idx = find(mdept>thres);
        exceed = mdept(idx)-thres;
        mdept(idx) = thres+(1-thres)*exceed;
        md    = min(mdept,ones(size(mdept)));

    case 1
        md    = min(mdept,ones(size(mdept)));
        md    = mdept-0.1*ones(size(mdept));
        md    = max(mdept,zeros(size(mdept)));
end

kp     = Ki(1,:).*Ki(2,:);
kpsign = (sign(kp));
kp     = abs(kp);

switch dataset
    case {0,90,99}
        fi = gzi.^p_g .* md.^p_m .* (kp.^p_k).*kpsign;
    case 1
        fi = gzi.^p_g*md.^p_m*kp.^p_k;
    otherwise
        error('Dataset does not include the calculation of fi')
end

mdept = md;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = il_Get_fluctuation_strength_params(N,fs)
% function params = il_Get_fluctuation_strength_params(N,fs)

params         = struct;
params.fs      = fs;
params.N       = N;
params.Chno    = 47;
params.debug   = 'none';

% dataset = 0; % 0 = Approved version
params.window_type = 'cosine';
params.filterbank = 'terhardt'; 
params.p_g     = 1; 
params.p_m     = 1.7;  
params.p_k     = 1.7; % warning('Temporal value') 
params.a0_in_time = 1;
params.a0_in_freq = ~params.a0_in_time;

params.cal     = 0.4980; % this value is twice 0.2490 on 15/06/2016
params.bIdle   = 1; % v5
%%%        

params.Hweight = Get_Hweight_fluctuation(fs);
params.Hweight = Get_Hweight_fluctuation(fs);
params.gzi     = il_Get_gzi_fluctuation(params.Chno);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_Do_cos_ramp(insig,fs,attack_ms,release_ms)
%  Applies a cosine ramp with attack and release times given in [ms]

sig_len = length(insig);
r =  cos_ramp(sig_len,fs,attack_ms,release_ms);
try
    outsig = transpose(r).*insig;
catch
    outsig = r.*insig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_PeripheralHearingSystem_t(insig,fs,struct_opt)
% function outsig = il_PeripheralHearingSystem_t(insig,fs,struct_opt)
% 
% Applies the effect of transmission from free field to the cochlea to a
% given signal. Time domain version.
% 
% Inputs:
% insig: The signal to process. insig has to be a row vector.
% fs: Sampling frequency,
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 2^12; % FIR filter order 

switch struct_opt.a0_type
    case 'fluctuationstrength_osses2016'
        B = calculate_a0(fs,K,'fluctuationstrength_osses2016');
    case 'fastl2007'
        B = calculate_a0(fs,K,'fastl2007');
    otherwise
        % Choosing the default:
        B = calculate_a0(fs,K,'fluctuationstrength_osses2016');
end

outsig = filter(B,1,[insig zeros(1,K/2)]);
outsig = outsig(K/2+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gzi = il_Get_gzi_fluctuation(Chno)
% function gzi = il_Get_gzi_fluctuation(Chno)
%
% Returns gzi parameters using the specified number of channels.

Chstep = 0.5;
    
% Hz:   100 250   519   717 926 1084 1255 1465 1571   1972 2730 4189   15550
g0 = [0,  1,  2.5,  4.9,6.5,  8,   9,  10,  11,  11.5,  13,  15,  17.5,   24;
      1,  1,  1  ,  1  ,1  ,  1,   1,   1,   1,   1  ,   1, 0.9,   0.7, 0.5];
g0 = transpose(g0);

gzi = interp1(g0(:,1),g0(:,2),(1:Chno)*Chstep);
gzi(isnan(gzi)) = g0(end,2); % 0
   
%**************************************************************************
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
%
%**************************************************************************