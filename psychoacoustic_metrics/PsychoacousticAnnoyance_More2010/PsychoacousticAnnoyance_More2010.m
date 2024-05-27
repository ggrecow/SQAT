function OUT = PsychoacousticAnnoyance_More2010(insig,fs,LoudnessField,time_skip,showPA,show)
% function OUT = PsychoacousticAnnoyance_More2010(insig,fs,LoudnessField,time_skip,showPA,show)
%
%   This function calculates the More's modified psychoacoustic annoyance 
%   model from an input acoustic signal
%
%   The modified psychoacoustic annoyance model is according to: (page 201)
%   [1] More, Shashikant. Aircraft noise characteristics and metrics. 
%       PhD Thesis, Purdue University, 2010
%
% - This metric combines 5 psychoacoustic metrics to quantitatively describe annoyance:
%
%    1) Loudness, N (sone) - calculated hereafter following ISO 532-1:2017
%       type <help Loudness_ISO532_1> for more info
%
%    2) Sharpness, S (acum) - calculated hereafter following DIN 45692:2009
%       NOTE: uses DIN 45692 weighting function by default, please change code if
%       the use of a different withgitng function is desired).
%       type <help Sharpness_DIN45692_from_loudness>
%
%    3) Roughness, R (asper) - calculated hereafter following Daniel & Weber model
%       type <help Roughness_Daniel1997> for more info
%
%    4) Fluctuation strength, FS (vacil) - calculated hereafter following 
%       Osses et al. model, type <help FluctuationStrength_Osses2016> for more info
%
%    5) Tonality, K (t.u.) - calculated hereafter following Aures' model
%       type <help Tonality_Aures1985> for more info
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   insig : array
%   acoustic signal [1,nTimeSteps], monophonic (Pa)
%
%   fs : integer
%   sampling frequency (Hz) - preferible 48 kHz or 44.1 kHz (default by the authors and takes less time to compute)
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations
%
%   LoudnessField : integer
%   chose field for loudness calculation; free field = 0; diffuse field = 1; (used in the loudness and tonality codes)
%   type <help Loudness_ISO532_1> for more info
%
%   show : logical(boolean)
%   optional parameter, display results of loudness, sharpness, roughness, fluctuation strength and tonality
%   'false' (disable, default value) or 'true' (enable).
%
%   showPA : logical(boolean)
%   optional parameter, display results of psychoacoustic annoyance
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUTS:
%   OUT: struct
%      * include results from the psychoacoustic annoyance
%             ** InstantaneousPA: instantaneous quantity (unity) vs time
%             ** ScalarPA : PA (scalar value) computed using the percentile values of each metric.
%                           NOTE: if the signal's length is smaller than 2s, this is the only output as no time-varying PA is calculated
%             ** time : time vector in seconds
%             ** wt : tonality and loudness weighting function (not squared)
%             ** wfr : fluctuation strength and roughness weighting function (not squared)
%             ** ws : sharpness and loudness weighting function (not squared)
%
%             ** Statistics
%               *** PAmean : mean value of psychoacoustic annoyance (unit)
%               *** PAstd : standard deviation of instantaneous psychoacoustic annoyance (unit)
%               *** PAmax : maximum of instantaneous psychoacoustic annoyance (unit)
%               *** PAmin : minimum of instantaneous psychoacoustic annoyance (unit)
%               *** PAx : x percentile of the PA metric exceeded during x percent of the time
%
%      * include structs with the results from the other metrics computed
%        **  L : struct with Loudness results, type <help Loudness_ISO532_1> for more info
%        **  S : struct with Sharpness, type <help Sharpness_DIN45692_from_loudness>
%        **  R : strcut with roughness results, type <help Roughness_Daniel1997> for more info
%        ** FS : struct with fluctuation strength results, type <help FluctuationStrength_Osses2016> for more info
%        **  K : struct with tonality results, type <help Tonality_Aures1985> for more info
%
%  NOTE: 1) Input signals should be in pascal values or calibrated .wav files
%
%        2) Fluctuation strength window has length of 2s. If the signal is 
%           less than 2s long, the FS calculation will be automatically
%           changed to stationary (i.e. uses a window with length equal to 
%           signal's size). in this case, no time-varying PA is available.
%
%        3) Be aware that, because of item 2), if the signal is more than 2s 
%           long, the last 2 seconds of the input signal are LOST !!!!
%
%        4) is a best practice to compute percentile values following a 
%           time_skip (s) after the signal's beginning to avoid misleading 
%           results caused by possible transient effects caused by digital filtering
%
%        5) because of item 2), the PA(t) outputs are also 2s smaller, but 
%           the percentile values are calculed inside each function before this cut
%
%        6) Loudness and sharpness have the same time vector, but roughness,
%           FS and tonality differ because of their window lengths.
%           Therefore, in order to have the same time vector, after each 
%           respective metric calculation, the outputs are interpolated 
%           with respect to the loudness time vector and all cutted in the 
%           end to the final time corresponding to the FS metric
%
% Author: Gil Felix Greco, Braunschweig 05.04.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help PsychoacousticAnnoyance_More2010;
    return;
end
if nargin < 6
    if nargout == 0
        show = 1;
    else
        show = 0;
    end
end
if nargin < 5
    if nargout == 0
        showPA = 1; 
    else
        showPA = 0;
    end
end

time_insig=(0 : length(insig)-1) ./ fs;  % time vector of the audio input, in seconds

if time_insig(end) < 2
    
    fprintf('\nWARNING: the signal''s length is smaller than 2 seconds.\nDue to the minimum window size used for the fluctuation strength, the computation of a time-varying psychoacoustic annoyance is not possible !!!\nOnly scalar psychoacoustic annoyance number will be calculated for this signal!\n');
    method_FS=0; % stationary method used for the fluctuation strength
else
    
    method_FS=1; % time-varying method used for the fluctuation strength
end

%% modified PA model constants (Ref. [1] pg. 204)

gamma_0=    -0.16;
gamma_1=    11.48;
gamma_2=    0.84;
gamma_3=    1.25;
gamma_4=    0.29;
gamma_5=    5.49;

%% Loudness (according to ISO 531-1:2017)

L = Loudness_ISO532_1( insig, fs,...   % input signal and sampling freq.
                   LoudnessField,...   % field; free field = 0; diffuse field = 1;
                               2,...   % method; stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2; 
                       time_skip,...   % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                               0);     % show results, 'false' (disable, default value) or 'true' (enable)
                           
OUT.L=L; % output loudness results
                                  
%% Sharpness (according to DIN 45692) from loudness input 

S = Sharpness_DIN45692_from_loudness(L.InstantaneousSpecificLoudness,...  % input (time-varying) specific loudness
                                                          'DIN45692',...  % type of weighting function used for sharpness calculation 
                                                              L.time,...  % time vector of the loudness calculation
                                                           time_skip,...  % time_skip (second) for statistics calculation
                                                                  0);     % show sharpness results; true or false
                                                              
OUT.S=S; % output sharpness results
                                           
%% Roughness (according to Daniel & Weber model)

R = Roughness_Daniel1997(insig,fs,...  % input signal and sampling freq.
                        time_skip,...  % time_skip, in seconds for statistical calculations
                                0);     % show results, 'false' (disable, default value) or 'true' (enable)  
                             
OUT.R=R; % output roughness results

%% Fluctuation strength (according to Osses et al. model)

% the output signal will be 2s smaller due to the windown length
FS = FluctuationStrength_Osses2016(insig,fs,...  % input signal and sampling freq.
                                  method_FS,...  % method, stationary analysis =0 - window size=length(insig), time_varying analysis - window size=2s
                                  time_skip,...  % time_skip, in seconds for statistical calculations
                                         0);     % show results, 'false' (disable, default value) or 'true' (enable)  
                                          
OUT.FS=FS; % output fluctuation strength results

%% Tonality (according to Aures' model)

K = Tonality_Aures1985(insig,fs,...  % input signal and sampling freq.
                  LoudnessField,...  % field for loudness calculation; free field = 0; diffuse field = 1;
                              0,...  % time_skip, in seconds for level (stationary signals) and statistics (stationary and time-varying signals) calculations
                             0);...  % show results, 'false' (disable, default value) or 'true' (enable)
                             
OUT.K=K; % output fluctuation strength results

%% for signal with length smaller than 2 s, only scalar psychoacoustic annoyance can be computed

if time_insig(end) < 2
    
    %% (scalar) psychoacoustic annoyance - computed directly from percentile values
    
    % sharpness influence
    if S.S5 > 1.75
        ws = (S.S5-1.75)*(log10(L.N5+10))/4; % in the Fastl&zwicker book, ln is used but it is not clear if it is natural log or log10, but most of subsequent literature uses log10
    else
        ws = 0;
    end
    
    ws( isinf(ws) | isnan(ws) ) = 0;  % replace inf and NaN with zeros
    
    % influence of roughness and fluctuation strength
    wfr = ( 2.18/(L.N5^(0.4)) ).*(0.4*FS.FS5 + 0.6*R.R5);
    
    wfr( isinf(wfr) | isnan(wfr) ) = 0;  % replace inf and NaN with zeros
    
    % Tonality influence
    wt = abs( ( 1-exp(-gamma_4*L.N5) )^2 * ( 1-exp(-gamma_5*K.K5) )^2 );
    
    wt( isinf(wt) | isnan(wt) ) = 0;  % replace inf and NaN with zeros
    
    % More's modified psychoacoustic annoyance
    PA_scalar = abs(L.N5*( 1 + sqrt( gamma_0 + (gamma_1*ws^2) + (gamma_2* wfr^2) + (gamma_3*wt) ) ));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   output struct for time-varying signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % main output results
    OUT.ScalarPA = PA_scalar;               % Annoyance calculated from the percentiles of each variable
    
else % for signals larger than 2 seconds
    
    %% interpolation due to different output lengths of the different metrics
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % step 1) find idx related to the last time step of output from the 
    %         fluctuation strength function (shorter output signal)
    % step 2) cut instaneous quantities - only 1st idx till index related 
    %         to the last time step of fluctuation strength remain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LastTime=FS.time(end);          % take last time of the fluctuation strength
    
    % loudness
    [~,idx_L] = min( abs(L.time-LastTime) ); % step 1) find idx
    
    L.time=L.time(1:idx_L); % step 2) cut signal's end according to idx from step 1)
    L.InstantaneousLoudness=L.InstantaneousLoudness(1:idx_L,1);  % step 2)
    
    % sharpness
    idx_S=idx_L;    % indice is the same as the loudness
    
    S.time=S.time(1:idx_S);  % step 2)
    S.InstantaneousSharpness=S.InstantaneousSharpness(1,1:idx_S);  % step 2)
    
    % roughness
    [~,idx_R] = min( abs(R.time-LastTime) ); % step 1) find idx
    
    R.time=R.time(1:idx_R,1);  % step 2)
    R.InstantaneousRoughness=R.InstantaneousRoughness(1:idx_R,1);  % step 2)
    
    R.time=transpose(R.time);  % step 2)
    R.InstantaneousRoughness=transpose(R.InstantaneousRoughness);  % step 2)
    
    % tonality
    [~,idx_K] = min( abs(K.time-LastTime) ); % step 1) find idx
    
    K.time=K.time(1:idx_K); % step 2) cut signal's end according to idx from step 1)
    K.InstantaneousTonality=K.InstantaneousTonality(1:idx_K,1);  % step 2)
    
    clear idx_R idx_L idx_S idx_K;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    roughness=interp1(R.time,R.InstantaneousRoughness,L.time,'spline'); % interpolation to have the same time vector as loudness metric
    
    fluctuation=interp1(FS.time,FS.InstantaneousFluctuationStrength,L.time,'spline'); % interpolation to have the same time vector as loudness metric
    
    tonality=interp1(K.time,K.InstantaneousTonality,L.time,'spline'); % interpolation to have the same time vector as loudness metric
    
    %% Time-varying psychoacoustic annoyance
    
    % declaring variables for pre allocating memory
    PA=zeros(1,length(L.time));
    ws=zeros(1,length(L.time));
    wfr=zeros(1,length(L.time));
    wt=zeros(1,length(L.time));
    
    for i=1:length(L.time)
        
        % sharpness influence
        if S.InstantaneousSharpness(i) > 1.75
            ws(i) = (S.InstantaneousSharpness(i)-1.75).*(log10(L.InstantaneousLoudness(i)+10))./4; % in the Fastl&zwicker book, ln is used but it is not clear if it is natural log or log10, but most of subsequent literature uses log10
        else
            ws(i) = 0;
        end
        
        ws( isinf(ws) | isnan(ws) ) = 0;  % replace inf and NaN with zeros
        
        % influence of roughness and fluctuation strength
        wfr(i) = ( 2.18./(L.InstantaneousLoudness(i).^(0.4)) ).*(0.4.*fluctuation(i)+0.6.*roughness(i));
        
        wfr( isinf(wfr) | isnan(wfr) ) = 0;  % replace inf and NaN with zeros
        
        % Tonality influence
        wt(i) = abs( ( 1-exp(-gamma_4.*L.InstantaneousLoudness(i)) ).^2 .*( 1-exp(-gamma_5.*tonality(i)) ).^2 );
        
        wt( isinf(wt) | isnan(wt) ) = 0;  % replace inf and NaN with zeros
        
        % More's modified psychoacoustic annoyance
        PA(i) = abs(L.InstantaneousLoudness(i).*( 1 + sqrt( gamma_0 + (gamma_1.*ws(i).^2) + (gamma_2.*wfr(i).^2) + (gamma_3.*wt(i)) ) ));
        
    end
    
    OUT.wt=sqrt(wt); % OUTPUT: tonality and loudness weighting function (not squared)
    OUT.wfr=wfr;     % OUTPUT: fluctuation strength and sharpness weighting function (not squared)
    OUT.ws=ws;       % OUTPUT: sharpness and loudness weighting function (not squared)
    
    %% (scalar) psychoacoustic annoyance - computed directly from percentile values
    
    % sharpness influence
    if S.S5 > 1.75
        ws = (S.S5-1.75)*(log10(L.N5+10))/4; % in the Fastl&zwicker book, ln is used but it is not clear if it is natural log or log10, but most of subsequent literature uses log10
    else
        ws = 0;
    end
    
    ws( isinf(ws) | isnan(ws) ) = 0;  % replace inf and NaN with zeros
    
    % influence of roughness and fluctuation strength
    wfr = ( 2.18/(L.N5^(0.4)) )*(0.4*FS.FS5 + 0.6*R.R5);
    
    wfr( isinf(wfr) | isnan(wfr) ) = 0;  % replace inf and NaN with zeros
    
    % Tonality influence
    wt = abs( ( 1-exp(-gamma_4*L.N5) )^2 * ( 1-exp(-gamma_5*K.K5) )^2 );
    
    wt( isinf(wt) | isnan(wt) ) = 0;  % replace inf and NaN with zeros
    
    % More's modified psychoacoustic annoyance
    PA_scalar = abs(L.N5*( 1 + sqrt( gamma_0 + (gamma_1*ws^2) + (gamma_2* wfr^2) + (gamma_3*wt) ) ));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output struct for time-varying signals
    
    % main output results
    OUT.InstantaneousPA = PA;               % instantaneous Annoyance
    OUT.ScalarPA = PA_scalar;               % Annoyance calculated from the percentiles of each variable
    OUT.time = L.time;                      % time vector
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Statistics from Time-varying PA
        
    [~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector
    
    OUT.PAmean = mean(PA(idx:end));
    OUT.PAstd = std(PA(idx:end));
    OUT.PAmax = max(PA(idx:end));
    OUT.PAmin = min(PA(idx:end));
    OUT.PA1 = get_percentile(PA(idx:end),1);
    OUT.PA2 = get_percentile(PA(idx:end),2);
    OUT.PA3 = get_percentile(PA(idx:end),3);
    OUT.PA4 = get_percentile(PA(idx:end),4);
    OUT.PA5 = get_percentile(PA(idx:end),5);
    OUT.PA10 = get_percentile(PA(idx:end),10);
    OUT.PA20 = get_percentile(PA(idx:end),20);
    OUT.PA30 = get_percentile(PA(idx:end),30);
    OUT.PA40 = get_percentile(PA(idx:end),40);
    OUT.PA50 = median(PA(idx:end));
    OUT.PA60 = get_percentile(PA(idx:end),60);
    OUT.PA70 = get_percentile(PA(idx:end),70);
    OUT.PA80 = get_percentile(PA(idx:end),80);
    OUT.PA90 = get_percentile(PA(idx:end),90);
    OUT.PA95 = get_percentile(PA(idx:end),95);
    
    %% plot
    
    if show == true
        
        il_plotter(OUT.time,L.InstantaneousLoudness,L.N5,'loudness');
        
        il_plotter(OUT.time,S.InstantaneousSharpness,S.S5,'sharpness');
        
        il_plotter(OUT.time,roughness,R.R5,'roughness');
        
        il_plotter(OUT.time,fluctuation,FS.FS5,'fluctuation');
        
        il_plotter(OUT.time,tonality,K.K5,'tonality');
        
    end
    
    if showPA == true
        
        il_plotter(OUT.time,OUT.InstantaneousPA,OUT.PA5,'annoyance')
        
    end
    
end

end % end PA function

%% function plotter, as an inline function:
function il_plotter(time,Instantaneous,percentile,variable)
        
x_axis='Time, $t$ (s)';  % string for the x-axis

switch variable
    case 'loudness'
        p='$N'; % string for percentile legend
        y_axis='Loudness, $N$ (sone)';  % string for the y-axis
        h  =figure('NAME','Loudness');

    case 'sharpness'
        p='$S'; % string for percentile legend
        y_axis='Sharpness, $S$ (acum)'; % string for the y-axis
        h  =figure('NAME','Sharpness');

    case 'roughness'
        p='$R'; % string for percentile legend
        y_axis='Roughness, $R$ (asper)'; % string for the y-axis
        h  =figure('NAME','Roughness');

    case 'fluctuation'
        p='FS$'; % string for percentile legend
        y_axis='Fluctuation strength, $F$ (vacil)'; % string for the y-axis
        h  =figure('NAME','Fluctuation strength');

    case 'tonality'
        p='K$'; % string for percentile legend
        y_axis='Aures tonality, $K$ (t.u.)'; % string for the y-axis
        h  =figure('NAME','Aures tonality');

    case 'annoyance'
        p='PA$'; % string for percentile legend
        y_axis='Psychoacoustic annoyance, PA (-)'; % string for the y-axis
        h  =figure('NAME','Psychoacoustic annoyance');
end

% begin plot

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

plot( time,Instantaneous,'k','Linewidth',0.5,'HandleVisibility','off'); hold on;
plot( time,percentile.*ones(length(time)),'r--','Linewidth',0.5);

legend(sprintf('%s_5=%.2f$',p,percentile),'Location','NorthEast','Interpreter','Latex');
legend boxoff

ylabel(sprintf('%s',y_axis),'Interpreter','Latex');
xlabel(sprintf('%s',x_axis),'Interpreter','Latex');
grid off

set(gcf,'color','w');

end % end plotter function

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
