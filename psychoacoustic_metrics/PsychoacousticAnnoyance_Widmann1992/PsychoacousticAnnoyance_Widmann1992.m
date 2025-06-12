function OUT = PsychoacousticAnnoyance_Widmann1992(insig,fs,LoudnessField,time_skip,showPA,show)
% function OUT = PsychoacousticAnnoyance_Widmann1992(insig,fs,LoudnessField,time_skip,showPA,show)
%
%   This function calculates Widmann's psychoacoustic annoyance model from an input acoustic signal
%
%   The psychoacoustic annoyance model is according to: (page 66) Widmann, U. (1992). Ein Modell der
%   Psychoakustischen Lästigkeit von Schallen und seine Anwendung in der Praxis der Lärmbeurteilung
%   (A model of the psychoacoustic annoyance of sounds and its application in noise assessment practice)
%   [Doctoral thesis, Technische Universität München (Technical University of Munich)].
%
%   As clarified by Lotinga, M. J. B. and A. J. Torija (2025) in
%   "Comment on "A study on calibration methods of noise annoyance data from listening tests"
%   [J. Acoust. Soc. Am. 156, 1877–1886 (2024)]." Journal of the Acoustical
%   Society of America 157(5): 3282–3285, this model is the same as that commonly
%   misattributed to (page 327) Zwicker, E. and Fastl, H. Second ed,
%   Psychoacoustics, Facts and Models, 2nd ed. Springer-Verlag, Berlin, 1999.
%
% - This metric combines 4 psychoacoustic metrics to quantitatively describe annoyance:
%
%    1) Loudness (sone) - calculated hereafter following ISO 532-1:2017
%       type <help Loudness_ISO532_1> for more info
%
%    2) Sharpness (acum) - calculated hereafter following DIN 45692:2009
%       NOTE: uses DIN 45692 weighting function by default, please change code if
%       the use of a different weighting function is desired. However, note
%       that the original PA model sharpness weighting is equal to the DIN
%       45692 weighting (i.e., Widmann weighting).
%       type <help Sharpness_DIN45692_from_loudness>
%
%    3) Roughness (asper) - calculated hereafter following Daniel & Weber model
%       type <help Roughness_Daniel1997> for more info
%
%    4) Fluctuation strength (vacil) - calculated hereafter following Osses et al. model
%       type <help FluctuationStrength_Osses2016> for more info
%
%   It should be noted however that the original model used metrics for
%   roughness and fluctuation strength that differ from those employed in
%   this implementation. The metrics employed in the original PA model
%   comprised Fastl's roughness and fluctuation strength (see Fastl &
%   Zwicker, 2007. Psychoacoustics: Facts and models.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   chose field for loudness calculation; free field = 0; diffuse field = 1;
%   type <help Loudness_ISO532_1> for more info
%
%   show : logical(boolean)
%   optional parameter, display results of loudness, sharpness, roughness and fluctuation strength
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
%             ** wfr : fluctuation strength and roughness weighting function (not squared)
%             ** ws : sharpness and loudness weighting function (not squared)
%
%             ** Statistics
%               *** PAmean : mean value of psychoacoustic annoyance (unit)
%               *** PAstd : standard deviation of instantaneous psychoacoustic annoyance (unit)
%               *** PAmax : maximum of instantaneous psychoacoustic annoyance (unit)
%               *** PAmin : minimum of instantaneous psychoacoustic annoyance (unit)
%               *** PAx : value exceeded x percent of the time
%
%      * include structs with the results from the other metrics computed
%        **  L : struct with Loudness results, type <help Loudness_ISO532_1> for more info
%        **  S : struct with Sharpness, type <help Sharpness_DIN45692_from_loudness>
%        **  R : strcut with roughness results, type <help Roughness_Daniel1997> for more info
%        ** FS : struct with fluctuation strength results, type <help FluctuationStrength_Osses2016> for more info
%
%
%  NOTE: 1) Input signals should be in pascal values or calibrated .wav files
%
%        2) Fluctuation strength window has length of 2s. If the signal is less than 2s long, the FS calculation will be automatically
%           changed to stationary (i.e. uses a window with length equal to signal's size) . in this case, no time-varying PA is available.
%
%        3) Be aware that, because of item 2), if the signal is more than 2s long, the last 2 seconds of the input signal are LOST !!!!
%
%        4) is a best practice to compute percentile values following a time_skip (s) after the signal's beginning to avoid misleading results caused by possible transient effects caused by digital filtering
%
%        5) because of item 2), the PA(t) outputs are also 2s smaller, but the percentile values are calculed inside each function before this cut
%
%        6) Loudness and sharpness have the same time vector, but roughness and FS differ because of their window lengths of 200ms and 2s, respectively.
%           Therefore, in order to have the same time vector, after each respective metric calculation, the outputs are interpolated with respect to the loudness time vector and all cutted in the end
%           to the final time corresponding to the FS metric
%
% Author: Gil Felix Greco, Braunschweig 04.03.2020 (updated 14.03.2023)
% Author: Gil Felix Greco, Braunschweig 16.02.2025 - introduced get_statistics function
% Modified: Mike Lotinga, 12.06.2025 - created from
% PsychoacousticAnnoyance_Zwicker1999.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help PsychoacousticAnnoyance_Widmann1992;
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

%% for signal with length smaller than 2 s, only scalar psychoacoustic annoyance can be computed

if time_insig(end) < 2
    
    %% (scalar) psychoacoustic annoyance - computed directly from percentile values
    
    % sharpness and loudness influence
    if S.S5 > 1.75
        % in Widmann (1992), log is used without specifying the base. In
        % Fastl&Zwicker (2007), lg is used and subsequent literature also uses log10
        ws = (S.S5-1.75).*(log10(L.N5+10))./4; 
    else
        ws = 0;
    end
    
    ws(isinf(ws)|isnan(ws)) = 0;  % replace inf and NaN with zeros
    
    % influence of roughness and fluctuation strength
    wfr = ( 2.18./(L.N5.^(0.4)) ).*(0.4.*FS.FS5 + 0.6.*R.R5);
    
    wfr(isinf(wfr)|isnan(wfr)) = 0;  % replace inf and NaN with zeros
    
    % psychoacoustic annoyance
    PA_scalar = L.N5.*(1 + sqrt (ws.^2 + wfr.^2));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   output struct for time-varying signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % main output results
    OUT.ScalarPA = PA_scalar;               % Annoyance calculated from the percentiles of each variable
    
else % for signals larger than 2 seconds
    
    %% interpolation due to different output lengths of the different metrics
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % step 1) find idx related to the last time step of output from the fluctuation strength function (shorter output signal)
    % step 2) cut instaneous quantities from 1st idx to index related to the last time step of fluctuation strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    clear idx_R idx_L idx_S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    roughness=interp1(R.time,R.InstantaneousRoughness,L.time,'spline'); % interpolation to have the same time vector as loudness metric
    
    fluctuation=interp1(FS.time,FS.InstantaneousFluctuationStrength,L.time,'spline'); % interpolation to have the same time vector as loudness metric
    
    %% Time-varying psychoacoustic annoyance
    
    % declaring variables for pre allocating memory
    PA=zeros(1,length(L.time));
    ws=zeros(1,length(L.time));
    wfr=zeros(1,length(L.time));
    
    for i=1:length(L.time)
        
        % sharpness influence
        if S.InstantaneousSharpness(i) > 1.75
            % in Widmann (1992), log is used without specifying the base. In
            % Fastl&Zwicker (2007), lg is used and subsequent literature also uses log10
            ws(i) = (S.InstantaneousSharpness(i)-1.75).*(log10(L.InstantaneousLoudness(i)+10))./4;
        else
            ws(i) = 0;
        end
        
        ws(isinf(ws)|isnan(ws)) = 0;  % replace inf and NaN with zeros
        
        % influence of roughness and fluctuation strength
        wfr(i) = ( 2.18./(L.InstantaneousLoudness(i).^(0.4)) ).*(0.4.*fluctuation(i)+0.6.*roughness(i));
        
        wfr(isinf(wfr)|isnan(wfr)) = 0;  % replace inf and NaN with zeros
        
        % psychoacoustic annoyance
        PA(i) = L.InstantaneousLoudness(i).*(1 + sqrt (ws(i).^2 + wfr(i).^2));
        
    end
    
    OUT.wfr=wfr;     % OUTPUT: fluctuation strength and sharpness weighting function (not squared)
    OUT.ws=ws;       % OUTPUT: sharpness and loudness weighting function (not squared)
    
    %% (scalar) psychoacoustic annoyance - computed directly from percentile values
    
    % sharpness influence
    if S.S5 > 1.75
        % in Widmann (1992), log is used without specifying the base. In
        % Fastl&Zwicker (2007), lg is used and subsequent literature also uses log10
        ws = (S.S5-1.75).*(log10(L.N5+10))./4;
    else
        ws = 0;
    end
    
    ws(isinf(ws)|isnan(ws)) = 0;  % replace inf and NaN with zeros
    
    % influence of roughness and fluctuation strength
    wfr = ( 2.18./(L.N5.^(0.4)) ).*(0.4.*FS.FS5 + 0.6.*R.R5);
    
    wfr(isinf(wfr)|isnan(wfr)) = 0;  % replace inf and NaN with zeros
    
    % psychoacoustic annoyance
    PA_scalar = L.N5.*(1 + sqrt (ws.^2 + wfr.^2));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   output struct for time-varying signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % main output results
    OUT.InstantaneousPA = PA;               % instantaneous Annoyance
    OUT.ScalarPA = PA_scalar;               % Annoyance calculated from the percentiles of each variable
    OUT.time = L.time;                      % time vector
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get statistics from Time-varying PA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector

    metric_statistics = 'PsychoacousticAnnoyance_Widmann1992';
    OUT_statistics = get_statistics( PA(idx:end), metric_statistics ); % get statistics

    % copy fields of <OUT_statistics> struct into the <OUT> struct
    fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

    for i = 1:numel(fields_OUT_statistics)
        fieldName = fields_OUT_statistics{i};
        if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
            OUT.(fieldName) = OUT_statistics.(fieldName);
        end
    end

    clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% plot
    
    if show == true
        
        il_plotter(OUT.time,L.InstantaneousLoudness,L.N5,'loudness');
        
        il_plotter(OUT.time,S.InstantaneousSharpness,S.S5,'sharpness');
        
        il_plotter(OUT.time,roughness,R.R5,'roughness');
        
        il_plotter(OUT.time,fluctuation,FS.FS5,'fluctuation');
        
    end
    
    if showPA == true
        
        il_plotter(OUT.time,OUT.InstantaneousPA,OUT.PA5,'annoyance')
        
    end
    
end

end % end PA function

%% function plotter
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
