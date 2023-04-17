function OUT = Sharpness_DIN45692(insig, fs, weight_type, field, method, time_skip, show_sharpness, show_loudness)

%% FUNCTION:
%   OUT = Sharpness_DIN45692(insig, fs, weight_type, field, method, time_skip, show_sharpness, show_loudness)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Stationary and time-varying sharpness calculation according to DIN 45692(2009)
%  from input signal (i.e. the loudness calculation is included within this code)
%
%  Loudness calculation is conducted according to ISO 532:1-2017
%  (type <help loudness_ISO532_1> for more info)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT ARGUMENTS
%   insig : [Nx1] array
%   calibrated audio signal (Pa), 1 channel only
%
%   fs : integer
%   sampling frequency (Hz). For method = 3, provide a dummy scalar
%
%   weight_type : string
%   sharpness calculation using weighting function according to:
%       - 'DIN45692'
%       - 'bismarck'
%       - 'aures' (dependent on the specific loudness level)
%
%   method : integer
%   method used for loudness calculation - method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0 (not accepted in this context); stationary = 1; time varying = 2;
%
%   field : integer
%   type of field used for loudness calculation; free field = 0; diffuse field = 1;
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics calculations (method=1 (time-varying) only)
%
%   show_loudness : logical(boolean)
%   optional parameter to display loudness results (only method=1)
%   'false' (disable, default value) or 'true' (enable).
%
%   show_sharpness : logical(boolean)
%   optional parameter to display sharpness results (only method=1)
%   'false' (disable, default value) or 'true' (enable).
%
% OUTPUTS (method==0; stationary)
%   OUT : struct containing the following fields
%
%       * Sharpness: sharpness (acum)
%
% OUTPUTS (method==1; time-varying)
%   OUT : struct containing the following fields
%
%       * loudness: output struct from loudness calculation (type <help loudness_ISO532_1> for more info)
%       * InstantaneousSharpness: instantaneous sharpness (acum) vs time
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousSharpness (acum)
%         ** Smean : mean value of InstantaneousSharpness (acum)
%         ** Sstd : standard deviation of InstantaneousSharpness (acum)
%         ** Smax : maximum of InstantaneousSharpness (acum)
%         ** Smin : minimum of InstantaneousSharpness (acum)
%         ** Sx : percentile sharpness exceeded during x percent of the signal (acum)
%           *** HINT: time-varying loudness calculation takes some time to
%                     have a steady-response (thus sharpness too!). Therefore, it is a good practice
%                     to consider a time_skip to compute the statistics
%
% Gil Felix Greco, Braunschweig 09.03.2023
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method==1 % stationary loudness calculation
    
    L = Loudness_ISO532_1(insig, fs,... % input signal and sampling freq.
                              field,... % free field = 0; diffuse field = 1;
                             method,... % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                          time_skip,... % time_skip
                      show_loudness);   % show loudness results
    
    n = size(L.SpecificLoudness,2);
    loudness_sones=zeros(size(L.SpecificLoudness,1),1); % pre allocate memory
    SpecificLoudness=L.SpecificLoudness;
    
    for i=1:size(L.SpecificLoudness,1)
        loudness_sones(i)=sum(L.SpecificLoudness(i,:),2).*0.10;
    end
    
elseif method==2 % time-varying loudness calculation
    
    L = Loudness_ISO532_1(insig, fs,... % input signal and sampling freq.
                              field,... % free field = 0; diffuse field = 1;
                             method,... % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                          time_skip,... % time_skip
                      show_loudness);   % show loudness results
    
    n = size(L.InstantaneousSpecificLoudness,2);
    loudness_sones=zeros(size(L.InstantaneousSpecificLoudness,1),1); % pre allocate memory
    SpecificLoudness=L.InstantaneousSpecificLoudness;
    
    for i=1:size(L.InstantaneousSpecificLoudness,1)
        loudness_sones(i)=sum(L.InstantaneousSpecificLoudness(i,:),2).*0.10;
    end
    
end

z=linspace(0.1,24,n); % create bark axis

%% Sharpness calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch weight_type
    case 'DIN45692'   % Widmann model
        
        g=sharpWeights(z,'standard',[]); % calculate sharpness weighting factors
        k=0.11; % adjusted to yield 1 acum using SQAT - DIN45692 allows 0.105<=k<=0.0115 for this weighting function
        
        for i=1:size(SpecificLoudness,1)
            s(i) = k * sum(SpecificLoudness(i,:).*g.*z.*0.10,2) ./ loudness_sones(i);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'aures'           % aures model
        
        for i=1:size(SpecificLoudness,1)
            g(i,:)=sharpWeights(z,'aures',loudness_sones(i)); % calculate sharpness weighting factor
            s(i) = 0.11 * sum(SpecificLoudness(i,:).*g(i,:).*z.*0.10,2) ./ loudness_sones(i);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'bismarck'        % von bismarck
        g=sharpWeights(z,'bismarck',[]); % calculate sharpness weighting factor
        
        for i=1:size(SpecificLoudness,1)
            s(i) = 0.11 * sum(SpecificLoudness(i,:).*g.*z.*0.10,2) ./ loudness_sones(i);
        end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   output struct for time-varying signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method==2 % (time-varying sharpness)
    
    OUT.InstantaneousSharpness = s;                       % instantaneous sharpness
    OUT.time = L.time;                                    % time vector
    OUT.loudness=L;                                       % output struct from the loudness calculation
    
    % statistics from Time-varying sharpness (acum)
    
    [~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector
    
    OUT.Smax = max(s(idx:end));
    OUT.Smin = min(s(idx:end));
    OUT.Smean = mean(s(idx:end));
    OUT.Sstd = std(s(idx:end));
    OUT.S1 = prctile(s(idx:end),99);
    OUT.S2 = prctile(s(idx:end),98);
    OUT.S3 = prctile(s(idx:end),97);
    OUT.S4 = prctile(s(idx:end),96);
    OUT.S5 = prctile(s(idx:end),95);
    OUT.S10 = prctile(s(idx:end),90);
    OUT.S20 = prctile(s(idx:end),80);
    OUT.S30 = prctile(s(idx:end),70);
    OUT.S40 = prctile(s(idx:end),60);
    OUT.S50 = median(s(idx:end));
    OUT.S60 = prctile(s(idx:end),40);
    OUT.S70 = prctile(s(idx:end),30);
    OUT.S80 = prctile(s(idx:end),20);
    OUT.S90 = prctile(s(idx:end),10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% show plots (time-varying)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if show_sharpness == true
        
        figure('NAME','Sharpness analysis (time-varying)');
        
        plot(L.time,OUT.S5*(ones(size(L.time))),'r--'); hold on;
        plot(L.time,s);
        
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('Sharpness, $S$ (acum)','Interpreter','Latex');
        
        legend( sprintf('$S_5$=%g',OUT.S5),'Location','best','Interpreter','Latex');
        legend boxoff
        
        set(gcf,'color','w')
        
    end
    
elseif method==1 % (stationary sharpness)
    
    OUT.Sharpness = s;                       % sharpness
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% embedded function (compute weighting functions according to required model type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function g=sharpWeights(z,type,N)
        
        g=zeros(1,length(z));
        
        switch type
            case 'standard' % Widmann model according to DIN 45692 (2009)
                g(z<15.8)=1;
                g(z>=15.8)=0.15.*exp( 0.42.*((z(z>=15.8))-15.8) ) + 0.85;
                
            case 'bismarck' % von bismark's model according to DIN 45692 (2009)
                g(z<15)=1;
                g(z>=15)=0.2.*exp( 0.308.*(z(z>=15)-15) ) + 0.8;
                
            case 'aures'    % Aure's model according to DIN 45692 (2009)
                for nt=1:length(N)
                    g(nt,:)=0.078.*( exp(0.171.*z)./z ).*( N(nt)./log(0.05.*N(nt)+1));
                end
                
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

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
