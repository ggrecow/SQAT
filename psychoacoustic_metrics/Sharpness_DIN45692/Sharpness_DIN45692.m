function OUT = Sharpness_DIN45692(insig, fs, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness)
% function OUT = Sharpness_DIN45692(insig, fs, weight_type, LoudnessField, LoudnessMethod, time_skip, show_sharpness, show_loudness)
%
%  Stationary and time-varying sharpness calculation according to DIN 45692
%    (2009) from an input signal. The loudness calculation, required as pre-
%    processing for sharpness, is included in this code.
%
%  Loudness calculation is conducted according to ISO 532:1-2017
%  (type <help Loudness_ISO532_1> for more info)
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
%   LoudnessField : integer
%   type of field used for loudness calculation; free field = 0; diffuse field = 1;
%
%   LoudnessMethod : integer
%   method used for loudness calculation - method used for loudness 
%       calculation: stationary (from input 1/3 octave unweighted SPL)=0 (not 
%       accepted in this context); stationary = 1; time varying = 2;
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics 
%        calculations (method=1 (time-varying) only)
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
%       * loudness: output struct from loudness calculation (type 
%         <help loudness_ISO532_1> for more info)
%       * InstantaneousSharpness: instantaneous sharpness (acum) vs time
%       * time : time vector in seconds
%       * Several statistics based on the InstantaneousSharpness (acum)
%         ** Smean : mean value of InstantaneousSharpness (acum)
%         ** Sstd : standard deviation of InstantaneousSharpness (acum)
%         ** Smax : maximum of InstantaneousSharpness (acum)
%         ** Smin : minimum of InstantaneousSharpness (acum)
%         ** Sx : sharpness value exceeded during x percent of the time (acum)
%
%           *** HINT: time-varying loudness calculation takes some time to
%                     have a steady-response (thus sharpness too!).
%                     Therefore, it is a good practice to consider a 
%                     time_skip to compute the statistics
%
% Author: Gil Felix Greco, Braunschweig 09.03.2023
% Author: Gil Felix Greco, Braunschweig 16.02.2025 - introduced get_statistics function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help Sharpness_DIN45692;
    return;
end
if nargin < 8
    if nargout == 0
        show_loudness = 1;
    else
        show_loudness = 0;
    end
end
if nargin < 7
    if nargout == 0
        show_sharpness = 1;
    else
        show_sharpness = 0;
    end
end

if LoudnessMethod==1 % stationary loudness calculation
    
    L = Loudness_ISO532_1(insig, fs,... % input signal and sampling freq.
                      LoudnessField,... % free field = 0; diffuse field = 1;
                     LoudnessMethod,... % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
                          time_skip,... % time_skip
                      show_loudness);   % show loudness results
    
    n = size(L.SpecificLoudness,2);
    loudness_sones=zeros(size(L.SpecificLoudness,1),1); % pre allocate memory
    SpecificLoudness=L.SpecificLoudness;
    
    for i=1:size(L.SpecificLoudness,1)
        loudness_sones(i)=sum(L.SpecificLoudness(i,:),2).*0.10;
    end
    
elseif LoudnessMethod==2 % time-varying loudness calculation
    
    L = Loudness_ISO532_1(insig, fs,... % input signal and sampling freq.
                      LoudnessField,... % free field = 0; diffuse field = 1;
                     LoudnessMethod,... % method used for loudness calculation: stationary (from input 1/3 octave unweighted SPL)=0; stationary = 1; time varying = 2;
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

%% Sharpness calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch weight_type
    case 'DIN45692'   % Widmann model
        
        g=il_sharpWeights(z,'standard',[]); % calculate sharpness weighting factors
        k=0.11; % adjusted to yield 1 acum using SQAT - DIN45692 allows 0.105<=k<=0.0115 for this weighting function
        
        for i=1:size(SpecificLoudness,1)
            s(i) = k * sum(SpecificLoudness(i,:).*g.*z.*0.10,2) ./ loudness_sones(i);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'aures' % Aures model
        
        for i=1:size(SpecificLoudness,1)
            g(i,:)=il_sharpWeights(z,'aures',loudness_sones(i)); % calculate sharpness weighting factor
            s(i) = 0.11 * sum(SpecificLoudness(i,:).*g(i,:).*z.*0.10,2) ./ loudness_sones(i);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'bismarck' % von Bismarck
        g=il_sharpWeights(z,'bismarck',[]); % calculate sharpness weighting factor
        
        for i=1:size(SpecificLoudness,1)
            s(i) = 0.11 * sum(SpecificLoudness(i,:).*g.*z.*0.10,2) ./ loudness_sones(i);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output struct for time-varying signals

if LoudnessMethod==2 % (time-varying sharpness)
    
    OUT.InstantaneousSharpness = s; % instantaneous sharpness
    OUT.time = L.time;              % time vector
    OUT.loudness=L;                 % output struct from the loudness calculation
       
   
    % get statistics from Time-varying sharpness (acum)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~,idx] = min( abs(OUT.time-time_skip) ); % find idx of time_skip on time vector

    metric_statistics = 'Sharpness_DIN45692';
    OUT_statistics = get_statistics( s(idx:end), metric_statistics ); % get statistics

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

      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show plots (time-varying)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
elseif LoudnessMethod==1 % (stationary sharpness)
    
    OUT.Sharpness = s;                       % sharpness
    
end

end % End of Sharpness_DIN45692

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Embedded function (compute weighting functions according to required model type)

function g = il_sharpWeights(z,type,N)

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

end % end of il_sharpWeights

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
