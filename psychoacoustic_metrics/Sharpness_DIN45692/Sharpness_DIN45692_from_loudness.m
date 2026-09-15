function OUT = Sharpness_DIN45692_from_loudness(SpecificLoudness, weight_type, time, time_skip, show_sharpness)
% function OUT = Sharpness_DIN45692_from_loudness(SpecificLoudness, weight_type, time, time_skip, show_sharpness)
%
%  Stationary and time-varying sharpness calculation according to DIN 45692(2009)
%  from input specific loudness (i.e. the loudness calculation is not included within this code)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT ARGUMENTS
%   SpecificLoudness : array
%   if method = 0 (stationary) - Specific loudness [1,sone/Bark]
%   if method = 1 (time-varying) - Instantaneous specific loudness [nTimeSteps,sone/Bark]
%
%   weight_type : string
%       weighting function used for sharpness calculation, according to:
%       - 'DIN45692'
%       - 'bismarck'
%       - 'aures' (dependent on the specific loudness level)
%
%   time : array
%       time vector of the specific loudness [1,nTimeSteps] - used only for
%       plot purposes if method = 1 (time-varying)
%
%   time_skip : integer
%   skip start of the signal in <time_skip> seconds for statistics 
%       calculations (method=1 (time-varying) only)
%
%   show : logical(boolean)
%   optional parameter for figures (results) display (only method=1)
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
% Modified: Mike Lotinga 29.05/2026 - vectorised for 95% reduction in
% compute time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file is part of the SQAT toolbox and is subject
% to the GPL-3.0 license, as detailed in <licenses/gpl-3.0.txt> in the SQAT
% repository root. Some files in SQAT carry a different license, always
% stated in their own header; where this file depends on them, the combined
% work remains governed by the GPL-3.0.
%
% As per the licensing information, this file is provided "as is", WITHOUT
% WARRANTY OF ANY KIND, express or implied, including but not limited to the
% warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    if nargout == 0
        show_sharpness = 1;
    else
        show_sharpness = 0;
    end
end

n = size(SpecificLoudness, 2);
z = linspace(0.1, 24, n);   % create bark axis

if size(SpecificLoudness, 1) == 1 % define method based on the size of the input specific loudness
    method = 0; % (stationary) - Specific loudness [1,sone/Bark]
else
    method = 1; % (time-varying) - Instantaneous specific loudness [nTimeSteps,sone/Bark]
end

loudness_sones = sum(SpecificLoudness, 2).*0.10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sharpness calculation

switch weight_type
    case 'DIN45692' % Widmann model
        
        g = il_sharpWeights(z, 'standard', []); % calculate sharpness weighting factors
        k = 0.11; % adjusted to yield 1 acum using SQAT - DIN45692 allows 0.105<=k<=0.0115 for this weighting function
        s = k * sum(SpecificLoudness.*g.*z.*0.10, 2) ./ loudness_sones;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'aures' % Aures model
        
        g = il_sharpWeights(z, 'aures', loudness_sones); % calculate sharpness weighting factor
        s = 0.11 * sum(SpecificLoudness.*g.*z.*0.10, 2) ./ loudness_sones;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'bismarck' % von Bismarck
        g = il_sharpWeights(z,'bismarck',[]); % calculate sharpness weighting factor
        s = 0.11 * sum(SpecificLoudness.*g.*z.*0.10, 2) ./ loudness_sones;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   output struct for time-varying signals

if method == 1 % (time-varying sharpness)
    
    OUT.InstantaneousSharpness = s; % instantaneous sharpness
    OUT.time = time;                % time vector
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show plots (time-varying)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if show_sharpness == true
        
        figure('name','Sharpness analysis (time-varying)')
        
        plot(time,OUT.S5*(ones(size(time))),'r--'); hold on;
        plot(time,s);
        
        xlabel('Time, $t$ (s)','Interpreter','Latex');
        ylabel('Sharpness, $S$ (acum)','Interpreter','Latex');
        
        legend( sprintf('$S_5$=%g',OUT.S5),'Location','best','Interpreter','Latex');
        legend boxoff
        
        set(gcf,'color','w')
        
    end
    
elseif method == 0 % (stationary sharpness)
    
    OUT.Sharpness = s;                       % sharpness
    
end
end % end of function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Embedded function (compute weighting functions according to required model type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = il_sharpWeights(z, type, N)
    switch type
        case 'standard' % Widmann model according to DIN 45692 (2009)
            g = zeros(1,length(z));
            g(z < 15.8) = 1;
            g(z >= 15.8) = 0.15.*exp(0.42.*((z(z >= 15.8)) - 15.8)) + 0.85;

        case 'bismarck' % von bismark's model according to DIN 45692 (2009)
            g = zeros(1,length(z));
            g(z < 15) = 1;
            g(z >= 15) = 0.2.*exp(0.308.*(z(z >= 15) - 15)) + 0.8;

        case 'aures'    % Aures' model according to DIN 45692 (2009)
            g = zeros(length(N), length(z));
            g = 0.078.*(exp(0.171.*z)./z ).*(N./log(0.05.*N + 1));
    end
end % end of il_sharpWeights