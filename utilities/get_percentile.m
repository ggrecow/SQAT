function OUT = get_percentile(input, PercentValue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the value exceeded  during <PercentValue> of the time 
%
% INPUT:
%   input: vector [nSamples x 1] - contain values to be evaluated (typically a time vector)
%   PercentValue: desired % value exceeded during nSamples
%
% OUTPUT:
%   OUT: value exceeded  during <PercentValue> of the time
%
% Standalone example:%
% compute the percentiles from a gaussian function, where x-axis is considered to be increasing time  
%   x = 0:0.1:10;
%   y = gaussmf(x,[2 5]);
%   plot(x,y); hold on;
%   P5=get_percentile(y,5); % value exceed 5% of the time 
%   yline(P5);
%   P90=get_percentile(y,90); % value exceed 90% of the time 
%   yline(P90);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_index = floor( (100-PercentValue)/100 * length(input) );

if X_index==0
    X_index=1;
else
end

sort_input = sort( input );

OUT = sort_input( X_index );

end