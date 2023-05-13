function OUT = get_percentile(input,PercentileValue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the value exceeded in a desired percentile % from the entire
% input signal
%
% INPUT:
%   PercentileValue: desired percentile value exceeded during nSamples
%   insig: vector containing values to be evaluated
%
% OUTPUT:
%   OUT: percentile value
%
% Standalone example:%
% compute the percentiles from a gaussian function 
%   x = 0:0.1:10;
%   y = gaussmf(x,[2 5]);
%   plot(x,y); hold on;
%   P5=get_percentile(y,5); % 5% percentile
%   yline(p5);
%   p90=get_percentile(y,90); % 90% percentile
%   yline(p90);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_index = floor( (100-PercentileValue)/100 * length(input) );

if X_index==0
    X_index=1;
else
end

sort_input = sort( input );

OUT = sort_input( X_index );

end