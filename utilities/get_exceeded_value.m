function OUT = get_exceeded_value(input, PercentValue)
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
%   P5=get_exceeded_value(y,5); % value exceed 5% of the time
%   yline(P5);
%   P90=get_exceeded_value(y,90); % value exceed 90% of the time
%   yline(P90);
%
% Last checked: Gil Felix Greco, 24.04.2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright statement: This file and code is part of the SQAT toolbox and
% is subject to the MIT license in its entirety, as detailed in the license
% text reproduced at the end of this file (see also <licenses/mit-license.txt> 
% file in the SQAT repository root.
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input dimension (only [Nx1], [Nx2] or [Nx3] are valid)
if  size(input,1) > 3 & size(input,2) > 3 % insig has more than 3 channels
    error('Error: Input signal has more than 3 channels. ')
elseif  size(input, 2) > 3  % insig is [1xN], [2xN] or [3xN]
    input = input';
    % fprintf('\nWarning: Input signal is not [Nx1], [Nx2], or or [Nx3] and was transposed.\n');
end

X_index = floor( (100-PercentValue)/100 * size(input,1) );

if X_index==0
    X_index=1;
else
end

OUT = zeros(1,size(input,2)); % declare output variable

for nColumns = 1:size(input,2) % operate on the columns (time vector must be a column vector)

    sort_input = sort( input(:,nColumns) );

    OUT(1,nColumns) = sort_input( X_index );

end

end

%**************************************************************************
%
% copyright © 2025 Gil Felix Greco
% 
% This software is licenced under the MIT license:
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the “Software”), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, a
% nd/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%
%**************************************************************************
