function [Bark,Bark_raw] = Get_Bark(N,qb,freqs)
% function [Bark,Bark_raw] = Get_Bark(N,qb,freqs)
%
% Extracted from FluctuationStrength_Osses2016.m on 12/05/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Bark_raw = [
    0   0       50      0.5
    1   100     150     1.5
    2   200     250     2.5
    3   300     350     3.5
    4   400     450     4.5
    5   510     570     5.5
    6   630     700     6.5
    7   770     840     7.5
    8   920     1000	8.5
    9   1080	1170	9.5
    10  1270	1370	10.5
    11  1480	1600	11.5
    12  1720	1850	12.5
    13  2000	2150	13.5
    14  2320	2500	14.5
    15  2700	2900	15.5
    16  3150	3400	16.5
    17  3700	4000	17.5
    18  4400	4800	18.5
    19  5300	5800	19.5
    20  6400	7000	20.5
    21  7700	8500	21.5
    22  9500	10500	22.5
    23  12000	13500	23.5
    24  15500   20000   24.5
]; 


Bark_sorted = [  sort([Bark_raw(:,2);Bark_raw(:,3)]),... % frequencies
                 sort([Bark_raw(:,1);Bark_raw(:,4)])];   % Bark

Bark     = zeros(1,round(N/2+1));
Bark(qb) = interp1(Bark_sorted(:,1),Bark_sorted(:,2),freqs);
