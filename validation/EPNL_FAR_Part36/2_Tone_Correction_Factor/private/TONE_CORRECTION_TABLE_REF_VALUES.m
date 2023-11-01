%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% TONE_CORRECTION_TABLE_REF_VALUES.m
%
% Table 3.7 Example of tone correction calculation for a turbofan engine
%
%  SOURCE: International Civil Aviation Organization (2015) Doc 9501, Environmental Technical Manual
%  Volume I, Procedures for the Noise Certification of Aircraft, Second Edition - ISBN 978-92-9249-721-7
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1/3-oct band(i), freq(Hz), SPL dB, S dB (Step 1), delS dB (Step 2), SPLP dB (Step 4), SP dB (Step 5), SB dB (Step 6), SPLPP dB (Step 7), F dB (Step 8), C dB (Step 9) 
1 50 0 inf inf inf inf inf inf 0 0
2 63 0 inf inf inf inf inf inf 0 0
3 80 70 inf inf 70 -8 -2.333333 70 0 0
4 100 62 -8 inf 62 -8 3.333333 67.6667 0 0
5 125 70 8 16 71 9 6.6667 71 0 0
6 160 80 10 2 80 9 2.6667 77.6667 2.333333 0.29
7 200 82 2 8 82 2 -1.333333 80.333333 1.6667 0.06
8 250 83 1 1 79 -3 -1.333333 79 4 0.61
9 315 76 -7 8 76 -3 0.333333 77.6667 0 0
10 400 80 4 11 78 2 1 78 2 0.17
11 500 80 0 4 80 2 0 79 0 0
12 630 79 -1 1 79 -1 0 79 0 0
13 800 78 -1 0 78 -1 -0.333333 79 0 0
14 1000 80 2 3 80 2 -0.6667 78.6667 0 0
15 1250 78 -2 4 78 -2 -0.333333 78 0 0
16 1600 76 -2 0 76 -2 0.333333 77.6667 0 0
17 2000 79 3 5 79 3 1 78 0 0
18 2500 85 6 3 79 0 -0.333333 79 6 2
19 3150 79 -6 12 79 0 -2.6667 78.6667 0 0
20 4000 78 -1 5 78 -1 -6.333333 76 2 0.33
21 5000 71 -7 6 71 -7 -8 69.6667 0 0
22 6300 60 -11 4 60 -11 -8.6667 61.6667 0 0
23 8000 54 -6 5 54 -6 -8 53 0 0
24 10000 45 -9 3 45 -9 inf 45 0 0