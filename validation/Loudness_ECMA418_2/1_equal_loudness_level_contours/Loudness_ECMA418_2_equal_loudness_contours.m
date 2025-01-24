% Script Loudness_ECMA418_2_equal_loudness_contours
%
% Recreate Fig. A.1 from ECMA-418-2:2024. Compute Loudness (ECMA 418-2:2024) using the 
% Loudness_ECMA418_2 implementation in SQAT and 
% compare with the equal-loudness-level contours from ISO 226:2003
%
% Author: Gil Felix Greco, Braunschweig 24.01.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% compute reference curves from ISO 226:2003

% loudness level (phon)
Ln = [80 60 40 20] ;  

% equal-loudness-level contours from ISO 226:2003 for each Ln
reference(1,:) = equal_loudness_contours(Ln (1:length(Ln) ) );  

%% plot 

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% reference values 
semilogx( reference.frequency, reference.Lp ); hold all;

xlabel( 'Frequency (Hz)','Interpreter','Latex' );
ylabel( 'Sound pressure level (dB SPL)','Interpreter','Latex' );

ylim( [-15 105] );
xlim( [80 11000] );

legend( sprintf('%g phon', Ln(1) ), ...
             sprintf('%g phon', Ln(2) ), ...
             sprintf('%g phon', Ln(3) ), ...
             sprintf('%g phon', Ln(4) ), 'Location','southwest' );

legend boxoff
set(gcf,'color','w');

%% reference results from ISO 226:2003 - values from Table B1

function OUT = equal_loudness_contours(Ln)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given a desired phon level, Ln, compute equal-loudness-levels contours using 
% the formulation from ISO 226:2003, Section 4.1. The equal-loudness-levels contours
% are given in sound pressure levels, Lp (dB SPL), as a function of the
% frequency, in Hz (29 values from 20 Hz - 12.5 kHz)
%
% INPUTS: 
%           Ln : scalar
%                  Desired phon level
%
% OUTPUT:
%           OUT : struct
%                      contains the following fields
%
%           Lp : [Nx1] column vector
%                  sound pressure level, in dB SPL  
%
%            frequency : [Nx1] column vector
%                             frequency, in Hz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency vector
frequency = [ 20; 25; 31.5; 40; 50; 63; 80; 100; 125; 160; ... 
                      200; 250; 315; 400; 500; 630; 800; 1000; 1250; 1600; ...
                      2000; 2500; 3150; 4000; 5000; 6300; 8000; 10000; 12500 ];

% Coefficients from Table 1
af = [ 0.532; 0.506; 0.480; 0.455; 0.432; 0.409; 0.387; 0.367; 0.349; 0.330; ...
         0.315; 0.301; 0.288; 0.276; 0.267; 0.259; 0.253; 0.250; 0.246; 0.244; 
         0.243; 0.243; 0.243; 0.242; 0.242; 0.245; 0.254; 0.271; 0.301 ]; 

Lu = [ -31.6; -27.2; -23.0; -19.1; -15.9; -13.0; -10.3; -8.1; -6.2; -4.5; -3.1; ...
          -2.0; -1.1; -0.4; 0.0; 0.3; 0.5; 0.0; -2.7; -4.1; -1.0; ...
           1.7; 2.5; 1.2; -2.1; -7.1; -11.2; -10.7; -3.1];

Tf = [ 78.5; 68.7; 59.5; 51.1; 44.0; 37.5; 31.5; 26.5; 22.1; 17.9; ...
          14.4; 11.4; 8.6; 6.2; 4.4; 3.0; 2.2; 2.4; 3.5; 1.7; -1.3; ...
          -4.2; -6.0; -5.4; -1.5; 6.0; 12.6; 13.9; 12.3 ]; 
 
% ISO 226:2003 - 4.1. Derive sound pressure level (Lp = dB SPL) from loudness level (Ln = phone)
Af = 4.47e-3 .* ( (10.^(0.025.*Ln) ) - 1.15) + ( (0.4 .* ( 10.^( ( (Tf + Lu) ./ 10 ) - 9 ) ) ).^af ) ;

Lp = ( (10 ./ af) .* log10( Af ) ) - Lu + 94; % Eq. (1)

OUT.frequency = frequency; 
OUT.Lp = Lp;

end
 