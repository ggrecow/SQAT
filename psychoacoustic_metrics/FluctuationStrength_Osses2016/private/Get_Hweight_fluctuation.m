function Hweight = Get_Hweight_fluctuation(fs)
% function Hweight = Get_Hweight_fluctuation(fs)
% 
% Returns the Hweight filter.
% 
% Inputs:
% params: Struct specifying filter characteristics.
% fs: Sampling frequency.
% 
% Outputs:
% Hweight: The digital filter.
% 
% Author: Alejandro Osses/Rodrigo Garcia
% Original file name: Get_Hweight_fluctuation2014
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

try
    load(sprintf('Hweight-%.0f-Hz-LP.mat',fs));
    load(sprintf('Hweight-%.0f-Hz-HP.mat',fs));
    Hweight = [Hweight_HP; Hweight_LP];

catch
    
    warning('Make sure you are using MATLAB 2014a or later...')
    
	% Design parameters of band-pass filter
    sf1 = 0.5; % 0.5
    pf1 = 3.1; % 2
    pf2 = 12; % Hz % 8
    sf2 = 20;
    passAtt1 = 17.5;
    passAtt2 = 14;

    Hweight_lp = designfilt(   'lowpassiir', ...
                            'PassbandFrequency'  , pf2, ...
                            'StopbandFrequency'  , sf2, ...
                            'PassbandRipple'      , 3, ...
                            'StopbandAttenuation', passAtt2, ... % 100
                            'SampleRate'          , fs);

    Hweight_hp = designfilt(   'highpassiir', ...
                            'StopbandFrequency'  , sf1, ...
                            'PassbandFrequency'  , pf1, ...
                            'StopbandAttenuation', passAtt1, ...
                            'PassbandRipple'      , 3, ...
                            'SampleRate'          , fs);
a='a';
    Hweight_HP = Hweight_hp.Coefficients; % second-order sections
    Hweight_LP = Hweight_lp.Coefficients; % second-order sections
    dirout = [Get_TUe_paths('MATLAB') 'Psychoacoustics' delim 'FluctuationStrength_TUe' delim 'private' delim];
    save(sprintf('%sHweight-%.0f-Hz-LP.mat',a,fs),'Hweight_LP');
    save(sprintf('%sHweight-%.0f-Hz-HP.mat',a,fs),'Hweight_HP');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
