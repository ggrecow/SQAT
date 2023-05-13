% Script run_verification_modulated_white_noise
%
% Compute roughness of:
%       1) unmodulated white noise (rms 70 dBSPL)
%
% Roughness computed using:
%   OUT = Roughness_Daniel1997(insig,fs,time_skip,show) 
%   type <help Roughness_Daniel1997> for more info
%
% Author: Gil Felix Greco, Braunschweig 13.05.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; % no=0; yes=1

dir_out = [fileparts(mfilename('fullpath')) filesep];

%% Specify signal temporal characteristics

L=10;                            % Length (seconds)
fs=48000;                        % Sampling frequency
dt=1/fs;                         % Time step
t = 0:dt:L;                      % Time vector

%% generate unmodulated white noise with rms=70dBSPL

levelOut = 70;                   % desired signal level (dBSPL)

s=randn(1,length(t));            % white noise

levelIn = 20*log10(rms(s)/2e-5); % compute level from generated signal

WhiteNoise = s * 10^((levelOut-levelIn)/20); % correct level according to levelOut

SPL_WhiteNoise=20*log10(rms(WhiteNoise)/2e-5);  % verify SPL of the signal

%% compute roughness

R_WhiteNoise=Roughness_Daniel1997(WhiteNoise,fs,...  % input signal and sampling freq.
                                                   0,...  % time_skip, in seconds for statistical calculations
                                                   0);    % show results, 'false' (disable, default value) or 'true' (enable)  

%% Plot white noise and roughness result

h  =figure;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% pressure 
yyaxis left
plot(t,WhiteNoise);
% axis([0 .11 -0.1 0.1]);
ylabel('Acoustic pressure, $p$ (Pa)','Interpreter','Latex');
xlabel('Time, $t$ (s)','Interpreter','Latex'); 

% roughness
yyaxis right
plot(R_WhiteNoise.time,R_WhiteNoise.InstantaneousRoughness);hold on;
a=plot(R_WhiteNoise.time,R_WhiteNoise.Rmean*ones(length(R_WhiteNoise.time)),'k--');
legend(a,sprintf('$R_{\\mathrm{mean}}=%.3g$ (asper)',R_WhiteNoise.Rmean),'Location','NorthEast','Interpreter','Latex'); %legend boxoff 

axis([0 L 0 0.3]);
ylabel('Roughness, $R$ (asper)','Interpreter','Latex');

set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'roughness_unmodulated_white_noise_70dBSPL';
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
