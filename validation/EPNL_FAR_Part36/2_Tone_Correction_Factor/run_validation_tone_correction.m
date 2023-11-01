% Script run_validation_tone_correction
%
% Verification of the tone correction factor calculation used
% in the EPNL code 
%   
%   - The Effective Perceived Noise Level (EPNL) is a metric used mainly in
%   the context of environmental aircraft noise, mainly for aircraft noise
%   certification. To the best of knowledge of the Author of this code, a
%   reference signal to validate/verify the implementation of the EPNL
%   calculation does not exist. However, it is possible to verify
%   particular parts of the EPNL calculation.
%
%   In this code, we use a table providing reference SPL spectra
%        values and the respective values that shall be computed in each step of the 
%        tone-correction factor calculation to verify the <get_PNLT> function implemented in the
%        <EPNL_FAR_Part36> code provided in SQAT. The table with the reference values is provided by:  
%
%   International Civil Aviation Organization (2015) Doc 9501, Environmental Technical Manual
%   Volume I, Procedures for the Noise Certification of Aircraft, Second
%   Edition - ISBN 978-92-9249-721-7 (see Table 3.7)
%
%  HOW TO RUN THIS CODE: this is a standalone code. Therefore, no additional steps are
%  necessary to run this code.
%
% Author: Gil Felix Greco, Braunschweig, 30.10.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

save_figs = 0; %% save figs flag

%% path settings 

dir_out = [fileparts(mfilename('fullpath')) filesep];
 
% Figure where the figures (and the results) will be stored:
figures_dir = [dir_out 'figs' filesep];
if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

%% Read Table A1-3. Example of tone correction calculation for a turbofan engine, which provides reference values for the tone correction calculation

tone_correction_tab = load('TONE_CORRECTION_TABLE_REF_VALUES.m');
% Band = tone_correction_tab(:,1);
freq_bands = tone_correction_tab(:,2);
SPL_ref = tone_correction_tab(:,3);
S_ref = tone_correction_tab(:,4);
delS_ref = tone_correction_tab(:,5);
SPLP_ref = tone_correction_tab(:,6);
SP_ref = tone_correction_tab(:,7);
SB_ref = tone_correction_tab(:,8);
SPLPP_ref = tone_correction_tab(:,9);
F_ref = tone_correction_tab(:,10);
C_ref = tone_correction_tab(:,11);
 
%% compute tone-correction using SQAT

PNL = 0; % dummy value for the PNL

[~, ~, ~, OUT] = get_PNLT( SPL_ref', freq_bands', PNL );

%% plot results 

figure('name','EPNL calculation - verification of tone-correction implementation part 1',...
    'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
   
Line_ref = 1;
Line_sqat= 1;

% plot input SPL levels per frequency band
subplot( 2, 6, [1,2] )
semilogx( freq_bands, SPL_ref );
xlabel('Frequency, $f$ (Hz)', 'Interpreter','Latex');
ylabel('Sound pressure level, $L_{\mathrm{p}}$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end-1)]);
title('Input SPL spectrum','Interpreter','Latex');

% 'Step 1
subplot( 2, 6, [3,4] )
semilogx( freq_bands, S_ref, 'b', 'Linewidth', Line_ref ); hold on;
semilogx( freq_bands, OUT.S, 'r--', 'Linewidth', Line_sqat ); hold on;
legend('Reference', 'SQAT' ,'Location', 'NE', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('$S$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(S_ref) max(S_ref)]);
title('Step 1','Interpreter','Latex');

% 'Step 2
subplot( 2, 6, [5,6] )
semilogx( freq_bands, delS_ref, 'b', 'Linewidth',Line_ref ); hold on;
semilogx( freq_bands, abs(OUT.diff), 'r--', 'Linewidth', Line_sqat ); hold on;
legend('Reference', 'SQAT' ,'Location', 'NE', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('$|\Delta S|$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(delS_ref) max(delS_ref)]);
title('Step 2','Interpreter','Latex');

% 'Step 4
subplot( 2, 6, [7,8] )
semilogx( freq_bands, SPLP_ref, 'b', 'Linewidth',Line_ref ); hold on;
semilogx( freq_bands, OUT.SPLP ,'r--', 'Linewidth', Line_sqat ); hold on;
legend('Reference', 'SQAT' ,'Location', 'SE', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('SPL$^{\prime}$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(SPLP_ref) max(SPLP_ref)]);
title('Step 4','Interpreter','Latex');

% 'Step 5
subplot( 2, 6, [9,10] )
semilogx( freq_bands, SP_ref, 'b', 'Linewidth', Line_ref ); hold on;
semilogx( freq_bands, OUT.SP(1:end-1), 'r--', 'Linewidth', Line_sqat ); hold on;
legend('Reference', 'SQAT' ,'Location', 'NE', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('S$^{\prime}$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(SP_ref) max(SP_ref)]);
title('Step 5','Interpreter','Latex');

% 'Step 6
subplot( 2, 6, [11,12] )
semilogx( freq_bands, SB_ref, 'b', 'Linewidth', Line_ref ); hold on;
semilogx( freq_bands, OUT.SB(1:end-1), 'r--', 'Linewidth', Line_sqat ); hold on;
legend('Reference', 'SQAT' ,'Location', 'NE', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('$\overline{\mathrm{S}}$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end-1)]);
ylim([min(SP_ref) max(SP_ref)]);
title('Step 6','Interpreter','Latex');

set(gcf,'color','w');
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'validation_Tone_Correction_1';
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot results 

figure('name','EPNL calculation - verification of tone-correction implementation part 2',...
    'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
   
% Step 7
subplot( 2, 6, [1,2] )
semilogx( freq_bands, SPLPP_ref,'b' ); hold on;
semilogx( freq_bands, OUT.SPLPP,'r--' );
legend('Reference', 'SQAT' ,'Location', 'SE', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('SPL$^{\prime\prime}$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(SPLPP_ref) max(SPLPP_ref)]);
title('Step 7','Interpreter','Latex');

% 'Step 8
subplot( 2, 6, [3,4] )
semilogx( freq_bands, F_ref,'b' ); hold on;
semilogx( freq_bands, OUT.F,'r--' );
legend('Reference', 'SQAT' ,'Location', 'NW', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('$F$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(F_ref) max(F_ref)]);
title('Step 8','Interpreter','Latex');

% 'Step 9
subplot( 2, 6, [5,6] )
semilogx( freq_bands, C_ref,'b' ); hold on;
semilogx( freq_bands, OUT.C,'r--' );
legend('Reference', 'SQAT' ,'Location', 'NW', 'Interpreter', 'Latex'); %legend boxoff;
xlabel('Frequency, $f$ (Hz)','Interpreter','Latex');
ylabel('$C$ (dB SPL)','Interpreter','Latex'); grid on;
xlim([freq_bands(1) freq_bands(end)]);
ylim([min(C_ref) max(C_ref)]);
title('Step 9','Interpreter','Latex');

set(gcf,'color','w');
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figs==1
    
    % Figure where the figures (and the results) will be stored:
    figures_dir = [dir_out 'figs' filesep];
    if ~exist(figures_dir,'dir')
        mkdir(figures_dir);
    end
    
    figname_short = 'validation_Tone_Correction_2';
    figname_out = [figures_dir figname_short];
    
    %     saveas(gcf,figname_out, 'fig');
    %     saveas(gcf,figname_out, 'pdf');
    saveas(gcf,figname_out, 'png');
    
    fprintf('%s.m: figure %s was saved on disk\n\t(full name: %s)\n',mfilename,figname_short,figname_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

