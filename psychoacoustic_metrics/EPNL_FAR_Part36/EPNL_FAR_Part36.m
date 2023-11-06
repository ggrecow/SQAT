function OUT = EPNL_FAR_Part36( insig, fs, method, dt, threshold, show )
% function OUT = EPNL_FAR_Part36( insig, fs, method, dt, threshold, show )
%
%   This function calculates the EFFECTIVE PERCEIVED NOISE LEVEL
%   based on the procedure from:
%
%   [1] Federal Aviation Regulations, 14 CFR Parts 36 and 91,
%       Docket No. FAA-2003-16526; Amendment No. 36-26, 91-288, (2005).
%       Website: https://www.ecfr.gov/current/title-14/appendix-Appendix%20A%20to%20Part%2036
%       (Last viewed 19 Oct 2023)
%
%   Other relevant sources describing the EPNL calculation procedure are:
%
%   [2] Annex 16 to the Convention on International Civil Aviation,
%        Environmental Protection, Volume I - Aircraft Noise, Eighth Edition,
%        July 2017, Internation Civil Aviation Organization
%
%   [3] International Civil Aviation Organization (2015) Doc 9501, Environmental Technical Manual
%        Volume I, Procedures for the Noise Certification of Aircraft, Second Edition - ISBN 978-92-9249-721-7
%
%   PLEASE NOTE - Requirements from the FAA regulations are:
%       1) The third octave frequency analysis is limited from 50 Hz to 10 kHz
%       2) The time step in which the EPNL calculation is conducted from the third octave SPLs is 0.5s
%       3) The threshold for the calculation of the duration correction is
%       a 10 TPNdB decay from PNLTM, i.e. the PNLT curve shall be integrated
%       considering the nearest points where the value (PNLTM-threshold) is observed, being
%       threshold = 10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% INPUT ARGUMENTS
%
%   insig :
%   for method = 0, insig is a SPL[nTime,nFreq] matrix with nFreq=24 columns containing unweighted SPL values
%   for each third-octave band from 50 Hz to 10 kHz, and nTime rows corresponding to the time intervals on which the SPL are provided
%        - in this case, insig is used as is (i.e. without any pre-processing before the EPNL calculation).
%          This means that it will be assumed that each nFreq column of the input matrix corresponds to third-octave bands
%          from 50 Hz to 10 kHz and that each nTime row are given in dt = 0.5 s intervals.
%          BE CAREFUL: These assumptions will be used for the EPNL calculation. If your SPL data has a different dt,
%          it will need to be pre-processed before used here.
%
%   for method = 1, insig is a [nTime,1] array corresponding to a calibrated audio signal (Pa)
%        - in this case, insig will be filtered to get the third octave
%          level from 50 Hz to 10 kHz. The third-octave filters conform with
%          the ones prescribed by ISO 532-1:2017, which comply with IEC
%          61260-1:2014. PLEASE NOTE:  the function <Do_OB13_ISO532_1> used to
%          filter the insig in 1/3 octave bands is hard-coded to work considering a sampling frequency fs=48 kHz.
%          Thus, insig will be resampled to this fs if necessary. After filtering, the prms per freq band is averaged
%          in dt time intervals. dt can be provided as input (see below), but
%          the default value is dt = 0.5 s.
%
%   fs : integer
%   sampling frequency (Hz). For method = 0, provide a dummy scalar
%
%   method : integer
%   for method = 0 - calculates EPNL from an input [nTime,nFreq] matrix with nFreq=24 columns containing unweighted SPL values
%         for each third-octave band from 50 Hz to 10 kHz, and nTime rows corresponding to the time intervals on which the SPL are provided
%
%   for method = 1 - calculates EPNL from an input calibrated audio file, in Pascal unit
%
%   dt : integer
%   time-step, in seconds, in which the third-octave SPLs are averaged to (when method = 1).
%   When method = 0, this parameter only affects the calculation of the duration correction
%   applied to compute the EPNL from the PNLT curve.
%   The default value is dt = 0.5.
%
%   threshold : integer
%   threshold value, in TPNdB, used to calculate the PNLT decay from PNLTM during the calculation of the
%   duration correction.
%   The default value is threshold = 10
%
%   show : logical(boolean)
%   optional parameter for figures (results) display
%   'false' (disable, default value) or 'true' (enable).
%
%%% OUTPUTS
%
%   OUT : struct containing the following fields
%            when method = 0, only quantities with (*) are provided
%            when method = 1, quantities with (*) and (**) are provided
%
%   - SPL quantities considering the original insig
%       ** InstantaneousSPL_insig - overall sound pressure level over time from the original insig (sound file), in dB SPL (from third-octave bands)
%       ** time_insig - time vector of the original input insig (sound file), in seconds
%
%   - SPL quantities averaged in dt time-steps
%       * InstantaneousSPL - overall sound pressure level over time, in dB SPL (from third-octave bands)
%       * time - time vector (sec). this vector is provided considering dt time-steps
%       ** SPL_TOB_spectra - [nTime,nFreq] matrix containing SPL over nTime rows for each nFreq third octave band
%       * TOB_freq - nominal third-octave central frequencies
%
%   - EPNL-related quantities
%       * PN - Perceived Noisiness vs. time, in Noys
%       * PNL - Perceived Noise Level vs. time, in PNdB
%       * PNLM - max. value of the Perceived Noise Level, in PNdB
%       * PNLT - Tone-Corrected Perceived Noise Level vs. time, in TPNdB
%       * PNLTM - max. value of the Tone-Corrected Perceived Noise Level, in TPNdB
%       * EPNL - Effective Perceived Noise Level, in EPNdB
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Source: This code is based on the one provided by:
% Shashikant R. More, Aircraft noise characteristics and metrics, Doctoral thesis, Purdue University, 2010 (permanent link: https://docs.lib.purdue.edu/dissertations/AAI3453255/)
%
% Author: Roberto Merino-Martinez, Delft University (2018) - MATLAB implementation, verification, adaptation
% Author: Gil Felix Greco, Braunschweig 27.10.2023 - adaptation/verification/tested for SQAT.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0 || nargin < 2
    help EPNL_FAR_Part36;
    return;
elseif nargin < 3 % default situation where insig is a sound file, dt and threshold are set to default values, and plots are not shown
    method = 1;
    dt = 0.5;
    threshold = 10;
    show = 0;
elseif ( method==0 && nargin < 4 ) % default situation for method == 0, where dt and threshold are set to default values, and plots are not shown
    dt = 0.5;
    threshold = 10;
    show = 0;
end

if method==0
    if size(insig,2)~=24 % insig matrix needs to have nFreq=24 columns
        warning('For method=0, the insig matrix should have nFreq=24 columns, which corresponds to 1/3 oct. bands from 50 Hz to 10 kHz. Please check the input matrix for the correct dimension!!!');
        return;
    end
end

%%  insig pre-processing stage

fc_TOB = [  50 63 80, 100 125 160, 200 250 315, 400 500 630, ...
    800 1000 1250, 1600 2000 2500, 3150 4000 5000, 6300 8000 10000 ];  % nominal center freq - preferred for freq. labeling (check Tabel E.1. of IEC 61260-1:2014)

num_freqs = length(fc_TOB); % number of freq bands = nFreq

switch method
    
    case 0 % insig is a [nTime,nFreq] matrix containing nFreq=24 columns containing unweighted SPL values for each third octave band from 50 Hz to 10 kHz
        
        SPL_TOB_spectra = insig;
        num_times = size(SPL_TOB_spectra,1); % number of time steps
        
        if num_freqs~=size(SPL_TOB_spectra,2) % insig matrix needs to have nFreq=24 columns
            warning('For method=0, the insig matrix should have nFreq=24 columns, which corresponds to 1/3 oct. bands from 50 Hz to 10 kHz. Please check the input matrix for the correct dimension!!!');
            return;
        end
        
        InstantaneousSPL = 10.*log10(sum(10.^(SPL_TOB_spectra./10),2));  % assumes insig is given in SPL
        
        time = linspace(0,dt*num_times,num_times); % time vector, in dt time-steps
        
        % OUTPUT
        OUT.InstantaneousSPL = InstantaneousSPL;
        OUT.time = time;
        OUT.TOB_freq = fc_TOB;
        
    case 1 % insig is a [nTime,1] array corresponding to a calibrated audio signal (Pa)
        
        if size(insig,2)~=1 % if the insig is not a [Nx1] array
            insig=insig';   % correct the dimension of the insig
        end
        
        % resample to 48 kHz if necessary (this is necessary because the function
        % <Do_OB13_ISO532_1> generates a 1/3 octave fitler bank which is hard-coded
        % to work on fs=48 kHz)
        if fs ~= 48000
            insig = resample(insig,48000,fs);
            fs = 48000;
            fprintf('\n%s.m: The 1/3 octave band filter bank used in this script has only been validated at a sampling frequency fs=48 kHz, resampling to this fs value\n',mfilename);
        end
        
        len_insig = size(insig,1);  % length of the (resample) input vector
        I_REF = 4e-10; % ref. pressure^2
        TINY_VALUE = 1e-12; % small value to avoid inf SPL values
        
        % filter insig to get 1/3-OB
        fmin = 50; % min freq of 1/3-OB is 50 Hz
        fmax = 10000; % max freq of 1/3-OB is 10 kHz
        
        [insig_P_TOB, ~] = Do_OB13_ISO532_1(insig, fs, fmin, fmax); % get 1/3-OB spectra from insig - output is p [nTime,nFreq]
        
        insig_Psquared_TOB = insig_P_TOB.^2;  % squaring the filtered signal - output is p^2 [nTime,nFreq]
        
        InstantaneousSPL_insig = 10*log10( (sum(insig_Psquared_TOB,2) + TINY_VALUE)/I_REF ); % sum energetically all 1/3-OB for each time step to get the overall SPL(t)
        time_insig = ( 1:len_insig)/fs; % time vector of insig
        
        % calculate SPL in dt steps
        Nbins = round(fs*dt); % define dt in N bins
        num_times = ceil(len_insig/Nbins); %  number of time steps of the signal in N blocks
        
        Psquared_TOB = zeros(num_times,num_freqs); % declare variable for memory preallocation
        for i = 1:num_freqs % for each i-th freq band, calculate  p^2 averaged along Nbins blocks related to dt
            
            Psquared_TOB(:,i) = mean( buffer( insig_Psquared_TOB(:,i),Nbins ) ,1 ); % output is p^2[nTime*,nFreq] , where nTime*=round(length(insig)/N)
            
        end
        
        SPL_TOB_spectra = 10*log10( (Psquared_TOB+TINY_VALUE)/I_REF ); % main SPL[nTime*,nFreq] matrix used for the EPNL calculation
        InstantaneousSPL = 10*log10( (sum(Psquared_TOB,2) + TINY_VALUE)/I_REF ); % overall SPL vs. time, in dt time-steps
        
        time = time_insig(1):dt:time_insig(end); % time vector, in dt time-steps
        
        % OUTPUT - quantities from the original insig
        OUT.InstantaneousSPL_insig = InstantaneousSPL_insig;
        OUT.time_insig = time_insig;
        
        % OUTPUT  - quantities averaged in dt time steps
        OUT.InstantaneousSPL = InstantaneousSPL;
        OUT.time = time;
        OUT.SPL_TOB_spectra = SPL_TOB_spectra;
        OUT.TOB_freq = fc_TOB;
        
        clear I_REF fmin fmax insig_Psquared_toct Nbins Psquared_toct
        
end

%% Calculate EPNL

% Convert SPL to Perceived Noisiness (PN) and compute Perceived Noisiness Level (PNL)
[PN, PNL, PNLM, PNLM_idx] = get_PNL(SPL_TOB_spectra);

% Calculate tone-correction and Tone-Corrected Perceived Noise Level (PNLT)
[PNLT, PNLTM, PNLTM_idx, ~] = get_PNLT(SPL_TOB_spectra, fc_TOB, PNL);

% Calculate duration correction factor
[D, idx_t1, idx_t2] = get_Duration_Correction( PNLT, PNLTM, PNLTM_idx, dt, threshold );

% Calculate Effective Perceived Noise Level, unit is EPNdB
OUT.EPNL = PNLTM + D;

% Print calculated EPNL value
fprintf( '\nThe calculated EPNL is %.4g (EPNdB)\n',OUT.EPNL );

% OUTPUTS
OUT.PN = PN; % PERCEIVED NOISINESS, unit is Noys
OUT.PNL = PNL;  % PERCEIVED NOISE LEVEL, unit is PNdB
OUT.PNLM = PNLM; % MAXIMUM PERCEIVED NOISE LEVEL, unit is PNdB
OUT.PNLT = PNLT; % TONE-CORRECTED PERCEIVED NOISE LEVEL, unit is TPNdB
OUT.PNLTM = PNLTM; % MAXIMUM TONE-CORRECTED PERCEIVED NOISE LEVEL (PNLTM)

%%  Show plots

if show == true
    
    xmax = time(end); % used to define the x-axis on the plots
    
    switch method
        
        case 0
            
            figure('name','EPNL calculation based on an input SPL matrix',...
                'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
            
            % plot instantaneous sound pressure level (dBSPL) from original signal and time-averaged over a given dt value
            subplot( 2, 6, [1,2] )
            plot( time, InstantaneousSPL,'Linewidth',2 );
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
            title('Instantaneous overall SPL (1/3 oct. bands)','Interpreter','Latex');
            
            % plot spectrogram (1/3 octave bands in dt time steps)
            subplot( 2, 6, [3,4,5,6] )
            fnom = fc_TOB./1000; % convert center freq to kHz to plot ?
            [xx,yy] = meshgrid( time,fnom );
            pcolor( xx,yy,SPL_TOB_spectra' );
            shading interp; colorbar; axis tight;
            colormap jet;
            
            % freq labels
            ax=gca;
            %     set(ax,'YScale', 'log');
            %     set(ax,'YTick', fnom);
            yticks( [fnom(1) fnom(14) fnom(17:24)] );
            ylabel('Center frequency, $f$ (kHz)','Interpreter','Latex');
            ax.YAxis.MinorTick = 'on';
            ax.YAxis.MinorTickValues = fnom;
            
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel(colorbar,'SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex');
            caxis([0 max(max(SPL_TOB_spectra))]);
            set(ax,'layer','top');
            box on
            title(sprintf('Spectrogram (1/3 oct. bands, dt=%.g sec)',dt),'Interpreter','Latex');
            
            % plot perceived noisiness (noys vs. time)
            subplot( 2, 6, [7,8] )
            plot( time, PN );
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('PN (noys)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
            title('Perceived noisiness','Interpreter','Latex');
            
            % plot perceived noise level (PNdB vs. time)
            subplot( 2, 6, [9,10] )
            plot( time, PNL ); hold on;
            a=plot( time(PNLM_idx), PNLM, 'ro','MarkerSize', 8 );
            legend(a,sprintf('PNLM=%.4g (PNdB)',PNLM));
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('PNL (PNdB)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
            title('Perceived noise level','Interpreter','Latex');
            
            % plot tone-corrected perceived noise level (TPNdB vs. time)
            subplot( 2, 6, [11,12] )
            plot( time, PNLT ); hold on;
            a = plot( time(PNLTM_idx), PNLTM, 'ro','MarkerSize',8);
            b = yline(PNLTM - threshold,'r-');
            c = plot(time( idx_t1 ), PNLT( idx_t1 ),'r*','MarkerSize',10);
            plot(time( idx_t2 ), PNLT( idx_t2 ),'r*','MarkerSize',10);
            
            legend([a,b,c],sprintf('PNLTM=%.4g (TPNdB)',PNLM),...
                sprintf('PNLTM-%.2g=%.4g (TPNdB)',threshold,PNLM-threshold),...
                'PNLT(t1) and PNLT(t2)','Location','SW');
            
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('PNLT (TPNdB)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.05]);
            box on;
            title(sprintf('Tone-corrected perceived noise level - EPNL=%.4g (EPNdB)',OUT.EPNL),'Interpreter','Latex');
            
            set(gcf,'color','w');
            
        case 1
            
            figure('name','EPNL calculation based on an input sound file',...
                'units','normalized','outerposition',[0 0 1 1]); % plot fig in full screen
            
            % plot input signal
            subplot( 2, 6, [1,2] )
            plot( time_insig, insig );
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('Sound pressure, $p$ (Pa)','Interpreter','Latex'); %grid on;
            ax = axis; axis([0 xmax max(insig)*-2 max(insig)*2]);
            title('Input signal','Interpreter','Latex');
            
            % plot instantaneous sound pressure level (dBSPL) from original signal and time-averaged over a given dt value
            subplot( 2, 6, [3,4] )
            plot( time_insig, InstantaneousSPL_insig ); hold on;
            plot( time, InstantaneousSPL,'Linewidth',2 );
            legend(sprintf('dt=%g sec',1/fs),sprintf('dt=%g sec',dt),'Location','SouthWest','Interpreter','Latex'); %legend boxoff;
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
            title('Instantaneous overall SPL (1/3 oct. bands)','Interpreter','Latex');
            
            % plot spectrogram (1/3 octave bands in dt time steps)
            
            subplot( 2, 6, [5,6] )
            fnom = fc_TOB./1000; % convert center freq to kHz to plot ?
            [xx,yy] = meshgrid( time,fnom );
            pcolor( xx,yy,SPL_TOB_spectra' );
            shading interp; colorbar; axis tight;
            colormap jet;
            
            % freq labels
            ax=gca;
            %         set(ax,'YScale', 'log');
            %     set(ax,'YTick', fnom);
            yticks( [fnom(1) fnom(14) fnom(17:24)] );
            ylabel('Center frequency, $f$ (kHz)','Interpreter','Latex');
            ax.YAxis.MinorTick = 'on';
            ax.YAxis.MinorTickValues = fnom;
            
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel(colorbar,'SPL, $L_{\mathrm{p}}$ (dB re 20~$\mu$Pa)','Interpreter','Latex');
            caxis([0 max(max(SPL_TOB_spectra))]);
            set(ax,'layer','top');
            box on
            title(sprintf('Spectrogram (1/3 oct. bands, dt=%.g sec)',dt),'Interpreter','Latex');
            
            % plot perceived noisiness (noys vs. time)
            subplot( 2, 6, [7,8] )
            plot( time, PN );
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('PN (noys)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
            title('Perceived noisiness','Interpreter','Latex');
            
            % plot perceived noise level (PNdB vs. time)
            subplot( 2, 6, [9,10] )
            plot( time, PNL ); hold on;
            a=plot( time(PNLM_idx), PNLM, 'ro','MarkerSize', 8 );
            legend(a,sprintf('PNLM=%.4g (PNdB)',PNLM));
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('PNL (PNdB)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.1]);
            title('Perceived noise level','Interpreter','Latex');
            
            % plot tone-corrected perceived noise level (TPNdB vs. time)
            subplot( 2, 6, [11,12] )
            plot( time, PNLT ); hold on;
            a = plot( time(PNLTM_idx), PNLTM, 'ro','MarkerSize',8);
            b = yline(PNLTM - threshold,'r-');
            c = plot(time( idx_t1 ), PNLT( idx_t1 ),'r*','MarkerSize',10);
            plot(time( idx_t2 ), PNLT( idx_t2 ),'r*','MarkerSize',10);
            
            legend([a,b,c],sprintf('PNLTM=%.4g (TPNdB)',PNLM),...
                sprintf('PNLTM-%.2g=%.4g (TPNdB)',threshold,PNLM-threshold),...
                'PNLT(t1) and PNLT(t2)','Location','SW');
            
            xlabel('Time, $t$ (s)','Interpreter','Latex');
            ylabel('PNLT (TPNdB)','Interpreter','Latex'); grid on;
            ax = axis; axis([0 xmax ax(3) ax(4)*1.05]);
            box on;
            title(sprintf('Tone-corrected perceived noise level - EPNL=%.4g (EPNdB)',OUT.EPNL),'Interpreter','Latex');
            
            set(gcf,'color','w');
            
    end
    
else % if show == 0, dont plot anything
end

end

%**************************************************************************
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%  * Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%  * Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  * Neither the name of the <ORGANISATION> nor the names of its contributors
%    may be used to endorse or promote products derived from this software
%    without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%**************************************************************************
