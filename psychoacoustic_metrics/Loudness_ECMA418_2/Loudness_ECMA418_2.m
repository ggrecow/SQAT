function OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
% OUT = Loudness_ECMA418_2(insig, fs, fieldtype, time_skip, show)
%
% Returns loudness values according to ECMA-418-2:2025 (using the Sottek
% Hearing Model) for an input calibrated single mono or single stereo
% audio (sound pressure) time-series signal, insig. For stereo signals,  
% Loudness is calculated for each channel [left ear, right ear],
% and also for the combination of both, denominated as 
% "combined binaural" (see Section 8.1.15 ECMA-418-2:2025). 
% 
% According to ECMA-418-2:2025 (Section 8.1.4), the representative 
% single value to express the overall loudness is obtained by a power
% average of the time-dependent loudness, provided here as the
% <loudnessPowAvg> output variable.
%
% NOTE: according to ECMA-418-2:2025 (Section 8.1.4), the loudness
% values corresponding to the first 300 ms of the input signal must be
% discarded due to the transient responses of the digital filters.
% Considering the temporal resolution of the loudness model, this means any
% values below 304 ms must be discarded.
% Therefore, <time_skip> must be greater or equal to 304 ms to compute any
% time-aggregated quantity. 
%
% Reference signal: pure tone with center frequency of 1 kHz and RMS value
% of 40 dBSPL equals 1 sone_HMS
%
% Inputs
% ------
% insig : column vector [Nx1] mono or [Nx2] binaural
%            the input signal as single mono or stereo audio (sound
%            pressure) signals
%
% fs : integer
%       the sample rate (frequency) of the input signal(s)
%
% fieldtype : keyword string (default: 'free-frontal')
%                 determines whether the 'free-frontal' or 'diffuse' field stages
%                 are applied in the outer-middle ear filter
%
% time_skip : integer (default: 304 ms seconds - see ECMA-418-2:2025, Section 8.1.4)
%                   skip start of the signal in <time_skip> seconds so that
%                   the transient response of the digital filters is avoided.
%                   Best-practice: <time_skip> must be equal or higher than
%                   default value
%
% show : Boolean true/false (default: false)
%             flag indicating whether to generate a figure from the output
%
% Returns
% -------
%
% OUT : structure
%            contains the following fields:
%
% specLoudness : matrix
%                time-dependent specific loudness for each (half) critical
%                band
%                arranged as [time, bands(, channels)]
%
% specloudnessPowAvg : matrix
%                      time-power-averaged specific loudness for each
%                      (half) critical band
%                      arranged as [bands(, channels)]
%                      OBS: takes <time_skip> into consideration
%
% loudnessTDep : vector or matrix
%                 time-dependent overall loudness
%                 arranged as [time(, channels)]
% 
% loudnessPowAvg : number or vector
%                  time-power-averaged overall loudness
%                  arranged as [loudness(, channels)]
%                  OBS: takes <time_skip> into consideration
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each (half)
%                   critical band rate scale width
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% timeInsig : vector
%           time (seconds) of insig
%
% soundField : string
%              identifies the soundfield type applied (the input argument
%              fieldtype)
%
% Several statistics based on loudnessTDep
%         ** Nmean : mean value of instantaneous loudness (sone_HMS)
%         ** Nstd : standard deviation of instantaneous loudness (sone_HMS)
%         ** Nmax : maximum of instantaneous loudness (sone_HMS)
%         ** Nmin : minimum of instantaneous loudness (sone_HMS)
%         ** Nx : loudness value exceeded during x percent of the time (sone_HMS)
%             in case of binaural input, Nx(1,3), being 1st, 2nd and 3rd column 
%             corresponding to [left ear, right ear, comb. binaural] 
%             OBS: all quantities here take <time_skip> into consideration
%
% In case of stereo inputs, the following additional fields are provided
% separately for the "comb. binaural" case (i.e. combination of left and right ears)  
%
% specLoudnessBin : matrix
%                 time-dependent specific loudness for each (half)
%                 critical band
%                 arranged as [time, bands]
%
% specLoudnessPowAvgBin : matrix
%                      time-power-averaged specific loudness for each
%                      (half) critical band
%                      arranged as [bands]
%                      OBS: takes <time_skip> into consideration
%
% loudnessTDepBin : vector or matrix
%                 time-dependent overall loudness
%                 arranged as [time]
%
% loudnessPowAvgBin : number or vector
%                  time-power-averaged overall loudness
%                  arranged as [time]
%                  OBS: takes <time_skip> into consideration
%
% If show==true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall loudness, with the latter also indicating the time-aggregated
% value. In case of stereo signals, a set of plots is returned for each input channel, 
% with another set for the combined binaural loudness. For the latter, the
% indicated time-averaged A-weighted sound level corresponds with the channel with 
% the highest sound level.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 22/09/2023
% Date last modified: 27/06/2025
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to GPL-3.0 license,
% as detailed in the original code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics). 
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% This code calls sub-component file 'cmap_viridis.txt'. The contents of
% the file includes a copy of data obtained from the repository 
% https://github.com/BIDS/colormap, and is CC0 1.0 licensed for modified
% use, see https://creativecommons.org/publicdomain/zero/1.0 for
% information.
%
% Checked by: Gil Felix Greco
% Date last checked: 16.02.2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Arguments validation
    arguments (Input) % Matlab R2018b or newer
        insig (:, :) double {mustBeReal}
        fs (1, 1) double {mustBePositive, mustBeInteger}
        fieldtype (1, :) string {mustBeMember(fieldtype,...
                                                       {'free-frontal',...
                                                        'diffuse'})} = 'free-frontal'
        time_skip (1, 1) double {mustBeReal} = 304e-3
        show {mustBeNumericOrLogical} = false
    end

%% Input checks

% define time threshold value from which all values before must be dropped.
t_threshold =  304e-3;

% check insig dimension (only [Nx1] or [Nx2] are valid)
if  size(insig,1) > 2 & size(insig,2) > 2 % insig has more than 2 channels
    error('Error: Input signal has more than 2 channels. ')
elseif  size(insig, 2) > 2  % insig is [1xN] or [2xN]
    insig = insig';
    fprintf('\nWarning: Input signal is not [Nx1] or [Nx2] and was transposed.\n');
end

% Check the length of the input data (must be at least 304 ms)
if size(insig, 1) <  t_threshold*fs
    error('Error: Input signal is too short along the specified axis to calculate loudness (must be at least 304 ms)')
end

% Check the channel number of the input data
if size(insig, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    chansIn = size(insig, 2);
    if chansIn > 1
        chans = ["Stereo left"; "Stereo right"];
        binaural = true;
    else
        chans = "Mono";
        binaural = false;
    end
end

if time_skip<t_threshold
        warning("Time_skip must be at least 304 ms to avoid transient responses of the digital filters (see ECMA-418-2:2025 (Section 8.1.4)). Setting time_skip to 304 ms!!!")
        time_skip = t_threshold;
end

%% Define constants

sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2025 [r_s]
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
c = 0.1618;  % Half-overlapping Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2025

dz = 0.5;  % critical band overlap [deltaz]
halfBark = 0.5:dz:26.5;  % half-overlapping band rate scale [z]
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025 [F(z)]

% Section 8.1.1 ECMA-418-2:2025
weight_n = 0.5331;  % Equations 113 & 114 ECMA-418-2:2025 [w_n]
% Table 12 ECMA-418-2:2025
a = 0.2918;
b = 0.5459;

% Output sample rate based on tonality hop sizes (Section 6.2.6
% ECMA-418-2:2025) [r_sd]
sampleRate1875 = sampleRate48k/256;

%% Signal processing

% get time vector of input signal
timeInsig = (0 : length(insig(:,1))-1) ./ fs;

% Input pre-processing
% --------------------
if fs ~= sampleRate48k  % Resample signal
    [p_re, ~] = shmResample(insig, fs);
else  % don't resample
    p_re = insig;
end

% Calculate specific loudnesses for tonal and noise components
% ------------------------------------------------------------

% Obtain tonal and noise component specific loudnesses from Sections 5 & 6 ECMA-418-2:2025
tonalitySHM = Tonality_ECMA418_2(p_re, sampleRate48k, fieldtype, time_skip, false);

specTonalLoudness = tonalitySHM.specTonalLoudness;  % [N'_tonal(l,z)]
specNoiseLoudness = tonalitySHM.specNoiseLoudness;  % [N'_noise(l,z)]

% Section 8.1.1 ECMA-418-2:2025
% Weight and combine component specific loudnesses
for chan = chansIn:-1:1

    % Equation 114 ECMA-418-2:2025 [e(z)]
    maxLoudnessFuncel = a./(max(specTonalLoudness(:, :, chan)...
                                + specNoiseLoudness(:, :, chan), [],...
                                2, "omitnan") + 1e-12) + b;

    % Equation 113 ECMA-418-2:2025 [N'(l,z)]
    specLoudness(:, :, chan) = (specTonalLoudness(:, :, chan).^maxLoudnessFuncel...
                                + abs((weight_n.*specNoiseLoudness(:, :, chan)).^maxLoudnessFuncel)).^(1./maxLoudnessFuncel);
end

if chansIn == 2 && binaural
    % Binaural loudness
    % Section 8.1.5 ECMA-418-2:2025 Equation 118 [N'_B(l,z)]
    specLoudness(:, :, 3) = sqrt(sum(specLoudness.^2, 3)/2);
    specTonalLoudness(:, :, 3) = sqrt(sum(specTonalLoudness.^2, 3)/2);
    specNoiseLoudness(:, :, 3) = sqrt(sum(specNoiseLoudness.^2, 3)/2);
    chansOut = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Combined binaural"];
else
    chansOut = chansIn;  % assign number of output channels
end

% time (s) corresponding with results output [t]
timeOut = (0:(size(specLoudness, 1) - 1))/sampleRate1875;

[~, time_skip_idx] = min( abs(timeOut-time_skip) ); % find idx of time_skip on timeOut
[~, idx_insig] = min( abs(timeInsig - time_skip) ); % find idx of time_skip on timeInsig

% Section 8.1.2 ECMA-418-2:2025
% Time-averaged specific loudness Equation 115 [N'(z)]
specLoudnessPowAvg = (sum(specLoudness(time_skip_idx:end, :, :).^(1/log10(2)), 1)./size(specLoudness(time_skip_idx:end, :, :), 1)).^log10(2); %<--- time index takes <time_skip> into consideration

% Section 8.1.3 ECMA-418-2:2025
% Time-dependent loudness Equation 116 [N(l)]
% Discard singleton dimensions
if chansOut == 1
    loudnessTDep = sum(specLoudness.*dz, 2);
    specLoudnessPowAvg = transpose(specLoudnessPowAvg);
else
    loudnessTDep = squeeze(sum(specLoudness.*dz, 2));
    specLoudnessPowAvg = squeeze(specLoudnessPowAvg);
end

% Section 8.1.4 ECMA-418-2:2025
% Overall loudness Equation 117 [N]
loudnessPowAvg = (sum(loudnessTDep(time_skip_idx:end, :).^(1/log10(2)), 1)./size(loudnessTDep(time_skip_idx:end, :), 1)).^log10(2); %<--- time index takes <time_skip> into consideration

%% Output assignment

% Assign outputs to structure
if chansOut == 3 % stereo case ["Stereo left"; "Stereo right"; "Combined binaural"];

    % outputs only with ["Stereo left"; "Stereo right"]
    OUT.specLoudness = specLoudness(:, :, 1:2);
    OUT.specTonalLoudness = specTonalLoudness(:, :, 1:2);
    OUT.specNoiseLoudness = specNoiseLoudness(:, :, 1:2);
    OUT.specLoudnessPowAvg = specLoudnessPowAvg(:, 1:2);
    OUT.loudnessTDep = loudnessTDep(:, 1:2);
    OUT.loudnessPowAvg = loudnessPowAvg(1:2);

    % outputs only with  "combined binaural"
    OUT.specLoudnessBin = specLoudness(:, :, 3);
    OUT.specTonalLoudnessBin = specTonalLoudness(:, :, 3);
    OUT.specNoiseLoudnessBin = specNoiseLoudness(:, :, 3);
    OUT.specLoudnessPowAvgBin = specLoudnessPowAvg(:, 3);
    OUT.loudnessTDepBin = loudnessTDep(:, 3);
    OUT.loudnessPowAvgBin = loudnessPowAvg(:, 3);

    % general outputs
    OUT.bandCentreFreqs = bandCentreFreqs;
    OUT.timeOut = timeOut;
    OUT.timeInsig = timeInsig;
    OUT.soundField = fieldtype;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loudness statistics based on loudnessTDep ["Stereo left"; "Stereo right"; "Combined binaural"];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    metric_statistics = 'Loudness_ISO532_1';
    OUT_statistics = get_statistics( loudnessTDep(time_skip_idx:end,1:chansOut), metric_statistics ); % get statistics

    % copy fields of <OUT_statistics> struct into the <OUT> struct
    fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

    for i = 1:numel(fields_OUT_statistics)
        fieldName = fields_OUT_statistics{i};
        if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
            OUT.(fieldName) = OUT_statistics.(fieldName);
        end
    end

    clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else % mono case

    OUT.specLoudness = specLoudness;
    OUT.specLoudnessPowAvg = specLoudnessPowAvg;
    OUT.loudnessTDep = loudnessTDep;
    OUT.loudnessPowAvg = loudnessPowAvg;

    OUT.bandCentreFreqs = bandCentreFreqs;
    OUT.timeOut = timeOut;
    OUT.timeInsig = timeInsig;
    OUT.soundField = fieldtype;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loudness statistics based on loudnessTDep (mono case)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    metric_statistics = 'Loudness_ECMA418_2';
    OUT_statistics = get_statistics( loudnessTDep(time_skip_idx:end), metric_statistics ); % get statistics

    % copy fields of <OUT_statistics> struct into the <OUT> struct
    fields_OUT_statistics = fieldnames(OUT_statistics);  % Get all field names in OUT_statistics

    for i = 1:numel(fields_OUT_statistics)
        fieldName = fields_OUT_statistics{i};
        if ~isfield(OUT, fieldName) % Only copy if OUT does NOT already have this field
            OUT.(fieldName) = OUT_statistics.(fieldName);
        end
    end

    clear OUT_statistics metric_statistics fields_OUT_statistics fieldName;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% Output plotting

if show

    % colormap
    cmap_viridis = load('cmap_viridis.txt');

    % generate A-weighting filter for LAeq calculation
    [b, a] = Gen_weighting_filters(fs, 'A');
    insig_A = filter(b, a, insig);  % filter signal
    LAeq_all = 20*log10(rms(insig_A(idx_insig:end, :))./2e-5);  % calculate LAeq

    for chan = chansOut:-1:1
        % Plot results
        fig = figure('name', sprintf( 'Loudness analysis - ECMA-418-2 (%s signal)', chans(chan) ) );
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specLoudness(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        %ax1.CLim = [0, ceil(max(specificLoudness(:, :, chan), [], 'all')*10)/10];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"];
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency (Hz)';
        ax1.XLabel.String = 'Time (s)';
        ax1.FontName =  'Times';
        ax1.FontSize = 11;
        colormap(cmap_viridis);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific Loudness,'; '(sone_{HMS}/Bark_{HMS})'});        

        
        if chan == 3 % the binaural channel

            % take the higher channel level as representative (PD ISO/TS 12913-3:2019 Annex D)
            [LAeq, LR] = max(LAeq_all);

            % if branch to identify which channel is higher
            if LR == 1
                whichEar = 'left ear';
            else
                whichEar = 'right ear';
            end  % end of if branch

        elseif chan == 2 % Stereo right

            LAeq = LAeq_all(chan);
            whichEar = 'right ear';

        elseif chan == 1 % Stereo left or mono

            LAeq = LAeq_all(chan);
            if chansOut~=1
                whichEar = 'left ear';
            else
                whichEar = 'mono';
            end
        end
         
        titleString = sprintf('%s signal, $L_{\\textrm{Aeq,%s}} =$ %.3g (dB SPL)', chans(chan), whichEar, LAeq);

        title(titleString, 'Interpreter','Latex' );

        ax2 = nexttile(2);
        plot(ax2, timeOut, loudnessTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold on;
        plot(ax2, timeOut, loudnessPowAvg(1, chan)*ones(size(timeOut)), '--', 'color',...
             cmap_viridis(34, :), 'LineWidth', 1, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];

        if max(loudnessTDep(:, chan)) > 0
            ax2.YLim = [0, 1.1*ceil(max(loudnessTDep(:, chan))*10)/10];
        end

        ax2.XLabel.String = 'Time (s)';
        ax2.YLabel.String = 'Loudness (sone_{HMS})';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Times';
        ax2.FontSize = 11;
        legend('Location', 'eastoutside', 'FontSize', 8);
        set(gcf,'color','w');
        clear LA; 
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting if outplot true

end % end of function
