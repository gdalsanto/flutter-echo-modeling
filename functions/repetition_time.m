function [repetitionTime, flutteryRange, lateRevTime] = repetition_time( ...
    filename, octaveBand, expectedRepetitionTime)
% REPETITION_TIME return the repetition time of the flutter echo.
% This function uses the energy decay curve to find the time spanned by the
% series of pulses covering a round trip path between the reflective 
% parallel surfaces [1].
% Requires ita_toolbox [2].
% 
% Args: 
%   filename (string)           : filename of room impulse response 
%   octaveBand (int)            : octave band from which the repetition
%                               time will be estimated
%   expectedRepetitionTime (float)  : expected repetition time in sec
% Returns:
%   repetitionTime (float)      : repetition time estimated from the EDC in sec
%   flutteryRange (int)         : vector of octave band central frequnecies
%                               containing flutter
%   lateRevTime (float)         : late reverberation time in sec
%
% Reference: 
% [1] G. Dal Santo, K. Prawda, and V. Välimäki, "Flutter Echo 
% Modeling", in Proc. 25rd Int. Conf. Digital Audio Effects (DAFx20in22), 
% Vienna, Austria, Sept. 6–10
% [2] Dietrich, M. Guski, J. Klein, M. Müller-Trapet, M. Pollow,
% R. Scharrer, and M. Vorlaender, "Measurements and room
% acoustic analysis with the ITA-Toolbox for MATLAB," Jan. 2013.

    % compute EDC using ITA Toolbox
    itaResults = ita_read(filename);
    raResultsEDC = ita_roomacoustics(itaResults, 'EDC', 'freqRange', ...
        [100 20000], 'bandsPerOctave', 1);
    sampleRate = itaResults.samplingRate;

    % extract octave band of interest
    timeDataBand = raResultsEDC.EDC.timeData(:,octaveBand);
    nSamplesEDC = length(timeDataBand);     % number of samples
    timeAxis = (0:nSamplesEDC-1)/sampleRate;
    
    % find late reveberation time as defined in [1]
    dBVal = -30;    
    for n = 1:size(raResultsEDC.EDC.time,2)
        EDCurve = 10*log10(raResultsEDC.EDC.timeData(:,n));
        % find the time where the curve goes below dBVal
        freqString = raResultsEDC.EDC.channelNames{n};
        freqBands(n) = str2double(freqString(1:end-3));
        idx = find(EDCurve < dBVal);
        if isempty(idx)
            endFrame = nFrames; % last frame
        else 
            endFrame = idx(1); % first frame idx after dBfit dB of decay
        end
        % estimate the energy decay slope
        segmentEDC = EDCurve(1:endFrame);
        coeffPoly = polyfit(timeAxis(1:endFrame), segmentEDC, 1);
        slopes(n) = coeffPoly(1);
    end

    th = mean(slopes);
    flutteryRange = freqBands(slopes>th);
    
    % find time where the fluttery octaves bands decay by -30dB
    energyDCs = 10*log10(raResultsEDC.EDC.timeData(:,slopes>th)); 
    for n = 1:sum(slopes>th)
        tLs(n) = find(energyDCs(:,n)>-30,1,'last');
    end
    % late reverberation point
    lateRevTime = mean(tLs)/sampleRate;
    
    Dvec = gradient(10*log10(timeDataBand));
    [~, locs] = findpeaks(-Dvec(~isnan(Dvec)), ...
        'MinPeakDistance',(1-0.25)*expectedRepetitionTime*sampleRate);

    % take average step length 
    locsRange = locs(locs >= lateRevTime*sampleRate & locs <= ...
        (lateRevTime+10*expectedRepetitionTime)*sampleRate);
    if (length(locsRange)/2)
        locsRange(end) = [];
    end
    
    repetitionTime = mean(diff(locsRange))/sampleRate;

    % ------ plot ------- %
    figure('Renderer', 'painters', 'Position', [10 10 900 300]); clf
    
    subplot(1,2,1);
    plot(timeAxis, -Dvec,"Color", [0, 0, 0, 0.75],'LineWidth',0.1)
    title('Estimated derivative')
    xlabel('Time (s)'); ylabel('Amplitude')
    xline(timeAxis(locs),'--', "LineWidth", 0.5,"Color", [0.3, 0.3, 0.3], "Alpha", 0.5)
    line([timeAxis(locs(13)) timeAxis(locs(14))], [0.01 0.01 ],"Color", [0, 0, 0], "Marker",'|',"LineWidth",1.5)
    text((timeAxis(locs(14))+timeAxis(locs(13)))/2-0.005,0.011,'$t_\textrm{r}$','interpreter','latex')
    axis([0.53 (0.53+10*0.0458) 0 0.015])
    set(gca,'Fontsize',12,'fontname','Times')
    
    subplot(1,2,2);
    plot(timeAxis, 10*log10(timeDataBand),"Color", [0, 0, 0, 1])
    title('Enerrgy decay curve')
    xline(timeAxis(locs),'--', "LineWidth", 0.5,"Color", [0.3, 0.3, 0.3], "Alpha", 0.5)
    xlabel('Time (s)'); ylabel('Energy (dB)')
    axis([0.53 (0.53+10*0.0458) -40 -25])
    set(gca,'Fontsize',12,'fontname','Times')
    
end

