function [fig] = plot_rir(timeData, sampleRate, nfft, plotTitle, xLim)
% PLOT_RIR plot the magnitude (dB) of timeData and its spectrogram 
% Args: 
%   timeData (1D float vector)  : input sequence in time domain 
%   sampleRate (int)            : sampling rate of timeData
%   nfft (int)                  : number of fft points
%   plotTitile (string)         : title of the plot
% Returns:
%   fig (Figure)                : handle of the current figure

    fig = figure('Renderer', 'painters', 'Position', [10 10 900 300]); clf
    sgtitle(plotTitle,'Fontsize',14,'fontname','Times'); 
    timeAxis = [0:1:length(timeData)-1]./sampleRate;
    % compute spectrogram
    [stftCoeff, freqAxis, timeAxisSTFT] = spectrogram(timeData, ...
    nfft, round(nfft*0.75), [], sampleRate, 'yaxis');

    % plot magnitude in dB
    subplot(1, 2, 1); 
    plot(timeAxis, 20*log10(abs(timeData/max(abs(timeData)))),"Color", [0, 0, 0, 0.2])
    title('Room Impulse Response')
    set(gca,'Fontsize',12,'fontname','Times')
    xlabel('Time (s)'); ylabel('Magnitude (dB)')
    axis([0 xLim -80 0]) 

    % plot spectrogram
    subplot(1, 2, 2); 
    pcol =  pcolor(timeAxisSTFT, freqAxis, 10*log10(stftCoeff.*conj(stftCoeff)/max(max(stftCoeff.*conj(stftCoeff)))));
    title('Spectrogram')
    pcol.EdgeColor = 'none';
    axis([0 xLim 50 16000]) 
    set(gca, 'YDir', 'normal','Yscale', 'log')
    xlabel('Time (s)'); ylabel('Frequency (kHz)')
    set(gca,'Fontsize',12,'YTick', [100 300 1000 3000 10000], 'YTicklabel',{'0.1' '0.3' '1' '3' '10'},'fontname','Times')

    colormap(flipud(gray))
    caxis([-60 0])
    clb = colorbar;
    box on
    set(gca, 'layer','top')
 
    set(get(clb,'Title'),'String','Energy (dB)', 'Fontsize',12,'fontname','Times')
end

