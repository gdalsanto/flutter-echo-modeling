function [initPulseWT] = pulse_shape(timeData, repetitionTime, ...
    lateTime, initSilence, onsetTime, slopeRT, sampleRate, plotPulse)
% PULSE_SHAPE return the estimated wavelet tranform of the flutter initial 
% pulse shape at time onsetTime
% 
% Args: 
%   timeData (vector)           : time series of the reference RIR
%   repetitionTime (float)      : repetition time in sec
%   lateTime (float)            : late reveberation time in sec
%   initSilence (int)           : sample index of the ref RIR onset 
%   onsetTime (float)           : target onset time of the flutter echo in sec
%   slopeRT (vector)            : decay slope siimulated by the attenuation 
%                               filter
%   sampleRate (int)            : sampling rate of timeData
% Returns:
%   repetitionTime (Figure)     : repetition time estimated from the EDC 
%   flutteryRange ()
%   lateRevTime (float)         : late reverberation time in sec
    
    samplesTR = round(repetitionTime*sampleRate);   % tR in samples

    % number of instances of flutter from tL and onset time
    periodsTR = floor(lateTime/repetitionTime);
    % start and end time of one period of flutter in the late reverberation
    startSample = initSilence + periodsTR*samplesTR;
    pulsesIndx = [startSample (startSample+samplesTR-1)]; % init and end time of out pulse 


    % define time range in samples of the rir to be analyse
    startSample = initSilence + periodsTR*samplesTR;
    timeRange = 0:5*samplesTR;
    timeRange = timeRange + startSample;
    % take a window of 5 tR seconds 
    lateData = timeData(timeRange);
    
    % MULTIRESOLUTION ANALYSIS 
    % use Daubechies 12 (db12) wavelet
    lev = 4;    % transform level
    levSel = 4; % chosen transform level
    % compute maximal overlap discrete wavelet transform MODWT 
    waveT = modwt(lateData,'db11',lev);
    mra = modwtmra(waveT,'db11');
    
    % normalize amplitude by removing the exponential envelope
    [envUp, envLow] = envelope(mra(levSel, :),samplesTR-10,'peak');
    negEnvUp = 1./(envUp./max(abs(envUp)));
    pulseTime = negEnvUp(1:samplesTR)'.*mra(levSel, pulsesIndx(1)-timeRange(1)+1:pulsesIndx(2)-timeRange(1)+1)';
    % compute fft of MRA 
    nfft = 2*samplesTR;
    fftRIR = fft(pulseTime, nfft);
    freqAxis = 0:sampleRate/nfft:sampleRate-sampleRate/nfft; 
    timeAxis = 0:1/sampleRate:repetitionTime;

    % to avoid that wrong that nosiy slope values may affect the estimation 
    % remove frequency coefficents that are out of the frequency range of 
    % the chosen wavelet
    minFreq = 1000;
    maxFreq = 3200;
    freqStart = find(freqAxis>=minFreq,1,'first'); 
    freqEnd = find(freqAxis>=maxFreq,1,'first');
    validRange = freqStart:freqEnd; 
    
    % derive the freq magnitude of the pulse at the beginning 
    initMag = 10*log10(fftRIR.*conj(fftRIR));
    
    % generate window to avoid abrupt chances at the edges
    temp = hamming(101);
    win = [temp(1:50); ones(length(validRange)-101,1); temp(51:end)];

    initMag(validRange) = initMag(validRange)+(-slopeRT(validRange)).* ...
        (pulsesIndx(1) - (initSilence+onsetTime))/sampleRate.*win; % hamming(length(validRange));

    symValidRange = length(fftRIR) - validRange;

    initMag(symValidRange) = initMag(validRange);
    initPhase = atan2(imag(fftRIR), real(fftRIR));
    initPulseWT = real(ifft(10.^(initMag./20).*exp(1i*initPhase), nfft));
    initPulseWT = initPulseWT(1:samplesTR);
    % shift the pulse 
    fractTR = rem(pulsesIndx(1) - initSilence,samplesTR);
    
    initPulseWT = circshift(initPulseWT,fractTR);

    % ------ plot ------- %
    if plotPulse
        figure('Position', [0 0 500 300]); clf
        
        plot(timeAxis+repetitionTime, initPulseWT/max(abs(initPulseWT)),"Color", [0, 0, 0, 1],'LineWidth',0.1)
        title('Pulse shape - MRA')
        xlabel('Time (s)'); ylabel('Amplitude')
        set(gca,'Fontsize',12,'fontname','Times')
        axis([repetitionTime 2*repetitionTime -1.1 1.1])
    end
end

