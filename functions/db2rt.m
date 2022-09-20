function [RT] = db2rt(gdB, len, sampleRate)
% Convert comb filter gains to RT values
%   gdB:    Comb filter gains
%   len:    Length of delay line 
%   sampleRate: sampling rate in Hz

    gdec = 10.^(gdB/20);             % db to linear 
    gains = gdec.^(sampleRate/len); % frequency dependent gain
    decdB = -20*log10(gains);       % decay rate in dB/s
    RT = 60./decdB ;                % reverberation time
end

