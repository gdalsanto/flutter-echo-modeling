function [averagePower] = compute_power(timeData, sampleRange, isWhitened)
% COMPUTE_POWER return the average power of a segment of the input time
% series after applying a whitening filter 
% 
% Args: 
%   timeData (vector)           : time series of the RIR
%   sampleRange (vector)        : indexes of the samples from which the 
%                               power is computed 
%   isWhitened (bool)           : if true whiten the 
% Returns:
%   averagePower (float)        : average power of timeData(sampleRange) 

    % LP whitening filter
    lpOrder = 10;
    [lpcCoeff,~] = lpc(timeData(sampleRange), lpOrder);
    if isWhitened
        % Whiten the signal
        whitenedData = filter(lpcCoeff, 1, timeData(sampleRange));
        % compute power
        averagePower = mean(whitenedData.^2);
    else 
        averagePower = mean(timeData(sampleRange).^2);
    end
end

