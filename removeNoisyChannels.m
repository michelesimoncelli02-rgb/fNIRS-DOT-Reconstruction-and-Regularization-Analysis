function remCh = removeNoisyChannels(nirsData,dRange,SNRrange)
    
    % Compute mean of each channel and SNR
    meanValue = mean(nirsData);
    SNRValue = mean(nirsData)./std(nirsData,[],1);
    
    % Find channels within dRange and > SNR
    remCh = zeros(size(nirsData,2),1);
    remCh(meanValue > dRange(1) & meanValue < dRange(2) & SNRValue > SNRrange) = 1;
    
    % Channels should be removed even if only 1 wavelength is bad quality
    tmp = [remCh(1:end/2) remCh(end/2 + 1: end)];
    remCh = zeros(size(nirsData,2),1);
    remCh(sum(tmp,2)==2) = 1; % Keep only channels that have both wavelengths at 1
    remCh(end/2+1:end) = remCh(1:end/2); % Copy the decision (0 or 1) to the second wavelength

end