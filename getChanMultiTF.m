function [trialTF, outTim, frexOut] = getChanMultiTF(trialDat, frex, srate, times, timeStepms)

%trialDat:            timepoints X trials data points



%initialize indexing
timewinms = 300; %time window in milliseconds
% timeStepms = 25; %steps in terms of milliseconds independent of sampling rate 
timewin = round(timewinms / (1000/srate)); %timewin length in data point steps
hz = linspace(0, srate, timewin);
numfrex = sum(hz>min(frex) & hz<max(frex));
L = size(trialDat,1); 
ii = 1; %index variable for output
winStart = L - round(timewinms/(1000/srate)/2); %starting data index
[curTim, timi] = min(times); %starting time and time index

%pad the data to prevent edge artifacts
padDat = mirrorPad(trialDat); 

%initialize output variables
trialTF = zeros([size(trialDat), numfrex]);
outTim = zeros(size(times)); 
frexOut = hz(hz>min(frex) & hz<max(frex));


%prepare the tapers
curDat = padDat(winStart:winStart+timewin-1, :);
tapers    = dpss(size(curDat,1),5);
tapers = tapers(:,1:end-1); %discard the last taper


check = true; 

while check
    curDat = padDat(winStart:winStart+timewin-1, :);
    tapSpect = zeros(size(tapers,2),timewin, size(curDat,2)); 
    for tapi = 1:size(tapers,2)
        test = abs(fft(curDat.*tapers(:,tapi))).^2;
        tapSpect(tapi,:,:) = test(1:timewin,:); 
    end
    
    %store to the outputs
    spect = squeeze(mean(tapSpect(:,hz>min(frex) & hz<max(frex),:), 1));
    trialTF(ii,:,:) = spect';
    outTim(ii) = curTim; 

    %update index variables
    ii = ii+1; 
    timStep = find(times(timi:end)>=curTim+timeStepms,1)-1;
    if  isempty(timStep) %out beyond the end of time
        check = false; 
    else
        timi = timi+timStep; 
        curTim = times(timi);
        winStart = winStart+timStep; 
    end
 
end

%trim down outTim and trialTF
remove = find(outTim==0);
if length(remove)>1
if remove(1) - remove(2) < -1 %if the outTim values include precisely time zero, don't lose that!
    remove(1) = []; 
end

outTim(remove) = []; 
trialTF(remove, :, :) = []; 

end


% 
% 
% 
% winStarts = [1:timeStepms:size(zeroPadTD,1)-timewin];
% %time resolution  = 5ms (.05 seconds)
% 
% %time X trials X freq
% trialTF = zeros([length(winStarts), size(trialDat,2), numfrex]); 
% outTim = times([1:timeStepms:length(times)]);
% 
%  
% 
% 
% 
% for tt = 1:length(winStarts)
%     curDat = zeroPadTD(winStarts(tt):winStarts(tt)+timewin-1, :);
%     curDat = detrend(curDat); 
%     tapers    = dpss(size(curDat,1),5);
%     tapers = tapers(:,1:end-1); %discard the last taper
%     tapSpect = zeros(size(tapers,2),timewin, size(curDat,2)); 
%     for tapi = 1:size(tapers,2)
%         test = abs(fft(curDat.*tapers(:,tapi))).^2;
%         tapSpect(tapi,:,:) = test(1:timewin,:); 
%     end
%     
%     spect = squeeze(mean(tapSpect(:,hz>min(frex) & hz<max(frex),:), 1));
%     trialTF(tt,:,:) = spect'; 
% 
%     
% 
% end
%     
% 
% 





end