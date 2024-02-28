function [allDat] = makeAllDat(trials, times, sampRate)


%error here: probably need to be upsampling the time as well. IMPLEMENTED
if sampRate ~= 1000
    %upsample the data
    for tt = 1:length(trials)
            [n,d] = rat(1000 / sampRate);
            test = resample(trials{tt}', n,d);
            upTim = resample(times{tt}, n,d); 
            trials{tt} = test'; 
            times{tt} = upTim; 

    end

end

%channels X time X trials
allDat = nan(size(trials{1},1),  4501, length(trials));

for tt=1:length(trials)
    curTrial = trials{tt}; 

    padDat = mirrorPad(curTrial);

    curTime = round(times{tt}*1000); 
    onset = find(curTime>=0,1);
    L = size(curTrial,2);

    try
        allDat(:,:,tt) = padDat(:,L+onset-1000:L+onset+3500); %t = 1001 is stim onset! 
    catch %likely error is that even padded the data are not long enough! 
        padDat2 = mirrorPad(padDat); 
        allDat(:,:,tt) = padDat2(:,(L*3)+L+onset-1000:(L*3)+L+onset+3500);

    end


end



end
