function [allDat] = makeAllDat(trials, times)


%channels X time X trials
allDat = nan(size(trials{1},1),  4501, length(trials));

for tt=1:length(trials)
    curTrial = trials{tt}; 

    padDat = mirrorPad(curTrial);

    curTime = round(times{tt}*1000); 
    onset = find(curTime==0);
    L = size(curTrial,2);

    try
        allDat(:,:,tt) = padDat(:,L+onset-1000:L+onset+3500); %t = 1001 is stim onset! 
    catch %likely error is that even padded the data are not long enough! 
        padDat2 = mirrorPad(padDat); 
        allDat(:,:,tt) = padDat2(:,(L*3)+onset-1000:(L*3)+onset+3500);

    end


end



end