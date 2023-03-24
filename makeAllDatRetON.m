function [allDat] = makeAllDatRetON(trials, times)


%channels X time X trials
allDat = zeros(size(trials{1},1), 4001, length(trials));

for tt = 1:length(trials)
    curTrial = trials{tt}; 

    padDat = mirrorPad(curTrial);

    curTime = times{tt}; 
    onset = find(curTime==0);
    L = size(curTrial,2);
    try
        allDat(:,:,tt) = padDat(:,L+onset-1000:L+onset+3000); %t = 1001 is stim onset! 
    catch %likely error is that even padded the data are not long enough! 
        padDat2 = mirrorPad(padDat); 
        allDat(:,:,tt) = padDat2(:,(L*3)+onset-1000:(L*3)+onset+3000);

    end

end









end