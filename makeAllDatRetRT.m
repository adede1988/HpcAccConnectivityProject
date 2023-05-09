function [allDat] = makeAllDatRetRT(trials, times, RT, errorTrials, sampRate)



if sampRate ~= 1000
    %upsample the data
    for tt = 1:length(trials)
            [n,d] = rat(1000 / sampRate);
            test = resample(trials{tt}', n,d);
            L = size(test,1); 
            if L >2501
                extra = L-2501;
                test(1:floor(extra/2), :) = []; 
                test(end-ceil(extra/2)+1:end, :) = []; 
            end
            trials{tt} = test'; 

    end

end



%channels X time X trials
allDat = zeros(size(trials{1},1), 2501, length(trials));

for tt = 1:length(trials)
    curTrial = trials{tt}; 

    padDat = mirrorPad(curTrial);

    curTime = round(times{tt}*1000); 
    onset = find(curTime==RT(tt));
    L = size(curTrial,2);

    if isempty(onset) %the experimenter advanced the trial before the subject responded? I don't get it nan it out
        allDat(:,:,tt) = nan; 

    else
        
        try
            allDat(:,:,tt) = padDat(:,L+onset-2000:L+onset+500); 
        catch %likely error is that even padded the data are not long enough! 
            padDat2 = mirrorPad(padDat); 
            allDat(:,:,tt) = padDat2(:,(L*3)+onset-2000:(L*3)+onset+500);
        end


    end

end



allDat(:,:,errorTrials) = []; 





end