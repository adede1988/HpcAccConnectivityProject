function [curDatSum] = HFBSummary(curDat, updown)


curDatSum = struct; 
ci = 1; 

%hardCode the time periods for threshold testing
threshTestWin = [[-50 2500]; [-1500, 500]; [-50 2000]; [-1500, 500]];



HFBConditions = fieldnames(curDat(1).HFB); 

sizeTest = reshape(cell2mat(arrayfun(...
    @(x) size(curDat(1).HFB.(HFBConditions{x})), 1:length(HFBConditions) ,...
    'UniformOutput',false)), [2,length(HFBConditions)])';

%the first short Idx in each block is the time index
shortIdx = find(sizeTest(:,1)==1);
timIdx = shortIdx(find(diff(shortIdx)==1));
timIdx = [1; timIdx]; 



for cond = 1:length(HFBConditions)

    if sizeTest(cond,1)>1 %this is data! 
        curDatSum(ci).condition = HFBConditions{cond}; 
        curDatSum(ci).HFB = zeros(size(curDat(1).HFB.(HFBConditions{cond}), 1 ), 10000); 
        curDatSum(ci).RT = zeros(10000,1); 
        curDatSum(ci).subID = cell(10000,1); 
        curDatSum(ci).chi = zeros(10000,1); 
        curDatSum(ci).time = curDat(1).HFB.(HFBConditions{timIdx(find(timIdx>cond,1))}); 
        curDatSum(ci).peakLat = zeros(10000,1);
        curDatSum(ci).peakAmp = zeros(10000,1);
        curDatSum(ci).threshCross = zeros(10000,1); 
        curDatSum(ci).maxAmp = zeros(10000,1); 
        curDatSum(ci).amp = zeros(10000,1); 
        curDatSum(ci).dur = zeros(10000,1); 
        curDatSum(ci).centerOfMass = zeros(10000,1); 
%         'here'
        ti = 1; %trial summary iterator
        %now loop over channels and actually get the data! 
        for ii = 1:length(curDat)
%             ii
            comboHFB = []; 
            %hard code hack: 
            add = 0; 
            if ci>2
                add = 2; 
            end
            for comboi = timIdx(find(timIdx>cond,1)-1)+add : timIdx(find(timIdx>cond,1))-1
                comboHFB = [comboHFB curDat(ii).HFB.(HFBConditions{comboi})];
            end
            test = mean(comboHFB,2); 
          
    
            if curDat(ii).HFBenc==updown || curDat(ii).HFBretOn ==updown || curDat(ii).HFBretRT==updown || curDat(ii).HFBencRT==updown
                test = true;
            else
                test = false; 
            end
            if test %only add in channels that have good threshold crossings! 
                
                L = size(curDat(ii).HFB.(HFBConditions{cond}),2);

                curDatSum(ci).HFB(:,ti:ti+L-1) = curDat(ii).HFB.(HFBConditions{cond});
                
                %get the RT
                if ~isempty(strfind(HFBConditions{cond}, 'sub')) %this is encoding data! 
                    if strfind(HFBConditions{cond}, 'Miss')>0 %miss trials
                        curDatSum(ci).RT(ti:ti+L-1) = curDat(ii).encInfo(curDat(ii).encInfo(:,1)==1 &...
                                                                         curDat(ii).encInfo(:,2)==2, 4);
                    else %hit trials
                        curDatSum(ci).RT(ti:ti+L-1) = curDat(ii).encInfo(curDat(ii).encInfo(:,1)==1 &...
                                                                         curDat(ii).encInfo(:,2)==1, 4);
                    end
                else %this is response data
                    if strfind(HFBConditions{cond}, 'hit')>0 %hit trials
                        curDatSum(ci).RT(ti:ti+L-1) = curDat(ii).retInfo(curDat(ii).retInfo(:,1)==1, 3);
                    elseif strfind(HFBConditions{cond}, 'miss')>0 %miss trials
                        curDatSum(ci).RT(ti:ti+L-1) = curDat(ii).retInfo(curDat(ii).retInfo(:,1)==2, 3);
                    elseif strfind(HFBConditions{cond}, 'cr')>0 %cr trials
                        curDatSum(ci).RT(ti:ti+L-1) = curDat(ii).retInfo(curDat(ii).retInfo(:,1)==3, 3);
                    elseif strfind(HFBConditions{cond}, 'fa')>0 %fa trials
                        curDatSum(ci).RT(ti:ti+L-1) = curDat(ii).retInfo(curDat(ii).retInfo(:,1)==4, 3);
                    end
                end
                
                %loop over trials to get the summary stats
                for tt = 1:L
%                     tt
                    curTrial = curDatSum(ci).HFB(:,ti+tt-1);
                    %record subID
                    curDatSum(ci).subID{ti+tt-1} = curDat(ii).subID; 
                    curDatSum(ci).chi(ti+tt-1) = curDat(ii).chi; 
                    winVals = threshTestWin(find(timIdx>cond,1)-1,:); 

                    %peakLat
                    if ci == 1 || ci ==2 || ci == 5 || ci ==6 || ci ==7 || ci == 8
                        if curDatSum(ci).RT(ti+tt-1) > max(curDatSum(ci).time)
                            curDatSum(ci).RT(ti+tt-1) = max(curDatSum(ci).time);
                        end
                        endTim = curDatSum(ci).RT(ti+tt-1);
                    else
                        endTim = winVals(2); 
                    end
                    [maxVal, peakLat] = max(curTrial(find(curDatSum(ci).time>=winVals(1),1):find(curDatSum(ci).time>=endTim,1) ) );
                    peakLat = peakLat + find(curDatSum(ci).time>=winVals(1),1) - 1; 
                    curDatSum(ci).peakLat(ti+tt-1) = peakLat; 

                    %peakamp

                    curDatSum(ci).peakAmp(ti+tt-1) = curTrial(peakLat); 
                    
                    %center of mass
%                     if ci == 7
%                         'here'
%                     end
                    test = curTrial( find(curDatSum(ci).time>=winVals(1),1) : find(curDatSum(ci).time>=endTim,1) );
                    testTim = curDatSum(ci).time(find(curDatSum(ci).time>=winVals(1),1) : find(curDatSum(ci).time>=endTim,1));
                    test = (test - min(test)); 
                    test = test./ max(test);
                    test(test<.8) = 0; 
                    testMean = wmean([1:length(testTim)]', test); 
                    curDatSum(ci).centerOfMass(ti+tt-1) = testTim(round(testMean));

                    curDatSum(ci).maxAmp(ti+tt-1) = maxVal; 
                    if maxVal>1.96
                        
    
                        %find first thresh cross
                        curDatSum(ci).threshCross(ti+tt-1) =  find(curTrial(find(curDatSum(ci).time>=winVals(1),1):...
                            find(curDatSum(ci).time>=endTim,1)) > 1.96, 1) + find(curDatSum(ci).time>=winVals(1),1);
    
                        %find amp
                        start = [flip(curTrial(1:peakLat-1)); 0];
                        startHalf = sum(start(1:find(start<1.96,1)));
                        back = [curTrial(peakLat:end); 0];
                        back = back(1:find(back<1/96,1)); 
                        backHalf = sum(back); 
                        curDatSum(ci).amp(ti+tt-1) = startHalf + backHalf; 
                        
                        %find duration
                        startHalf = find(start<1.96,1);
                        backHalf = find(back<1/96,1);
                        if isempty(startHalf)
                            curDatSum(ci).dur(ti+tt-1) = backHalf; 
                        elseif isempty(backHalf)
                            curDatSum(ci).dur(ti+tt-1) = startHalf;
                        else
                            curDatSum(ci).dur(ti+tt-1) = startHalf + backHalf; 
                        end
                    else
                        curDatSum(ci).amp(ti+tt-1) = nan;
                        curDatSum(ci).threshCross(ti+tt-1) = nan;
                        curDatSum(ci).dur(ti+tt-1) = nan; 

                    end

                end
                ti = ti+L; 

            end
        end

        %trim zeros
%         outnames = fieldnames(curDatSum); 
        test = curDatSum(ci).HFB(1,:); 
        cutPoint = 10000 - find(flip(test)~=0,1) + 2;
        
        curDatSum(ci).HFB(:,cutPoint:end) = []; 
        curDatSum(ci).RT(cutPoint:end) = []; 
        curDatSum(ci).subID(cutPoint:end) = []; 
        curDatSum(ci).chi(cutPoint:end) = []; 
        curDatSum(ci).peakLat(cutPoint:end) = []; 
        curDatSum(ci).peakAmp(cutPoint:end) = [];
        curDatSum(ci).threshCross(cutPoint:end) = []; 
        curDatSum(ci).maxAmp(cutPoint:end) = []; 
        curDatSum(ci).amp(cutPoint:end) = []; 
        curDatSum(ci).dur(cutPoint:end) = []; 
        curDatSum(ci).centerOfMass(cutPoint:end) = [];


     
        disp(ci)


        ci = ci+1; 

    end

end


















end