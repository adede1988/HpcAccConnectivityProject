function [errorLog] = conditionFileAudit(fileInfo, RTerrors)

cndDat = load([fileInfo.folder '/' fileInfo.name]).cndDat;
% frex = logspace(log10(2),log10(80),100);
% numTim = size(cndDat.datForPerm,1); 
% numFrex = size(cndDat.datForPerm,2); 
%check in allDat for missing data
load("C:\Users\dtf8829\Documents\QuestConnect\allDat.mat")
errorLog = struct; 
ei = 1; 

varName = cndDat.targVar; 
pairVar = cndDat.pairVar; 

IDs = {cndDat.metaDat.subID};

IDs_uni = unique(IDs);

%each subject should be observed twice since there are two conditions...not
%necessarily. It's okay for a subject to only have one of two conditions.
for sub = 1:length(IDs_uni)
   
    if sum(cellfun(@(x) strcmp(IDs_uni(sub), x), IDs))==1
        errorLog(ei).subID = IDs_uni(sub);
        errorLog(ei).subID = errorLog(ei).subID{1};

        border = find(diff(cndDat.condCode));
        %is it missing the first or second variable? 
        if find(cellfun(@(x) strcmp(IDs_uni(sub), x), IDs)) > border
            errorLog(ei).missingVar = varName; 
        else
            errorLog(ei).missingVar = pairVar; 
        end
        subi = find(cellfun(@(x) strcmp(IDs_uni(sub), x), IDs), 1);
        %check the ROI matrix to see if the issue is caused there
        if sum(sum(cndDat.metaDat(subi).roiNote))==0
            roiMat = cndDat.metaDat(subi).roimni; 
            errorLog(ei).anatUsed = 'MNI'; 
        else
            roiMat = cndDat.metaDat(subi).roiNote; 
            errorLog(ei).anatUsed = 'Note'; 
        end

        %which ROIs are involved in the missing data? 
        splitMiss = split(errorLog(ei).missingVar, '_');
        rois = {'dlPFC', 'hip', 'phg', 'acc'}; 
        curTarg = [0,0,0,0]; 
        anat = find(sum(cell2mat(cellfun(@(x) contains(rois, x), splitMiss, 'uniformoutput', false)))>0);
        
        if sum(cell2mat(cellfun(@(y) sum(cell2mat(cellfun(@(x) strcmp(x, y), splitMiss, 'uniformOutput', false))), rois, 'uniformOutput', false))) ==1
            errorLog(ei).type = "power"; 
        else
            errorLog(ei).type = "connect"; 
        end

        errorLog(ei).curVarRois = anat; 

        allIDs = {allDat.subID}; 
        %check if data are missing in allDat
        alli = find(cellfun(@(x) strcmp(IDs_uni(sub), x), allIDs), 1);
        errorLog(ei).allDati = alli; 
        test = allDat(alli).(errorLog(ei).missingVar); 
        RTi = find(cellfun(@(x) strcmp(IDs_uni(sub), x), {RTerrors.subID}), 1);



        if isempty(test)
            errorLog(ei).errorType = "missing"; 
            if contains(errorLog(ei).missingVar, 'fa_') && sum(allDat(alli).retInfo(:,1)==4)<=5
                errorLog(ei).output = "insufficient FA trials"; 
            elseif contains(errorLog(ei).missingVar, 'cr_') && sum(allDat(alli).retInfo(:,1)==3)<=5
                errorLog(ei).output = "insufficient CR trials"; 
            elseif contains(errorLog(ei).missingVar, 'hit_') && sum(allDat(alli).retInfo(:,1)==1)<=5
                errorLog(ei).output = "insufficient HIT trials"; 
            elseif contains(errorLog(ei).missingVar, 'miss_') && sum(allDat(alli).retInfo(:,1)==2)<=5
                errorLog(ei).output = "insufficient MISS trials"; 
            elseif contains(errorLog(ei).missingVar, 'Miss_') && sum(allDat(alli).retInfo(:,1)==2)<=5
                errorLog(ei).output = "insufficient MISS trials"; 
            elseif contains(errorLog(ei).missingVar, 'Hit_') && sum(allDat(alli).retInfo(:,1)==1)<=5
                errorLog(ei).output = "insufficient HIT trials"; 
            end



        elseif isnan(test(1))
            errorLog(ei).errorType = "nanVals"; 
            if RTerrors(RTi).numTrials>0
                errorLog(ei).output = "RT mistake"; 
            elseif min(cndDat.metaDat(subi).retInfo(:,3)) < 0
                errorLog(ei).output = "negative RT"; 
            end
        elseif test(1) == 0 
            errorLog(ei).errorType = "zero";
            if strcmp(errorLog(ei).anatUsed, 'MNI') && strcmp(errorLog(ei).type, "connect")
                errorLog(ei).output = "doubleCount";
            end
        elseif min(cndDat.metaDat(subi).retInfo(:,3)) < 0
            errorLog(ei).errorType = "value available"; 
            errorLog(ei).output = "negative RT"; 
        end
        
        if ~isfield(errorLog, 'output')
            errorLog(ei).output = "NOT KNOWN!!!!"; 
        elseif isempty(errorLog(ei).output)
            errorLog(ei).output = "NOT KNOWN!!!!"; 
        end
     

        ei = ei+1; 
    end

end









end


