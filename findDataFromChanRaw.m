function [curOut] = findDataFromChanRaw(chanFiles, chan, IDs)


    curChan = load([chanFiles(chan).folder '/' chanFiles(chan).name]).chanDat; 

    if sum(cellfun(@(x) strcmp(x, curChan.subID), IDs)) == 0

        if strcmp('NO NOTES', curChan.labels{1,1})
%             disp('skip')
            curOut = []; 
    
    
        else
    
            curOut = struct; 
    
            curOut.site = curChan.site; 
            curOut.subID = curChan.subID;
            
            %get other metadata and behavioral data from curChan
    
            phgCheck = cellfun(@(x) strcmp(x, 'phg'), curChan.labels(:, 3)); 
            
    %         curLabs = curChan.labels(:,3); 
    %         phgCheck = zeros(length(curChan.labels(:,3)),1); 
    %         for ii = 1:length(phgCheck)
    %             phgCheck(ii) = strcmp(curLabs{ii}, 'phg');
    %         end
    
            test = split(curChan.dataDir, ':');
            rawDat = load(['Z:' test{2} '/' curChan.encDatFn]).data; 
            
            curOut.phgChans = rawDat.elec.label(phgCheck==1);
            
    
    
        end


    else %we've already seen this one

        curOut = []; 


    end

%     curOut = 'hello'; 


end