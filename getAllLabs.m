function allBrod = getAllLabs(allDat)


    allBrod = struct; 
    allBrod(1).lab = 'mtl'; 
    allBrod(1).count = 0; 
    allBrod(1).subN = 0; 
    allBrod(1).subi = []; 
    
    for ii = 1:length(allDat)
        if ~isempty(allDat{ii})
            curDat = allDat{ii}; 
            curMeet = curDat.meetLabs(:,3); 
            for chi = 1:length(curMeet)
                if isempty(curMeet{chi})
                    curMeet{chi} = 'ZZZ'; 
                end
            end
            knownLabs = {allBrod.lab}; 
            curUni = unique(curMeet); 
            for Li = 1:length(curUni)
                prevComp = cellfun(@(x) strcmp(curUni{Li}, x), knownLabs); 
                if sum(prevComp)>0 %this label is already in the data
                    areaCount = cellfun(@(x) strcmp(curUni{Li}, x), curMeet); 
                    allBrod(prevComp).count = allBrod(prevComp).count + sum(areaCount); 
                    allBrod(prevComp).subi = [allBrod(prevComp).subi, ii];
                    allBrod(prevComp).subN = allBrod(prevComp).subN + 1; 

                else
                    curL = length(allBrod)+1; 
                    allBrod(curL).lab = curUni{Li}; 
                    areaCount = cellfun(@(x) strcmp(curUni{Li}, x), curMeet); 
                    allBrod(curL).count = sum(areaCount); 
                    allBrod(curL).subi = [ii];
                    allBrod(curL).subN = 1; 

                end


            end
    
        end
    end














end