function [datForPerm, condCode] = getPairedDat(varNames, allDat, ii, pairi)


    curHave = arrayfun(@(x) ~isempty(allDat(x).(varNames{ii})),  1:length(allDat));
    pairHave = arrayfun(@(x) ~isempty(allDat(x).(varNames{pairi})),  1:length(allDat));
    
   
    condCode = [ones(sum(curHave),1); zeros(sum(pairHave),1)];
    
    curidx = find(curHave); 
    pairidx = find(pairHave); 
    test = nan; 
    ti = 1; %test index to find the first non nan set to initialize with
    while isnan(test)
        datForPerm = allDat(curidx(ti)).(varNames{ii}); 
        test = datForPerm(1); 
        if isnan(test)
            condCode(ti) = 999; 
        end
        ti = ti+1; 
    end
    dimVal = length(size(datForPerm)); 
    %loop forward grabbing condition 1 data 
    for curi = ti:length(curidx)
        temp = allDat(curidx(curi)).(varNames{ii});

        if isnan(temp(1,1))
            
            condCode(curi) = 999;
        else
            datForPerm = cat(dimVal+1, datForPerm, temp);
        end
    end
    for pri = 1:length(pairidx)
        temp = allDat(pairidx(pri)).(varNames{pairi});
        if isnan(temp(1,1))
            condCode(curi+pri) = 999; 
        else
            datForPerm = cat(dimVal+1, datForPerm, temp);
        end
    end

    condCode(condCode==999) = []; 


end