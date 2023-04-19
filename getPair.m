function [pairi, pairName] = getPair(curVar, varNames, ii)


  if contains(curVar, 'sub')
        %subsequent memory
        if contains(curVar, 'Miss')
            splitVar = split(curVar, 'Miss');
            temp = varNames; 
            temp{ii} = 'not this one'; 
            idx = find(cellfun(@(x) contains(x, splitVar{1}) & contains(x, splitVar{2}), temp )); 
            [~, idxi] = min(abs(idx-ii)); 
            pairi = idx(idxi);

        elseif contains(curVar, 'Hit')
            splitVar = split(curVar, 'Hit');
            temp = varNames; 
            temp{ii} = 'not this one'; 
            idx = find(cellfun(@(x) contains(x, splitVar{1}) & contains(x, splitVar{2}), temp )); 
            [~, idxi] = min(abs(idx-ii)); 
            pairi = idx(idxi);

        end
  else
        %onset locked OR rt locked
        if contains(curVar, 'miss')
            splitVar = split(curVar, 'miss');
            temp = varNames; 
            temp{ii} = 'not this one'; 
            idx = find(cellfun(@(x) contains(x, 'hit') & contains(x, splitVar{2}), temp )); 
            [~, idxi] = min(abs(idx-ii)); 
            pairi = idx(idxi);

        elseif contains(curVar, 'hit')
            splitVar = split(curVar, 'hit');
            temp = varNames; 
            temp{ii} = 'not this one'; 
            idx = find(cellfun(@(x) contains(x, 'miss') & contains(x, splitVar{2}), temp )); 
            [~, idxi] = min(abs(idx-ii)); 
            pairi = idx(idxi);

        elseif contains(curVar, 'cr')
            splitVar = split(curVar, 'cr'); 
            temp = varNames; 
            temp{ii} = 'not this one'; 
            idx = find(cellfun(@(x) contains(x, 'fa') & contains(x, splitVar{2}), temp )); 
            [~, idxi] = min(abs(idx-ii)); 
            pairi = idx(idxi);
        else
            splitVar = split(curVar, 'fa'); 
            temp = varNames; 
            temp{ii} = 'not this one'; 
            idx = find(cellfun(@(x) contains(x, 'cr') & contains(x, splitVar{2}), temp )); 
            [~, idxi] = min(abs(idx-ii)); 
            pairi = idx(idxi);

        end

  end



pairName  = varNames{pairi}; 








end