function [check] = checkForThreshold(test, timeVals, timeLims)

    check  = false; 
    test = test>1.96;
    testidx = find(test([find(timeVals>=timeLims(1),1):find(timeVals>=timeLims(2),1)]));
    if ~isempty(testidx)
        breakPoints = [1 find(diff(testidx)>1)']; 
        if length(breakPoints)==1 && length(testidx)>10
            check = true; 
        else
            breakPoints = [breakPoints, length(testidx)]; 
            for bb = 1:length(breakPoints)-1
                if breakPoints(bb+1) - breakPoints(bb) > 10
                    check = true; 
                end
            end

        end
    end


end