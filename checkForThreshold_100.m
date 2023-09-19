function [check] = checkForThreshold_100(test, timeVals, timeLims)
    tmp = test; 
    check  = 0; 
    test = test>1.96;
    testidx = find(test([find(timeVals>=timeLims(1),1):find(timeVals>=timeLims(2),1)]));
    if ~isempty(testidx)
        breakPoints = [1 find(diff(testidx)>1)']; 
        if length(breakPoints)==1 && length(testidx)>4
            check = 1; 
        else
            breakPoints = [breakPoints, length(testidx)]; 
            for bb = 1:length(breakPoints)-1
                if breakPoints(bb+1) - breakPoints(bb) > 4
                    check = 1; 
                end
            end

        end
    end

    if check == 0 
        test = tmp<-1.96;
        testidx = find(test([find(timeVals>=timeLims(1),1):find(timeVals>=timeLims(2),1)]));
        if ~isempty(testidx)
            breakPoints = [1 find(diff(testidx)>1)']; 
            if length(breakPoints)==1 && length(testidx)>4
                check = -1; 
            else
                breakPoints = [breakPoints, length(testidx)]; 
                for bb = 1:length(breakPoints)-1
                    if breakPoints(bb+1) - breakPoints(bb) > 4
                        check = -1; 
                    end
                end
    
            end
        end

    end




end