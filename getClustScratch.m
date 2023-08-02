function [errorChans, hipClustUp, hipClustDown, notClustUp, notClustDown] = getClustScratch(chanFiles, sub)



try
    chanDat = load([chanFiles(sub).folder '/' chanFiles(sub).name]).chanDat; 
    
    %get the overlapped clusters
    subMemUp = zeros(length(chanDat.leadLag.encTim), 301); 
    subMemDown = subMemUp; 
    tim = chanDat.leadLag.encTim; 
    subMem = chanDat.leadLag.subMem; 
    LLtim = -150:150; 
    for chan2 = 1:size(subMem,1)
    for chan = 1:size(subMem, 1)
        for cc = 1:30
            if ~isnan(subMem(chan2, chan,1,cc,1)) %there's an up cluster!
                subMemUp(tim>=subMem(chan2, chan,1,cc,3) & tim<= subMem(chan2,chan,1,cc,4), LLtim>=subMem(chan2,chan,1,cc,6) & LLtim<=subMem(chan2,chan,1,cc,7)) = ...
                subMemUp(tim>=subMem(chan2,chan,1,cc,3) & tim<= subMem(chan2,chan,1,cc,4), LLtim>=subMem(chan2,chan,1,cc,6) & LLtim<=subMem(chan2,chan,1,cc,7)) + 1; 
            end
             if ~isnan(subMem(chan2,chan,2,cc,1)) %there's a down cluster!
                
                subMemDown(tim>=subMem(chan2,chan,2,cc,3) & tim<= subMem(chan2,chan,2,cc,4), LLtim>=subMem(chan2,chan,2,cc,6) & LLtim<=subMem(chan2,chan,2,cc,7)) = ...
                subMemDown(tim>=subMem(chan2,chan,2,cc,3) & tim<= subMem(chan2,chan,2,cc,4), LLtim>=subMem(chan2,chan,2,cc,6) & LLtim<=subMem(chan2,chan,2,cc,7)) + 1; 
            end 

        

        end
    end
    end
    if chanDat.hip == 1 
        hipClustUp = subMemUp; 
        hipClustDown = subMemDown; 
        notClustUp = []; 
        notClustDown = []; 
    else
        hipClustUp = []; 
        hipClustDown = []; 
        notClustUp = subMemUp; 
        notClustDown = subMemDown;
    end


    errorChans = 0; 

catch

    hipClustUp = []; 
        hipClustDown = []; 
        notClustUp = []; 
        notClustDown = [];

    errorChans = 1; 

end

















end