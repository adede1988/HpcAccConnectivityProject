function [HFB] = HFBdown(HFB)


% downsample for size
    HFB_names = fieldnames(HFB); 
    
    datidx = {[1,2], [5,6], [9,10,11,12], [15,16,17,18]}; 
    timidx = [3,7,13,19]; 

    for dati = 1:4
        curtim = HFB.(HFB_names{timidx(dati)});
        HFB.(HFB_names{timidx(dati)}) = curtim(1:5:end);
        curDi = datidx{dati}; 
        for fi = 1:length(curDi)
            cur = HFB.(HFB_names{curDi(fi)});

            HFB.(HFB_names{curDi(fi)}) = cur(1:5:end,:);
        end
    end


end