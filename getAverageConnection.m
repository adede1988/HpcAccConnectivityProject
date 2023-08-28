function [allConSub, allConRet, allConN] = getAverageConnection(targBrod, allDat)


%area X area X time X offset X up/down
allConSub = zeros([length(targBrod), length(targBrod), 2, 301, length(allDat{1}.leadLag.encTim)]);
allConRet = zeros([length(targBrod), length(targBrod), 2, 301, length(allDat{1}.leadLag.retTim)]);
allConN = zeros(length(targBrod));

for sub = 1:length(allDat)
    sub
    if ~isempty(allDat{sub})
        c = allDat{sub}.leadLag; 
        b = allDat{sub}.brodmann; 
        for chan1 = 1:size(c.subMem,1)
            chan1i = find(cellfun(@(y) sum(y), cellfun(@(x) strcmp(b{chan1}, x), {targBrod.lab}, 'uniformoutput', false)));
            if sum(chan1i)>0 %this is one of the target channels! 
            for chan2 = 1:size(c.subMem,1)
                if chan2 ~= chan1
                    chan2i = find(cellfun(@(y) sum(y), cellfun(@(x) strcmp(b{chan2}, x), {targBrod.lab}, 'uniformoutput', false)));
              
                if sum(chan2i)>0 %this is one of the target channels! 
                    allConN(chan1i, chan2i) = allConN(chan1i, chan2i) + 1; 
 
                    %increment subsequent memory
                    allConSub(chan1i, chan2i, :, :, :) = allConSub(chan1i, chan2i, :,:,:) + c.subMem(chan1, chan2, :, :, :);  
                    %increment retrieval data
                    allConRet(chan1i, chan2i, :, :, :) = allConRet(chan1i, chan2i, :,:,:) + c.retMem(chan1, chan2, :, :, :);

                       
                   
                end
                end
            end
            end
        end
    




    end


end









end