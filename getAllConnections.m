function [allCon, allConN] = getAllConnections(targBrod, allDat)

targBrod(cellfun(@(x) strcmp('ERROR', x), {targBrod.lab})) = []; 
%area X area X time X offset X up/down
allCon = zeros([length(targBrod), length(targBrod), length(allDat{1}.leadLag.encTim), 301, 2]);
allConN = zeros(length(targBrod));
tim = allDat{1}.leadLag.encTim; 
LLtim = -150:150; 
for sub = 1:length(allDat)
    if ~isempty(allDat{sub})
        c = allDat{sub}.leadLag; 
        b = allDat{sub}.brodmann; 
        for chan1 = 1:size(c.subMem,1)
            chan1i = cellfun(@(x) strcmp(b{chan1}, x), {targBrod.lab});
            if sum(chan1i)>0 %this is one of the target channels! 
            for chan2 = 1:size(c.subMem,1)
                if chan2 ~= chan1
                chan2i = cellfun(@(x) strcmp(b{chan2}, x), {targBrod.lab});
                if sum(chan2i)>0 %this is one of the target channels! 
                    allConN(chan1i, chan2i) = allConN(chan1i, chan2i) + 1; 
                    for cc = 1:size(c.subMem, 4)
                        if ~isnan(c.subMem(chan1, chan2, 1, cc, 1)) %&& length(find(tim>=c.subMem(chan1, chan2,1,cc,3) & tim<= c.subMem(chan1,chan2,1,cc,4)))>1
                            if c.subMem(chan1, chan2, 1, cc, 1) > 0
                          

                            allCon(chan1i, chan2i, tim>=c.subMem(chan1, chan2,1,cc,3) & tim<= c.subMem(chan1,chan2,1,cc,4),...
                                                   LLtim>=c.subMem(chan1, chan2,1,cc,6) & LLtim<=c.subMem(chan1, chan2,1,cc,7),...
                                                   1) = allCon(chan1i, chan2i, ...
                                                   tim>=c.subMem(chan1, chan2,1,cc,3) & tim<= c.subMem(chan1,chan2,1,cc,4),...
                                                   LLtim>=c.subMem(chan1, chan2,1,cc,6) & LLtim<=c.subMem(chan1, chan2,1,cc,7),...
                                                   1) + 1;  
                            end
    
                        end

                        if ~isnan(c.subMem(chan1, chan2, 2, cc, 2)) %&& length(find(tim>=c.subMem(chan1, chan2,2,cc,3) & tim<= c.subMem(chan1,chan2,2,cc,4)))>1
                           if c.subMem(chan1, chan2, 2, cc, 1) > 0
                            allCon(chan1i, chan2i, tim>=c.subMem(chan1, chan2,2,cc,3) & tim<= c.subMem(chan1,chan2,2,cc,4),...
                                                   LLtim>=c.subMem(chan1, chan2,2,cc,6) & LLtim<=c.subMem(chan1, chan2,2,cc,7),...
                                                   2) = allCon(chan1i, chan2i,...
                                                   tim>=c.subMem(chan1, chan2,2,cc,3) & tim<= c.subMem(chan1,chan2,2,cc,4),...
                                                   LLtim>=c.subMem(chan1, chan2,2,cc,6) & LLtim<=c.subMem(chan1, chan2,2,cc,7),...
                                                   2) + 1; 
                           end
    
                        end
                    end
                end
                end
            end
            end
        end
    




    end


end









end