function [] = makeLLBox(LLdat, hitVals, missVals, tim, clusLoc)

    Xidx = repmat(tim, [size(hitVals,2),1]);
    Yidx = repmat(-150:150, [size(hitVals,3), 1])';
    
    chanMeans = zeros(LLdat.n_pair*2,1); %hits, misses
    hmSort = chanMeans; 
    ti = 1; 
    for pi = 1:LLdat.n_pair
        tmp = squeeze(hitVals(pi, :, :)); 
        chanMeans(ti) = mean(tmp(clusLoc), 'all'); 
        hmSort(ti) = 1; 
        ti = ti+1; 
        tmp = squeeze(missVals(pi, :, :)); 
        chanMeans(ti) = mean(tmp(clusLoc), 'all'); 
        hmSort(ti) = 0; 
        ti = ti+1; 
       
    end

    hold off

    b = boxchart(hmSort, chanMeans);
    b.MarkerStyle = 'none'; 
    xticks([0,1])
    xticklabels({'Miss', 'Hit'})
    meanX = round(mean(Xidx(clusLoc), 'all')); 
    meanY = round(mean(Yidx(clusLoc), 'all')); 
%             LLvals = -150:150; 
    title(['time: ' num2str(round(meanX)) ', offset: ' num2str(round(meanY)) ])
    PLH = chanMeans(hmSort==1); 
    PLM = chanMeans(hmSort==0); 
    randVals = (rand(length(PLH),1)-.5)*.5;
    hold on 
    scatter(randVals, PLM, 10,  'blue')
    scatter(randVals+1, PLH, 10, 'blue')
    
    for pi = 1:length(PLH)
        plot([0+randVals(pi),1+randVals(pi)], [PLM(pi),PLH(pi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
        

    end

    %make subject means
    subHits = zeros(LLdat.n_sub,1); 
    subMisses = zeros(LLdat.n_sub,1); 
    uniqueSubs = unique(LLdat.regSubIDs); 
    regSubs = LLdat.regSubIDs; 
    for sub = 1:length(subHits)
        subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
        subHits(sub) = mean(PLH(subidx));
        subMisses(sub) = mean(PLM(subidx));
        plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

    end
    scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
    scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
    
    h = ttest(PLH - PLM);
    if ~isnan(h)
    if h
        text(-.8, max([PLH;PLM])*.5, "E: *")
    end
    end
    ylim([min([PLH;PLM]), max([PLH; PLM])*1.2])
    xlim([-1,1.5])
    
    h = ttest(subHits - subMisses);
    if ~isnan(h)
    if h
        text(-.8, max([PLH;PLM])*.2, "P: *")
    end
    end
    ylabel("LL correlation")




end