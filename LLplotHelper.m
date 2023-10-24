function [] = LLplotHelper(LLdat, ii)


  hitVals = permute(squeeze(LLdat.regRes(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(LLdat.regRes(2,:,:,:)), [3,1,2]);
    tim = LLdat.encTim; 
    roiLabs = {LLdat.aggTargs{LLdat.reg1}, LLdat.aggTargs{LLdat.reg2}};
    f = figure('visible', false);
    f.Position = [0 0 1500 700];
    subplot 251
    imagesc(tim, -150:150, squeeze(mean(hitVals)))
    allmin = min(min(mean(hitVals), [], 'all'), min(mean(missVals),[], 'all')); 
    allmax = max(max(mean(hitVals), [], 'all'), max(mean(missVals),[], 'all')); 
    allmax = max([abs(allmin), allmax]); 
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['encode Hit '  num2str(LLdat.n_sub) ' subjects'])
    yline(0)
    xline(0)
    caxis([.0, allmax])
%     colorbar

    subplot 252
    imagesc(tim, -150:150, squeeze(mean(missVals)))
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['encode Miss ' num2str(LLdat.n_pair) ' pairs'])
    yline(0)
    xline(0)
    caxis([.0, allmax])
%     caxis([-allmax, allmax])
%     colorbar


    subplot 253
    hold off
    imagesc(LLdat.tVals_sub)
    allmax = max([abs(min(LLdat.tVals_sub, [], 'all')), max(LLdat.tVals_sub, [], 'all') ]);
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    caxis([-allmax*.7, allmax*.7])
%     colorbar
    if min(LLdat.p_sub, [], 'all')<.05
        
        addRedOutline(LLdat.p_sub', .05, 'red');
        title('t-value; p<.05')

    else
        minp = min(LLdat.p_sub, [], 'all');
        addRedOutline(LLdat.p_sub', minp+.01, 'green');
        title(['t-value; p=' num2str(round(minp, 2))])
    end



    xticks([21,61,101,141])
    xticklabels(tim([21,61,101, 141]))
    yticks([1:50:301])
    yticklabels([-150:50:150])
    yline(151)
    xline(21)
    hold on 


    subplot 254
    
    hold off
    plot(tim, squeeze(mean(LLdat.regResOff(1,:,:),3)),'linewidth', 3, 'color', 'blue')
    hold on 
    plot(tim, squeeze(mean(LLdat.regResOff(2,:,:),3)),'linewidth', 3, 'color', 'red')
    set(gca, 'YDir', 'reverse');

    sdHit = std(squeeze(LLdat.regResOff(1,:,:)), [], 2) ./ sqrt(LLdat.n_pair);
    sdMiss = std(squeeze(LLdat.regResOff(2,:,:)), [], 2) ./ sqrt(LLdat.n_pair); 
    
    x = [tim, flip(tim)];
    y = [squeeze(mean(LLdat.regResOff(1,:,:),3))' - sdHit; flip(squeeze(mean(LLdat.regResOff(1,:,:),3))' + sdHit)]; 
    fill(flip(x), flip(y'), 'blue', 'FaceAlpha', .2)
    x = [tim, flip(tim)];
    y = [squeeze(mean(LLdat.regResOff(2,:,:),3))' - sdHit; flip(squeeze(mean(LLdat.regResOff(2,:,:),3))' + sdHit)];
    fill(flip(x), flip(y'), 'red', 'FaceAlpha', .2)
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    ttestRes = zeros(length(tim),1); 
    for ti = 1:length(tim)
        ttestRes(ti) = ttest(LLdat.regResOff(1,ti,:) - LLdat.regResOff(2,ti,:));
    end
    scatter(tim(ttestRes==1), ttestRes(ttestRes==1)*0, 30, 'k', 'filled')
    yline(0)
    xline(0)
    xlim([-500, 3000])
    title('mean offset')

    subplot(2, 5, 5)
    clusLoc = LLdat.p_sub == min(LLdat.p_sub, [], 'all');
    makeLLBox(LLdat, hitVals, missVals, tim, clusLoc)
  


%% retrieval! 

    hitVals = permute(squeeze(LLdat.regRes2(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(LLdat.regRes2(2,:,:,:)), [3,1,2]);
    tim = LLdat.retTim; 
    roiLabs = {LLdat.aggTargs{LLdat.reg1}, LLdat.aggTargs{LLdat.reg2}};
   
    subplot 256
    imagesc(tim, -150:150, squeeze(mean(hitVals)))
    allmin = min(min(mean(hitVals), [], 'all'), min(mean(missVals),[], 'all')); 
    allmax = max(max(mean(hitVals), [], 'all'), max(mean(missVals),[], 'all')); 
    allmax = max([abs(allmin), allmax]); 
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['retrieve Hit'])
    yline(0)
    xline(0)
    caxis([.0, allmax])
%     colorbar

    subplot 257
    imagesc(tim, -150:150, squeeze(mean(missVals)))
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['retrieve Miss'])
    yline(0)
    xline(0)
    caxis([.0, allmax])
%     colorbar


    subplot 258
    hold off
    imagesc(LLdat.tVals_ret)
    allmax = max([abs(min(LLdat.tVals_ret, [], 'all')), max(LLdat.tVals_ret, [], 'all') ]);
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    caxis([-allmax*.7, allmax*.7])
    title('t-value')
%     colorbar

    if min(LLdat.p_ret, [], 'all')<.05
        
        addRedOutline(LLdat.p_ret', .05, 'red');
        title('t-value; p<.05')

    else
        minp = min(LLdat.p_ret, [], 'all');
        addRedOutline(LLdat.p_ret', minp+.01, 'green');
        title(['t-value; p=' num2str(round(minp, 2))])
    end

    xticks([21,61,101])
    xticklabels(tim([21,61,101]))
    yticks([1:50:301])
    yticklabels([-150:50:150])
    yline(151)
    xline(21)


    subplot 259
    
    hold off
    plot(tim, squeeze(mean(LLdat.regRes2Off(1,:,:),3)),'linewidth', 3, 'color', 'blue')
    hold on 
    plot(tim, squeeze(mean(LLdat.regRes2Off(2,:,:),3)),'linewidth', 3, 'color', 'red')
    set(gca, 'YDir', 'reverse');

    sdHit = std(squeeze(LLdat.regRes2Off(1,:,:)), [], 2) ./ sqrt(LLdat.n_pair);
    sdMiss = std(squeeze(LLdat.regRes2Off(2,:,:)), [], 2) ./ sqrt(LLdat.n_pair); 
    
    x = [tim, flip(tim)];
    y = [squeeze(mean(LLdat.regRes2Off(1,:,:),3))' - sdHit; flip(squeeze(mean(LLdat.regRes2Off(1,:,:),3))' + sdHit)]; 
    fill(flip(x), flip(y'), 'blue', 'FaceAlpha', .2)
    x = [tim, flip(tim)];
    y = [squeeze(mean(LLdat.regRes2Off(2,:,:),3))' - sdHit; flip(squeeze(mean(LLdat.regRes2Off(2,:,:),3))' + sdHit)];
    fill(flip(x), flip(y'), 'red', 'FaceAlpha', .2)
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])

    ttestRes = zeros(length(tim),1); 
    for ti = 1:length(tim)
        ttestRes(ti) = ttest(LLdat.regRes2Off(1,ti,:) - LLdat.regRes2Off(2,ti,:));
    end
    scatter(tim(ttestRes==1), ttestRes(ttestRes==1)*0, 30, 'k', 'filled')
    title('mean offset')
    xlim([-500, 2500])

    subplot(2, 5, 10)
    clusLoc = LLdat.p_ret == min(LLdat.p_ret, [], 'all'); 
    makeLLBox(LLdat, hitVals, missVals, tim, clusLoc)
%     Xidx = repmat(tim, [size(hitVals,2),1]);
%     Yidx = repmat(-150:150, [size(hitVals,3), 1])';
%     clusLoc = LLdat.p_ret == min(LLdat.p_ret, [], 'all'); 
%     chanMeans = zeros(LLdat.n_pair*2,1); %hits, misses
%     hmSort = chanMeans; 
%     ti = 1; 
%     for pi = 1:LLdat.n_pair
%         tmp = squeeze(hitVals(pi, :, :)); 
%         chanMeans(ti) = mean(tmp(clusLoc), 'all'); 
%         hmSort(ti) = 1; 
%         ti = ti+1; 
%         tmp = squeeze(missVals(pi, :, :)); 
%         chanMeans(ti) = mean(tmp(clusLoc), 'all'); 
%         hmSort(ti) = 0; 
%         ti = ti+1; 
%        
%     end
% 
%     hold off
% 
%     b = boxchart(hmSort, chanMeans);
%     b.MarkerStyle = 'none'; 
%     xticks([0,1])
%     xticklabels({'Miss', 'Hit'})
%     meanX = round(mean(Xidx(clusLoc), 'all')); 
%     meanY = round(mean(Yidx(clusLoc), 'all')); 
% %             LLvals = -150:150; 
%     title(['time: ' num2str(round(meanX)) ', offset: ' num2str(round(meanY)) ])
%     PLH = chanMeans(hmSort==1); 
%     PLM = chanMeans(hmSort==0); 
%     randVals = (rand(length(PLH),1)-.5)*.5;
%     hold on 
%     scatter(randVals, PLM, 10,  'blue')
%     scatter(randVals+1, PLH, 10, 'blue')
%     
%     for pi = 1:length(PLH)
%         plot([0+randVals(pi),1+randVals(pi)], [PLM(pi),PLH(pi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
%         
% 
%     end
% 
%     %make subject means
%     subHits = zeros(LLdat.n_sub,1); 
%     subMisses = zeros(LLdat.n_sub,1); 
%     uniqueSubs = unique(LLdat.regSubIDs); 
%     regSubs = LLdat.regSubIDs; 
%     for sub = 1:length(subHits)
%         subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
%         subHits(sub) = mean(PLH(subidx));
%         subMisses(sub) = mean(PLM(subidx));
%         plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')
% 
%     end
%     scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
%     scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
%     
%     h = ttest(PLH - PLM);
%     if h
%         text(-.8, max([PLH;PLM])*.5, "E: *")
%     end
%     ylim([min([PLH;PLM]), max([PLH; PLM])*1.2])
%     xlim([-1,1.5])
%     h = ttest(subHits - subMisses);
%     if h
%         text(-.8, max([PLH;PLM])*.2, "P: *")
%     end
%     ylabel("LL correlation")



    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\LL_finalizedFigs\' roiLabs{1} '_' roiLabs{2} '.jpg'], '-r300')

    end