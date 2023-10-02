

% LL stats plotting

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_LL_KEY_STATS';

statFiles = dir(path); 

statFiles(1:2) = []; 
LLdat = load([statFiles(1).folder '/' statFiles(1).name]).LLdat; 

allSigEnc = zeros(length(LLdat.encTim), 11); 
allSigRet = zeros(length(LLdat.retTim), 11); 
regNames = cell(11,1); 

%hit/miss X region X early/late X encode/retrieve
%early = 0-1000ms; late = 1000-2000ms
allMeanEarlyLate = zeros(2, 11, 2, 2); 


parfor ii = 1:length(statFiles)
ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    hitVals = permute(squeeze(LLdat.regRes(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(LLdat.regRes(2,:,:,:)), [3,1,2]);
    tim = LLdat.encTim; 
    roiLabs = {LLdat.aggTargs(LLdat.reg1).ROI, LLdat.aggTargs(LLdat.reg2).ROI};
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
    caxis([.02, allmax])
%     colorbar

    subplot 252
    imagesc(tim, -150:150, squeeze(mean(missVals)))
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['encode Miss ' num2str(LLdat.n_pair) ' pairs'])
    yline(0)
    xline(0)
    caxis([.02, allmax])
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
    roiLabs = {LLdat.aggTargs(LLdat.reg1).ROI, LLdat.aggTargs(LLdat.reg2).ROI};
   
    subplot 256
    imagesc(tim, -150:150, squeeze(mean(hitVals)))
    allmin = min(min(mean(hitVals), [], 'all'), min(mean(missVals),[], 'all')); 
    allmax = max(max(mean(hitVals), [], 'all'), max(mean(missVals),[], 'all')); 
    allmax = max([abs(allmin), allmax]); 
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['retrieve Hit'])
    yline(0)
    xline(0)
    caxis([.02, allmax])
%     colorbar

    subplot 257
    imagesc(tim, -150:150, squeeze(mean(missVals)))
    ylabel([roiLabs{2} ' leads      ' roiLabs{1} ' leads'])
    title(['retrieve Miss'])
    yline(0)
    xline(0)
    caxis([.02, allmax])
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


%% figure of all encoding latencies 
ROInames = {LLdat.aggTargs.ROI};
figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
for ii = 1:length(statFiles)
    ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    roiLabs = {LLdat.aggTargs(LLdat.reg1).ROI, LLdat.aggTargs(LLdat.reg2).ROI};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))

    curDat = [LLdat.regResLat(1,:), LLdat.regResLat(2,:)];
    hmSort = [ones(length(LLdat.regResLat(1,:)),1); zeros(length(LLdat.regResLat(1,:)),1)];
    modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
                'VariableNames', {'latency', 'hitMiss', 'sub'}); 
    lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
    tval = [0,0];
    tval(1) = lme.Coefficients(2,4); 
    tval(2) = lme.Coefficients(2,6); 
    
    tval(2) = round(tval(2),2);

    if tval(2)<.2
        tval(2) = tval(2)*5; 
    else 
        tval(2) = 1; 
    end
    [~, coli] = min(abs(tval(2) - findMap));
    
  
    
    hold off
    fill([0,2000,2000,0], [0,0,.8,.8], colormap(coli,:), 'FaceAlpha', .8)
    xlim([0,2000])
%     ylim([0, .8])
    hold on
    h = histogram(LLdat.regResLat(1,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(LLdat.regResLat(2,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(LLdat.regResLat(1,:), [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(LLdat.regResLat(2,:), [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


   
%     title(num2str(round(tval(2),2)))
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
      roiLabs = {LLdat.aggTargs(LLdat.reg1).ROI, LLdat.aggTargs(LLdat.reg2).ROI};
%     title(['encoding latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

%     subplot 122
%     hold off
%     histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability')
%     hold on
%     histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability')
%     
%     curDat = [LLdat.regRes2Lat(1,:), LLdat.regRes2Lat(2,:)];
%     hmSort = [ones(length(LLdat.regResLat(1,:)),1); zeros(length(LLdat.regResLat(1,:)),1)];
% 
% 
%     modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
%                 'VariableNames', {'latency', 'hitMiss', 'sub'}); 
%     lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
%     tval = [0,0];
%     tval(1) = lme.Coefficients(2,4); 
%     tval(2) = lme.Coefficients(2,6); 
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
%       
%     title(['retrieval latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

    
  

end

export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\LL_finalizedFigs\' 'LL_latency_allEncode' '.jpg'],''), '-r300')


figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
for ii = 1:length(statFiles)
    ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    roiLabs = {LLdat.aggTargs(LLdat.reg1).ROI, LLdat.aggTargs(LLdat.reg2).ROI};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))

    curDat = [LLdat.regRes2Lat(1,:), LLdat.regRes2Lat(2,:)];
    hmSort = [ones(length(LLdat.regRes2Lat(1,:)),1); zeros(length(LLdat.regRes2Lat(1,:)),1)];
    modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
                'VariableNames', {'latency', 'hitMiss', 'sub'}); 
    lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
    tval = [0,0];
    tval(1) = lme.Coefficients(2,4); 
    tval(2) = lme.Coefficients(2,6); 
    
    tval(2) = round(tval(2),2);

    if tval(2)<.2
        tval(2) = tval(2)*5; 
    else 
        tval(2) = 1; 
    end
    [~, coli] = min(abs(tval(2) - findMap));
    
  
    
    hold off
    fill([0,2000,2000,0], [0,0,.8,.8], colormap(coli,:), 'FaceAlpha', .8)
    xlim([0,2000])
%     ylim([0, .8])
    hold on
    h = histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


   
%     title(num2str(round(tval(2),2)))
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
      roiLabs = {LLdat.aggTargs(LLdat.reg1).ROI, LLdat.aggTargs(LLdat.reg2).ROI};
%     title(['encoding latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

%     subplot 122
%     hold off
%     histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability')
%     hold on
%     histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability')
%     
%     curDat = [LLdat.regRes2Lat(1,:), LLdat.regRes2Lat(2,:)];
%     hmSort = [ones(length(LLdat.regResLat(1,:)),1); zeros(length(LLdat.regResLat(1,:)),1)];
% 
% 
%     modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
%                 'VariableNames', {'latency', 'hitMiss', 'sub'}); 
%     lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
%     tval = [0,0];
%     tval(1) = lme.Coefficients(2,4); 
%     tval(2) = lme.Coefficients(2,6); 
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
%       
%     title(['retrieval latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

    
  

end

export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\LL_finalizedFigs\' 'LL_latency_allRetrieve' '.jpg'],''), '-r300')

%% count pairwise electrodes
pairwiseTrodes = zeros(11); 
pairwiseSubs = zeros(11); 
for ii = 1:length(statFiles)
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    pairwiseTrodes(LLdat.reg1, LLdat.reg2) = LLdat.n_pair; 
    pairwiseSubs(LLdat.reg1, LLdat.reg2) = LLdat.n_sub; 
end





%% looking at all latencies across all regions


%enc/ret X hit/miss X pair
totalLatency = zeros(2, 2, 20098); 
ti = 1; 
for ii = 1:length(statFiles)
ii
 LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
 
 N = LLdat.n_pair; 

 totalLatency(1,1,ti:ti+N-1) = LLdat.regResLat(1,:) ./ LLdat.subhitRT'; 
 totalLatency(1,2,ti:ti+N-1) = LLdat.regResLat(2,:)./ LLdat.submissRT'; 
 totalLatency(2,1,ti:ti+N-1) = LLdat.regRes2Lat(1,:)./ LLdat.rethitRT'; 
 totalLatency(2,2,ti:ti+N-1) = LLdat.regRes2Lat(2,:)./ LLdat.retmissRT';

ti = ti+N;
end



figure
subplot 121
histogram(totalLatency(1,1,:), [0:.03:1], 'normalization', 'probability')
hold on 
histogram(totalLatency(1,2,:), [0:.03:1], 'normalization', 'probability')
title('encoding latency across all ROI pairs')
subplot 122
histogram(totalLatency(1,1,:) - totalLatency(1,2,:), [-1:.03:1], 'normalization', 'probability')
title('within electrode pair hit - miss difference')
xline(0)


figure
subplot 121
histogram(totalLatency(2,1,:), [0:.03:1], 'normalization', 'probability')
hold on 
histogram(totalLatency(2,2,:), [0:.03:1], 'normalization', 'probability')
title('retrieval latency across all ROI pairs')
subplot 122
histogram(totalLatency(2,1,:) - totalLatency(2,2,:), [-1:.03:1], 'normalization', 'probability')
title('within electrode pair hit - miss difference')
xline(0)



scatter(squeeze(totalLatency(2,1,:)), squeeze(totalLatency(2,2,:)),200, 'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.02,'MarkerEdgeAlpha',.02)
hold on 
plot([200,2000], [200,2000])

