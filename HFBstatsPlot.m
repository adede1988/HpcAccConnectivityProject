

%HFBstatsPlot script to take in HFB summary stats and construct plots

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS';

statFiles = dir(path); 

statFiles(1:2) = []; 
HFBdat = load([statFiles(1).folder '/' statFiles(1).name]).HFBdat; 

allSigEnc = zeros(length(HFBdat.encTim), 11); 
allSigRet = zeros(length(HFBdat.retTim), 11); 
regNames = cell(11,1); 

%hit/miss X region X early/late X encode/retrieve
%early = 0-1000ms; late = 1000-2000ms
allMeanEarlyLate = zeros(2, 11, 2, 2); 

for ii = 1:length(statFiles)

    HFBdat = load([statFiles(ii).folder '/' statFiles(ii).name]).HFBdat; 

    allSigEnc(HFBdat.p_sub<.05 & HFBdat.tVals_sub>0, ii) = 1; 
    allSigEnc(HFBdat.p_sub<.05 & HFBdat.tVals_sub<0, ii) = -1;
    allSigRet(HFBdat.p_ret<.05 & HFBdat.tVals_ret<0, ii) = -1; 
    allSigRet(HFBdat.p_ret<.05 & HFBdat.tVals_ret>0, ii) = 1; 
    regNames{ii} = HFBdat.aggTargs(HFBdat.reg1).ROI;
    
    %get all mean early late
    wins = [0,1000,2000]; 

    for hm = 1:2 %hit/miss loop
        for w = 1:2 %window
            allMeanEarlyLate(hm,ii,w, 1) = median(HFBdat.regRes(hm,HFBdat.encTim>=wins(w) & HFBdat.encTim<=wins(w+1),:), 'all');
            allMeanEarlyLate(hm,ii,w, 2) = median(HFBdat.regRes2(hm,HFBdat.retTim>=wins(w) & HFBdat.retTim<=wins(w+1),:), 'all');
        end
    end

    figure('visible', true, 'Position', [0,0,1000,500])
    
    subplot(3, 4, [1:2])
    hold off
    plot(HFBdat.encTim, HFBdat.hitVals_sub, 'linewidth', 3, 'color', 'blue')
    hold on 
    plot(HFBdat.encTim, HFBdat.missVals_sub, 'linewidth', 3, 'color', 'red')
    xlim([HFBdat.encTim(1), HFBdat.encTim(end)])

    sdHit = std(squeeze(HFBdat.regRes(1,:,:)), [], 2) ./ sqrt(HFBdat.n_pair);
    sdMiss = std(squeeze(HFBdat.regRes(2,:,:)), [], 2) ./ sqrt(HFBdat.n_pair); 
    
    x = [HFBdat.encTim, flip(HFBdat.encTim)];
    y = [HFBdat.hitVals_sub' - sdHit; flip(HFBdat.hitVals_sub' + sdHit)]; 
    fill(flip(x), flip(y'), 'blue', 'FaceAlpha', .2)
    x = [HFBdat.encTim, flip(HFBdat.encTim)];
    y = [HFBdat.missVals_sub' - sdMiss; flip(HFBdat.missVals_sub' + sdMiss)]; 
    fill(flip(x), flip(y'), 'red', 'FaceAlpha', .2)

%     maxidx = [1:length(HFBdat.encTim)]; 
%     maxidx(HFBdat.p_sub>.05) = []; 
%     if isempty(maxidx)
%         maxidx = [1:length(HFBdat.encTim)];
%     end
%     [~, maxLoc] = max(abs(HFBdat.tVals_sub(maxidx)));
%     maxLoc = maxidx(maxLoc); 
% %     maxLoc = 41;
%     xline(HFBdat.encTim(maxLoc), 'linewidth', 2, 'linestyle', '--', 'color', 'red')

    plot(HFBdat.encTim, HFBdat.p_sub<.05, 'color', 'k', 'linewidth', 2)
    ylabel("HFB amplitude (z-score)")
    title(['HFB ' HFBdat.aggTargs(HFBdat.reg1).ROI ' encoding'])


    subplot(3, 4, [3:4])
    hold off
    plot(HFBdat.retTim, HFBdat.hitVals_ret, 'linewidth', 3, 'color', 'blue')
    hold on 
    plot(HFBdat.retTim, HFBdat.missVals_ret, 'linewidth', 3, 'color', 'red')
    xlim([HFBdat.retTim(1), HFBdat.retTim(end)])

    sdHit = std(squeeze(HFBdat.regRes2(1,:,:)), [], 2) ./ sqrt(HFBdat.n_pair);
    sdMiss = std(squeeze(HFBdat.regRes2(2,:,:)), [], 2) ./ sqrt(HFBdat.n_pair); 
    
    x = [HFBdat.retTim, flip(HFBdat.retTim)];
    y = [HFBdat.hitVals_ret' - sdHit; flip(HFBdat.hitVals_ret' + sdHit)]; 
    fill(flip(x), flip(y'), 'blue', 'FaceAlpha', .2)
    x = [HFBdat.retTim, flip(HFBdat.retTim)];
    y = [HFBdat.missVals_ret' - sdMiss; flip(HFBdat.missVals_ret' + sdMiss)]; 
    fill(flip(x), flip(y'), 'red', 'FaceAlpha', .2)
%     maxidx = [1:length(HFBdat.retTim)]; 
%     maxidx(HFBdat.p_ret>.05) = []; 
%     if isempty(maxidx)
%         maxidx = [1:length(HFBdat.retTim)];
%     end
%     [~, maxLoc2] = max(abs(HFBdat.tVals_ret(maxidx)));
%     maxLoc2 = maxidx(maxLoc2); 
% %     maxLoc2 = 51; 
%     xline(HFBdat.retTim(maxLoc2), 'linewidth', 2, 'linestyle', '--', 'color', 'red')
    ylabel("HFB amplitude (z-score)")

    plot(HFBdat.retTim, HFBdat.p_ret<.05, 'color', 'k', 'linewidth', 2)

    title(['HFB ' HFBdat.aggTargs(HFBdat.reg1).ROI ' retrieval'])


    subplot(3, 4, [5:6])

    combo = [squeeze(HFBdat.regRes(1,:,:)), squeeze(HFBdat.regRes(2,:,:))];
    imagesc(HFBdat.encTim, [], combo')
    yline(HFBdat.n_pair+.5, 'linewidth', 2, 'color', 'red')
%     xline(HFBdat.encTim(maxLoc), 'linewidth', 2, 'linestyle', '--', 'color', 'red')
    allMax = max(combo, [], 'all');
    caxis([0, 5])
    yticklabels('')
    ylabel('Misses           Hits')
    title('electrode means')

    subplot(3, 4, [7:8])

    combo = [squeeze(HFBdat.regRes2(1,:,:)), squeeze(HFBdat.regRes2(2,:,:))];
    imagesc(HFBdat.retTim, [], combo')
    yline(HFBdat.n_pair+.5, 'linewidth', 2, 'color', 'red')
%     xline(HFBdat.retTim(maxLoc2), 'linewidth', 2, 'linestyle', '--', 'color', 'red')
    allMax = max(combo, [], 'all');
    caxis([0, 5])
    yticklabels('')
    ylabel('Misses           Hits')
    title('electrode means')



    %plot the point of maximum difference in more detail
    subplot(3, 4, 9)
    chanMeans = zeros(HFBdat.n_pair*2,1); %hits, misses
    hmSort = chanMeans; 
    ti = 1; 
    for pi = 1:HFBdat.n_pair
        chanMeans(ti) = mean(HFBdat.regRes(1,HFBdat.encTim>=0 & HFBdat.encTim<=1000,pi), 'all'); 
        hmSort(ti) = 1; 
        ti = ti+1; 
        chanMeans(ti) = mean(HFBdat.regRes(2,HFBdat.encTim>=0 & HFBdat.encTim<=1000,pi), 'all'); 
        hmSort(ti) = 0; 
        ti = ti+1; 
       
    end

    hold off

    b = boxchart(hmSort, chanMeans);
    b.MarkerStyle = 'none'; 
    xticks([0,1])
    xticklabels({'Miss', 'Hit'})
%             meanX = round(mean(Xidx(tmp), 'all')); 
%             meanY = round(mean(Yidx(tmp), 'all')); 
%             LLvals = -150:150; 
    title(['time:  0-1000 ms' ])
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
    subHits = zeros(HFBdat.n_sub,1); 
    subMisses = zeros(HFBdat.n_sub,1); 
    uniqueSubs = unique(HFBdat.regSubIDs); 
    regSubs = HFBdat.regSubIDs; 
    for sub = 1:length(subHits)
        subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
        subHits(sub) = mean(PLH(subidx));
        subMisses(sub) = mean(PLM(subidx));
        plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

    end
    scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
    scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
    
    h = ttest(PLH - PLM);
    if h
        text(-.8, max([PLH;PLM])*.5, "E: *")
    end
    ylim([min([PLH;PLM]), max([PLH; PLM])*1.2])
    xlim([-1,1.5])
    h = ttest(subHits - subMisses);
    if h
        text(-.8, max([PLH;PLM])*.2, "P: *")
    end
    ylabel("HFB amplitude (z-score)")



     subplot(3, 4, 10)
    chanMeans = zeros(HFBdat.n_pair*2,1); %hits, misses
    hmSort = chanMeans; 
    ti = 1; 
    for pi = 1:HFBdat.n_pair
        chanMeans(ti) = mean(HFBdat.regRes(1,HFBdat.encTim>=1000 & HFBdat.encTim<=2000,pi), 'all'); 
        hmSort(ti) = 1; 
        ti = ti+1; 
        chanMeans(ti) = mean(HFBdat.regRes(2,HFBdat.encTim>=1000 & HFBdat.encTim<=2000,pi), 'all'); 
        hmSort(ti) = 0; 
        ti = ti+1; 
       
    end

    hold off

    b = boxchart(hmSort, chanMeans);
    b.MarkerStyle = 'none'; 
    xticks([0,1])
    xticklabels({'Miss', 'Hit'})
%             meanX = round(mean(Xidx(tmp), 'all')); 
%             meanY = round(mean(Yidx(tmp), 'all')); 
%             LLvals = -150:150; 
    title(['time:  1000-2000 ms' ])
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
    subHits = zeros(HFBdat.n_sub,1); 
    subMisses = zeros(HFBdat.n_sub,1); 
    uniqueSubs = unique(HFBdat.regSubIDs); 
    regSubs = HFBdat.regSubIDs; 
    for sub = 1:length(subHits)
        subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
        subHits(sub) = mean(PLH(subidx));
        subMisses(sub) = mean(PLM(subidx));
        plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

    end
    scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
    scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
    
    h = ttest(PLH - PLM);
    if h
        text(-.8, max([PLH;PLM])*.5, "E: *")
    end
    ylim([min([PLH;PLM]), max([PLH; PLM])*1.2])
    xlim([-1,1.5])
    h = ttest(subHits - subMisses);
    if h
        text(-.8, max([PLH;PLM])*.2, "P: *")
    end
    ylabel("HFB amplitude (z-score)")


      subplot(3, 4, 11)
     
    chanMeans = zeros(HFBdat.n_pair*2,1); %hits, misses
    hmSort = chanMeans; 
    ti = 1; 
    for pi = 1:HFBdat.n_pair
        chanMeans(ti) = mean(HFBdat.regRes2(1,HFBdat.retTim>=0 & HFBdat.retTim<=1000,pi), 'all'); 
        hmSort(ti) = 1; 
        ti = ti+1; 
        chanMeans(ti) = mean(HFBdat.regRes2(2,HFBdat.retTim>=0 & HFBdat.retTim<=1000,pi), 'all'); 
        hmSort(ti) = 0; 
        ti = ti+1; 
       
    end

    hold off

    b = boxchart(hmSort, chanMeans);
    b.MarkerStyle = 'none'; 
    xticks([0,1])
    xticklabels({'Miss', 'Hit'})
%             meanX = round(mean(Xidx(tmp), 'all')); 
%             meanY = round(mean(Yidx(tmp), 'all')); 
%             LLvals = -150:150; 
    title(['time:  0-1000 ms' ])
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
    subHits = zeros(HFBdat.n_sub,1); 
    subMisses = zeros(HFBdat.n_sub,1); 
    uniqueSubs = unique(HFBdat.regSubIDs); 
    regSubs = HFBdat.regSubIDs; 
    for sub = 1:length(subHits)
        subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
        subHits(sub) = mean(PLH(subidx));
        subMisses(sub) = mean(PLM(subidx));
        plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

    end
    scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
    scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
    
    h = ttest(PLH - PLM);
    if h
        text(-.8, max([PLH;PLM])*.5, "E: *")
    end
    ylim([min([PLH;PLM]), max([PLH; PLM])*1.2])
    xlim([-1,1.5])
    h = ttest(subHits - subMisses);
    if h
        text(-.8, max([PLH;PLM])*.2, "P: *")
    end
    ylabel("HFB amplitude (z-score)")



     subplot(3, 4, 12)
    chanMeans = zeros(HFBdat.n_pair*2,1); %hits, misses
    hmSort = chanMeans; 
    ti = 1; 
    for pi = 1:HFBdat.n_pair
        chanMeans(ti) = mean(HFBdat.regRes2(1,HFBdat.retTim>=1000 & HFBdat.retTim<=2000,pi), 'all'); 
        hmSort(ti) = 1; 
        ti = ti+1; 
        chanMeans(ti) = mean(HFBdat.regRes2(2,HFBdat.retTim>=1000 & HFBdat.retTim<=2000,pi), 'all'); 
        hmSort(ti) = 0; 
        ti = ti+1; 
       
    end

    hold off

    b = boxchart(hmSort, chanMeans);
    b.MarkerStyle = 'none'; 
    xticks([0,1])
    xticklabels({'Miss', 'Hit'})
%             meanX = round(mean(Xidx(tmp), 'all')); 
%             meanY = round(mean(Yidx(tmp), 'all')); 
%             LLvals = -150:150; 
    title(['time:  1000-2000 ms' ])
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
    subHits = zeros(HFBdat.n_sub,1); 
    subMisses = zeros(HFBdat.n_sub,1); 
    uniqueSubs = unique(HFBdat.regSubIDs); 
    regSubs = HFBdat.regSubIDs; 
    for sub = 1:length(subHits)
        subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
        subHits(sub) = mean(PLH(subidx));
        subMisses(sub) = mean(PLM(subidx));
        plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')

    end
    scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
    scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
    
    h = ttest(PLH - PLM);
    if h
        text(-.8, max([PLH;PLM])*.5, "E: *")
    end
    ylim([min([PLH;PLM]), max([PLH; PLM])*1.2])
    xlim([-1,1.5])
    h = ttest(subHits - subMisses);
    if h
        text(-.8, max([PLH;PLM])*.2, "P: *")
    end
    ylabel("HFB amplitude (z-score)")



    export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\' 'HFB_' HFBdat.aggTargs(HFBdat.reg1).ROI '.jpg'],''), '-r300')



end


figure('visible', true, 'Position', [0,0,500,500])
imagesc(HFBdat.encTim, [], allSigEnc')
yticks(1:11)
yticklabels(regNames)
xline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'red')
title("encoding significant hit/miss effects")
 export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\' 'HFB_allReg_encoding.jpg'],''), '-r300')

figure('visible', true, 'Position', [0,0,500,500])
imagesc(HFBdat.retTim, [], allSigRet')
yticks(1:11)
yticklabels(regNames)
xline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'red')
xline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'red')
title("retrieval significant hit/miss effects")
 export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\' 'HFB_allReg_retrieval.jpg'],''), '-r300')



%% latency 

%chan X stats
%col 1: subID
%col 2: chi
%col 3: center of mass
%col 4: peak latency
%col 5: peak value
%col 6: encode v. retrieve
%col 7: condition
%col 8: region
%col 9: mean RT
%col 10: mean ( centerOfMass / RT)
aovDat = table;
aovDat.subID = repmat("askj", 1000,1); 
aovDat.chi = zeros(1000,1); 
aovDat.centerOfMass = zeros(1000,1); 
aovDat.peakLat = zeros(1000,1); 
aovDat.peakVal = zeros(1000,1); 
aovDat.encRet = repmat("askj", 1000,1); 
aovDat.cond = repmat("askj", 1000,1); 
aovDat.reg = repmat("askj", 1000,1); 
aovDat.RT = zeros(1000,1); 
aovDat.adjTime = zeros(1000,1); 
aovi = 1; 

for ii = 1:length(statFiles)

    HFBdat = load([statFiles(ii).folder '/' statFiles(ii).name]).HFBdat; 
% 
%      tim = HFBdat.encTim; 
%     hitVals = permute(squeeze(HFBdat.regRes(1,:,:)), [2,1]);
%     missVals = permute(squeeze(HFBdat.regRes(2,:,:)), [2,1]);
%     
%    
%     hitMeanTim = arrayfun(@(x) wmean(tim(21:81), hitVals(x,21:81) - min(hitVals(x,21:81)) ), 1:size(hitVals,1));
%     missMeanTim = arrayfun(@(x) wmean(tim(21:81), missVals(x,21:81) - min(missVals(x,21:81)) ), 1:size(missVals,1));
%     
% 
%     HFBdat.hitTim_sub = hitMeanTim; 
%     HFBdat.missTim_sub = missMeanTim; 
% 
%      tim = HFBdat.retTim; 
%     hitVals = permute(squeeze(HFBdat.regRes2(1,:,:)), [2,1]);
%     missVals = permute(squeeze(HFBdat.regRes2(2,:,:)), [2,1]);
%     
%    
%     hitMeanTim = arrayfun(@(x) wmean(tim(21:81), hitVals(x,21:81) - min(hitVals(x,21:81)) ), 1:size(hitVals,1));
%     missMeanTim = arrayfun(@(x) wmean(tim(21:81), missVals(x,21:81) - min(missVals(x,21:81)) ), 1:size(missVals,1));
%     
% 
%     HFBdat.hitTim_ret = hitMeanTim; 
%     HFBdat.missTim_ret = missMeanTim;
    hmSort = zeros(HFBdat.n_pair*2,1);; 
    ti = 1; 
    for pi = 1:HFBdat.n_pair
        hmSort(ti) = 1; 
        ti = ti+1; 
        hmSort(ti) = 0; 
        ti = ti+1; 
       
    end


    figure('position', [100,100,500, 300])
    subplot 121
    hold off
    histogram(HFBdat.hitTim_sub, [0:100:1500])
    hold on
    histogram(HFBdat.missTim_sub, [0:100:1500])
    
    curDat = [HFBdat.hitTim_sub, HFBdat.missTim_sub];

    modDat = table(curDat', hmSort, [HFBdat.realID; HFBdat.realID], ...
                'VariableNames', {'latency', 'hitMiss', 'sub'}); 
    lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
    tval = [0,0];
    tval(1) = lme.Coefficients(2,4); 
    tval(2) = lme.Coefficients(2,6); 
    text(0, HFBdat.n_pair*.1, ['t-val: ' num2str(round(tval(1),1)) ])
    text(0, HFBdat.n_pair*.15, [' p-val: '  num2str(round(tval(2),2))])
    ylabel('time (ms)')
    title(['encoding latency ' HFBdat.aggTargs(HFBdat.reg1).ROI ])


    subplot 122
    hold off
    histogram(HFBdat.hitTim_ret, [0:100:1500])
    hold on
    histogram(HFBdat.missTim_ret, [0:100:1500])

    curDat = [HFBdat.hitTim_ret, HFBdat.missTim_ret];

    modDat = table(curDat', hmSort, [HFBdat.realID; HFBdat.realID], ...
                'VariableNames', {'latency', 'hitMiss', 'sub'}); 
    lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
    tval = [0,0];
    tval(1) = lme.Coefficients(2,4); 
    tval(2) = lme.Coefficients(2,6); 
    text(0, HFBdat.n_pair*.1, ['t-val: ' num2str(round(tval(1),1)) ])
    text(0, HFBdat.n_pair*.15, [' p-val: '  num2str(round(tval(2),2))])
    title(['retrieval latency ' HFBdat.aggTargs(HFBdat.reg1).ROI ])

    
    export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\' 'HFB_latency_' HFBdat.aggTargs(HFBdat.reg1).ROI '.jpg'],''), '-r300')


    cond = { 'subHit', 'subMiss', 'hit_on', 'miss_on'}; 
    grabFrom = {'hitTim_sub', 'missTim_sub', 'hitTim_ret', 'missTim_ret'}; 
    RTgrab = {'subhitRT', 'submissRT', 'rethitRT', 'retmissRT'}; 


    for eleci = 1:HFBdat.n_pair

        for cndi = 1:2
            for hm = 1:2 %loop conditions
            
                %pull data for stats: 
                %col 1: subID
                aovDat(aovi, 1) = {string(HFBdat.regSubIDs{eleci})}; 
                %col 2: chi
                aovDat(aovi, 2) = {HFBdat.chani(eleci)}; 

                %col 3: center of mass
                aovDat(aovi, 3) = {HFBdat.(grabFrom{(cndi-1)*2+hm})(eleci)}; 
                %col 4: peak latency
%                 aovDat(aovi, 4) = {mean(curDatSum(targConditions(ci)).time(curDatSum(targConditions(ci)).peakLat(curChanMask) ))}; 
                %col 5: peak value
%                 aovDat(aovi, 5) = {mean(curDatSum(targConditions(ci)).peakAmp(curChanMask))}; 
                %col 6: encode v. retrieve
        %         aovDat(aovi, 6) = {string(curDatSum(targConditions(ci)).condition)};
                %col 7: condition
                aovDat(aovi, 7) = {cond{(cndi-1)*2+hm}};
                %col 8: region
                aovDat(aovi, 8) = {string(HFBdat.aggTargs(HFBdat.reg1).ROI)}; 
                %col 9: RT
                aovDat(aovi, 9) = {HFBdat.(RTgrab{(cndi-1)*2+hm})(eleci)}; 
                %col 10: mean ( centerOfMass / RT)
%                 aovDat(aovi, 10)= {mean( curDatSum(targConditions(ci)).centerOfMass(curChanMask)./...
%                                          curDatSum(targConditions(ci)).RT(curChanMask) ) };
                aovi = aovi + 1; 
        
            end
        end
    end

end



writetable(aovDat, join(['R:\MSS\Johnson_Lab\dtf8829\' 'latencyAOVdatNEW.csv'],''))



