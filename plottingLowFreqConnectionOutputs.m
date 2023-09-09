%% plot results for connectivity analyses


%local paths: 

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% single parameter model 

datFolder = [datPre 'LowFreqConDat_memOnly']; 
cndFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({cndFiles.name}, '.mat'));
cndFiles = cndFiles(test); 



for ii = 1:length(cndFiles)

    cd = load([cndFiles(ii).folder '/' cndFiles(ii).name]).connectionDat;
    if isfield(cd, 'reg2')
        regVals = {cd.reg1, cd.reg2}; 
    else
        regVals = {cd.reg, 'all'}; 
    end

    
    figure('visible', false, 'Position', [0 0 1000 2000])
    tim = cd.tim; 

    subplot 421
    lowDat = cd.lowBand(cd.hmSort,:,2); 
    yyaxis left
    plot(tim, mean(lowDat), 'linewidth', 2, 'color', '#540b0e')
    ylabel('mean over electrode pairs')
    ylim([-.02, .05])
    yyaxis right
    plot(tim, var(lowDat), 'linewidth', 2, 'color', '#e09f3e')
    ylabel('var over electrode pairs')
    title(['PPC at 3 Hz between ' regVals{1} ' and ' regVals{2}])
   

    subplot 422
    highDat = cd.highBand(cd.hmSort,:,2); 
    yyaxis left
    plot(tim, mean(highDat), 'linewidth', 2, 'color', '#540b0e')
    ylabel('mean over electrode pairs')
    ylim([-.02, .05])
    yyaxis right
    plot(tim, var(highDat), 'linewidth', 2, 'color', '#e09f3e')
    ylabel('var over electrode pairs')
    title(['PPC at 8-9 Hz between ' regVals{1} ' and ' regVals{2}])


    subplot 423
    hold off
    lowDat = cd.lowBand(cd.hmSort,:,2); 
    [vals, order] = sort(cd.d); 
    imagesc(tim, [], lowDat(order, :))
    caxis([0,.1])
    subBorders = find(diff(vals)>0); 
    for si = 1:length(subBorders)
        yline(subBorders(si)+.5, 'linewidth', 2, 'color', [1,.55,0,.5], 'linestyle', '--')
    end
    yticklabels('')
    ylabel('good memory    bad memory')
    title(['PPC at 3 Hz between ' regVals{1} ' and ' regVals{2}])

    subplot 424
    hold off
    highDat = cd.highBand(cd.hmSort,:,2); 
    [vals, order] = sort(cd.d); 
    imagesc(tim, [], highDat(order, :))
    caxis([0,.1])
    subBorders = find(diff(vals)>0); 
    for si = 1:length(subBorders)
        yline(subBorders(si)+.5, 'linewidth', 2, 'color', [1,.55,0,.5], 'linestyle', '--')
    end
    yticklabels('')
    ylabel('good memory    bad memory')
    title(['PPC at 8-9 Hz between ' regVals{1} ' and ' regVals{2}])

    subplot 425
    hold off
    plot(tim, cd.lowtVals, 'linewidth', 3, 'color', 'k')
    hold on 
    plot(tim, cd.low975, 'linewidth', 1, 'color', 'red', 'linestyle', '--')
    plot(tim, cd.low025, 'linewidth', 1, 'color', 'red', 'linestyle', '--')
    sigPlot = zeros(size(tim));
    sigPlot(cd.lowpRaw<.05) = max([cd.low975; cd.lowtVals]);
    sigX = tim(sigPlot>0); 
    sigPlot(sigPlot==0) = []; 
    scatter(sigX, sigPlot, 'green', 'filled')
    sigPlot = zeros(size(tim)); 
    sigPlot(cd.lowp<.05) = max([cd.low975; cd.lowtVals])*.9;
    sigX = tim(sigPlot>0); 
    sigPlot(sigPlot==0) = []; 
    scatter(sigX, sigPlot, 'magenta', 'filled')
    title('rank(ppc) ~ rank(memory) + (1|sub)')
    ylabel('t-value for memory')

    subplot 426
    hold off
    plot(tim, cd.hightVals, 'linewidth', 3, 'color', 'k')
    hold on 
    plot(tim, cd.high975, 'linewidth', 1, 'color', 'red', 'linestyle', '--')
    plot(tim, cd.high025, 'linewidth', 1, 'color', 'red', 'linestyle', '--')
    sigPlot = zeros(size(tim));
    sigPlot(cd.highpRaw<.05) = max([cd.high975; cd.hightVals]);
    sigX = tim(sigPlot>0); 
    sigPlot(sigPlot==0) = []; 
    scatter(sigX, sigPlot, 'green', 'filled')
    sigPlot = zeros(size(tim)); 
    sigPlot(cd.highp<.05) = max([cd.high975; cd.hightVals])*.9;
    sigX = tim(sigPlot>0); 
    sigPlot(sigPlot==0) = []; 
    scatter(sigX, sigPlot, 'magenta', 'filled')

    subplot 427
    hold off
    [vals, order] = sort(abs(cd.lowtVals)); 
    jitter = .1*(rand(length(cd.dOrig),1)-.5)';
    scatter(cd.dOrig+jitter, lowDat(:, order(end)), 20, [.9492, .7344, .1797], 'filled')
    subMeans = zeros(length(cd.uniqueSubs),1);
    for si = 1:length(subMeans)
        subMeans(si) = mean(lowDat(cellfun(@(x) strcmp(cd.uniqueSubs{si}, x), cd.allSubs), order(end)));
    end
    hold on 
    scatter(cd.subDOrig, subMeans, 40, [156,39,6]./256, 'filled')
    ylabel('PPC')
    xlabel("memory performance (d')")
    title(['max t: ' num2str(round(vals(end),2)) ' at time: ' num2str(tim(order(end))) ])

    subplot 428
    hold off
    [vals, order] = sort(abs(cd.hightVals)); 
    jitter = .1*(rand(length(cd.dOrig),1)-.5)';
    scatter(cd.dOrig+jitter, highDat(:, order(end)), 20, [.9492, .7344, .1797], 'filled')
    subMeans = zeros(length(cd.uniqueSubs),1);
    for si = 1:length(subMeans)
        subMeans(si) = mean(highDat(cellfun(@(x) strcmp(cd.uniqueSubs{si}, x), cd.allSubs), order(end)));
    end
    hold on 
    scatter(cd.subDOrig, subMeans, 40, [156,39,6]./256, 'filled')
    ylabel('PPC')
    xlabel("memory performance (d')")
    title(['max t: ' num2str(round(vals(end),2)) ' at time: ' num2str(tim(order(end))) ])


    set(gcf,'color','w');
    
    
    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\ISPC_regional_stats\' 'singleVar_' regVals{1} '_' regVals{2} '.jpg'])









end



%% 



datFolder = [datPre 'LowFreqConDat_mem_HM']; 
cndFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({cndFiles.name}, '.mat'));
cndFiles = cndFiles(test); 

%get a matrix of all significant connections across time
allCon = zeros([11,11,length(tim)]);

for ii = 1:length(cndFiles)

    cd = load([cndFiles(ii).folder '/' cndFiles(ii).name]).connectionDat;
    if isfield(cd, 'reg2')
        regVals = {cd.reg1, cd.reg2}; 
        regidx = [cd.reg1i, cd.reg2i]; 
    else
        regVals = {cd.reg, 'all'}; 
        regidx = [cd.regi, cd.regi]; 
    end

    
%     figure('visible', false, 'Position', [0 0 1000 2000])
%     tim = cd.tim; 
%     
%     subplot 421
%     
%     lowDat = cd.lowBand(cd.hmSort,:,2); 
%     hold off
%     plot(tim, mean(lowDat), 'linewidth', 2, 'color', '#283618')
%     hold on 
%     lowDat = cd.lowBand(~cd.hmSort,:,2); 
%     plot(tim, mean(lowDat), 'linewidth', 2, 'color', '#bc6c25', 'linestyle', '-')
%     ylim([-.02, .05])
%     ylabel('mean over electrode pairs')
% 
%     title(['PPC at 3 Hz between ' regVals{1} ' and ' regVals{2}])
%    
% 
%     subplot 422
%     
%     highDat = cd.highBand(cd.hmSort,:,2); 
%     hold off
%     plot(tim, mean(highDat), 'linewidth', 2, 'color', '#283618')
%     hold on 
%     highDat = cd.highBand(~cd.hmSort,:,2); 
%     plot(tim, mean(highDat), 'linewidth', 2, 'color', '#bc6c25', 'linestyle', '-')
%     ylim([-.02, .05])
%     ylabel('mean over electrode pairs')
% 
%     title(['PPC at 8 Hz between ' regVals{1} ' and ' regVals{2}])
% 
% 
%     subplot 423
%     hold off
    lowDat = cd.lowBand(cd.hmSort,:,2) - cd.lowBand(~cd.hmSort, :,2); 
%     [vals, order] = sort(cd.d(1:length(cd.dOrig))); 
%     imagesc(tim, [], lowDat(order, :))
%     caxis([-.1,.1])
%     subBorders = find(diff(vals)>0); 
%     for si = 1:length(subBorders)
%         yline(subBorders(si)+.5, 'linewidth', 2, 'color', '#E11584', 'linestyle', '--')
%     end
%     yticklabels('')
%     ylabel('good memory    bad memory')
%     title(['Hit - Miss PPC at 3 Hz between ' regVals{1} ' and ' regVals{2}])
% 
%     subplot 424
%     hold off
%     highDat = cd.highBand(cd.hmSort,:,2) - cd.highBand(~cd.hmSort, :,2); 
%     [vals, order] = sort(cd.d(1:length(cd.dOrig))); 
%     imagesc(tim, [], highDat(order, :))
%     caxis([-.1,.1])
%     subBorders = find(diff(vals)>0); 
%     for si = 1:length(subBorders)
%         yline(subBorders(si)+.5, 'linewidth', 2, 'color', '#E11584', 'linestyle', '--')
%     end
%     yticklabels('')
%     ylabel('good memory    bad memory')
%     title(['Hit - Miss PPC at 8-9 Hz between ' regVals{1} ' and ' regVals{2}])

%     subplot 425
%     hold off
%     plot(tim, cd.lowtVals(:,1), 'linewidth', 3, 'color', 'k')
%     hold on 
%     plot(tim, cd.low975(:,1), 'linewidth', 1, 'color', 'red', 'linestyle', '--')
%     plot(tim, cd.low025(:,1), 'linewidth', 1, 'color', 'red', 'linestyle', '--')
%     sigPlot = zeros(size(tim));
%     sigPlot(cd.lowpRaw(:,1)<.05) = max([cd.low975(:,1); cd.lowtVals(:,1)]);
%     sigX = tim(sigPlot>0); 
%     sigPlot(sigPlot==0) = []; 
%     scatter(sigX, sigPlot, 'green', 'filled')
    sigPlot = zeros(size(tim)); 
    test = sum(cd.lowCrossP(:,:,1)<.05,1) == size(cd.lowCrossP,1); 
    sigPlot(test) = max([cd.low975(:,1); cd.lowtVals(:,1)])*.9;
    sigX = tim(sigPlot>0); 
    sigPlot(sigPlot==0) = []; 
%     scatter(sigX, sigPlot, 'magenta', 'filled')
%     title('rank(ppc) ~ rank(memory)* hit/miss + (1|sub)')
%     ylabel('t-value for hit/miss')
    timidx = arrayfun(@(x) find(tim==x), sigX);
    meanDif = mean(lowDat); 
    allCon(regidx(1), regidx(2), timidx) = meanDif(timidx); 
    
% 
%     subplot 426
%     hold off
%     plot(tim, cd.hightVals(:,1), 'linewidth', 3, 'color', 'k')
%     hold on 
%     plot(tim, cd.high975(:,1), 'linewidth', 1, 'color', 'red', 'linestyle', '--')
%     plot(tim, cd.high025(:,1), 'linewidth', 1, 'color', 'red', 'linestyle', '--')
% %     sigPlot = zeros(size(tim));
% %     sigPlot(cd.highpRaw<.05) = max([cd.high975; cd.hightVals]);
% %     sigX = tim(sigPlot>0); 
% %     sigPlot(sigPlot==0) = []; 
% %     scatter(sigX, sigPlot, 'green', 'filled')
%     sigPlot = zeros(size(tim)); 
%     sigPlot(cd.highp(:,1)<.05) = max([cd.high975(:,1); cd.hightVals(:,1)])*.9;
%     sigX = tim(sigPlot>0); 
%     sigPlot(sigPlot==0) = []; 
%     scatter(sigX, sigPlot, 'magenta', 'filled')
%     title('rank(ppc) ~ rank(memory)* hit/miss + (1|sub)')
%     ylabel('t-value for hit/miss')
% 
%     subplot 427
%     hold off
%     idx = find(cd.lowp(:,1)<.2); 
%     if isempty(idx)
%         idx = 1:length(cd.highp(:,1)); 
%     end
%     [vals, order] = sort(abs(cd.lowtVals(idx,1))); 
%     plotVals = cd.lowBand(:,idx(order(end)), 2); 
% 
%     b = boxchart(cd.hmSort*1, plotVals);
%     b.MarkerStyle = 'none'; 
%     xticks([0,1])
%     xticklabels({'subMiss', 'subHit'})
%     title(['PPC 3 Hz at ' num2str(tim(idx(order(end)))) ' ms'])
%     PLH = plotVals(cd.hmSort); 
%     PLM = plotVals(~cd.hmSort); 
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
%     ylabel('PPC')
% 
% 
%     subplot 428
%     hold off
%     idx = find(cd.highp(:,1)<.2); 
%     if isempty(idx)
%         idx = 1:length(cd.highp(:,1)); 
%     end
%     [vals, order] = sort(abs(cd.hightVals(idx,1))); 
%     plotVals = cd.highBand(:,idx(order(end)), 2); 
% 
%     b = boxchart(cd.hmSort*1, plotVals);
%     b.MarkerStyle = 'none'; 
%     xticks([0,1])
%     xticklabels({'subMiss', 'subHit'})
%     title(['PPC 8 Hz at ' num2str(tim(idx(order(end)))) ' ms'])
%     PLH = plotVals(cd.hmSort); 
%     PLM = plotVals(~cd.hmSort); 
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
%     ylabel('PPC')
% 
% 
% 
%     set(gcf,'color','w');
%     
%     
%     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\ISPC_regional_stats\' 'doubleVar_' regVals{1} '_' regVals{2} '.jpg'])


    









end


%make a video of the hit > miss 3 Hz significant connections
hitOverMiss = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\ISPC_regional_stats\' 'doubleVar_hitOverMiss_sig.avi'],''));
open(hitOverMiss); 
f = figure;
f.Position = [100 100 1000 600];
for tt = 1:141
    makeTimePointPlot2(allCon, tt, tim)
    frame = getframe(gcf);
    writeVideo(hitOverMiss, frame); 
end

close(hitOverMiss)



