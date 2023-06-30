function [aovDat, aovi] = plotConditionCompare(curDatSum, targConditions, regName, histBins, aovDat, aovi)

    b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)']]/255;
    b2w2r(129:end, 1) = 1; 
    b2w2r(1:128, 3) = 1; 

    figure('visible', 'off', 'Position', [1,1,800,1200])
    %hardCode the time periods for threshold testing
    threshTestWin = [[-50 2500]; [-50 2000]; [-1500, 500]];
    colors = {[75, 122, 71]./255, [236, 146, 72]./255};

    for ci = 1:2
        %overall heatmap of all trials
    %eliminate long RT trials
    badTrials = curDatSum(targConditions(ci)).RT>3000; 
    sum(badTrials)
    curDatSum(targConditions(ci)).HFB(:,badTrials) = []; 
    curDatSum(targConditions(ci)).RT(badTrials) = []; 
    curDatSum(targConditions(ci)).subID(badTrials) = []; 
    curDatSum(targConditions(ci)).chi(badTrials) = []; 
    curDatSum(targConditions(ci)).peakLat(badTrials) = []; 
    curDatSum(targConditions(ci)).peakAmp(badTrials) = []; 
    curDatSum(targConditions(ci)).threshCross(badTrials) = []; 
    curDatSum(targConditions(ci)).maxAmp(badTrials) = []; 
    curDatSum(targConditions(ci)).amp(badTrials) = []; 
    curDatSum(targConditions(ci)).dur(badTrials) = []; 
    curDatSum(targConditions(ci)).centerOfMass(badTrials) = []; 



   
    [sortedRT, order] = sort(curDatSum(targConditions(ci)).RT ); 
   

    subplot(4,2, [ci])
    imagesc(curDatSum(targConditions(ci)).HFB(:,order)')
    caxis([-7,7])
    hold on 
    plot(sortedRT/25 + find(curDatSum(targConditions(ci)).time>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
    xticks([9:16:size(curDatSum(targConditions(ci)).HFB,1)])
    xticklabels(curDatSum(targConditions(ci)).time([9:16:size(curDatSum(targConditions(ci)).HFB,1)]))
    xline(find(curDatSum(targConditions(ci)).time>=0,1), '--', 'linewidth', 4, 'color', 'green')
    title([regName ' ' curDatSum(targConditions(ci)).condition], 'interpreter', 'none')


    %time course of mean channels
    subplot(4,2, 3)
    hold on 
    subIDs = unique(curDatSum(targConditions(ci)).subID); 
    chanMeans = []; 
    chanMeanAmp = []; 
    for subi = 1:length(subIDs)
        tempDati = cellfun(@(x) strcmp(subIDs{subi}, x), curDatSum(targConditions(ci)).subID);
        tempGood = tempDati; 
        chanNums = unique(curDatSum(targConditions(ci)).chi(tempDati)); 
        for chani = 1:length(chanNums)
            curChanMask = arrayfun(@(x) x==chanNums(chani), curDatSum(targConditions(ci)).chi);
            curChanMask = curChanMask==1 & tempDati==1; 
            curMeanTrial = mean(curDatSum(targConditions(ci)).HFB(:,curChanMask), 2); 
            curMeanAmp = mean(curDatSum(targConditions(ci)).peakAmp(curChanMask)); 
            chanMeanAmp = [chanMeanAmp curMeanAmp]; 
            chanMeans = [chanMeans, curMeanTrial];

            %pull data for stats: 
            %col 1: subID
            aovDat(aovi, 1) = {string(subIDs{subi})}; 
            %col 2: chi
            aovDat(aovi, 2) = {chanNums(chani)}; 
            %col 3: center of mass
            aovDat(aovi, 3) = {mean(curDatSum(targConditions(ci)).centerOfMass(curChanMask))}; 
            %col 4: peak latency
            aovDat(aovi, 4) = {mean(curDatSum(targConditions(ci)).time(curDatSum(targConditions(ci)).peakLat(curChanMask) ))}; 
            %col 5: peak value
            aovDat(aovi, 5) = {mean(curDatSum(targConditions(ci)).peakAmp(curChanMask))}; 
            %col 6: encode v. retrieve
            aovDat(aovi, 6) = {string(curDatSum(targConditions(ci)).condition)};
            %col 7: condition
            aovDat(aovi, 7) = {string(curDatSum(targConditions(ci)).condition)};
            %col 8: region
            aovDat(aovi, 8) = {string(regName)}; 
            %col 9: RT
            aovDat(aovi, 9) = {mean(curDatSum(targConditions(ci)).RT(curChanMask))}; 
            %col 10: mean ( centerOfMass / RT)
            aovDat(aovi, 10)= {mean( curDatSum(targConditions(ci)).centerOfMass(curChanMask)./...
                                     curDatSum(targConditions(ci)).RT(curChanMask) ) };
            aovi = aovi + 1; 



        end
    end
    cndMean = mean(chanMeans,2);
    cndStd = std(chanMeans,[],2) ./ sqrt(size(chanMeans,2)); 
    cndLow = cndMean - cndStd; %prctile(allHit, 2.5, 2); 
    cndHigh = cndMean + cndStd; %prctile(allHit, 97.5, 2); 
    plot(cndMean, 'color', colors{ci}, 'linewidth', 2)
    hold on 
    xticks([101:100:size(curDatSum(targConditions(ci)).HFB,1)])
    xticklabels(curDatSum(targConditions(ci)).time([101:100:size(curDatSum(targConditions(ci)).HFB,1)]))
    xline(find(curDatSum(targConditions(ci)).time>=0,1), '--', 'linewidth', 4, 'color', 'green')
    x = [1:length(cndMean)]; 
    x = [x flip(x)]; 
    y = [cndLow' flip(cndHigh')];
    h = fill(x,y,colors{ci},'LineStyle','none'); 
    set(h, 'facealpha', .5)
    xlim([5,length(cndMean)-5])
    title('mean of channel means')
    xlabel('time (ms)')


    subplot(4,2,4)
    hold on 
    histogram(chanMeanAmp, [0:5:50], ...
        'FaceColor', colors{ci}, 'normalization', 'probability')
    title('channel mean peak amplitude')
    xlabel('HFB (z-score)')

    subplot(4,2,5)
    hold on 
    cndMean = mean(curDatSum(targConditions(ci)).HFB,2);
    cndStd = std(curDatSum(targConditions(ci)).HFB,[],2) ./ sqrt(size(curDatSum(targConditions(ci)).HFB,2)); 
    cndLow = cndMean - cndStd; %prctile(allHit, 2.5, 2); 
    cndHigh = cndMean + cndStd; %prctile(allHit, 97.5, 2); 

    plot(cndMean, 'color', colors{ci}, 'linewidth', 2)
    hold on 
    xticks([101:100:size(curDatSum(targConditions(ci)).HFB,1)])
    xticklabels(curDatSum(targConditions(ci)).time([101:100:size(curDatSum(targConditions(ci)).HFB,1)]))
    xline(find(curDatSum(targConditions(ci)).time>=0,1), '--', 'linewidth', 4, 'color', 'green')
    x = [1:length(cndMean)]; 
    x = [x flip(x)]; 
    y = [cndLow' flip(cndHigh')];
    h = fill(x,y,colors{ci},'LineStyle','none'); 
    set(h, 'facealpha', .5)
    xlim([5,length(cndMean)-5])
    title('mean of individual trials')
    xlabel('time (ms)')

    subplot(4,2,6)
    hold on 
    histogram(curDatSum(targConditions(ci)).peakAmp, [0:1:50], ...
        'FaceColor', colors{ci}, 'normalization', 'probability')
    title('peak amplitude')
    xlabel('HFB (z-score)')
    

    subplot(4,2,[7])
    hold on 
    histogram(curDatSum(targConditions(ci)).time(curDatSum(targConditions(ci)).peakLat),histBins, ...
        'FaceColor', colors{ci}, 'normalization', 'probability')
    title('peak latency')
    xlabel('time (ms)')

    subplot(4,2,[8])
    hold on 

    crossTimes = curDatSum(targConditions(ci)).centerOfMass;
    crossTimes(isnan(crossTimes)) = []; 

    histogram(crossTimes, histBins, ...
        'FaceColor', colors{ci}, 'normalization', 'probability')
    title('center of mass')
    xlabel('time (ms)')

    
    



    end
    subplot(4,2,8)
    legend({curDatSum(targConditions).condition})
    set(gcf,'color','w');
   
    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\HFBsummary_' regName '_' ...
        curDatSum(targConditions(1)).condition '_' ...
        curDatSum(targConditions(2)).condition '.jpg'])

%     subplot(3,2, [2,4])
%     [sortedRT, order] = sort(curDatSum(targConditions(1)).RT ); 
%     subplot(3,2, [1,3])
%     imagesc(curDatSum(targConditions(1)).HFB(:,order)')
%     caxis([-7,7])
%     hold on 
%     plot(sortedRT/5 + find(curDatSum(targConditions(1)).time>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
%     xticks([101:100:size(curDatSum(targConditions(1)).HFB,1)])
%     xticklabels(curDatSum(targConditions(1)).time([101:100:size(curDatSum(targConditions(1)).HFB,1)]))
%     xline(find(curDatSum(targConditions(1)).time>=0,1), '--', 'linewidth', 4, 'color', 'green')
%     title([regName ' ' curDatSum(targConditions(1)).condition])
% 







end