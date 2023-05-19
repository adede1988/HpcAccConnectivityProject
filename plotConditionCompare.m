function [] = plotConditionCompare(curDatSum, targConditions, regName, histBins)

    b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)']]/255;
    b2w2r(129:end, 1) = 1; 
    b2w2r(1:128, 3) = 1; 

    figure
    %hardCode the time periods for threshold testing
    threshTestWin = [[-50 2500]; [-50 2000]; [-1500, 500]];
    colors = {[75, 122, 71]./255, [236, 146, 72]./255};

    for ci = 1:2
    [sortedRT, order] = sort(curDatSum(targConditions(ci)).RT ); 
    subplot(4,2, [ci,ci+2])
    imagesc(curDatSum(targConditions(ci)).HFB(:,order)')
    caxis([-7,7])
    hold on 
    plot(sortedRT/5 + find(curDatSum(targConditions(ci)).time>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
    xticks([101:100:size(curDatSum(targConditions(ci)).HFB,1)])
    xticklabels(curDatSum(targConditions(ci)).time([101:100:size(curDatSum(targConditions(ci)).HFB,1)]))
    xline(find(curDatSum(targConditions(ci)).time>=0,1), '--', 'linewidth', 4, 'color', 'green')
    title([regName ' ' curDatSum(targConditions(ci)).condition], 'interpreter', 'none')


    subplot(4,2,[5,6])
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


    

    subplot(4,2,[7])
    hold on 
    histogram(curDatSum(targConditions(ci)).time(curDatSum(targConditions(ci)).peakLat),histBins, ...
        'FaceColor', colors{ci}, 'normalization', 'probability')

    subplot(4,2,[8])
    hold on 

    crossTimes = curDatSum(targConditions(ci)).centerOfMass;
    crossTimes(isnan(crossTimes)) = []; 

    histogram(crossTimes, histBins, ...
        'FaceColor', colors{ci}, 'normalization', 'probability')







    end


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