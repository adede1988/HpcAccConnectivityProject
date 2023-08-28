function [p] = makeScaledHeat(plotMat, xScale, yScale, adjFact)

    plotMat = squeeze(plotMat);
    plotMat = 10.^plotMat; 
    plotMat = plotMat + adjFact - .01; 
    imagesc(squeeze(mean(plotMat))')
%     tVals = myArrayT_oneSamp(nullVals); 
    
   

%     [h, p, clusterinfo] = cluster_test(tVals, conditionNull); 
% 
%     addRedOutline(p, alpha);
    
    yticks([1:2:length(yScale)])
    yticklabels(round(yScale(1:2:length(yScale))))
    xticks([1:20:length(xScale)])
    xticklabels(xScale([1:20:length(xScale)]))
    set(gca, 'YDir','normal')
    colorbar

p=plotMat;


end