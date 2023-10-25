function [maxV] = quickPPCplotHelp(PPCvals, tit, regions, reg1, reg2, tim, frex)
    plotMat = squeeze(mean(PPCvals));
    imagesc(tim, [], plotMat')
    yticks([2:2:20])
    yticklabels(round(frex(2:2:20)))
    set(gca, 'ydir', 'normal')
    colorbar
    caxis([0,.1])
    title([tit ' ' regions{reg1} ' ' regions{reg2}])
    ylabel('frequency (Hz)')
    xlabel('time (ms)')
    maxV = max(max(plotMat)); 





end