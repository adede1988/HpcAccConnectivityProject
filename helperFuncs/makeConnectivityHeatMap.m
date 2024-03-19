function [] = makeConnectivityHeatMap(inMat, p, ss, dat, regions, ...
    regii, regjj, phases, pp, statName)


imagesc(inMat')
colorbar
set(gca, 'ydir', 'normal')
addRedOutline(p, .05, 'red'); 
yticks([2:4:20])
yticklabels(round(dat.frex([2:4:20])))
if ss == 1
    xticks([1:5:41])
    xticklabels([-500:125:500])
    xline(21, 'linewidth', 2, 'color', 'green', 'linestyle', '--')
else
    xticks([9:16:161])
    xticklabels(dat.tim([9:16:161]))
    xline(41, 'linewidth', 2, 'color', 'green', 'linestyle', '--')
end
        
% xlabel(['time aligned to ' regions{regii} ' ' statName{ss}])

ylabel('frequency')



end