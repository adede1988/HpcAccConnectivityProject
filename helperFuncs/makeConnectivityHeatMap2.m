function [] = makeConnectivityHeatMap2(inMat, p, ss, ...
    frex, tim)

%called from publicationFigs.m

imagesc(inMat')
colorbar
if max(inMat(:) > 1)
    caxis([-4, 4])
% else
%     caxis([.02, .1])
end
set(gca, 'ydir', 'normal')
addRedOutline(p, .05, 'white'); 
yticks([2:6:20]) %hard code that there are 20 frequencies 
yticklabels(round(frex([2:6:20])))
if ss == 1 %hard coded time for HFB locked
    xticks([1:10:41])
    xticklabels([-500:250:500])
    xline(21, 'linewidth', 2, 'color', 'green', 'linestyle', '--')
else
    xticks([9:20:length(tim)])
    xticklabels(tim([9:20:length(tim)]))
    xline(19, 'linewidth', 2, 'color', 'green', 'linestyle', '--')
end
        
% xlabel(['time aligned to ' regions{regii} ' ' statName{ss}])
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
ylabel('frequency (Hz)')
 set(gca, 'FontSize', 24);


end