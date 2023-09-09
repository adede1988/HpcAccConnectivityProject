function [] = makeTimePointPlot2(allCon, tt, tim)
subplot(1,5,[1:4])
locsX = [9 ,8.5,10,5.5,3.5,2,1,4,5,6.5,3]; 
locsY = [10,8  ,6 ,6.5,2.5,9,4,1,9,1.5,6];

scatter(locsX, locsY, 50, 'k', 'filled')
hold on 
regNames = {'dlPFC', 'mlPFC', 'piPFC', 'ACC', 'lTemp', 'Par', 'Vis', 'iTemp', 'motor', 'MTL', 'PCC'}; 
xlim([0,11])
ylim([0,11])

% colorMap = [logspace(log10(1),log10(10),301)./10; zeros(301,1)'; logspace(log10(10),log10(1),301)./10];
for ii = 1:length(locsX)
    for jj = 1:length(locsX)
%         if ii <=jj
        connection = squeeze(allCon(ii, jj, tt)); 
       if connection>0
        px = [locsX(ii),locsX(jj)]; 
        py = [locsY(ii), locsY(jj)]; 
               
                 
        plot(px, py, 'linewidth', (connection*100)^2, 'color', 'red')

       elseif connection < 0 
        px = [locsX(ii),locsX(jj)]; 
        py = [locsY(ii), locsY(jj)]; 
               
                 
        plot(px, py, 'linewidth', (connection*100)^2, 'color', 'blue')

       end

        
    end
end
title([num2str(tim(tt)) ' ms'])
scatter(locsX, locsY)
for ai = 1:length(locsX)
    text(locsX(ai), locsY(ai), regNames{ai})
end
axis off
set(gcf, 'color', 'w')
hold off

subplot(1,5,5)
hold on 
vals = [.005: .005: .05]; 
for ii = 1:10

    plot([0,1], [ii,ii], 'linewidth', (vals(ii)*100)^2, 'color', 'k')



end
ylim([0,11])
xlim([-.5,1])
yticks([1:10])
yticklabels(vals)
title('PPC diff hit - miss')
xticklabels([])
text(ones(10,1).*-.5, [1:10], arrayfun(@(x) num2str(x), vals, 'UniformOutput',false))
axis off

% subplot(1,5,5)
% hold on 
% vals = [.1: .1: 1]; 
% for ii = 1:10
% 
%     scatter(1, ii, (vals(ii)*5)^5,'k','filled'  )
% 
% 
% 
% end
% ylim([0,11])
% xlim([-.5,2])
% yticks([1:10])
% yticklabels(vals)
% title('z-score HFB diff hit - miss')
% xticklabels([])
% text(ones(10,1).*-.5, [1:10], arrayfun(@(x) num2str(x), vals, 'UniformOutput',false))
% axis off





end