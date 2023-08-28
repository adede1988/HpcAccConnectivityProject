function [] = makeTimePointPlot(allCon, aggTargs, sigHFBSub, tt, tim)
subplot(1,5,[1:3])
locsX = [9 ,8.5,10,5.5,3,2,1,4,5,6.5,3]; 
locsY = [10,8  ,6 ,6.5,4,9,4,1,9,1.5,6];
allVals = []; 
scatter(locsX, locsY)
hold on 

xlim([0,11])
ylim([0,11])

% colorMap = [logspace(log10(1),log10(10),301)./10; zeros(301,1)'; logspace(log10(10),log10(1),301)./10];
for ii = 1:length(locsX)
    for jj = 1:length(locsX)
%         if ii <=jj
        connection = squeeze(allCon(ii, jj, 2, :, tt)); 
        p = squeeze(allCon(ii,jj,3,:,tt)); 
        sigi = find(p<.15);
        if ~isempty(sigi)
             if tt>90
                    if (ii == 1 || jj ==1) && (ii==3 || jj == 3)
                        x= 5; 
                    end 
                    end
            offSet = mean(sigi);
            [~, coni] = max(abs(connection(sigi)));
            val = connection(sigi(coni));
            allVals = [allVals val];
            offSet = sigi(coni); 
            if val>0 %&& offSet>150 %only leading hit > miss connections
                p1 = [locsX(ii),locsY(ii)]; 
                p2 = [locsX(jj), locsY(jj)]; 
                if offSet>150
                    tmp = p1; 
                    p1 = p2; 
                    p2 = tmp; 
                end
                
                dp = (p2 - p1); 
                quiver(p1(1), p1(2), dp(1), dp(2), 0,'linewidth', (val*100)^2, 'color', 'red')



            elseif val<0 %&& offSet>150
              
                p1 = [locsX(ii),locsY(ii)]; 
                p2 = [locsX(jj), locsY(jj)]; 
                if offSet>150
                    tmp = p1; 
                    p1 = p2; 
                    p2 = tmp; 
                end
                dp = p2 - p1; 
                quiver(p1(1), p1(2), dp(1), dp(2), 0,'linewidth', (val*100)^2, 'color', 'blue')




            end



        end
        

        
    end
end
title([num2str(tim(tt)) ' ms'])
scatter(locsX, locsY)
for ai = 1:length(locsX)
    text(locsX(ai), locsY(ai), aggTargs(ai).ROI)
    if sigHFBSub(ai,3,tt)<.15
        if sigHFBSub(ai,2,tt)>0
            scatter(locsX(ai), locsY(ai), (abs(sigHFBSub(ai,2,tt))*5)^5, 'green', 'filled')
        else
            scatter(locsX(ai), locsY(ai), (abs(sigHFBSub(ai,2,tt))*5)^5, 'magenta', 'filled')
        end
    end
end
axis off
set(gcf, 'color', 'w')
hold off

subplot(1,5,4)
hold on 
vals = [.005: .005: .05]; 
for ii = 1:10

    plot([0,1], [ii,ii], 'linewidth', (vals(ii)*100)^2, 'color', 'k')



end
ylim([0,11])
xlim([-.5,1])
yticks([1:10])
yticklabels(vals)
title('cor diff hit - miss')
xticklabels([])
text(ones(10,1).*-.5, [1:10], arrayfun(@(x) num2str(x), vals, 'UniformOutput',false))
axis off

subplot(1,5,5)
hold on 
vals = [.1: .1: 1]; 
for ii = 1:10

    scatter(1, ii, (vals(ii)*5)^5,'k','filled'  )



end
ylim([0,11])
xlim([-.5,2])
yticks([1:10])
yticklabels(vals)
title('z-score HFB diff hit - miss')
xticklabels([])
text(ones(10,1).*-.5, [1:10], arrayfun(@(x) num2str(x), vals, 'UniformOutput',false))
axis off





end