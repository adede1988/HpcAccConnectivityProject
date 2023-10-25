function latency = visualizeHFBsingleTrialDat_RTsort(allRes, reg, regions, phase)

    %create vector to identify the data for the target region 
    test = cellfun(@(x) strcmp(x, regions{reg}), {allRes{:,5}});
    curReg = allRes(test,:);
    hits = [];
    misses = []; 
    hitRT = [];
    missRT = []; 
    for ii = 1:size(curReg,1)
        hits = [hits, curReg{ii,1}];
        hitRT = [hitRT; curReg{ii, 3}]; 
        misses = [misses, curReg{ii,2}];
        missRT = [missRT; curReg{ii, 4}]; 
    end
    tim = curReg{1,6}; 

    eliminate = []; 
    latency = []; 
    
    %scale within the 0 to RT time period to make 0 to 1
    for tt = 1:size(hits,2)
        if hitRT(tt)<100 || hitRT(tt)>3000
            eliminate = [eliminate, tt]; 
        
        else
            trial = hits(:,tt); 
           
            test = gausswin(11);
            test = test ./ sum(test); 
            trial = [zeros(5,1); trial; zeros(5,1)];
            smoothT = conv(trial, test, 'valid'); 
            [~, idx] = max(smoothT(find(tim==0):find(tim>hitRT(tt),1))); 
            testTim = tim(find(tim==0):find(tim>hitRT(tt),1));
           
            latency = [latency, testTim(idx)];



        end
    end


    eliminate2 = []; 
    latency2 = []; 
    for tt = 1:length(missRT)
        if missRT(tt)<100 || missRT(tt)>3000
            eliminate2 = [eliminate2, tt]; 
        else
            trial = misses(:,tt); 
            test = gausswin(11);
            test = test ./ sum(test); 
            trial = [zeros(5,1); trial; zeros(5,1)];
            smoothT = conv(trial, test, 'valid'); 
            [~, idx] = max(smoothT(find(tim==0):find(tim>missRT(tt),1))); 
            testTim = tim(find(tim==0):find(tim>missRT(tt),1));
            latency2 = [latency2,  testTim(idx)];



        end
    end

    
    hits(:,eliminate) = [];
    hitRT(eliminate) = []; 
    misses(:,eliminate2) = []; 
    missRT(eliminate2) = []; 
    
    
    

   
   %plotting subsequent all trials
      %plotting subsequent all trials
    figure('visible', true, 'position', [0,0,600,1000])

    subplot(3,2,[5,6])
    x = tim(10:end-10)'; 
    y = mean(hits(10:end-10,:),2);
    e = std(hits(10:end-10,:),[],2) ./ sqrt(length(latency));
    plot(x, y, 'color', 'blue', 'linewidth', 2); 
    hold on 
    fill([x; flip(x)], [y-e; flip(y+e)], 'blue', 'facealpha', .2, 'edgealpha', .2)

    y = mean(misses(10:end-10,:),2);
    plot(x,y, 'color', 'red', 'linewidth', 2); 
    fill([x; flip(x)], [y-e; flip(y+e)], 'red', 'facealpha', .2, 'edgealpha', .2)

    xlim([tim(10), tim(end-10)])


    difVals = mean(hits(10:end-10,:),2) - mean(misses(10:end-10,:),2);
    [~, sortIDX] = max(movmean(abs(difVals), 5));

    subplot(3, 2, [1,3])
    [~, order] = sort(hitRT); 
%     [~, order] = sort(mean(encHit((10+sortIDX)-3:(10+sortIDX)+3,:)), 'descend'); 
    imagesc(tim(10:end-10), [], hits(10:end-10,order)')
    caxis([-5,10])
    hold on 
    plot(hitRT(order), [1:length(hitRT)], 'color', 'red', 'linewidth', 1)
    xline(0, 'color', 'green', 'linewidth', 3)

%     plot(latency(order), [1:length(encHitRT)], 'color', 'green')
    title([regions{reg} phase ' Hit'])

    

    subplot(3, 2, [2,4])


    
%     [~, order] = sort(mean(misses((10+sortIDX)-3:(10+sortIDX)+3,:)), 'descend'); 
     [~, order] = sort(missRT); 
%     [~, order] = sort(latency2);
    imagesc(tim(10:end-10), [], misses(10:end-10,order)')
       caxis([-5,10])
    hold on 
    plot(missRT(order), [1:length(missRT)], 'color', 'red', 'linewidth', 1)
    xline(0, 'color', 'green', 'linewidth', 3)
%     plot(latency2(order), [1:length(encMissRT)], 'color', 'green')
    title([regions{reg} phase ' Miss'])

%      subplot(4,2,6)
%      hold off
%      w = 20; 
% 
%      hitLatIdx = arrayfun(@(x) find(x==tim), latency); 
%      hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
%      peakHit = arrayfun(@(x) hits(hitLatIdx(x)-w:...
%                                          hitLatIdx(x)+w,x), [1:length(hitLatIdx)], ...
%                                          'uniformoutput', false);
%      peakHit = cell2mat(peakHit);
%      y = mean(peakHit, 2); 
%      x = [-w*25:25:w*25]';
%      plot(x, y, 'color', 'blue', 'linewidth', 2); 
%      hold on 
% 
%      e = std(peakHit,[],2) ./ sqrt(length(latency));
%     
%      fill([x; flip(x)], [y-e; flip(y+e)], 'blue', 'facealpha', .2, 'edgealpha', .2)
% 
%     
% 
% 
%      missLatIdx = arrayfun(@(x) find(x==tim), latency2); 
%      missLatIdx(missLatIdx>length(tim)-w) = length(tim)-w; 
%      peakMiss = arrayfun(@(x) misses(missLatIdx(x)-w:...
%                                          missLatIdx(x)+w,x), [1:length(missLatIdx)], ...
%                                          'uniformoutput', false);
%      peakMiss = cell2mat(peakMiss);
% 
%      y2 = mean(peakMiss, 2); 
%      x = [-w*25:25:w*25]';
%      plot(x, y2, 'color', 'red', 'linewidth', 2); 
%      hold on 
% 
%      e = std(peakMiss,[],2) ./ sqrt(length(latency));
%     
%      fill([x; flip(x)], [y2-e; flip(y2+e)], 'red', 'facealpha', .2, 'edgealpha', .2)
% 
%     subplot(4,2,8)
%     hold off
%     plot(x, y-y2, 'color', 'k', 'linewidth', 2)
% 
%      %do a quick t-test across time 
%      tVals = zeros(length(y), 1); 
%      for ti = 1:w*2+1
%         tVals(ti) =ttest2(peakMiss(ti,:), peakHit(ti,:));
%      end
% 
%      subplot(4,2,7)
%      hold off
%      latency = latency';% ./ hitRT; 
%      latency2 = latency2';% ./ missRT; 
% 
%      latPrc = arrayfun(@(x) prctile(latency, x), [0:.5:100]);
%      latPrc2 = arrayfun(@(x) prctile(latency2, x), [0:.5:100]);
% 
% 
%      plot( latPrc,[0:.5:100], 'color', 'blue', 'linewidth', 2)
%      hold on 
% 
%      plot( latPrc2,[0:.5:100], 'color', 'red', 'linewidth', 2)
%     set(gca, 'YDir','reverse')
%      histogram(latency'./hitRT, [0:.05:1], 'normalization', 'probability')
%      hold on 
%      histogram(latency2'./missRT, [0:.05:1], 'normalization', 'probability')
set(gcf, 'color', 'w')
     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\' 'singleTrial_' ...
         phase '_RTsort_' regions{reg} '.jpg'], '-r300')


    latency = {latency ./ hitRT, latency2 ./ missRT};  



















end