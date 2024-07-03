function latency = visualizeHFBsingleTrialDat(allRes, reg, regions, phase)
    statInfo = struct; %output for stats
    statInfo.reg = regions{reg}; 
    statInfo.regi = reg; 
    statInfo.phase = phase; 

    %create vector to identify the data for the target region 
    test = cellfun(@(x) strcmp(x, regions{reg}), {allRes{:,5}});
    curReg = allRes(test,:);
    hits = [];
    misses = []; 
    hitRT = [];
    missRT = []; 
    hitChi = []; 
    missChi = []; 
    hitTriali = []; 
    missTriali = []; 
    hitSub = cell(5000, 1); 
    missSub = cell(5000,1); 
    hitAge = []; 
    hitAcc = []; 
    hitd = []; 
    missAge = [];
    missAcc = [];
    missd = [];
    hi = 1; 
    mi = 1; 


    for ii = 1:size(curReg,1)
        hits = [hits, curReg{ii,1}];
        hitRT = [hitRT; curReg{ii, 3}]; 
        misses = [misses, curReg{ii,2}];
        missRT = [missRT; curReg{ii, 4}]; 
        

        L = size(curReg{ii,1},2);
        hitChi = [hitChi; repmat(curReg{ii, 8}, L, 1) ]; 
        hitAge = [hitAge; repmat(curReg{ii, 21}, L, 1) ]; 
        hitAcc = [hitAcc; repmat(curReg{ii, 20}, L, 1) ]; 
        hitd = [hitd; repmat(curReg{ii, 19}, L, 1) ]; 
        hitTriali = [hitTriali; [1:L]' ]; 
        for jj = 0:L-1
            hitSub{hi+jj} = curReg{ii, 7}; 
        end
        hi = hi + L;
        L1 = L; 

        L = size(curReg{ii,2},2);
        missChi = [missChi; repmat(curReg{ii, 8}, L, 1) ];
        missAge = [missAge; repmat(curReg{ii, 21}, L, 1) ];
        missAcc = [missAcc; repmat(curReg{ii, 20}, L, 1) ];
        missd = [missd; repmat(curReg{ii, 19}, L, 1) ];
        missTriali = [missTriali; [L1+1:L1+L]' ];
        for jj = 0:L-1
            missSub{mi+jj} = curReg{ii, 7}; 
        end
        mi = mi + L; 
%         missSub = [missSub; repmat(curReg{ii, 7}, size(curReg{ii,2},2), 1) ];
    end
    hitSub(hi:end) = []; 
    missSub(mi:end) = []; 

    tim = curReg{1,6}; 
   

    %calcualte the latency values
    latency = gausLat(hits, tim, hitRT, 1, 5);
    latency2 = gausLat(misses, tim, missRT, 1, 5); 



    eliminate = latency == -1;
    eliminate2 = latency2 == -1; 

    latency(eliminate) = []; 
    latency2(eliminate2) = []; 
    hits(:,eliminate) = [];
    hitRT(eliminate) = []; 
    misses(:,eliminate2) = []; 
    missRT(eliminate2) = []; 
    hitChi(eliminate) = []; 
    missChi(eliminate2) = []; 
    hitSub(eliminate) = []; 
    missSub(eliminate2) = []; 
    hitTriali(eliminate) = []; 
    missTriali(eliminate2) = []; 
    hitAge(eliminate) = []; 
    hitAcc(eliminate) = []; 
    hitd(eliminate) = []; 
    missAge(eliminate2) = [];
    missAcc(eliminate2) = [];
    missd(eliminate2) = [];



    statInfo.hits = hits; 
    statInfo.misses = misses; 
    statInfo.tim = tim; 
    statInfo.hitRT = hitRT; 
    statInfo.missRT = missRT;
    statInfo.hitLat = latency; 
    statInfo.missLat = latency2; 
    statInfo.hitChi = hitChi; 
    statInfo.missChi = missChi; 
    statInfo.hitSub = hitSub; 
    statInfo.missSub = missSub; 
    statInfo.hitTriali = hitTriali; 
    statInfo.missTriali = missTriali; 
    statInfo.hitAge = hitAge; 
    statInfo.hitAcc = hitAcc ;
    statInfo.hitd = hitd ;
    statInfo.missAge = missAge;
    statInfo.missAcc = missAcc;
    statInfo.missd = missd;

    save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_singleTrial\' ...
        statInfo.reg '_' statInfo.phase '.mat'], "statInfo")
    

    
   
   %plotting subsequent all trials
      %plotting subsequent all trials
    figure('visible', true, 'position', [0,0,600,1000])

    subplot(4,2,5)
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

    subplot(4, 2, [1,3])
    [~, order] = sort(latency); 
%     [~, order] = sort(mean(encHit((10+sortIDX)-3:(10+sortIDX)+3,:)), 'descend'); 
    imagesc(tim(10:end-10), [], hits(10:end-10,order)')
    caxis([-10,10])
    hold on 
%     plot(hitRT(order), [1:length(hitRT)], 'color', 'red', 'linewidth', 1)
    xline(0, 'color', 'green', 'linewidth', 3)

%     plot(latency(order), [1:length(encHitRT)], 'color', 'green')
    title([regions{reg} phase ' Hit'])

    

    subplot(4, 2, [2,4])


    
    [~, order] = sort(mean(misses((10+sortIDX)-3:(10+sortIDX)+3,:)), 'descend'); 
%      [~, order] = sort(encMissRT); 
    [~, order] = sort(latency2);
    imagesc(tim(10:end-10), [], misses(10:end-10,order)')
       caxis([-10,10])
%     hold on 
%     plot(encMissRT(order), [1:length(encMissRT)], 'color', 'red', 'linewidth', 3)
    xline(0, 'color', 'green', 'linewidth', 3)
%     plot(latency2(order), [1:length(encMissRT)], 'color', 'green')
    title([regions{reg} phase ' Miss'])

     subplot(4,2,6)
     hold off
     w = 20; 

     hitLatIdx = arrayfun(@(x) find(x==tim), latency); 
     hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
     peakHit = arrayfun(@(x) hits(hitLatIdx(x)-w:...
                                         hitLatIdx(x)+w,x), [1:length(hitLatIdx)], ...
                                         'uniformoutput', false);
     peakHit = cell2mat(peakHit);
     y = mean(peakHit, 2); 
     x = [-w*25:25:w*25]';
     plot(x, y, 'color', 'blue', 'linewidth', 2); 
     hold on 

     e = std(peakHit,[],2) ./ sqrt(length(latency));
    
     fill([x; flip(x)], [y-e; flip(y+e)], 'blue', 'facealpha', .2, 'edgealpha', .2)

    


     missLatIdx = arrayfun(@(x) find(x==tim), latency2); 
     missLatIdx(missLatIdx>length(tim)-w) = length(tim)-w; 
     peakMiss = arrayfun(@(x) misses(missLatIdx(x)-w:...
                                         missLatIdx(x)+w,x), [1:length(missLatIdx)], ...
                                         'uniformoutput', false);
     peakMiss = cell2mat(peakMiss);

     y2 = mean(peakMiss, 2); 
     x = [-w*25:25:w*25]';
     plot(x, y2, 'color', 'red', 'linewidth', 2); 
     hold on 

     e = std(peakMiss,[],2) ./ sqrt(length(latency));
    
     fill([x; flip(x)], [y2-e; flip(y2+e)], 'red', 'facealpha', .2, 'edgealpha', .2)

    subplot(4,2,8)
    hold off
    plot(x, y-y2, 'color', 'k', 'linewidth', 2)

     %do a quick t-test across time 
     tVals = zeros(length(y), 1); 
     for ti = 1:w*2+1
        tVals(ti) =ttest2(peakMiss(ti,:), peakHit(ti,:));
     end

     subplot(4,2,7)
     hold off
     latency = latency';% ./ hitRT; 
     latency2 = latency2';% ./ missRT; 

     latPrc = arrayfun(@(x) prctile(latency, x), [0:.5:100]);
     latPrc2 = arrayfun(@(x) prctile(latency2, x), [0:.5:100]);


     plot( latPrc,[0:.5:100], 'color', 'blue', 'linewidth', 2)
     hold on 

     plot( latPrc2,[0:.5:100], 'color', 'red', 'linewidth', 2)
    set(gca, 'YDir','reverse')
%      histogram(latency'./hitRT, [0:.05:1], 'normalization', 'probability')
%      hold on 
%      histogram(latency2'./missRT, [0:.05:1], 'normalization', 'probability')
     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\RTfix\'...
         'FinalizedHFB\' 'singleTrial_' ...
         phase '_new_' regions{reg} '.jpg'], '-r300')


    latency = {latency ./ hitRT, latency2 ./ missRT};  



















end