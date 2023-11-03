function [] = quickPPCtrialPlot_TFHelper(hits, misses, hitLat, missLat, tim, ...
    phase, regNam, fi, hitRT, missRT, freqLab)

figure('visible', true, 'position', [0,0,1000,1000])
  
        subplot(4,2,5)
        hold off
        x = tim(10:end-10)'; 
        y = mean(hits(:, 10:end-10,fi),[1,3]);
        e = std(hits(:, 10:end-10,fi),[],[1,3]) ./ sqrt(length(hitLat));
        plot(x, y, 'color', 'blue', 'linewidth', 2); 
        hold on 
        fill([x; flip(x)], [y'-e'; flip(y'+e')], 'blue', 'facealpha', .2, 'edgealpha', .2)

        y2 = mean(misses(:, 10:end-10,fi),[1,3]);
        plot(x,y2, 'color', 'red', 'linewidth', 2); 
        e = std(misses(:, 10:end-10,fi),[],[1,3]) ./ sqrt(length(missLat));
        fill([x; flip(x)], [y2'-e'; flip(y2'+e')], 'red', 'facealpha', .2, 'edgealpha', .2)

        xlim([tim(10), tim(end-10)])

        %get the hit TF latencies to compare to the HFB latencies
        hitTFLat = gausLat(squeeze(mean(hits(:,:,fi),3))', tim, hitRT);
        missTFLat = gausLat(squeeze(mean(misses(:,:,fi),3))', tim, missRT);
        subplot(4,2,7)
        hold off
        histogram(hitTFLat, [0:75:3000], ...
            'edgecolor', 'blue', 'facecolor', 'blue',...
            'edgealpha', .4, 'facealpha', .4, ...
            'normalization', 'probability')
        hold on 
        histogram(missTFLat , [0:75:3000], ...
            'edgecolor', 'blue', 'facecolor', 'red',...
            'edgealpha', .4, 'facealpha', .4, ...
            'normalization', 'probability')
       xlabel('Theta Lat (ms)')
        subplot(4,2,8)
        hold off
        histogram(hitTFLat - hitLat, [-2000:75:2000], ...
            'edgecolor', 'blue', 'facecolor', 'blue',...
            'edgealpha', .4, 'facealpha', .4, ...
            'normalization', 'probability')
        hold on 
        histogram(missTFLat - missLat, [-2000:75:2000], ...
            'edgecolor', 'blue', 'facecolor', 'red',...
            'edgealpha', .4, 'facealpha', .4, ...
            'normalization', 'probability')
       xline(0, 'linewidth', 2, 'linestyle', '--')
       xlabel('Theta Lat - HFB Lat (ms)')


        subplot(4, 2, [1,3])
        hold off
        [~, order] = sort(hitLat(:,1)); 
        imagesc(tim(10:end-10), [], squeeze(mean(hits(order, 10:end-10,fi),3))) 
           
        caxis([1,10])
        hold on 
        plot(hitLat(order), [1:length(order)], 'color', 'red', 'linewidth', 1)
        xline(0, 'color', 'green', 'linewidth', 3)
        title([char(regNam) ' ' char(phase) ' Hit'])


        subplot(4, 2, [2,4])
        hold off
        [~, order] = sort(missLat(:,1)); 
        imagesc(tim(10:end-10), [], squeeze(mean(misses(order, 10:end-10,fi),3))) 
           
        caxis([1,10])
        hold on 
        plot(missLat(order), [1:length(order)], 'color', 'red', 'linewidth', 1)
        xline(0, 'color', 'green', 'linewidth', 3)
        title([char(regNam) ' ' char(phase) ' miss'])


         subplot(4,2,6)
     hold off
     w = 20; 

     hitLatIdx = arrayfun(@(x) find(x==tim), hitLat); 
     hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
     peakHit = arrayfun(@(x) mean(hits(x, hitLatIdx(x)-w:...
                                         hitLatIdx(x)+w,fi),3), [1:length(hitLatIdx)], ...
                                         'uniformoutput', false);
     peakHit = cell2mat(peakHit);
     peakHit = reshape(peakHit, w*2+1,[]);
     y = mean(peakHit, 2); 
     x = [-w*25:25:w*25]';
     plot(x, y, 'color', 'blue', 'linewidth', 2); 
     hold on 

     e = std(peakHit,[],2) ./ sqrt(length(hitLat));
    
     fill([x; flip(x)], [y-e; flip(y+e)], 'blue', 'facealpha', .2, 'edgealpha', .2)

    


     missLatIdx = arrayfun(@(x) find(x==tim), missLat); 
     missLatIdx(missLatIdx>length(tim)-w) = length(tim)-w; 
     peakMiss = arrayfun(@(x) mean(misses(x, missLatIdx(x)-w:...
                                         missLatIdx(x)+w,fi),3), [1:length(missLatIdx)], ...
                                         'uniformoutput', false);
     peakMiss = cell2mat(peakMiss);
     peakMiss = reshape(peakMiss, w*2+1,[]);
     y2 = mean(peakMiss, 2); 
     x = [-w*25:25:w*25]';
     plot(x, y2, 'color', 'red', 'linewidth', 2); 
     hold on 

     e = std(peakMiss,[],2) ./ sqrt(length(missLat));
    
     fill([x; flip(x)], [y2-e; flip(y2+e)], 'red', 'facealpha', .2, 'edgealpha', .2)

     xline(0, 'linewidth', 2, 'linestyle', '--')

%      subplot(4,2,8)
%         plot(x, y - y2, 'color', 'k', 'linewidth', 2)
   
%     if fi(1)==1
%         freqLab = 'low';
%     else
%         freqLab = 'high';
%     end
     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\pairedOnly\'...
         regNam '_' char(phase) '_' freqLab '.jpg'], '-r300')




end