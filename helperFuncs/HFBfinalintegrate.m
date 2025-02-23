function [] = HFBfinalintegrate(reg, headFiles, ...
    outStatFiles, regions, phase)


%down select files to the current target region and phase of experiment
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, regions{reg}));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, phase));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'stat0'));
HFBperms= outStatFiles(test); 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'stat1'));
Imageperms= outStatFiles(test); 
%do the same for the head files, there will only be 1
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, regions{reg}));
headFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, phase));
headFiles(~test) = []; 


%load in the main data
dat = load([headFiles.folder '/' headFiles.name]).statInfo;
%stats for Image locked
perm = load([Imageperms(1).folder '/' Imageperms(1).name]).outDat;
%aggregate the perms
nullTs = zeros([size(perm.nulls, [1]), length(Imageperms)*50]);
for ii = 1:length(Imageperms)
    perm = load([Imageperms(ii).folder '/' Imageperms(ii).name]).outDat;
    nullTs(:,(ii-1)*50+1:ii*50) = perm.nulls; 

end


[h, p, clusterinfo] = cluster_test(perm.tVals, nullTs);

%for publication figures, I want to grab only the key data that I'm going
%to plot, so I'm flagging what I know from the summary figures is significant: 

outDat = struct; 
tim = dat.tim; 
outDat.reg = regions{reg}; 
outDat.hits = dat.hits(tim>=-450 & tim<=3000,:); 
outDat.misses = dat.misses(tim>=-450 & tim<=3000,:); 
outDat.hitRT = dat.hitRT; 
outDat.missRT = dat.missRT;
outDat.hitLat = dat.hitLat; 
outDat.missLat = dat.missLat; 
outDat.hitSub = dat.hitSub; 
outDat.missSub = dat.missSub; 
outDat.hitChi = dat.hitChi; 
outDat.missChi = dat.missChi; 
outDat.phase = phase; 
tim(tim<-450 | tim>3000) = []; 
outDat.tim = tim; 
save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure1/' ...
    '/HFB_singleTrial' ...
    regions{reg} '_' phase '_image.mat'], "outDat")


realID = arrayfun(@(i) [dat.hitSub{i} '_' num2str(dat.hitChi(i))],...
    1:length(dat.hitChi), 'UniformOutput', false);
uniID = unique(realID); 




outDat = perm; 
outDat.realID = uniID; 
outDat.clusterinfo = clusterinfo; 
outDat.h = h; 
outDat.p = p; 
outDat.tim = tim; 
outDat.reg = regions{reg}; 
outDat.phase = phase; 
outDat.meanHitRT = mean(dat.hitRT); 
outDat.meanMissRT = mean(dat.missRT); 
save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure1/HFB' ...
    regions{reg} '_' phase '_image.mat'], "outDat")

%repeat for HFB locked
perm2 = load([HFBperms(1).folder '/' HFBperms(1).name]).outDat;
%aggregate the perms
nullTs = zeros([size(perm2.nulls, [1]), length(HFBperms)*50]);
for ii = 1:length(HFBperms)
    perm2 = load([HFBperms(ii).folder '/' HFBperms(ii).name]).outDat;
    nullTs(:,(ii-1)*50+1:ii*50) = perm2.nulls; 

end

[h, p2, clusterinfo] = cluster_test(perm2.tVals, nullTs); 

if ~isempty(perm.eliminate)
    keep = 1:size(perm.hitVals,1);
    keep(perm.eliminate) = []; 
    perm.missVals = perm.missVals(keep,:); 
    perm2.missVals = perm2.missVals(keep,:); 
    perm.hitVals = perm.hitVals(keep,:); 
    perm2.hitVals = perm2.hitVals(keep,:); 
end


outDat = perm2; 
outDat.realID = uniID; 
outDat.clusterinfo = clusterinfo; 
outDat.h = h; 
outDat.p = p2; 
outDat.tim = [-500:25:500]; 
outDat.reg = regions{reg}; 
outDat.phase = phase; 
outDat.meanHitRT = mean(dat.hitRT); 
outDat.meanMissRT = mean(dat.missRT); 
save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure2/HFB' ...
    regions{reg} '_' phase '_hfb.mat'], "outDat")


figure('visible', false, 'position', [0,0,600, 1000])

%hit miss Image aligned plot
subplot(4, 2, 1)
hold off
%hits
x = tim; 
y = mean(perm.hitVals);
y(isnan(y)) = 0; 
se = std(perm.hitVals) ./ sqrt(size(perm.hitVals,1));
se(isnan(se)) = 0; 
plot(x, y, 'color', 'blue', 'linewidth', 2)
hold on 
x = [x, flip(x)]; 
y = [y+se, flip(y-se)]; 
fill(x, y, [0 0.4470 0.7410], 'FaceAlpha', .2)
%misses
x = tim; 
y = mean(perm.missVals);
y(isnan(y)) = 0;  
se = std(perm.missVals) ./ sqrt(size(perm.missVals,1));
se(isnan(se)) = 0; 
plot(x, y, 'color', 'red', 'linewidth', 2)
hold on 
x = [x, flip(x)]; 
y = [y+se, flip(y-se)]; 
fill(x, y, 'red', 'FaceAlpha', .2)
xticks(tim([3:16:139]))
xticklabels(tim([3:16:139]))
title([regions{reg} ' Image locked ' phase])

%hit miss Image aligned difference and significance
subplot(4,2,3)
hold off
tim = dat.tim; 
tim(tim<-450 | tim>3000) = []; 
%hits
x = tim; 
y = mean(perm.hitVals) - mean(perm.missVals);
y(isnan(y)) = 0; 
allDifs = perm.hitVals - perm.missVals; 
se = std(allDifs) ./ sqrt(size(allDifs,1));
se(isnan(se)) = 0; 
plot(x, y, 'color', 'k', 'linewidth', 2)
hold on 
x = [x, flip(x)]; 
y = [y+se, flip(y-se)]; 
y(isnan(y)) = 0; 
fill(x, y, 'k', 'FaceAlpha', .2)
yline(0, '--')
xticks(tim([3:16:139]))
xticklabels(tim([3:16:139]))
x = tim; 
x(p>.05) = []; 
breaks = find(diff(x)>25); 
breaks = [0, breaks, length(x)];
for ii = 1:length(breaks)-1
    plot(x(breaks(ii)+1:breaks(ii+1)), zeros(breaks(ii+1) - breaks(ii),1),...
        'color', 'red', 'linewidth', 2)
end
title([regions{reg} ' Image locked ' phase ' difference'])


%hit miss HFB aligned plot
subplot(4, 2, 2)
hold off
tim = [-500:25:500]; 
%hits
x = tim; 
y = mean(perm2.hitVals); 
se = std(perm2.hitVals) ./ sqrt(size(perm2.hitVals,1));
plot(x, y, 'color', 'blue', 'linewidth', 2)
hold on 
x = [x, flip(x)]; 
y = [y+se, flip(y-se)]; 
fill(x, y, 'blue', 'FaceAlpha', .2)
%misses
x = tim; 
y = mean(perm2.missVals); 
se = std(perm2.missVals) ./ sqrt(size(perm2.missVals,1));
plot(x, y, 'color', 'red', 'linewidth', 2)
hold on 
x = [x, flip(x)]; 
y = [y+se, flip(y-se)]; 
fill(x, y, 'red', 'FaceAlpha', .2)
title([regions{reg} ' HFB locked ' phase])


%hit miss HFB aligned difference and significance
subplot(4,2,4)
hold off
tim = [-500:25:500]; 
%hits
x = tim; 
y = mean(perm2.hitVals) - mean(perm2.missVals);
allDifs = perm2.hitVals - perm2.missVals; 
se = std(allDifs) ./ sqrt(size(allDifs,1));
plot(x, y, 'color', 'k', 'linewidth', 2)
hold on 
x = [x, flip(x)]; 
y = [y+se, flip(y-se)]; 
fill(x, y, 'k', 'FaceAlpha', .2)
yline(0, '--')
x = tim; 
x(p2>.05) = []; 
breaks = find(diff(x)>25); 
breaks = [0, breaks, length(x)];
for ii = 1:length(breaks)-1
    plot(x(breaks(ii)+1:breaks(ii+1)), zeros(breaks(ii+1) - breaks(ii),1),...
        'color', 'red', 'linewidth', 2)
end
title([regions{reg} ' HFB locked ' phase ' difference'])


hitVals = dat.hits'; 
missVals = dat.misses';
tim = dat.tim; 

hitVals(:, tim<-450 | tim>3000) = []; 
missVals(:, tim<-450 | tim>3000) = []; 
tim(tim<-450 | tim>3000) = []; 
hitVals(:, 139) = 0; 
missVals(:, 139) = 0; 

subplot(4,2,5)
[sortedLat, order] = sort(dat.hitLat); 
imagesc(tim, [], hitVals(order,:)) 
caxis([-10, 20])
yticks([])
xticks(tim([3:16:139]))
xticklabels(tim([3:16:139]))
xline(0, '--', 'color', 'green', 'linewidth', 2)
title('single trial hits')

subplot(4,2,6)
[sortedLat, order] = sort(dat.missLat); 
imagesc(tim, [], missVals(order,:)) 
caxis([-10, 20])
yticks([])
xticks(tim([3:16:139]))
xticklabels(tim([3:16:139]))
xline(0, '--', 'color', 'green', 'linewidth', 2)
title('single trial misses')

subplot(4,2,7)
hold off
y = [1:.5:100]; 
x = prctile(dat.hitLat, y);
plot(x, y, 'color', 'blue', 'linewidth', 2)
set(gca, 'ydir', 'reverse')
hold on 
x = prctile(dat.missLat, y); 
plot(x, y, 'color', 'red', 'linewidth', 2)
[h,p,ks2stat] = kstest2(dat.hitLat, dat.missLat)
text(1200, 25, ['KS test: p=' num2str(round(p,2))])
title('latency of HFB')
ylabel('percentile')

subplot(4,2,8)
hold off
y = [1:.5:100]; 
x = prctile(dat.hitLat ./ dat.hitRT, y);
plot(x, y, 'color', 'blue', 'linewidth', 2)
set(gca, 'ydir', 'reverse')
hold on 
x = prctile(dat.missLat./ dat.missRT, y); 
plot(x, y, 'color', 'red', 'linewidth', 2)
[h,p,ks2stat] = kstest2(dat.hitLat./ dat.hitRT, dat.missLat./ dat.missRT)
xlim([0,1])
text(.5, 25, ['KS test: p=' num2str(round(p,2))])
title('latency of HFB adjusted for reaction time')
ylabel('percentile')

export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\FinalizedHFB\fieldtripBased'...
    '/' regions{reg} '_' phase '.jpg'], '-r300')

% 
% 
% 
% 
% 
% 
% 
% 
% 
% %miss Image plot
% subplot(4, 2, 3)
% imagesc(squeeze(perm.missVals)')
% set(gca, 'ydir', 'normal')
% caxis([-2, 2])
% colorbar
% xticks([3:16:160])
% xticklabels(tim([3:16:139]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' miss'])
% 
% 
% 
% %hit - miss Image plot
% subplot(4, 2, 5)
% imagesc(squeeze(perm.hitVals)' - ...
%         squeeze(perm.missVals)')
% addRedOutline(p, .05, 'red');
% set(gca, 'ydir', 'normal')
% caxis([-2, 2])
% colorbar
% xticks([3:16:160])
% xticklabels(tim([3:16:139]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit - miss'])
% 
% %tval plot
% subplot(4, 2, 7)
% hold off
% imagesc(perm.tVals')
% addRedOutline(p, .05, 'red');
% set(gca, 'ydir', 'normal')
% caxis([-4, 4])
% colorbar
% xticks([3:16:160])
% xticklabels(tim([3:16:139]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit - miss tVals'])
% 
% 
% %hit HFB plot 
% subplot(4, 2, 2) 
% timCut = [-500:25:500]; 
% imagesc(squeeze(perm2.hitVals)')
% set(gca, 'ydir', 'normal')
% caxis([-2, 10])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% xline(21, 'linewidth', 2, 'linestyle', '--')
% title([regions{reg} ' ' phase ' hit'])
% 
% %miss Image plot
% subplot(4, 2, 4) 
% imagesc(squeeze(perm2.missVals)')
% set(gca, 'ydir', 'normal')
% caxis([-2, 10])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% xline(21, 'linewidth', 2, 'linestyle', '--')
% title([regions{reg} ' ' phase ' miss'])
% 
% 
% 
% %hit - miss Image plot
% subplot(4, 2, 6)
% imagesc(squeeze(perm2.hitVals)' - ...
%         squeeze(perm2.missVals)')
% addRedOutline(p2, .05, 'red');
% set(gca, 'ydir', 'normal')
% caxis([-2, 2])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit - miss'])
% 
% %tval plot
% subplot(4, 2, 8)
% imagesc(perm2.tVals')
% addRedOutline(p2, .05, 'red');
% set(gca, 'ydir', 'normal')
% caxis([-4, 4])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit - miss tVals'])
% 
% 
% 
% export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
%     '/' regions{reg} '_' phase '.jpg'], '-r300')
% 
% 
% figure('visible', false, 'position', [0,0,600,1000])
% 
% subplot(5,2,1)
% hold off
% tim = dat.tim; 
% hitPhase = squeeze(dat.hits_p(:,:,5)); 
% missPhase = squeeze(dat.misses_p(:,:,5)); 
% %time point to align at will be determined by peak ITPC
% test = arrayfun(@(y) arrayfun(@(x) size(hitPhase,2)*...
%     abs(mean(exp(1i*squeeze(dat.hits_p(x,:,y)) ))).^2,...
%     1:size(hitPhase,1) ), 1:100, 'uniformoutput', false);
% test = cat(1, test{:}); 
% imagesc(test)
% set(gca, 'ydir', 'normal')
% caxis([2, 10])
% colorbar
% xticks([1:20:161])
% xticklabels(tim([1:20:131]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit phase reset'])
% 
% subplot(5,2,2)
% hold off
% hitPeakPhase = arrayfun(@(fi)...
%     getPeaks(dat.hits_p(:,:,fi), dat.hitLat, tim),...
%     1:100, 'uniformoutput', false); 
% hitPeakPhase = cat(3, hitPeakPhase{:}); 
% [tL, trialL, frexL] = size(hitPeakPhase);
% %time point to align at will be determined by peak ITPC
% testHFB = arrayfun(@(fi) arrayfun(@(ti) trialL*...
%     abs(mean(exp(1i*squeeze(hitPeakPhase(ti,:,fi)) ))).^2,...
%     1:tL ), 1:frexL, 'uniformoutput', false);
% testHFB = cat(1, testHFB{:}); 
% imagesc(testHFB)
% set(gca, 'ydir', 'normal')
% caxis([2, 10])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit phase reset'])
% 
% 
% subplot(5,2,3)
% hold off
% hitPhase = squeeze(dat.hits_p(:,:,5)); 
% missPhase = squeeze(dat.misses_p(:,:,5)); 
% %time point to align at will be determined by peak ITPC
% test2 = arrayfun(@(y) arrayfun(@(x) size(missPhase,2)*...
%     abs(mean(exp(1i*squeeze(dat.misses_p(x,:,y)) ))).^2,...
%     1:size(missPhase,1) ), 1:100, 'uniformoutput', false);
% test2 = cat(1, test2{:}); 
% imagesc(test2)
% set(gca, 'ydir', 'normal')
% caxis([2, 10])
% colorbar
% xticks([1:20:161])
% xticklabels(tim([1:20:161]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' miss phase reset'])
% 
% subplot(5,2,4)
% hold off
% missPeakPhase = arrayfun(@(fi)...
%     getPeaks(dat.misses_p(:,:,fi), dat.missLat, tim),...
%     1:100, 'uniformoutput', false); 
% missPeakPhase = cat(3, missPeakPhase{:}); 
% [tL, trialL, frexL] = size(missPeakPhase);
% %time point to align at will be determined by peak ITPC
% testHFB2 = arrayfun(@(fi) arrayfun(@(ti) trialL*...
%     abs(mean(exp(1i*squeeze(missPeakPhase(ti,:,fi)) ))).^2,...
%     1:tL ), 1:frexL, 'uniformoutput', false);
% testHFB2 = cat(1, testHFB2{:}); 
% imagesc(testHFB2)
% set(gca, 'ydir', 'normal')
% caxis([2, 10])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit phase reset'])
% 
% subplot(5,2,5)
% imagesc(abs(test - test2))
% set(gca, 'ydir', 'normal')
% caxis([4, 10])
% colorbar
% xticks([1:20:161])
% xticklabels(tim([1:20:161]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit - miss phase reset'])
% 
% subplot(5,2,6)
% imagesc(abs(testHFB - testHFB2))
% set(gca, 'ydir', 'normal')
% caxis([4, 10])
% colorbar
% xticks([1:10:41])
% xticklabels(timCut([1:10:41]))
% yticks([1:19:100])
% yticklabels(round(dat.frex([1:19:100])))
% title([regions{reg} ' ' phase ' hit - miss phase reset'])
% 
% 
% subplot(5,2,7)
% plot(squeeze(perm2.hitVals(1,21,:)), 'color', 'blue',...
%     'linewidth', 2)
% hold on 
% plot(squeeze(perm2.missVals(1,21,:)), 'color', 'red',...
%     'linewidth', 2)
% xticks([1:19:100])
% xticklabels(round(dat.frex([1:19:100])))
% title('broadband power at HFB peak')
% 
% subplot(5,2,8)
% hold off
% plot(squeeze(testHFB(:,21)), 'color', 'blue',...
%     'linewidth', 2)
% hold on 
% plot(squeeze(testHFB2(:,21)), 'color', 'red',...
%     'linewidth', 2)
% xticks([1:19:100])
% xticklabels(round(dat.frex([1:19:100])))
% title('phase reseting at HFB peak')
% 
% % 
% % 
% % subplot(4,2,2)
% % hold off
% % tim =dat.tim; 
% % peakIdx = arrayfun(@(x) find(x==tim), dat.hitLat);
% % vals = arrayfun(@(x) dat.hits_p(peakIdx(x), x, peakToUse),...
% %     [1:length(peakIdx)]); 
% % histogram(vals, [-pi:pi/16:pi], 'normalization', 'probability')
% % hold on 
% % peakIdx = arrayfun(@(x) find(x==tim), dat.missLat);
% % vals = arrayfun(@(x) dat.misses_p(peakIdx(x), x, peakToUse),...
% %     [1:length(peakIdx)]); 
% % histogram(vals, [-pi:pi/16:pi], 'normalization', 'probability')
% % yline(1/32, 'linewidth', 2, 'linestyle', '--')
% % title('phase of peak oscillation at HFB peak')
% % ylim([.01, .05])
% % 
% % 
% % subplot(4,2,3)
% % 
% % 
% % 
% % 
% % subplot(4, 2, 4)
% % 
% % test2 = arrayfun(@(x) size(missPhase,2)*...
% %     abs(mean(exp(1i*missPhase(x,:) ))).^2,...
% %     1:size(missPhase,1) );
% % tim = dat.tim; 
% % plot(tim, test, 'color', 'blue',...
% %     'linewidth', 2)
% % hold on 
% % plot(tim, test2, 'color', 'red',...
% %     'linewidth', 2)
% % xlim([-450, 2500])
% % title("ITPC Rayleigh's z")
% % ylim([0,10])
% % 
% % %get the latencies of the peak frequencies
% % hit_lowLat = gausLat(squeeze(dat.hits(:,:,peakToUse)), tim, dat.hitRT);
% % miss_lowLat = gausLat(squeeze(dat.misses(:,:,peakToUse)), tim, dat.missRT);
% % subplot(4,2,4)
% % hold off
% % histogram(dat.hitLat - hit_lowLat, [-500:25:500],...
% %     'normalization', 'probability')
% % hold on 
% % histogram(dat.missLat -miss_lowLat , [-500:25:500],...
% %     'normalization', 'probability')
% % title('relative timing of HFB peak - low freq peak')
% % xlabel('HFB first               low first')
% % xline(0, 'linewidth', 2)
% % ylim([0,.15])
% % 
% % %plot latencies on their own
% % subplot(4,2,[5,7])
% % hold off
% % [sortedLat,order] = sort(dat.hitLat); 
% % lowPow = squeeze(dat.hits(:,:,peakToUse)); 
% % imagesc(tim, [], lowPow(:,order)')
% % caxis([-5,10])
% % hold on 
% % plot(sortedLat, [1:length(sortedLat)], 'color', 'red', 'linewidth',1)
% % xline(0, 'color', 'green', 'linewidth', 2)
% % title(['low freq pow ' phase ' hits'])
% % xlim([-450, 2500])
% % 
% % 
% % subplot(4,2,[6,8])
% % hold off
% % [sortedLat,order] = sort(dat.missLat); 
% % lowPow = squeeze(dat.misses(:,:,peakToUse)); 
% % imagesc(tim, [], lowPow(:,order)')
% % caxis([-5,10])
% % hold on 
% % plot(sortedLat, [1:length(sortedLat)], 'color', 'red', 'linewidth',1)
% % xline(0, 'color', 'green', 'linewidth', 2)
% % title(['low freq pow ' phase ' misses'])
% % xlim([-450, 2500])
% 
% 
% 
% export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
%     '/' regions{reg} '_' phase '_2.jpg'], '-r300')










%aligning the HFB by the phase of the underlying rhythm didn't seem helpful

% 
% %find alignment phase and time point (based on hits only for now)
% [~, alignPoint] = max(test(10:end-10)); 
% alignPoint = alignPoint+9; 
% angles = hitPhase(alignPoint); 
% sumSin = sum(sin(angles));
% sumCos = sum(cos(angles));
% meanAngle = atan2(sumSin, sumCos);
% % Ensure the result is in the [0, 2*pi) range
% alignPhase = mod(meanAngle + 2*pi, 2*pi);
% 
% %bring in the HFB
% HFBpre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_singleTrial/';
% HFBdat = split(headFiles.name, '_all'); 
% HFBdat = load([HFBpre HFBdat{1} HFBdat{2}]).statInfo;
% 
% hitPhase = squeeze(dat.hits_p(:,:,peakToUse)); 
% hitHFB = HFBdat.hits; 
% alignedHFB = zeros(size(hitPhase)); 
% 
% tMax = max(HFBdat.tim); 
% tMin = min(HFBdat.tim); 
% %how many 25ms timesteps in one cycle of target frequency? 
% cycleLen = round(1000/dat.frex(peakToUse));
% phaseTime = linspace(-pi, pi, cycleLen);
% totalTim = tMax - tMin;
% cycles = round(totalTim / cycleLen); 
% testPhase = repmat(phaseTime, [1,cycles*5]);
% testTim = [tMin:tMax]; 
% phasei = find(testPhase(length(testTim):end) > alignPhase, 1)-1;
% testTim(alignPoint*25:4501) = testPhase(length(testTim)+phasei:...
%                length(testTim)*2+phasei-alignPoint*25); 
% testTim(1:alignPoint*25-1) = ...
%     testPhase(length(testTim)+phasei-alignPoint*25+1:...
%                length(testTim)+phasei-1);
% testTim = testTim([1:25:end]);
% 
% 
% %now, to align data around the time point specified by
% %alignPoint and the phase angle specified by meanAngle
% 
% circularDifference = @(angle1, angle2) mod(angle1 - angle2 + pi, 2*pi) - pi;
% searchDist = round(cycleLen/2/25); 
% 
% for triali = 1:size(hitHFB,2)
%     curHFB = hitHFB(:,triali);
%     curPhase = hitPhase(:,triali);
%     alignedDat = alignedHFB(:, triali); 
%     alignedN = alignedHFB(:, triali); 
% 
%     for tt = 1:length(curHFB) 
%         if tt<=searchDist
%             difVals = abs(circularDifference(curPhase(tt), ...
%                 testTim(1:tt+searchDist)));
%             [~,idx] = min(difVals); 
%         elseif tt>length(curPhase)-searchDist 
%             difVals = abs(circularDifference(curPhase(tt), ...
%                 testTim(tt-searchDist:end)));
%             [~,idx] = min(difVals); 
%             curi = [tt-searchDist:length(curPhase)];
%             idx = curi(idx); 
%         else
%             difVals = abs(circularDifference(curPhase(tt), ...
%                 testTim(tt-searchDist:tt+searchDist)));
%             [~,idx] = min(difVals); 
%             curi = [tt-searchDist:tt+searchDist];
%             idx = curi(idx); 
%         end
%         %for each time point, find the nearest testTim to the curPhase
%         alignedDat(idx) = curHFB(tt); 
%         alignedN(idx) = alignedN(idx)+1; 
%           
% 
%        
% 
%     end
% 
%    
%  nonZeroIndices = alignedN > 0;
% 
% % Calculate mean values at each time point with observations
% meanValues = alignedDat(nonZeroIndices) ./ alignedN(nonZeroIndices);
% 
% % Interpolate mean values where there are no observations
% interpolatedValues = interp1(tim(nonZeroIndices), meanValues, tim,...
%     'linear', 'extrap');
% alignedHFB(:,triali) = interpolatedValues; 
% 
% end
% 
% 
% 
% subplot(4,2,4)
% hold off
% histogram(alignedHFB(testTim>2.5 | testTim<-2.5, :), [-10:.5:30], ...
%     'normalization', 'probability')
% hold on 
% histogram(alignedHFB(testTim>-1 & testTim<1, :), [-10:.5:30], ...
%     'normalization', 'probability')
% 
% 
% 

















end