function [outDat2] = TFfinalintegrate(reg, headFiles, ...
    outStatFiles, regions, phase, outStatFilesPhase)
try

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

%down select the phase permutation files
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, regions{reg}));
outStatFilesPhase(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, phase));
outStatFilesPhase(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'stat0'));
HFBpermsPhase= outStatFilesPhase(test); 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'stat1'));
ImagepermsPhase= outStatFilesPhase(test); 

%do the same for the head files, there will only be 1
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, regions{reg}));
headFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, phase));
headFiles(~test) = []; 


%load in the main data
dat = load([headFiles.folder '/' headFiles.name]).statInfo;
perm = load([Imageperms(1).folder '/' Imageperms(1).name]).outDat;
%aggregate the perms
nullTs = zeros([size(perm.nulls, [1,2]), length(Imageperms)*50]);
for ii = 1:length(Imageperms)
    perm = load([Imageperms(ii).folder '/' Imageperms(ii).name]).outDat;
    nullTs(:,:,(ii-1)*50+1:ii*50) = perm.nulls; 

end

[h, p, clusterinfo] = cluster_test(perm.tVals, nullTs); 

perm2 = load([HFBperms(1).folder '/' HFBperms(1).name]).outDat;
%aggregate the perms
nullTs = zeros([size(perm2.nulls, [1,2]), length(HFBperms)*50]);
for ii = 1:length(HFBperms)
    perm2 = load([HFBperms(ii).folder '/' HFBperms(ii).name]).outDat;
    nullTs(:,:,(ii-1)*50+1:ii*50) = perm2.nulls; 

end

[h, p2, clusterinfo] = cluster_test(perm2.tVals, nullTs); 


%permutations on TF power took the mean, 
% but I want channel level data back
[perm, perm2] = fixMissing(dat, perm, perm2);


%% create output structures for later publication plotting
%power

%image locked power 
outDat = struct; 
tim = dat.tim; 
tim(tim<-450 | tim>3000) = []; 
outDat.tim = tim; 
outDat.reg = regions{reg}; 
outDat.frex = dat.frex; 
outDat.hits_image = perm.hitVals; 
outDat.misses_image = perm.missVals;  
outDat.p_image = p; 
outDat.t_image = perm.tVals; 
outDat.hitRT = dat.hitRT; 
outDat.missRT = dat.missRT;
outDat.hitLat = dat.hitLat; 
outDat.missLat = dat.missLat; 
outDat.hitChi = dat.hitChi; 
outDat.missChi = dat.missChi; 
outDat.hitSub = dat.hitSub; 
outDat.missSub = dat.missSub; 
outDat.phase = phase; 
save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure2/' ...
    '/TF_' ...
    regions{reg} '_' phase '_image.mat'], "outDat")


%HFB locked power
outDat = struct; 
tim = [-500:25:500];
outDat.tim = tim; 
outDat.reg = regions{reg}; 
outDat.frex = dat.frex; 
outDat.hits_hfb = perm2.hitVals; 
outDat.misses_hfb = perm2.missVals;  
outDat.p_hfb = p2; 
outDat.t_hfb = perm2.tVals; 
outDat.hitRT = dat.hitRT; 
outDat.missRT = dat.missRT;
outDat.hitLat = dat.hitLat; 
outDat.missLat = dat.missLat; 
outDat.hitChi = dat.hitChi; 
outDat.missChi = dat.missChi; 
outDat.hitSub = dat.hitSub; 
outDat.missSub = dat.missSub; 
outDat.phase = phase; 
save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure2/' ...
    '/TF_' ...
    regions{reg} '_' phase '_HFB.mat'], "outDat")
 



scatterFig = figure('visible', true, 'position', [0,0,1000,1000]);

outDat2 = cell(2,2,2); 
if strcmp(phase, 'sub')
    outi = 1; 
else
    outi = 2; 
end

set(0, 'currentfigure', scatterFig);
subplot 441
makeHFBImageScatter(perm,perm2, 1:23, 'power')
title([regions{reg} ' ' phase ' freq: ' num2str(round(dat.frex(1),1)) ':'...
       num2str(round(dat.frex(23),1)) 'Hz power'])

subplot 442
outDat2(1,1,:) = makeHFBImageHist(perm, perm2, 1:23);

subplot(4,4,5)
makeHFBTimeScatter(dat, perm2, 1:23, 'power')

subplot(4,4,6)
makeIndexTimeScatter(dat, perm, perm2, 1:23)

subplot(4,4,9)
makeHFBImageScatter(perm,perm2, 23:38, 'power')
title([regions{reg} ' ' phase ' freq: ' num2str(round(dat.frex(23),1)) ':'...
       num2str(round(dat.frex(38),1)) 'Hz power'])

subplot(4,4,10)
outDat2(2, 1,:) = makeHFBImageHist(perm, perm2, 23:38);

subplot(4,4,13)
makeHFBTimeScatter(dat, perm2, 23:38, 'power')

subplot(4,4,14)
makeIndexTimeScatter(dat, perm, perm2, 23:38)






figure('visible', true, 'position', [0,0,600, 1000])
%hit Image plot
subplot(4, 2, 1)
tim = dat.tim; 
tim(tim<-450 | tim>3000) = []; 
imagesc(squeeze(mean(perm.hitVals))')
set(gca, 'ydir', 'normal')
caxis([-2, 2])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit'])

%miss Image plot
subplot(4, 2, 3)
imagesc(squeeze(mean(perm.missVals))')
set(gca, 'ydir', 'normal')
caxis([-2, 2])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' miss'])



%hit - miss Image plot
subplot(4, 2, 5)
imagesc(squeeze(mean(perm.hitVals))' - ...
        squeeze(mean(perm.missVals))')
addRedOutline(p, .05, 'red');
set(gca, 'ydir', 'normal')
caxis([-2, 2])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit - miss'])

%tval plot
subplot(4, 2, 7)
hold off
imagesc(perm.tVals')
addRedOutline(p, .05, 'red');
set(gca, 'ydir', 'normal')
caxis([-4, 4])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit - miss tVals'])


%hit HFB plot 
subplot(4, 2, 2) 
timCut = [-500:25:500]; 
imagesc(squeeze(mean(perm2.hitVals))')
set(gca, 'ydir', 'normal')
caxis([-2, 10])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
xline(21, 'linewidth', 2, 'linestyle', '--')
title([regions{reg} ' ' phase ' hit'])

%miss Image plot
subplot(4, 2, 4) 
imagesc(squeeze(mean(perm2.missVals))')
set(gca, 'ydir', 'normal')
caxis([-2, 10])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
xline(21, 'linewidth', 2, 'linestyle', '--')
title([regions{reg} ' ' phase ' miss'])



%hit - miss Image plot
subplot(4, 2, 6)
imagesc(squeeze(mean(perm2.hitVals))' - ...
        squeeze(mean(perm2.missVals))')
addRedOutline(p2, .05, 'red');
set(gca, 'ydir', 'normal')
caxis([-2, 2])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit - miss'])

%tval plot
subplot(4, 2, 8)
imagesc(perm2.tVals')
addRedOutline(p2, .05, 'red');
set(gca, 'ydir', 'normal')
caxis([-4, 4])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit - miss tVals'])



export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' regions{reg} '_' phase '.jpg'], '-r300')


%% repeat the process for the phase locking statistics 

perm = load([ImagepermsPhase(1).folder '/'...
    ImagepermsPhase(1).name]).outDat;
%aggregate the perms
nullTs = zeros([size(perm.nulls, [1,2]), length(ImagepermsPhase)*50]);
for ii = 1:length(ImagepermsPhase)
    perm = load([ImagepermsPhase(ii).folder '/' ...
        ImagepermsPhase(ii).name]).outDat;
    nullTs(:,:,(ii-1)*50+1:ii*50) = perm.nulls; 

end

[h, p, clusterinfo] = cluster_test(perm.tVals, nullTs); 

perm2 = load([HFBpermsPhase(1).folder '/' HFBpermsPhase(1).name]).outDat;
%aggregate the perms
nullTs = zeros([size(perm2.nulls, [1,2]), length(HFBpermsPhase)*50]);
for ii = 1:length(HFBpermsPhase)
    perm2 = load([HFBpermsPhase(ii).folder '/' HFBpermsPhase(ii).name]).outDat;
    nullTs(:,:,(ii-1)*50+1:ii*50) = perm2.nulls; 

end

[h, p2, clusterinfo] = cluster_test(perm2.tVals, nullTs); 


if ~isempty(perm.eliminate)
    keep = 1:size(perm.hitVals,1);
    keep(perm.eliminate) = []; 
    perm.missVals = perm.missVals(keep,:,:); 
    perm2.missVals = perm2.missVals(keep,:,:); 
    perm.hitVals = perm.hitVals(keep,:,:); 
    perm2.hitVals = perm2.hitVals(keep,:,:); 
end

%% create output structures for later publication plotting
%phase

%image locked phase 
outDat = struct; 
tim = dat.tim; 
tim(tim<-450 | tim>3000) = []; 
outDat.tim = tim; 
outDat.reg = regions{reg}; 
outDat.frex = dat.frex; 
outDat.hits_image = perm.hitVals; 
outDat.misses_image = perm.missVals;  
outDat.p_image = p; 
outDat.t_image = perm.tVals;
outDat.hitRT = dat.hitRT; 
outDat.missRT = dat.missRT;
outDat.hitLat = dat.hitLat; 
outDat.missLat = dat.missLat; 
outDat.hitChi = dat.hitChi; 
outDat.missChi = dat.missChi; 
outDat.hitSub = dat.hitSub; 
outDat.missSub = dat.missSub; 
outDat.phase = phase; 

%get the mean phase angle for each channel
[~,idx] = max(sum(outDat.p_image<.05, 2)); %when is max significance
maxTim = outDat.tim(idx); %need to translate to uncut tim
outDat.maxTim = maxTim; 
uncutIdx = find(dat.tim==maxTim);
hitAngles = squeeze(dat.hits_p(uncutIdx, :, :)); 
hitrealID = arrayfun(@(x) [dat.hitSub{x} '_' num2str(dat.hitChi(x))], ...
    [1:length(dat.hitChi)], 'uniformoutput', false); 
uni_ID = unique(hitrealID); 
meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 12th frequency which is 3 Hz
    tmp = hitAngles(ismember(hitrealID, uni_ID{ui}), 12); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.hit_angles = meanAngles; 

meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 33rd frequency which is 6.5 Hz
    tmp = hitAngles(ismember(hitrealID, uni_ID{ui}), 33); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.hit_angles_high = meanAngles; 
%do it again for misses: 
missAngles = squeeze(dat.misses_p(uncutIdx, :, :)); 
missrealID = arrayfun(@(x) [dat.missSub{x} '_' num2str(dat.missChi(x))], ...
    [1:length(dat.missChi)], 'uniformoutput', false); 
uni_ID = unique(missrealID); 
meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 12th frequency which is 3 Hz
    tmp = missAngles(ismember(missrealID, uni_ID{ui}), 12); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.miss_angles = meanAngles; 

meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 33rd frequency which is 6.5 Hz
    tmp = missAngles(ismember(missrealID, uni_ID{ui}), 33); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.miss_angles_high = meanAngles; 

save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure2/' ...
    '/TFphase_' ...
    regions{reg} '_' phase '_image.mat'], "outDat")


%HFB locked phase
outDat = struct; 
tim = [-500:25:500];
outDat.tim = tim; 
outDat.reg = regions{reg}; 
outDat.frex = dat.frex; 
outDat.hits_hfb = perm2.hitVals; 

%get the mean phase angle for each channel
idx = arrayfun(@(x) find(dat.tim>=x,1), dat.hitLat); 
hitAngles = arrayfun(@(x,y) squeeze(dat.hits_p(x,y,:)),...
    idx', [1:length(idx)], 'uniformoutput', false);
hitAngles =cat(2, hitAngles{:});
hitrealID = arrayfun(@(x) [dat.hitSub{x} '_' num2str(dat.hitChi(x))], ...
    [1:length(dat.hitChi)], 'uniformoutput', false); 
uni_ID = unique(hitrealID); 
meanAngles = []; 
for ui = 1:length(uni_ID)
    tmp = hitAngles(12, ismember(hitrealID, uni_ID{ui})); %hard code of the 12th frequency which is 3 Hz
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.hit_angles = meanAngles; 

meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 33rd frequency which is 6.5 Hz
    tmp = hitAngles(33, ismember(hitrealID, uni_ID{ui})); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.hit_angles_high = meanAngles; 

outDat.misses_hfb = perm2.missVals;  
%get the mean phase angle for each channel
idx = arrayfun(@(x) find(dat.tim>=x,1), dat.missLat); 
missAngles = arrayfun(@(x,y) squeeze(dat.misses_p(x,y,:)),...
    idx', [1:length(idx)], 'uniformoutput', false);
missAngles =cat(2, missAngles{:});
missrealID = arrayfun(@(x) [dat.missSub{x} '_' num2str(dat.missChi(x))], ...
    [1:length(dat.missChi)], 'uniformoutput', false); 
uni_ID = unique(missrealID); 
meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 12th frequency which is 3 Hz
    tmp = missAngles(12, ismember(missrealID, uni_ID{ui})); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.miss_angles = meanAngles;

meanAngles = []; 
for ui = 1:length(uni_ID)
    %hard code of the 33rd frequency which is 6.5 Hz
    tmp = missAngles(33, ismember(missrealID, uni_ID{ui})); 
    meanAngles = [meanAngles, angle(mean(exp(1i * tmp)))];
%     figure
%     histogram(tmp, [-pi:pi/6:pi])
end
outDat.miss_angles_high = meanAngles;
outDat.p_hfb = p2; 
outDat.hitRT = dat.hitRT; 
outDat.missRT = dat.missRT;
outDat.hitLat = dat.hitLat; 
outDat.missLat = dat.missLat; 
outDat.hitChi = dat.hitChi; 
outDat.missChi = dat.missChi; 
outDat.hitSub = dat.hitSub; 
outDat.missSub = dat.missSub; 
outDat.phase = phase; 
save(['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData/Figure2/' ...
    '/TFphase_' ...
    regions{reg} '_' phase '_HFB.mat'], "outDat")







set(0, 'currentfigure', scatterFig);
subplot 443
makeHFBImageScatter(perm,perm2, 1:23, 'phase')
title([regions{reg} ' ' phase ' freq: ' num2str(round(dat.frex(1),1)) ':'...
       num2str(round(dat.frex(23),1)) 'Hz phase'])

subplot 444
outDat2(1, 2,:) = makeHFBImageHist(perm, perm2, 1:23);

subplot(4,4,7)
makeHFBTimeScatter(dat, perm2, 1:23, 'phase')

subplot(4,4,8)
makeIndexTimeScatter(dat, perm, perm2, 1:23)

subplot(4,4,11)
makeHFBImageScatter(perm,perm2, 23:38, 'phase')
title([regions{reg} ' ' phase ' freq: ' num2str(round(dat.frex(23),1)) ':'...
       num2str(round(dat.frex(38),1)) 'Hz phase'])

subplot(4,4,12)
outDat2(2, 2,:) = makeHFBImageHist(perm, perm2, 23:38);

subplot(4,4,15)
makeHFBTimeScatter(dat, perm2, 23:38, 'phase')

subplot(4,4,16)
makeIndexTimeScatter(dat, perm, perm2, 23:38)


export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' regions{reg} '_' phase '_3.jpg'], '-r300')



% % update significant frequencies
% pSum = sum(p<.05);
% pSumHFB = sum(p2<.05); 
% if strcmp(phase, 'sub')
%     sigFreqs(:,1,1,2) = sigFreqs(:,1,1,2) + pSum'; 
%     sigFreqs(:,1,2,2) = sigFreqs(:,1,2,2) + pSumHFB'; 
% else
%     sigFreqs(:,2,1,2) = sigFreqs(:,2,1,2) + pSum'; 
%     sigFreqs(:,2,2,2) = sigFreqs(:,2,2,2) + pSumHFB'; 
% 
% end

figure('visible', true, 'position', [0,0,600,1000])

subplot(4,2,1)
hold off
tim = dat.tim; 
tim(tim<-450 | tim>3000) = []; 
hitPhase = squeeze(dat.hits_p(:,:,5)); 
% missPhase = squeeze(dat.misses_p(:,:,5)); 
% %time point to align at will be determined by peak ITPC
% test = arrayfun(@(y) arrayfun(@(x) size(hitPhase,2)*...
%     abs(mean(exp(1i*squeeze(dat.hits_p(x,:,y)) ))).^2,...
%     1:size(hitPhase,1) ), 1:100, 'uniformoutput', false);
% test = cat(1, test{:}); 
test = squeeze(mean(perm.hitVals,1))';
imagesc(test)
set(gca, 'ydir', 'normal')
caxis([2, 5])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit phase reset'])

subplot(4,2,2)
hold off
timCut = [-500:25:500]; 
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
testHFB = squeeze(mean(perm2.hitVals,1))';
imagesc(testHFB)
set(gca, 'ydir', 'normal')
caxis([2, 5])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit phase reset'])


subplot(4,2,3)
hold off
% hitPhase = squeeze(dat.hits_p(:,:,5)); 
% missPhase = squeeze(dat.misses_p(:,:,5)); 
% %time point to align at will be determined by peak ITPC
% test2 = arrayfun(@(y) arrayfun(@(x) size(missPhase,2)*...
%     abs(mean(exp(1i*squeeze(dat.misses_p(x,:,y)) ))).^2,...
%     1:size(missPhase,1) ), 1:100, 'uniformoutput', false);
% test2 = cat(1, test2{:}); 
test2 = squeeze(mean(perm.missVals))';
imagesc(test2)
set(gca, 'ydir', 'normal')
caxis([2, 5])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' miss phase reset'])

subplot(4,2,4)
hold off
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
testHFB2 = squeeze(mean(perm2.missVals))';
imagesc(testHFB2)
set(gca, 'ydir', 'normal')
caxis([2, 5])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit phase reset'])

subplot(4,2,5)
imagesc(test - test2)
set(gca, 'ydir', 'normal')
addRedOutline(p, .05, 'red');
caxis([-3, 3])
colorbar
xticks([3:16:160])
xticklabels(tim([3:16:139]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit - miss phase reset'])

subplot(4,2,6)
hold off
imagesc(testHFB - testHFB2)
set(gca, 'ydir', 'normal')
addRedOutline(p2, .05, 'red');
caxis([-3, 3])
colorbar
xticks([1:10:41])
xticklabels(timCut([1:10:41]))
yticks([1:19:100])
yticklabels(round(dat.frex([1:19:100])))
title([regions{reg} ' ' phase ' hit - miss phase reset'])


test = mean(perm.tVals,2);
[~, timepoint] = max(test);
subplot(4,2,7)
hold off
plot(squeeze(mean(perm.hitVals(:,timepoint,:))), 'color', 'blue',...
    'linewidth', 2)
hold on 
plot(squeeze(mean(perm.missVals(:,timepoint,:))), 'color', 'red',...
    'linewidth', 2)
xticks([1:19:100])
xticklabels(round(dat.frex([1:19:100])))
title(['phase reset at peak significance ' num2str(tim(timepoint))])

subplot(4,2,8)
hold off
plot(squeeze(testHFB(:,21)), 'color', 'blue',...
    'linewidth', 2)
hold on 
plot(squeeze(testHFB2(:,21)), 'color', 'red',...
    'linewidth', 2)
xticks([1:19:100])
xticklabels(round(dat.frex([1:19:100])))
title('phase reseting at HFB peak')

% 
% 
% subplot(4,2,2)
% hold off
% tim =dat.tim; 
% peakIdx = arrayfun(@(x) find(x==tim), dat.hitLat);
% vals = arrayfun(@(x) dat.hits_p(peakIdx(x), x, peakToUse),...
%     [1:length(peakIdx)]); 
% histogram(vals, [-pi:pi/16:pi], 'normalization', 'probability')
% hold on 
% peakIdx = arrayfun(@(x) find(x==tim), dat.missLat);
% vals = arrayfun(@(x) dat.misses_p(peakIdx(x), x, peakToUse),...
%     [1:length(peakIdx)]); 
% histogram(vals, [-pi:pi/16:pi], 'normalization', 'probability')
% yline(1/32, 'linewidth', 2, 'linestyle', '--')
% title('phase of peak oscillation at HFB peak')
% ylim([.01, .05])
% 
% 
% subplot(4,2,3)
% 
% 
% 
% 
% subplot(4, 2, 4)
% 
% test2 = arrayfun(@(x) size(missPhase,2)*...
%     abs(mean(exp(1i*missPhase(x,:) ))).^2,...
%     1:size(missPhase,1) );
% tim = dat.tim; 
% plot(tim, test, 'color', 'blue',...
%     'linewidth', 2)
% hold on 
% plot(tim, test2, 'color', 'red',...
%     'linewidth', 2)
% xlim([-450, 2500])
% title("ITPC Rayleigh's z")
% ylim([0,10])
% 
% %get the latencies of the peak frequencies
% hit_lowLat = gausLat(squeeze(dat.hits(:,:,peakToUse)), tim, dat.hitRT);
% miss_lowLat = gausLat(squeeze(dat.misses(:,:,peakToUse)), tim, dat.missRT);
% subplot(4,2,4)
% hold off
% histogram(dat.hitLat - hit_lowLat, [-500:25:500],...
%     'normalization', 'probability')
% hold on 
% histogram(dat.missLat -miss_lowLat , [-500:25:500],...
%     'normalization', 'probability')
% title('relative timing of HFB peak - low freq peak')
% xlabel('HFB first               low first')
% xline(0, 'linewidth', 2)
% ylim([0,.15])
% 
% %plot latencies on their own
% subplot(4,2,[5,7])
% hold off
% [sortedLat,order] = sort(dat.hitLat); 
% lowPow = squeeze(dat.hits(:,:,peakToUse)); 
% imagesc(tim, [], lowPow(:,order)')
% caxis([-5,10])
% hold on 
% plot(sortedLat, [1:length(sortedLat)], 'color', 'red', 'linewidth',1)
% xline(0, 'color', 'green', 'linewidth', 2)
% title(['low freq pow ' phase ' hits'])
% xlim([-450, 2500])
% 
% 
% subplot(4,2,[6,8])
% hold off
% [sortedLat,order] = sort(dat.missLat); 
% lowPow = squeeze(dat.misses(:,:,peakToUse)); 
% imagesc(tim, [], lowPow(:,order)')
% caxis([-5,10])
% hold on 
% plot(sortedLat, [1:length(sortedLat)], 'color', 'red', 'linewidth',1)
% xline(0, 'color', 'green', 'linewidth', 2)
% title(['low freq pow ' phase ' misses'])
% xlim([-450, 2500])



export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' regions{reg} '_' phase '_2.jpg'], '-r300')










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





catch

disp(['failure on: ' regions{reg} ' ' phase])
end











end