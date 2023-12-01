%I want to make plots with hit/miss difference for regions both image
%locked and HFB locked


codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath(genpath([codePre 'mni2atlas']))
% addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\fieldtrip-20230118')
% ft_defaults
HFBdatFolder = [datPre 'HFB_singleTrial\out']; 
TFdatFolder = [datPre 'TF_singleTrial\out']; 
HFBdatFolder = dir(HFBdatFolder);
TFdatFolder = dir(TFdatFolder);
test = cellfun(@(x) length(x)>0, strfind({HFBdatFolder.name}, '.mat'));
HFBdatFolder = HFBdatFolder(test); 
test = cellfun(@(x) length(x)>0, strfind({TFdatFolder.name}, '.mat'));
TFdatFolder = TFdatFolder(test); 



for ii = 1: length(HFBdatFolder)
    fileBits = split(HFBdatFolder(ii).name, 'stat'); 
    HFBdatFolder(ii).stat = str2num(fileBits{2}(1));
    fileBits = split(fileBits{2}, '_'); 
    HFBdatFolder(ii).reg = fileBits{2}; 
    fileBits = split(fileBits{3}, '.'); 
    HFBdatFolder(ii).phase = fileBits{1}; 

end

for ii = 1: length(TFdatFolder)
    fileBits = split(TFdatFolder(ii).name, 'stat'); 
    TFdatFolder(ii).stat = str2num(fileBits{2}(1));
    fileBits = split(fileBits{2}, '_'); 
    TFdatFolder(ii).reg = fileBits{2}; 
    TFdatFolder(ii).freq = split(fileBits{4}, '.'); 
    TFdatFolder(ii).freq = TFdatFolder(ii).freq{1}; 
    TFdatFolder(ii).phase = fileBits{3}; 

end



regions = unique({HFBdatFolder.reg});
phases = {'sub', 'ret'}; 

%phase (tvals:1-2; pvals: 3-4) X time X regions;component
allResHFB = zeros(4, 41, length(regions)*3); 
allResImage = zeros(4, 139, length(regions)*3); 

for reg = 1:length(regions)

    curHFB = ismember({HFBdatFolder.reg}, regions{reg}); 
    curHFB = HFBdatFolder(curHFB); 
    curTF = ismember({TFdatFolder.reg}, regions{reg});
    curTF = TFdatFolder(curTF); 

%arrange subplots with 3 signals tall (HFB, TFlow, TFhigh)
%2 columns on the left are stat 1 encode/retrieve
%2 columns on the right are stat 0 encode/retrieve

figure
for phasei = 1:2
    HFBphase = curHFB(ismember({curHFB.phase}, phases{phasei})); 
    TFphase = curTF(ismember({curTF.phase}, phases{phasei})); 

    for stat = 0:1
        HFBstat = HFBphase([HFBphase.stat]==stat); 
        TFstat = TFphase([TFphase.stat]==stat); 
        statInfo = load([HFBstat.folder '/' HFBstat.name]).statInfo;
        if stat==0
            suf = 'HFB'; 
            mainTim = [-500:25:500]; 
        else
            suf = 'image';
            mainTim = statInfo.tim; 
            mainTim(mainTim<-450 | mainTim>3000) = []; 
        end

        

        %HFB peak
        subplot(9,4,[stat*2+1+phasei-1, stat*2+1+4+phasei-1])
        [hitVals, missVals] = getSingleTrialHitMiss(statInfo, stat);
        standardLinePlot(hitVals, missVals, mainTim); 
        if stat==1
            xlim([-500, 2500])
        end

        title([statInfo.reg ' ' statInfo.phase ' HFB peak difference'])

        %HFB peak t-value
        subplot(9,4,stat*2+1+8+phasei-1)
        tValplot(statInfo.(['tVals_' suf]), ...
            mainTim, statInfo.(['p_' suf]))
        if stat==1
            xlim([-500, 2500])
            allResImage(phasei, :, (reg-1)*3+1) = ...
                statInfo.(['tVals_' suf]);
            allResImage(phasei+2, :, (reg-1)*3+1) = ...
                statInfo.(['p_' suf]);
        else
            allResHFB(phasei, :, (reg-1)*3+1) = ...
                statInfo.(['tVals_' suf]);
            allResHFB(phasei+2, :, (reg-1)*3+1) = ...
                statInfo.(['p_' suf]);
        end
        title('t-value')

        %low frequency at HFB peak
        subplot(9,4,[stat*2+1+12+phasei-1, stat*2+1+16+phasei-1])
        statInfo = load([TFstat(1).folder '/stat' num2str(stat) ...
            '_' regions{reg} '_' TFstat(1).phase '_low.mat']).statInfo;
        [hitVals, missVals] = getSingleTrialHitMiss(statInfo, stat);
        standardLinePlot(hitVals, missVals, mainTim); 
        if stat==1
            xlim([-500, 2500])
        end
        title([statInfo.reg ' ' statInfo.phase ' low freq difference'])

        %low frequency at HFB peak t-value
        subplot(9,4,stat*2+1+20+phasei-1)
        tValplot(statInfo.(['tVals_' suf]), ...
            mainTim, statInfo.(['p_' suf]))
        
        if stat==1
            xlim([-500, 2500])
            allResImage(phasei, :, (reg-1)*3+2) = ...
                statInfo.(['tVals_' suf]);
            allResImage(phasei+2, :, (reg-1)*3+2) = ...
                statInfo.(['p_' suf]);
        else
            allResHFB(phasei, :, (reg-1)*3+2) = ...
                statInfo.(['tVals_' suf]);
            allResHFB(phasei+2, :, (reg-1)*3+2) = ...
                statInfo.(['p_' suf]);
        end
        title('t-value')

         %high frequency at HFB peak
        subplot(9,4,[stat*2+1+24+phasei-1, stat*2+1+28+phasei-1])
        statInfo = load([TFstat(1).folder '/stat' num2str(stat) ...
            '_' regions{reg} '_' TFstat(1).phase '_high.mat']).statInfo;
        [hitVals, missVals] = getSingleTrialHitMiss(statInfo, stat);
        standardLinePlot(hitVals, missVals, mainTim); 
        if stat==1
            xlim([-500, 2500])
        end
        title([statInfo.reg ' ' statInfo.phase ' high freq difference'])

        %high frequency at HFB peak t-value
        subplot(9,4,stat*2+1+32+phasei-1)
        tValplot(statInfo.(['tVals_' suf]), ...
            mainTim, statInfo.(['p_' suf]))
        if stat==1
            xlim([-500, 2500])
            allResImage(phasei, :, (reg-1)*3+3) = ...
                statInfo.(['tVals_' suf]);
            allResImage(phasei+2, :, (reg-1)*3+3) = ...
                statInfo.(['p_' suf]);
        else
            allResHFB(phasei, :, (reg-1)*3+3) = ...
                statInfo.(['tVals_' suf]);
            allResHFB(phasei+2, :, (reg-1)*3+3) = ...
                statInfo.(['p_' suf]);
        end
        title('t-value')

    end





end











end

climits = [-3, 3]; 
figure
subplot 221
imagesc(squeeze(allResHFB(1,:,:))')
addRedOutline(squeeze(allResHFB(3,:,:)), .1, 'green');
addRedOutline(squeeze(allResHFB(3,:,:)), .05, 'red');
caxis(climits)
for ii = [3.5:3:27]
    yline(ii, 'linewidth', 2', 'linestyle', '--')
end
xticks([1:10:41])
xticklabels([-500:250:500])
xline(21)
yticks([2:3:27])
yticklabels(regions)
colorbar
title("encoding HFB locked")

subplot 222
imagesc(squeeze(allResHFB(2,:,:))')
addRedOutline(squeeze(allResHFB(4,:,:)), .1, 'green');
addRedOutline(squeeze(allResHFB(4,:,:)), .05, 'red');
caxis(climits)
for ii = [3.5:3:27]
    yline(ii, 'linewidth', 2', 'linestyle', '--')
end
xticks([1:10:41])
xticklabels([-500:250:500])
xline(21)
yticks([2:3:27])
yticklabels(regions)
colorbar
title("retrieval HFB locked")

mainTim = statInfo.tim; 
mainTim(mainTim<-450 | mainTim>3000) = []; 

subplot 223
imagesc(squeeze(allResImage(1,:,:))')
addRedOutline(squeeze(allResImage(3,:,:)), .1, 'green');
addRedOutline(squeeze(allResImage(3,:,:)), .05, 'red');
caxis(climits)
for ii = [3.5:3:27]
    yline(ii, 'linewidth', 2', 'linestyle', '--')
end
xticks([1:18:140])
xticklabels(mainTim([1:18:140]))
xline(19)
yticks([2:3:27])
yticklabels(regions)
colorbar
title("encoding image locked")

subplot 224
imagesc(squeeze(allResImage(2,:,:))')
addRedOutline(squeeze(allResImage(4,:,:)), .1, 'green');
addRedOutline(squeeze(allResImage(4,:,:)), .05, 'red');
caxis(climits)
for ii = [3.5:3:27]
    yline(ii, 'linewidth', 2', 'linestyle', '--')
end
xticks([1:18:140])
xticklabels(mainTim([1:18:140]))
xline(19)
yticks([2:3:27])
yticklabels(regions)
colorbar
title("retrieval image locked")







