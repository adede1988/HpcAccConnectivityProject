%script to perform plotting of simultaneous pair recording lead lag
%information


%% set up the environment
codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';
addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse\export_fig_repo'])

savFolder = [datPre 'PPC/'];

chanFiles = load([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat']).chanFiles;
chanFiles(~[chanFiles.targPair]) = []; 
chanFiles(~[chanFiles.bothReact]) = []; 
regions = unique([chanFiles.chan1Reg]);


%% do an audit of the pair data files: 
% parfor ii = 1:length(chanFiles)
%     ii
%     try
%        
%         %check each file to see if it can be loaded
%         nameBits = split(chanFiles(ii).name, '_'); 
%         ch1 = split(nameBits{3}, '.mat');
%         ch1 = str2num(ch1{1});
%         nameBits2 = split(chanFiles(ii).name2, '_'); 
%         ch2 = split(nameBits2{3}, '.mat');
%         ch2 = str2num(ch2{1});  
%     
%         fn = [nameBits{2} '_' num2str(ch1) '_' num2str(ch2) '_'...
%             chanFiles(ii).chan1Reg{1} '_' chanFiles(ii).chan2Reg{1} '.mat'];
%         pairDat = load([savFolder fn]).pairDat; 
%         chanFiles(ii).error = 0; 
%         if ~isfield(pairDat, 'LL')
%             chanFiles(ii).error = 1; 
%         end
%          
% 
%     catch
% 
%         chanFiles(ii).error = 1; 
%     end
% 
% 
% end
% 
% chanFiles2 = load([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat']).chanFiles;
% chanFiles2(~[chanFiles2.targPair]) = []; 
% chanFiles2(~[chanFiles2.bothReact]) = []; 
% parfor ii = 1:length(chanFiles)
%     chanFiles2(ii).error = chanFiles(ii).error; 
% end
% chanFiles = chanFiles2; 
% save([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat'], 'chanFiles')


%% Get key information by region pair
%goal is to get this information out in a saved file that can be sent to
%Quest for stats processing

%Do for encoding and retrieval separately 

%get the key information out for all channel pairs, but leave empty if it's not
%the current target region pair

%1: HFB trials X time hit, Chan 1
%2: HFB trials X time miss, Chan 1
%3: HFB trials X time hit, Chan 2
%4: HFB trials X time miss, Chan 2
%5: TF trials X time LOW HIT, chan 1
%6: TF trials X time LOW MISS, chan 1
%7: TF trials X time HIGH HIT, chan 1
%8: TF trials X time HIGH MISS chan 1
%9: TF trials X time LOW HIT, chan 2
%10: TF trials X time LOW MISS, chan 2
%11: TF trials X time HIGH HIT, chan 2
%12: TF trials X time HIGH MISS chan 2
%13: LL trials X offset X time hit LOW 
%14: LL trials X offset X time miss LOW
%15: LL trials X offset X time hit HIGH 
%16: LL trials X offset X time miss HIGH
%17: trial RTs Hit
%18: trial RTs Miss
%19: region1
%20: region2
%21: time
%22: subID
%23: chani 1
%24: chani 2
phases = {'enc', 'ret'};
HFBmisses = {'subMiss', 'miss_on'}; 
HFBhits = {'subHit', 'hit_on'}; 
HFBtim = {'enc', 'on'}; 
LLphase = {'sub', 'ret'}; 
for phase_i = 1:2
for reg1 = 1:5
    for reg2 = 1:5
disp(['working on ' num2str(phase_i) ' ' regions{reg1} ' ' regions{reg2}])

allDat = cell(length(chanFiles), 24); 
% reg1 = 4; 
% reg2 = 3; 

% allHits = zeros(10000, 61, 141); 
% allMisses = allHits; 
% allHitLat = zeros(10000,1); 
% allMissLat = allHitLat; 
% shi = 1; 
% smi = 1; 
parfor ii = 1:length(chanFiles)
 
    if strcmp(chanFiles(ii).chan1Reg{1}, regions{reg1}) && ...
        strcmp(chanFiles(ii).chan2Reg{1}, regions{reg2})
        %select only the target region pair
    nameBits = split(chanFiles(ii).name, '_'); 
    ch1 = split(nameBits{3}, '.mat');
    ch1 = str2num(ch1{1});
    nameBits2 = split(chanFiles(ii).name2, '_'); 
    ch2 = split(nameBits2{3}, '.mat');
    ch2 = str2num(ch2{1});  

    fn = [nameBits{2} '_' num2str(ch1) '_' num2str(ch2) '_'...
        chanFiles(ii).chan1Reg{1} '_' chanFiles(ii).chan2Reg{1} '.mat'];
    pairDat = load([savFolder fn]).pairDat; 
    chanFiles(ii).error = 0; 
    
    chanDat1 = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\'...
        chanFiles(ii).name]).chanDat; 
    chanDat2 = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\'...
        chanFiles(ii).name2]).chanDat; 

    %save the outputs
    out = cell(24,1); 
    %1: HFB trials X time hit, Chan 1
    out{1} = chanDat1.HFB.(HFBhits{phase_i});  
    %2: HFB trials X time miss, Chan 1
    out{2} = chanDat1.HFB.(HFBmisses{phase_i}); 
    %3: HFB trials X time hit, Chan 2
    out{3} = chanDat2.HFB.(HFBhits{phase_i});  
    %4: HFB trials X time miss, Chan 2
    out{4} = chanDat2.HFB.(HFBmisses{phase_i}); 
    %5: TF trials X time LOW HIT, chan 1
    out{5} = squeeze(mean(pairDat.([HFBhits{phase_i} 'TF1'])(:,:,1:4), 3)); 
    %6: TF trials X time LOW MISS, chan 1
    out{6} = squeeze(mean(pairDat.([HFBmisses{phase_i} 'TF1'])(:,:,1:4), 3)); 
    %7: TF trials X time HIGH HIT, chan 1
    out{7} = squeeze(mean(pairDat.([HFBhits{phase_i} 'TF1'])(:,:,5:8), 3)); 
    %8: 8F trials X time HIGH MISS chan 1
    out{8} = squeeze(mean(pairDat.([HFBmisses{phase_i} 'TF1'])(:,:,5:8), 3)); 
    %9: TF trials X time LOW HIT, chan 2
    out{9} = squeeze(mean(pairDat.([HFBhits{phase_i} 'TF2'])(:,:,1:4), 3)); 
    %10: TF trials X time LOW MISS, chan 2
    out{10} = squeeze(mean(pairDat.([HFBmisses{phase_i} 'TF2'])(:,:,1:4), 3)); 
    %11: TF trials X time HIGH HIT, chan 2
    out{11} = squeeze(mean(pairDat.([HFBhits{phase_i} 'TF2'])(:,:,5:8), 3)); 
    %12: TF trials X time HIGH MISS chan 2
    out{12} = squeeze(mean(pairDat.([HFBmisses{phase_i} 'TF2'])(:,:,5:8), 3)); 
    %13: LL trials X offset X time hit LOW 
    out{13} = pairDat.LL.([LLphase{phase_i} 'Hit_low']); 
    %14: LL trials X offset X time miss LOW
    out{14} = pairDat.LL.([LLphase{phase_i} 'Miss_low']); 
    %15: LL trials X offset X time hit HIGH 
    out{15} = pairDat.LL.([LLphase{phase_i} 'Hit_high']); 
    %16: LL trials X offset X time miss HIGH
    out{16} = pairDat.LL.([LLphase{phase_i} 'Miss_high']); 
    if phase_i == 1
    %17: trial RTs Hit
    out{17} = chanDat1.encInfo(chanDat1.use & chanDat1.hits, 4); 
    %18: trial RTs Miss
    out{18} = chanDat1.encInfo(chanDat1.use & chanDat1.misses, 4); 
    else
    %17: trial RTs Hit
    out{17} = chanDat1.retInfo(chanDat1.retInfo(:,1)==1, 3); 
    %18: trial RTs Miss
    out{18} = chanDat1.retInfo(chanDat1.retInfo(:,1)==2, 3); 
    end
    %19: region1
    out{19} = regions{reg1}; 
    %20: region2
    out{20} = regions{reg2}; 
    %21: time (HFB and TF)
    out{21} = chanDat1.HFB.([HFBtim{phase_i} 'MulTim']);  
    %22: subID
    out{22} = chanDat1.subID; 
    %23: chani 1
    out{23} = ch1; 
    %24: chani 2
    out{24} = ch2; 
    

    allDat(ii,:) = out; 
        

    end



end

test = cellfun(@(x) isempty(x), allDat(:,1)); 
allDat(test,:) = []; 

save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\ThetaLL\'...
    regions{reg1} '_' regions{reg2} '_' phases{phase_i} '.mat'], 'allDat', '-v7.3')

disp('saved summary dat')
freqNames = {'LOW', 'HIGH'};
for freqi = 1:2

%get the key hit info
HFB = cat(2, allDat{:,1} ); 
RT = cat(1, allDat{:,17});
Lat = gausLat(HFB, allDat{1,21}, RT); 
LL = cat(1, allDat{:,13+(freqi-1)*2});
LL(Lat==-1, : ,: ) = []; 
Lat(Lat==-1) = []; 
LLtim = [-500:25:max(allDat{1,21})-500];

figure('visible', false, 'position', [0,0,1500, 2500])

subplot(5,2,5)
imagesc(LLtim, [-150:5:150], squeeze(mean(LL, 1)))
allMax = max(max(squeeze(mean(LL, 1))));
allMin = min(min(squeeze(mean(LL, 1))));
colorbar
yline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'green')
xline(0, 'linewidth', 2, 'color', 'k')
title(['hit image locked LL ' freqNames{freqi}])

subplot(5,2,6)
Lat(Lat>max(LLtim)-500) =max(LLtim)-500; 
 test = arrayfun(@(x) squeeze(LL(x,:, ...
        find(LLtim>=Lat(x),1)-20 : ...
        find(LLtim>=Lat(x),1)+20)), ...
        1:length(Lat), 'uniformoutput', false);
    test2 = reshape(cell2mat(test), 61, 41, []); 
imagesc([-500:250:500], [-150:5:150], squeeze(mean(test2, 3)))
allMax = [allMax, max(max(squeeze(mean(test2, 3))))];
allMin = [allMin, min(min(squeeze(mean(test2, 3))))];
yline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'green')
xline(0, 'linewidth', 2, 'color', 'k')
lagFinder = squeeze(mean(test2, 3)); 
[~, lagFinder] = max(sum(lagFinder,2));
lagVals = [-150:5:150]; 
yline(lagVals(lagFinder), 'color', 'red', 'linestyle', '--')
title(['hit HFB locked LL ' freqNames{freqi}])
colorbar

subplot(5,2,[1,3])
[sortedLat, order] = sort(Lat); 
plotMat = squeeze(LL(order, lagFinder, :));
test = gausswin(31);
test = test ./ sum(test); 
plotMat = arrayfun(@(x) conv(plotMat(:,x), test, 'same'), ...
    [1:size(plotMat,2)], 'uniformoutput', false);
plotMat = cat(2, plotMat{:});  
test = gausswin(5); 
test = test ./ sum(test); 
plotMat = arrayfun(@(x) conv(plotMat(x,:), test, 'same'), ...
    [1:size(plotMat,1)], 'uniformoutput', false);
plotMat = cat(1, plotMat{:}); 
imagesc(LLtim, [], plotMat)
caxis([-0, .2])
colorbar
hold on 
plot(Lat(order), [1:length(sortedLat)], 'color', 'red')
test = gausLat(squeeze(LL(:, lagFinder, :))', [LLtim, max(allDat{1,21})], RT);
tmpLat = Lat; 
tmpLat(test==-1) = []; 
test(test==-1) = []; 
title([phases{phase_i} ' hits ' regions{reg1} ...
    ' to ' regions{reg2} ' at: ' num2str(lagVals(lagFinder))])

subplot(5,2,10)
histogram(tmpLat - test, [-2500:75:2500], 'normalization', 'probability')%,...
%     'facecolor', 'blue', 'edgecolor', 'blue', 'facealpha', .5)
xline(0, 'color', 'k', 'linewidth', 2, 'linestyle', '--')
early = sum(tmpLat-test <0) / length(tmpLat);
text(-2500, .05, ['LL before HFB: ' num2str(round(early*100)) '%'], 'color', 'blue')
% hold on 
% histogram(randsample(tmpLat, length(tmpLat), false) - test, ...
%     [-2500:75:2500], 'normalization', 'probability', 'facecolor', 'blue', ...
%     'edgecolor', 'blue','facealpha', .05, 'edgealph', .9)

subplot(5,2,9)
plot(arrayfun(@(x) prctile(test, x), [1:100]), [1:100], 'color', 'blue', 'linewidth', 2)

set(gca, 'YDir','reverse')



%get the key miss info
HFB = cat(2, allDat{:,2} ); 
RT = cat(1, allDat{:,18});
Lat = gausLat(HFB, allDat{1,21}, RT); 
LL = cat(1, allDat{:,14+(freqi-1)*2});
LL(Lat==-1, : ,: ) = []; 
Lat(Lat==-1) = []; 

subplot(5,2,7)
imagesc(LLtim, [-150:5:150], squeeze(mean(LL, 1)))
allMax = [allMax, max(max(squeeze(mean(LL, 1))))];
allMin = [allMin, min(min(squeeze(mean(LL, 1))))];
yline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'green')
xline(0, 'linewidth', 2, 'color', 'k')
title('miss image locked LL')
colorbar

subplot(5,2,8)
Lat(Lat>max(LLtim)-500) =max(LLtim)-500; 
 test = arrayfun(@(x) squeeze(LL(x,:, ...
        find(LLtim>=Lat(x),1)-20 : ...
        find(LLtim>=Lat(x),1)+20)), ...
        1:length(Lat), 'uniformoutput', false);
    test2 = reshape(cell2mat(test), 61, 41, []); 
imagesc([-500:250:500], [-150:5:150], squeeze(mean(test2, 3)))
allMax = [allMax, max(max(squeeze(mean(test2, 3))))];
allMin = [allMin, min(min(squeeze(mean(test2, 3))))];
yline(0, 'linewidth', 2, 'linestyle', '--', 'color', 'green')
xline(0, 'linewidth', 2, 'color', 'k')
% lagFinder = squeeze(mean(test2, 3)); 
% [~, lagFinder] = max(sum(lagFinder,2));
lagVals = [-150:5:150]; 
yline(lagVals(lagFinder), 'color', 'red', 'linestyle', '--')
title('miss HFB locked LL')
colorbar


subplot(5,2,[2,4])
[sortedLat, order] = sort(Lat); 
plotMat = squeeze(LL(order, lagFinder, :));
test = gausswin(31);
test = test ./ sum(test); 
plotMat = arrayfun(@(x) conv(plotMat(:,x), test, 'same'), ...
    [1:size(plotMat,2)], 'uniformoutput', false);
plotMat = cat(2, plotMat{:});  
test = gausswin(5); 
test = test ./ sum(test); 
plotMat = arrayfun(@(x) conv(plotMat(x,:), test, 'same'), ...
    [1:size(plotMat,1)], 'uniformoutput', false);
plotMat = cat(1, plotMat{:}); 
imagesc(LLtim, [], plotMat)
caxis([-0, .2])
hold on 
plot(Lat(order), [1:length(sortedLat)], 'color', 'red')
test = gausLat(squeeze(LL(:, lagFinder, :))', [LLtim, max(allDat{1,21})], RT);
tmpLat = Lat; 
tmpLat(test==-1) = []; 
test(test==-1) = []; 
title([phases{phase_i} ' misses ' regions{reg1} ...
    ' to ' regions{reg2} ' at: ' num2str(lagVals(lagFinder))])
colorbar


subplot(5,2,10)
hold on 
histogram(tmpLat - test, [-2500:75:2500], 'normalization', 'probability')%,...
%     'facecolor', 'red', 'edgecolor', 'red', 'facealpha', .5)
xline(0, 'color', 'k', 'linewidth', 2, 'linestyle', '--')
ylim([0,.1])
early = sum(tmpLat-test <0) / length(tmpLat);
text(-2500, .04, ['LL before HFB: ' num2str(round(early*100)) '%'], 'color', 'red')
title('timing of LL peak relative to HFB peak')
xlabel('HFB peak - LL peak (ms)')

subplot(5,2,9)
hold on 
plot(arrayfun(@(x) prctile(test, x), [1:100]), [1:100],...
    'color', 'red', 'linewidth', 2)
title(['LL peak at offset: ' num2str(lagVals(lagFinder)) ' ms'])
xlabel('time after image (ms)')

subplot(5,2,5)
clim([min(allMin)+.01 max(allMax) - .01])
subplot(5,2,6)
clim([min(allMin)+.01 max(allMax) - .01])
subplot(5,2,7)
clim([min(allMin)+.01 max(allMax) - .01])
subplot(5,2,8)
clim([min(allMin)+.01 max(allMax) - .01])

export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\Theta_LL\'...
    regions{reg1} '_' regions{reg2} '_' phases{phase_i} '_' freqNames{freqi} '.jpg'], '-r300')
end %frequency loop

    end %inner region loop
end %outer region loop
end %encoding/retrieval loop


% 
% enctim = [-1000:3500];
% retOtim = [-1000:3000];
% 
% 
%     leadLagEncTim = enctim(501:25:end-500);
%     leadLagRetTim = retOtim(501:25:end-500); 
% 
% 
%     figure
% 
%     tmp = allHitLat;
%     tmp2 = allHits(tmp>-1,:,:); 
%     tmp(tmp==-1) = []; 
% %     tmp(tmp>=max(pairDat.encTim)-500) = max(pairDat.encTim)-500;
%     test = arrayfun(@(x) squeeze(tmp2(x,:, ...
%         find(leadLagRetTim>=tmp(x),1)-20 : ...
%         find(leadLagRetTim>=tmp(x),1)+20)), ...
%         1:length(tmp), 'uniformoutput', false);
%     test2 = reshape(cell2mat(test), 61, 41, []); 
% 
% 
% 
% 
% 
% subplot(4,2,5)
% imagesc(leadLagEncTim, -150:5:150, squeeze(mean(tmp2)))
% colorbar
% % caxis([.01, .08])
% yline(0)
% title('mean sub Hit')
% 
% subplot(4,2,6)
% imagesc(-500:25:500, -150:5:150, squeeze(mean(test2, 3)))
% colorbar
% % caxis([.01, .08])
% yline(0)
% xline(0)
% title('sub Hit HFB aligned')
% 
% subplot(4,2, [1,3])
%     [sortedLat, order] = sort(tmp); 
%     LLchoice = squeeze(mean(test2, 3)); 
%     [~, LLchoice] = max(sum(LLchoice, 2)); 
%     imagesc(leadLagEncTim, [], squeeze(tmp2(order, LLchoice, :)))
%     hold on 
%     plot(sortedLat, 1:length(sortedLat), 'color', 'red')
% %     caxis([-.1, .7])
%     title([regions{reg1} ' to ' regions{reg2} ' hit LL' ])
% 
% 
%     tmp = allMissLat;
%     tmp2 = allMisses(tmp>-1,:,:); 
%     tmp(tmp==-1) = []; 
% %     tmp(tmp>=max(pairDat.encTim)-500) = max(pairDat.encTim)-500;
%     test = arrayfun(@(x) squeeze(tmp2(x,:, ...
%         find(leadLagRetTim>=tmp(x),1)-20 : ...
%         find(leadLagRetTim>=tmp(x),1)+20)), ...
%         1:length(tmp), 'uniformoutput', false);
%     test2 = reshape(cell2mat(test), 61, 41, []); 
% 
% 
% 
% 
% 
% subplot(4,2,7)
% imagesc(leadLagEncTim, -150:5:150, squeeze(mean(tmp2)))
% colorbar
% % caxis([.01, .08])
% yline(0)
% title('mean sub Miss')
% 
% subplot(4,2,8)
% imagesc(-500:25:500, -150:5:150, squeeze(mean(test2, 3)))
% colorbar
% % caxis([.01, .08])
% yline(0)
% xline(0)
% title('sub Miss HFB aligned')
% 
% subplot(4,2, [2,4])
%     [sortedLat, order] = sort(tmp); 
% %     LLchoice = squeeze(mean(test2, 3)); 
% %     [~, LLchoice] = max(sum(LLchoice, 2)); 
%     imagesc(leadLagEncTim, [], squeeze(tmp2(order, LLchoice, :)))
%     hold on 
%     plot(sortedLat, 1:length(sortedLat), 'color', 'red')
% %     caxis([-.1, .7])
%     title([regions{reg1} ' to ' regions{reg2} ' miss LL' ])
% 
% 
%      tmp = allHitLat;
%     tmp2 = allHits(tmp>-1,:,:); 
%     tmp(tmp==-1) = []; 
% %     tmp(tmp>=max(pairDat.encTim)-500) = max(pairDat.encTim)-500;
%     test = arrayfun(@(x) squeeze(tmp2(x,:, ...
%         find(leadLagRetTim>=tmp(x),1)-20 : ...
%         find(leadLagRetTim>=tmp(x),1)+20)), ...
%         1:length(tmp), 'uniformoutput', false);
%     test2 = reshape(cell2mat(test), 61, 41, []); 
% 
% 
% 
