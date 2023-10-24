

% LL stats plotting

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_LL_KEY_STATS';

statFiles = dir(path); 

statFiles(1:2) = []; 
LLdat = load([statFiles(1).folder '/' statFiles(1).name]).LLdat; 

allSigEnc = zeros(length(LLdat.encTim), 11); 
allSigRet = zeros(length(LLdat.retTim), 11); 
regNames = cell(11,1); 

%hit/miss X region X early/late X encode/retrieve
%early = 0-1000ms; late = 1000-2000ms
allMeanEarlyLate = zeros(2, 11, 2, 2); 


parfor ii = 1:length(statFiles)
ii
try
   
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    LLplotHelper(LLdat, ii)
  
catch
end

end


%% figure of all encoding latencies 
ROInames = {LLdat.aggTargs.ROI};
figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
for ii = 1:length(statFiles)
    ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    roiLabs = {LLdat.aggTargs{LLdat.reg1}, LLdat.aggTargs{LLdat.reg2}};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))

    curDat = [LLdat.regResLat(1,:), LLdat.regResLat(2,:)];
    hmSort = [ones(length(LLdat.regResLat(1,:)),1); zeros(length(LLdat.regResLat(1,:)),1)];
    modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
                'VariableNames', {'latency', 'hitMiss', 'sub'}); 
    lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
    tval = [0,0];
    tval(1) = lme.Coefficients(2,4); 
    tval(2) = lme.Coefficients(2,6); 
    
    tval(2) = round(tval(2),2);

    if tval(2)<.2
        tval(2) = tval(2)*5; 
    else 
        tval(2) = 1; 
    end
    [~, coli] = min(abs(tval(2) - findMap));
    
  
    
    hold off
    fill([0,2000,2000,0], [0,0,.8,.8], colormap(coli,:), 'FaceAlpha', .8)
    xlim([0,2000])
%     ylim([0, .8])
    hold on
    h = histogram(LLdat.regResLat(1,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(LLdat.regResLat(2,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(LLdat.regResLat(1,:), [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(LLdat.regResLat(2,:), [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


   
%     title(num2str(round(tval(2),2)))
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
      roiLabs = {LLdat.aggTargs{LLdat.reg1}, LLdat.aggTargs{LLdat.reg2}};
%     title(['encoding latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

%     subplot 122
%     hold off
%     histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability')
%     hold on
%     histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability')
%     
%     curDat = [LLdat.regRes2Lat(1,:), LLdat.regRes2Lat(2,:)];
%     hmSort = [ones(length(LLdat.regResLat(1,:)),1); zeros(length(LLdat.regResLat(1,:)),1)];
% 
% 
%     modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
%                 'VariableNames', {'latency', 'hitMiss', 'sub'}); 
%     lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
%     tval = [0,0];
%     tval(1) = lme.Coefficients(2,4); 
%     tval(2) = lme.Coefficients(2,6); 
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
%       
%     title(['retrieval latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

    
  

end

export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\LL_finalizedFigs\' 'LL_latency_allEncode' '.jpg'],''), '-r300')


figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
for ii = 1:length(statFiles)
    ii
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    roiLabs = {LLdat.aggTargs{LLdat.reg1}, LLdat.aggTargs{LLdat.reg2}};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))

    curDat = [LLdat.regRes2Lat(1,:), LLdat.regRes2Lat(2,:)];
    hmSort = [ones(length(LLdat.regRes2Lat(1,:)),1); zeros(length(LLdat.regRes2Lat(1,:)),1)];
    modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
                'VariableNames', {'latency', 'hitMiss', 'sub'}); 
    lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
    tval = [0,0];
    tval(1) = lme.Coefficients(2,4); 
    tval(2) = lme.Coefficients(2,6); 
    
    tval(2) = round(tval(2),2);

    if tval(2)<.2
        tval(2) = tval(2)*5; 
    else 
        tval(2) = 1; 
    end
    [~, coli] = min(abs(tval(2) - findMap));
    
  
    
    hold off
    fill([0,2000,2000,0], [0,0,.8,.8], colormap(coli,:), 'FaceAlpha', .8)
    xlim([0,2000])
%     ylim([0, .8])
    hold on
    h = histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


   
%     title(num2str(round(tval(2),2)))
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
      roiLabs = {LLdat.aggTargs{LLdat.reg1}, LLdat.aggTargs{LLdat.reg2}};
%     title(['encoding latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

%     subplot 122
%     hold off
%     histogram(LLdat.regRes2Lat(1,:), [0:100:2000], 'normalization', 'probability')
%     hold on
%     histogram(LLdat.regRes2Lat(2,:), [0:100:2000], 'normalization', 'probability')
%     
%     curDat = [LLdat.regRes2Lat(1,:), LLdat.regRes2Lat(2,:)];
%     hmSort = [ones(length(LLdat.regResLat(1,:)),1); zeros(length(LLdat.regResLat(1,:)),1)];
% 
% 
%     modDat = table(curDat', hmSort, [LLdat.chani; LLdat.chani], ...
%                 'VariableNames', {'latency', 'hitMiss', 'sub'}); 
%     lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
%     tval = [0,0];
%     tval(1) = lme.Coefficients(2,4); 
%     tval(2) = lme.Coefficients(2,6); 
%     text(0, .7, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, .6, [' p-val: '  num2str(round(tval(2),2))])
%     ylabel('time (ms)')
%       
%     title(['retrieval latency ' roiLabs{1} ' to ' roiLabs{2} ])
%     ylim([0,.8])

    
  

end

export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\LL_finalizedFigs\' 'LL_latency_allRetrieve' '.jpg'],''), '-r300')

%% count pairwise electrodes
pairwiseTrodes = zeros(11); 
pairwiseSubs = zeros(11); 
for ii = 1:length(statFiles)
    LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
    pairwiseTrodes(LLdat.reg1, LLdat.reg2) = LLdat.n_pair; 
    pairwiseSubs(LLdat.reg1, LLdat.reg2) = LLdat.n_sub; 
end





%% looking at all latencies across all regions


%enc/ret X hit/miss X pair
totalLatency = zeros(2, 2, 20098); 
ti = 1; 
for ii = 1:length(statFiles)
ii
 LLdat = load([statFiles(ii).folder '/' statFiles(ii).name]).LLdat; 
 
 N = LLdat.n_pair; 

 totalLatency(1,1,ti:ti+N-1) = LLdat.regResLat(1,:) ./ LLdat.subhitRT'; 
 totalLatency(1,2,ti:ti+N-1) = LLdat.regResLat(2,:)./ LLdat.submissRT'; 
 totalLatency(2,1,ti:ti+N-1) = LLdat.regRes2Lat(1,:)./ LLdat.rethitRT'; 
 totalLatency(2,2,ti:ti+N-1) = LLdat.regRes2Lat(2,:)./ LLdat.retmissRT';

ti = ti+N;
end



figure
subplot 121
histogram(totalLatency(1,1,:), [0:.03:1], 'normalization', 'probability')
hold on 
histogram(totalLatency(1,2,:), [0:.03:1], 'normalization', 'probability')
title('encoding latency across all ROI pairs')
subplot 122
histogram(totalLatency(1,1,:) - totalLatency(1,2,:), [-1:.03:1], 'normalization', 'probability')
title('within electrode pair hit - miss difference')
xline(0)


figure
subplot 121
histogram(totalLatency(2,1,:), [0:.03:1], 'normalization', 'probability')
hold on 
histogram(totalLatency(2,2,:), [0:.03:1], 'normalization', 'probability')
title('retrieval latency across all ROI pairs')
subplot 122
histogram(totalLatency(2,1,:) - totalLatency(2,2,:), [-1:.03:1], 'normalization', 'probability')
title('within electrode pair hit - miss difference')
xline(0)



scatter(squeeze(totalLatency(2,1,:)), squeeze(totalLatency(2,2,:)),200, 'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.02,'MarkerEdgeAlpha',.02)
hold on 
plot([200,2000], [200,2000])

