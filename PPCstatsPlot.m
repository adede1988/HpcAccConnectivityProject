
%PPC stats plot

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\PCC_KEY_STATS_HM';
path2 = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\PCC_KEY_STATS_HM2';
statFiles = dir(path); 
statFiles2 = dir(path2); 

statFiles(1:2) = []; 
statFiles2(1:2) = []; 

% PPC = load([statFiles(1).folder '/' statFiles(1).name]).connectionDat; 
% 
% allSigEnc = zeros(length(PPC.encTim), 11); 
% allSigRet = zeros(length(PPC.retTim), 11); 
% regNames = cell(11,1); 

%low/high, hit/miss, ROI, ROI, timeWindow (baseline, early, late)
allConnections = zeros(2,2,11,11, 3); 
allConnections2 = zeros(2,2,11,11, 3); 

%low/high, hit/miss, connection, time
rawAllConnections = zeros(2, 2, 20000, 141); 
rawAllConnections2 = zeros(2, 2, 20000, 121); %retrieval

%low/high, hit/miss, connection, sumInfo
%sumInfo: 
%1: subIDnum
%2: chanIDnum 1
%3: chanIDnum 2
%4: ROIi1
%5: ROIi2

rawALLconnectionsSum = zeros(20000, 5);
ci = 1;
 

for ii = 1:length(statFiles)
    try
    ii
    PPC = load([statFiles(ii).folder '/' statFiles(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    roi_idx = [PPC.reg1i, PPC.reg2i]; 
  
    f = figure('visible', true);
    f.Position = [0 0 600 1000];

    hitVals = squeeze(PPC.lowBand(PPC.hmSort,:) );
    missVals = squeeze(PPC.lowBand(~PPC.hmSort,:) );
    tim = PPC.tim; 

    nChan = size(hitVals,1); 
    if roi_idx(1) >= roi_idx(2)
    rawAllConnections(1,1,ci:ci+nChan-1, :) = hitVals; 
    rawAllConnections(1,2,ci:ci+nChan-1, :) = missVals; 
    %1: subIDnum
    subnums =cellfun(@(y) find(cellfun(@(x) strcmp(x, y), PPC.uniqueSubs)), PPC.allSubs);
    rawAllConnectionsSum(ci:ci+nChan-1,1) = subnums;
    %2-3: chanIDnum
    chani = cellfun(@(x) split(x, "_"), PPC.chani, 'uniformoutput', false);
    chani1 = cellfun(@(x) str2num(x{1}), chani);
    chani2 = cellfun(@(x) str2num(x{2}), chani);
    rawAllConnectionsSum(ci:ci+nChan-1,2) = chani1;
    rawAllConnectionsSum(ci:ci+nChan-1,3) = chani2;
    %4: ROIi1
    rawAllConnectionsSum(ci:ci+nChan-1,4) = ones(nChan,1)*roi_idx(1);
    %5: ROIi2
    rawAllConnectionsSum(ci:ci+nChan-1,5) = ones(nChan,1)*roi_idx(2);
    end

    allConnections(1,1,roi_idx(1), roi_idx(2), 1) = mean(hitVals(:, tim<0), 'all');
    allConnections(1,1,roi_idx(1), roi_idx(2), 2) = mean(hitVals(:, tim>=0 & tim<=1000), 'all');
    allConnections(1,1,roi_idx(1), roi_idx(2), 3) = mean(hitVals(:, tim>=1000 & tim<=2000), 'all'); 
    allConnections(1,2,roi_idx(1), roi_idx(2), 1) = mean(missVals(:, tim<0), 'all');
    allConnections(1,2,roi_idx(1), roi_idx(2), 2) = mean(missVals(:, tim>=0 & tim<=1000), 'all');
    allConnections(1,2,roi_idx(1), roi_idx(2), 3) = mean(missVals(:, tim>=1000 & tim<=2000), 'all'); 



    subplot 421
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    plot(tim, mean(missVals,1), 'linewidth', 3, 'color', 'red')
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
    sdMiss = std(missVals, [], 1) ./ sqrt(PPC.N); 
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
    x = [tim, flip(tim)];
     y = [mean(missVals,1) - sdMiss, flip(mean(missVals,1)) + flip(sdMiss)];  
    fill(flip(x), flip(y), 'red', 'FaceAlpha', .2)


    scatter(tim(PPC.lowp_sub<.05), PPC.lowp_sub(PPC.lowp_sub<.05)*0, 30,  'k', 'filled')
    ylabel("3Hz connectivity")
    title(['PPC ' roiLabs{1} ' to ' roiLabs{2} ' encoding'])



    subplot 423
    hold off
    plot(tim, PPC.lowtVals, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.low975, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.low025, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value hit/miss")





    hitVals = squeeze(PPC.highBand(PPC.hmSort,:) );
    missVals = squeeze(PPC.highBand(~PPC.hmSort,:) );
    tim = PPC.tim; 

    if roi_idx(1) >= roi_idx(2)
    rawAllConnections(2,1,ci:ci+nChan-1, :) = hitVals; 
    rawAllConnections(2,2,ci:ci+nChan-1, :) = missVals; 
    end


    allConnections(2,1,roi_idx(1), roi_idx(2), 1) = mean(hitVals(:, tim<0), 'all');
    allConnections(2,1,roi_idx(1), roi_idx(2), 2) = mean(hitVals(:, tim>=0 & tim<=1000), 'all');
    allConnections(2,1,roi_idx(1), roi_idx(2), 3) = mean(hitVals(:, tim>=1000 & tim<=2000), 'all'); 
    allConnections(2,2,roi_idx(1), roi_idx(2), 1) = mean(missVals(:, tim<0), 'all');
    allConnections(2,2,roi_idx(1), roi_idx(2), 2) = mean(missVals(:, tim>=0 & tim<=1000), 'all');
    allConnections(2,2,roi_idx(1), roi_idx(2), 3) = mean(missVals(:, tim>=1000 & tim<=2000), 'all'); 


    subplot 425
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    plot(tim, mean(missVals,1), 'linewidth', 3, 'color', 'red')
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
    sdMiss = std(missVals, [], 1) ./ sqrt(PPC.N); 
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
    x = [tim, flip(tim)];
     y = [mean(missVals,1) - sdMiss, flip(mean(missVals,1)) + flip(sdMiss)];  
    fill(flip(x), flip(y), 'red', 'FaceAlpha', .2)


    scatter(tim(PPC.highp_sub<.05), PPC.highp_sub(PPC.highp_sub<.05)*0, 30,  'k', 'filled')
    ylabel("8Hz connectivity")
    title(['PPC ' roiLabs{1} ' to ' roiLabs{2} ' encoding'])


    subplot 427
    hold off
    plot(tim, PPC.hightVals, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.high975, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.high025, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value hit/miss")


    %% retrieval 
    PPC = load([statFiles2(ii).folder '/' statFiles2(ii).name]).connectionDat; 
    hitVals = squeeze(PPC.lowBand2(PPC.hmSort,:) );
    missVals = squeeze(PPC.lowBand2(~PPC.hmSort,:) );
    tim = PPC.tim2; 

    if roi_idx(1) >= roi_idx(2)
    rawAllConnections2(1,1,ci:ci+nChan-1, :) = hitVals; 
    rawAllConnections2(1,2,ci:ci+nChan-1, :) = missVals; 
    end

    allConnections2(1,1,roi_idx(1), roi_idx(2), 1) = mean(hitVals(:, tim<0), 'all');
    allConnections2(1,1,roi_idx(1), roi_idx(2), 2) = mean(hitVals(:, tim>=0 & tim<=1000), 'all');
    allConnections2(1,1,roi_idx(1), roi_idx(2), 3) = mean(hitVals(:, tim>=1000 & tim<=2000), 'all'); 
    allConnections2(1,2,roi_idx(1), roi_idx(2), 1) = mean(missVals(:, tim<0), 'all');
    allConnections2(1,2,roi_idx(1), roi_idx(2), 2) = mean(missVals(:, tim>=0 & tim<=1000), 'all');
    allConnections2(1,2,roi_idx(1), roi_idx(2), 3) = mean(missVals(:, tim>=1000 & tim<=2000), 'all'); 


    subplot 422
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    plot(tim, mean(missVals,1), 'linewidth', 3, 'color', 'red')
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
    sdMiss = std(missVals, [], 1) ./ sqrt(PPC.N); 
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
    x = [tim, flip(tim)];
     y = [mean(missVals,1) - sdMiss, flip(mean(missVals,1)) + flip(sdMiss)];  
    fill(flip(x), flip(y), 'red', 'FaceAlpha', .2)


    scatter(tim(PPC.lowp_ret<.05), PPC.lowp_ret(PPC.lowp_ret<.05)*0, 30,  'k', 'filled')
    ylabel("3Hz connectivity")
    title(['PPC ' roiLabs{1} ' to ' roiLabs{2} ' retrieval'])



    subplot 424
    hold off
    plot(tim, PPC.lowtVals_ret, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.low975_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.low025_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value hit/miss")





    hitVals = squeeze(PPC.highBand2(PPC.hmSort,:) );
    missVals = squeeze(PPC.highBand2(~PPC.hmSort,:) );
%   
    if roi_idx(1) >= roi_idx(2)
    rawAllConnections2(2,1,ci:ci+nChan-1, :) = hitVals; 
    rawAllConnections2(2,2,ci:ci+nChan-1, :) = missVals; 
    ci = ci+nChan;
    end
    allConnections2(2,1,roi_idx(1), roi_idx(2), 1) = mean(hitVals(:, tim<0), 'all');
    allConnections2(2,1,roi_idx(1), roi_idx(2), 2) = mean(hitVals(:, tim>=0 & tim<=1000), 'all');
    allConnections2(2,1,roi_idx(1), roi_idx(2), 3) = mean(hitVals(:, tim>=1000 & tim<=2000), 'all'); 
    allConnections2(2,2,roi_idx(1), roi_idx(2), 1) = mean(missVals(:, tim<0), 'all');
    allConnections2(2,2,roi_idx(1), roi_idx(2), 2) = mean(missVals(:, tim>=0 & tim<=1000), 'all');
    allConnections2(2,2,roi_idx(1), roi_idx(2), 3) = mean(missVals(:, tim>=1000 & tim<=2000), 'all');

    subplot 426
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    plot(tim, mean(missVals,1), 'linewidth', 3, 'color', 'red')
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
    sdMiss = std(missVals, [], 1) ./ sqrt(PPC.N); 
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
    x = [tim, flip(tim)];
     y = [mean(missVals,1) - sdMiss, flip(mean(missVals,1)) + flip(sdMiss)];  
    fill(flip(x), flip(y), 'red', 'FaceAlpha', .2)


    scatter(tim(PPC.highp_ret<.05), PPC.highp_ret(PPC.highp_ret<.05)*0, 30,  'k', 'filled')
    ylabel("8Hz connectivity")
    title(['PPC ' roiLabs{1} ' to ' roiLabs{2} ' retrieval'])


    subplot 428
    hold off
    plot(tim, PPC.hightVals_ret, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.high975_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.high025_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value hit/miss")


    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' roiLabs{1} '_' roiLabs{2} '.jpg'], '-r300')
    catch
        disp(['missing data for ' num2str(ii)])

    end
end




load("R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\ACC.mat")
ROInames = {HFBdat.aggTargs.lab};

clVals = [.01, .05];

figure('position', [0,0,1000,1000])
subplot 341
makeConnectionPlot(squeeze(allConnections(1, 1, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("encoding hit")

subplot 342
makeConnectionPlot(squeeze(allConnections(1, 2, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("encoding miss")
subplot 343
makeConnectionPlot(squeeze(allConnections2(1, 1, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("retrieval hit")
subplot 344
makeConnectionPlot(squeeze(allConnections2(1, 2, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("retrieval miss")



subplot 345
makeConnectionPlot(squeeze(allConnections(1, 1, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot 346
makeConnectionPlot(squeeze(allConnections(1, 2, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot 347
makeConnectionPlot(squeeze(allConnections2(1, 1, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot 348
makeConnectionPlot(squeeze(allConnections2(1, 2, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)


subplot 349
makeConnectionPlot(squeeze(allConnections(1, 1, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot(3,4,10)
makeConnectionPlot(squeeze(allConnections(1, 2, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot(3,4,11)
makeConnectionPlot(squeeze(allConnections2(1, 1, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot(3,4,12)
makeConnectionPlot(squeeze(allConnections2(1, 2, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)

export_fig("G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\allLowFreq.jpg", '-r300')


figure('position', [0,0,1000,1000])
subplot 341
makeConnectionPlot(squeeze(allConnections(2, 1, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("encoding hit")

subplot 342
makeConnectionPlot(squeeze(allConnections(2, 2, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("encoding miss")
subplot 343
makeConnectionPlot(squeeze(allConnections2(2, 1, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("retrieval hit")
subplot 344
makeConnectionPlot(squeeze(allConnections2(2, 2, :, :, 1)), clVals, ROInames, ROInames, 1:11, 1:11)
title("retrieval miss")



subplot 345
makeConnectionPlot(squeeze(allConnections(2, 1, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot 346
makeConnectionPlot(squeeze(allConnections(2, 2, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot 347
makeConnectionPlot(squeeze(allConnections2(2, 1, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot 348
makeConnectionPlot(squeeze(allConnections2(2, 2, :, :, 2)), clVals, ROInames, ROInames, 1:11, 1:11)


subplot 349
makeConnectionPlot(squeeze(allConnections(2, 1, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot(3,4,10)
makeConnectionPlot(squeeze(allConnections(2, 2, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot(3,4,11)
makeConnectionPlot(squeeze(allConnections2(2, 1, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)
subplot(3,4,12)
makeConnectionPlot(squeeze(allConnections2(2, 2, :, :, 3)), clVals, ROInames, ROInames, 1:11, 1:11)

export_fig("G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\allHighFreq.jpg", '-r300')











%% latency plots: 
% 
% for ii = 1:length(statFiles)
%     ii
%     PPC = load([statFiles(ii).folder '/' statFiles(ii).name]).connectionDat; 
%     
%     roiLabs = {PPC.reg1, PPC.reg2};
%     roi_idx = [PPC.reg1i, PPC.reg2i]; 
%   
%     f = figure('visible', false);
% %     f.Position = [0 0 600 1000];
% 
% 
%     %% encoding 
% 
%     hitVals = squeeze(PPC.lowBand(PPC.hmSort,:,2) );
%     missVals = squeeze(PPC.lowBand(~PPC.hmSort,:,2) );
%     tim = PPC.tim; 
% 
%     hitLat = zeros(size(hitVals,1),1); 
%     missLat = zeros(size(missVals,1),1); 
% 
%     for ci = 1:PPC.N
%         %hit
%         cur = hitVals(ci,:); 
%         cur = cur(tim>=0 & tim<=PPC.subhitRT(ci));
%         testTim = tim(tim>=0 & tim<=PPC.subhitRT(ci));
%         cur = cur - min(cur); 
%         cur = cur ./ max(cur); 
%         latidx = wmean(1:length(testTim), cur); 
%         hitLat(ci) = testTim(round(latidx)); 
%         %miss
%         cur = missVals(ci,:); 
%         cur = cur(tim>=0 & tim<=PPC.submissRT(ci));
%         testTim = tim(tim>=0 & tim<=PPC.submissRT(ci));
%         cur = cur - min(cur); 
%         cur = cur ./ max(cur); 
%         latidx = wmean(1:length(testTim), cur); 
%         missLat(ci) = testTim(round(latidx));
% 
% 
%     end
% 
%     subplot 121
%     hold off
%     histogram(hitLat, [0:100:2000])
%     hold on
%     histogram(missLat, [0:100:2000])
%     curDat = [hitLat; missLat];
%     hmSort = [ones(length(hitLat),1); zeros(length(missLat),1)];
%     modDat = table(curDat, hmSort, [PPC.allSubs; PPC.allSubs], ...
%                 'VariableNames', {'latency', 'hitMiss', 'sub'}); 
%     lme = fitlme(modDat, 'latency ~  hitMiss +  (1|sub)'); 
%     tval = [0,0];
%     tval(1) = lme.Coefficients(2,4); 
%     tval(2) = lme.Coefficients(2,6); 
%     text(0, PPC.N*.1, ['t-val: ' num2str(round(tval(1),1)) ])
%     text(0, PPC.N*.15, [' p-val: '  num2str(round(tval(2),2))])
%     
%     title(['encoding latency ' PPC.reg1 ' ' PPC.reg2 ])
%   
% 
%     %% retrieval 
%     PPC = load([statFiles2(ii).folder '/' statFiles2(ii).name]).connectionDat; 
%     hitVals = squeeze(PPC.lowBand2(PPC.hmSort,:,2) );
%     missVals = squeeze(PPC.lowBand2(~PPC.hmSort,:,2) );
%     tim = PPC.tim2; 
% 
%    
% 
%     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' roiLabs{1} '_' roiLabs{2} '.jpg'], '-r300')
% 
% end

%% latency plot encoding copying the LL method to put all on one

figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
load("R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\ACC.mat")
ROInames = {HFBdat.aggTargs.ROI};
for ii = 1:length(statFiles)
    ii
    PPC = load([statFiles(ii).folder '/' statFiles(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))


    hitVals = squeeze(PPC.lowBand(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.lowBand(~PPC.hmSort,:,2) );
    tim = PPC.tim; 

    hitLat = zeros(size(hitVals,1),1); 
    missLat = zeros(size(missVals,1),1); 

    for ci = 1:PPC.N
        %hit
        cur = hitVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.subhitRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.subhitRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        hitLat(ci) = testTim(round(latidx)); 
        %miss
        cur = missVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.submissRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.submissRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        missLat(ci) = testTim(round(latidx));


    end


    curDat = [hitLat; missLat];


    hmSort = [ones(length(hitLat),1); zeros(length(missLat),1)];
    modDat = table(curDat, hmSort, [PPC.allSubs; PPC.allSubs], ...
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
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


   


    
  

end


export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' 'PPC_latency_lowEncode' '.jpg'],''), '-r300')



%% latency plot retrieval copying the LL method to put all on one

figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
load("R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\ACC.mat")
ROInames = {HFBdat.aggTargs.ROI};
for ii = 1:length(statFiles)
    ii
    PPC = load([statFiles2(ii).folder '/' statFiles2(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))


    hitVals = squeeze(PPC.lowBand2(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.lowBand2(~PPC.hmSort,:,2) );
    tim = PPC.tim2; 

    hitLat = zeros(size(hitVals,1),1); 
    missLat = zeros(size(missVals,1),1); 

    for ci = 1:PPC.N
        %hit
        cur = hitVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.rethitRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.rethitRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        hitLat(ci) = testTim(round(latidx)); 
        %miss
        cur = missVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.retmissRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.retmissRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        missLat(ci) = testTim(round(latidx));


    end


    curDat = [hitLat; missLat];


    hmSort = [ones(length(hitLat),1); zeros(length(missLat),1)];
    modDat = table(curDat, hmSort, [PPC.allSubs; PPC.allSubs], ...
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
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


    
  

end


export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' 'PPC_latency_lowRetrieve' '.jpg'],''), '-r300')


%% latency plot encoding copying the LL method to put all on one HIGH FREQUENCY

figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
load("R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\ACC.mat")
ROInames = {HFBdat.aggTargs.ROI};
for ii = 1:length(statFiles)
    ii
    PPC = load([statFiles(ii).folder '/' statFiles(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))


    hitVals = squeeze(PPC.highBand(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.highBand(~PPC.hmSort,:,2) );
    tim = PPC.tim; 

    hitLat = zeros(size(hitVals,1),1); 
    missLat = zeros(size(missVals,1),1); 

    for ci = 1:PPC.N
        %hit
        cur = hitVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.subhitRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.subhitRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        hitLat(ci) = testTim(round(latidx)); 
        %miss
        cur = missVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.submissRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.submissRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        missLat(ci) = testTim(round(latidx));


    end


    curDat = [hitLat; missLat];


    hmSort = [ones(length(hitLat),1); zeros(length(missLat),1)];
    modDat = table(curDat, hmSort, [PPC.allSubs; PPC.allSubs], ...
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
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


   


    
  

end


export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' 'PPC_latency_highEncode' '.jpg'],''), '-r300')



%% latency plot retrieval copying the LL method to put all on one HIGH FREQUENCY

figure('visible', true, 'position', [100,100,1000, 1000])
colormap = [ flip([.01:.01:1])', ones(100,1), flip([.01:.01:1])'];
findMap = flip([.01:.01:1]); 
load("R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\ACC.mat")
ROInames = {HFBdat.aggTargs.ROI};
for ii = 1:length(statFiles)
    ii
    PPC = load([statFiles2(ii).folder '/' statFiles2(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    roi_idx = zeros(2,1); 
    roi_idx(1) = find(cellfun(@(x) strcmp(roiLabs{1}, x), ROInames)); 
    roi_idx(2) = find(cellfun(@(x) strcmp(roiLabs{2}, x), ROInames)); 
    subplot(11, 11, (roi_idx(1)-1)*11+roi_idx(2))


    hitVals = squeeze(PPC.highBand2(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.highBand2(~PPC.hmSort,:,2) );
    tim = PPC.tim2; 

    hitLat = zeros(size(hitVals,1),1); 
    missLat = zeros(size(missVals,1),1); 

    for ci = 1:PPC.N
        %hit
        cur = hitVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.rethitRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.rethitRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        hitLat(ci) = testTim(round(latidx)); 
        %miss
        cur = missVals(ci,:); 
        cur = cur(tim>=0 & tim<=PPC.retmissRT(ci));
        testTim = tim(tim>=0 & tim<=PPC.retmissRT(ci));
        cur = cur - min(cur); 
        cur = cur ./ max(cur); 
        latidx = wmean(1:length(testTim), cur); 
        missLat(ci) = testTim(round(latidx));


    end


    curDat = [hitLat; missLat];


    hmSort = [ones(length(hitLat),1); zeros(length(missLat),1)];
    modDat = table(curDat, hmSort, [PPC.allSubs; PPC.allSubs], ...
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
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability', 'facealpha', 1, 'facecolor', 'white');
    h = histogram(hitLat, [0:100:2000], 'normalization', 'probability', 'facecolor', 'blue');
    
    h2 = histogram(missLat, [0:100:2000], 'normalization', 'probability',  'facecolor', 'red');
    allMax = max(max(h.Values), max(h2.Values)); 
    ylim([0, allMax*1.1])


    xticklabels([]); 
    yticklabels([]); 


    
  

end


export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' 'PPC_latency_highRetrieve' '.jpg'],''), '-r300')






%% GED: doesn't work! needs longer vectors going into the covariance matrix
% 
% rawAllConnections(:,:,ci:20000,:) = []; 
% rawAllConnections2(:,:,ci:20000,:) = []; 
% rawAllConnectionsSum(ci:end,:) = []; 
% 
% 
% distsHit = squeeze(rawAllConnections(1, 1, :, :));
% distsHit = distsHit - min(distsHit,[],2); 
% distsHit = distsHit ./ max(distsHit,[], 2); 
% distsMiss = squeeze(rawAllConnections(1, 2, :, :)); 
% distsMiss = distsMiss - min(distsMiss,[],2); 
% distsMiss = distsMiss ./ max(distsMiss,[], 2); 
%  
% S = cov(distsHit'); 
% R = cov(distsMiss'); 
% 
% % regularize R
% gamma = .01;
% Rr = R *(1-gamma) + eye(length(R))*gamma*mean(eig(R));
% 
% % global variance normalize
% S = S / (std(S(:))/std(R(:)));
% 
% test = DBscanDynamicEpi(S, 3, "mapDist", 1, 1); 
% 
% 
% %eigen decomposition won't work because with only 141 timesteps going in
% %and 141<<12000 pairs, the covariance matrix is going to be singular
% [V, D] = eigs(S, R, 141); 
% [L, sidx] = sort(diag(D), 'descend'); 
% V = V(:,sidx); 
% 
% figure
% plot(diag(D))
% xlim([0,50])
% 
% distsHit = squeeze(rawAllConnections(1, 1, :, :));
% distsMiss = squeeze(rawAllConnections(1, 2, :, :)); 
% 
% for ii = 1:10
% figure
% subplot 211
% hit = (V(:,ii)' * distsHit); 
% hit = reshape(hit, size(allConSub, [3,4])); 
% imagesc(LLdat.encTim, [], hit)
% propVar = L(ii) / sum(L); 
% title(['GED componenent: ' num2str(ii) ' prop Var:' num2str(round(propVar,2))])
% colorbar
% 
% hit = V(:,ii)' * distsHitCov; 
% test = squeeze(mean(distsHit(hit>prctile(hit, 90), :),1)); 
% test = reshape(test, size(allConSub, [3,4])); 
% subplot 212
% imagesc(test)
% colorbar
% 
% export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\GEDFigs\cmp' num2str(ii) '.jpg'] )
% 
% 
% end


