
%PPC stats plot

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\';


addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

path = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\PCC_KEY_STATS';
path2 = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\PCC_KEY_STATS2';
statFiles = dir(path); 
statFiles2 = dir(path2); 

statFiles(1:2) = []; 
statFiles2(1:2) = []; 



parfor ii = 1:length(statFiles)
    ii
    PPC = load([statFiles(ii).folder '/' statFiles(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    roi_idx = [PPC.reg1i, PPC.reg2i]; 
  
    f = figure('visible', false);
    f.Position = [0 0 600 1000];

    hitVals = squeeze(PPC.lowBand(PPC.hmSort,:,2) );
%     missVals = squeeze(PPC.lowBand(~PPC.hmSort,:,2) );
    tim = PPC.tim; 

   

    subplot 621
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
   
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
 
    scatter(tim(PPC.lowp_sub<.05), PPC.lowp_sub(PPC.lowp_sub<.05)*0, 30,  'k', 'filled')
    ylabel("3Hz connectivity")
    title(['PPC ' roiLabs{1} ' to ' roiLabs{2} ' encoding'])



    subplot 623
    hold off
    plot(tim, PPC.lowtVals, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.low975, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.low025, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value mem")

    





   
   
    tim = PPC.tim; 

    subplot 625
    hold off
    tidx = 1:length(tim); 
    tidx = tidx(PPC.lowp_sub<.05);
    if isempty(tidx)
        tidx = 1:length(tim); 
    end
    [~, ti] = max(abs(PPC.lowtVals(tidx)));
    maxT = PPC.lowtVals(tidx(ti));
    ti = tidx(ti); 
    dVals = PPC.d;
    plotVals = hitVals(:, ti);

    scatter(dVals, plotVals)
    title(['time: ' num2str(tim(ti)) ' p: ' num2str(round(PPC.lowp_sub(ti),2)) ' t: ' num2str(round(maxT,2))])
    ylabel('PPC')
    xlabel('memory performance')
   
    hitVals = squeeze(PPC.highBand(PPC.hmSort,:,2) );

    subplot 627
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
    


    scatter(tim(PPC.highp_sub<.05), PPC.highp_sub(PPC.highp_sub<.05)*0, 30,  'k', 'filled')
    ylabel("8Hz connectivity")
    


    subplot 629
    hold off
    plot(tim, PPC.hightVals, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.high975, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.high025, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value mem")


    subplot(6,2,11)
    hold off
    tidx = 1:length(tim); 
    tidx = tidx(PPC.highp_sub<.05);
    if isempty(tidx)
        tidx = 1:length(tim); 
    end
   
    [~, ti] = max(abs(PPC.hightVals(tidx)));
    maxT = PPC.hightVals(tidx(ti));
    ti = tidx(ti); 
    dVals = PPC.d;
    plotVals = hitVals(:, ti);

    scatter(dVals, plotVals)
    title(['time: ' num2str(tim(ti)) ' p: ' num2str(round(PPC.highp_sub(ti),2)) ' t: ' num2str(round(maxT,2))])
    ylabel('PPC')
    xlabel('memory performance')

    %% retrieval 
    PPC = load([statFiles2(ii).folder '/' statFiles2(ii).name]).connectionDat; 
    hitVals = squeeze(PPC.lowBand2(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.lowBand2(~PPC.hmSort,:,2) );
    tim = PPC.tim2; 



       subplot 622
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
   
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
 
    scatter(tim(PPC.lowp_ret<.05), PPC.lowp_ret(PPC.lowp_ret<.05)*0, 30,  'k', 'filled')
    ylabel("3Hz connectivity")
    title(['PPC ' roiLabs{1} ' to ' roiLabs{2} ' encoding'])



    subplot 624
    hold off
    plot(tim, PPC.lowtVals_ret, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.low975_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.low025_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value mem")

    





   
   
 

    subplot 626
    hold off
    tidx = 1:length(tim); 
    tidx = tidx(PPC.lowp_ret<.05);
    if isempty(tidx)
        tidx = 1:length(tim); 
    end
    [~, ti] = max(abs(PPC.lowtVals_ret(tidx)));
    maxT = PPC.lowtVals_ret(tidx(ti));
    ti = tidx(ti); 
    dVals = PPC.d;
    plotVals = hitVals(:, ti);

    scatter(dVals, plotVals)
    title(['time: ' num2str(tim(ti)) ' p: ' num2str(round(PPC.lowp_ret(ti),2)) ' t: ' num2str(round(maxT,2))])
    ylabel('PPC')
    xlabel('memory performance')
   
    hitVals = squeeze(PPC.highBand2(PPC.hmSort,:,2) );

    subplot 628
    hold off
    plot(tim, mean(hitVals,1), 'linewidth', 3, 'color', 'blue')
    hold on 
    xlim([tim(1), tim(end)])


    sdHit = std(hitVals, [], 1) ./ sqrt(PPC.N);
    
    x = [tim, flip(tim)];
    y = [mean(hitVals,1) - sdHit, flip(mean(hitVals,1)) + flip(sdHit)]; 
    fill(flip(x), flip(y), 'blue', 'FaceAlpha', .2)
    


    scatter(tim(PPC.highp_ret<.05), PPC.highp_ret(PPC.highp_ret<.05)*0, 30,  'k', 'filled')
    ylabel("8Hz connectivity")
    


    subplot(6,2,10)
    hold off
    plot(tim, PPC.hightVals_ret, 'linewidth', 3, 'color', 'green')
    hold on 
    plot(tim, PPC.high975_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    plot(tim, PPC.high025_ret, 'linewidth', 1, 'color', 'k', 'linestyle' ,'--')
    xlim([tim(1), tim(end)])
    ylabel("t-value mem")


    subplot(6,2,12)
    hold off
    tidx = 1:length(tim); 
    tidx = tidx(PPC.highp_ret<.05);
    if isempty(tidx)
        tidx = 1:length(tim); 
    end
  
    [~, ti] = max(abs(PPC.hightVals_ret(tidx)));
    maxT = PPC.hightVals_ret(tidx(ti));
    ti = tidx(ti); 
    dVals = PPC.d;
    plotVals = hitVals(:, ti);

    scatter(dVals, plotVals)
    title(['time: ' num2str(tim(ti)) ' p: ' num2str(round(PPC.highp_ret(ti),2)) ' t: ' num2str(round(maxT,2))])
    ylabel('PPC')
    xlabel('memory performance')


    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\PPC_finalizedFigs\' roiLabs{1} '_' roiLabs{2} '_mem' '.jpg'], '-r300')

end





load("R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\ACC.mat")
ROInames = {HFBdat.aggTargs.ROI};

clVals = [.005, .02];

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

