
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



for ii = 1:length(statFiles)

    PPC = load([statFiles(ii).folder '/' statFiles(ii).name]).connectionDat; 
    
    roiLabs = {PPC.reg1, PPC.reg2};
    f = figure('visible', false);
    f.Position = [0 0 600 1000];

    hitVals = squeeze(PPC.lowBand(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.lowBand(~PPC.hmSort,:,2) );
    tim = PPC.tim; 


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





    hitVals = squeeze(PPC.highBand(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.highBand(~PPC.hmSort,:,2) );
    tim = PPC.tim; 
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
    hitVals = squeeze(PPC.lowBand2(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.lowBand2(~PPC.hmSort,:,2) );
    tim = PPC.tim2; 


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





    hitVals = squeeze(PPC.highBand2(PPC.hmSort,:,2) );
    missVals = squeeze(PPC.highBand2(~PPC.hmSort,:,2) );
    
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

end


