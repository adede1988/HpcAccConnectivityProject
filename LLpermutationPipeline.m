function [] = LLpermutationPipeline(paths)
    tic
    %load data
    LL = load([paths.folder '/' paths.name]).LLout;
    disp('data loaded')
    %make output struct
    LLstats = rmfield(LL, 'LL');
    LLstats.encTim = LL(1).LL.encTim;
    LLstats.subMiss = 1; 
    LLstats.subHit = 1; 
    LLstats.contrastEnc = 1; 
    LLstats.retTim = LL(1).LL.retTim; 
    LLstats.missRet = 1; 
    LLstats.hitRet = 1; 
    LLstats.contrastRet = 1; 

    %encoding conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tim = LL.LL.encTim;
    LLstats.encTim = tim; 
    %subMiss
    allObs = LL.LL.subMiss; 
    LLstats.subMiss = conditionCluTest(allObs, tim); 
    disp(['subMiss finished; time: ' num2str(round(toc/60) ) ])

    %subHit
    allObs = LL.LL.subHit; 
    LLstats.subHit = conditionCluTest(allObs, tim); 
    disp(['subHit finished; time: ' num2str(round(toc/60) ) ])

    %contrastEnc
    allObs = LL.LL.subHit - LL.LL.subMiss;
    LLstats.contrastEnc = conditionCluTest(allObs, tim); 
    disp(['encContrast finished; time: ' num2str(round(toc/60) ) ])
    

    %retrieval conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tim = LL.LL.retTim;

    LLstats.retTim = tim; 
    %retMiss
    allObs = LL.LL.missRet; 
    LLstats.missRet = conditionCluTest(allObs, tim); 
    disp(['missRet finished; time: ' num2str(round(toc/60) ) ])

    %retHit
    allObs = LL.LL.hitRet; 
    LLstats.hitRet = conditionCluTest(allObs, tim); 
    disp(['hitRet finished; time: ' num2str(round(toc/60) ) ])

    %contrastEnc
    allObs = LL.LL.hitRet - LL.LL.missRet;
    LLstats.contrastRet = conditionCluTest(allObs, tim); 
    disp(['retContrast finished; time: ' num2str(round(toc/60) ) ] )


    %save output
    outname = split(paths.name, '.'); 

    save([paths.folder '/' outname{1} '_stats.mat'], "LLstats", '-v7.3')








end