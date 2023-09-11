function [] = LLstatspipeline(statFiles, idx, figPath)

        LLdat = load([statFiles(idx).folder '/' statFiles(idx).name]).LLdat; 

        
        regRes = LLdat.regRes; 
        regRes2 =LLdat.regRes2; 
        regSubs = LLdat.regSubs; 
        regSubIDs = LLdat.regSubIDs; 
        aggTargs = LLdat.aggTargs; 
        reg1 = LLdat.reg1; 
        reg2 = LLdat.reg2; 
        regd = LLdat.d; 
        regacc = LLdat.acc; 


        %% encoding mean difference

        hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
        missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
        tVals = myArrayT(hitVals, missVals,1);
        perms = 1000; 

        nullTs = zeros([size(tVals), perms]); 
        for ii = 1:perms

            nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
        end

        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        %store results 
        LLdat.tVals_sub = tVals; 
        LLdat.hitVals_sub = squeeze(mean(hitVals,1));
        LLdat.missVals_sub = squeeze(mean(missVals,1)); 
        LLdat.p_sub = p; 

        tim = LLdat.encTim; 
        roiLabs = {aggTargs(reg2).ROI, aggTargs(reg1).ROI};
        
        makeLLplot(tim, hitVals, missVals, tVals,...
            clusterinfo, regSubs, figPath , roiLabs, 'sub')


         



        %% retrieval mean difference
        

        hitVals = permute(squeeze(regRes2(1,:,:,:)), [3,1,2]);
        missVals = permute(squeeze(regRes2(2,:,:,:)), [3,1,2]);
        tVals = myArrayT(hitVals, missVals,1);
        perms = 1000; 

        nullTs = zeros([size(tVals), perms]); 
        for ii = 1:perms

            nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
        end

        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        %store results    
        LLdat.tVals_ret = tVals; 
        LLdat.hitVals_ret = squeeze(mean(hitVals,1));
        LLdat.missVals_ret = squeeze(mean(missVals,1)); 
        LLdat.p_ret = p; 

        tim = LLdat.retTim; 
        
        makeLLplot(tim, hitVals, missVals, tVals,...
            clusterinfo, regSubs, figPath , roiLabs, 'ret')
        
       
         


        %% encoding latency

        tim = LLdat.encTim; 
        LLtim = -150:150; 
        hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
        missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
        hitMeanTim = arrayfun(@(x) wmean(tim', squeeze(mean(hitVals(x,:,:),2)) - min(squeeze(mean(hitVals(x,:,:),2))) ), 1:size(hitVals,1));
        
        missMeanTim = arrayfun(@(x) wmean(tim', squeeze(mean(missVals(x,:,:),2)) - min(squeeze(mean(missVals(x,:,:),2))) ), 1:size(missVals,1));
        

        timBorders = [-500, 0, 500, 1000, 1500, 2000]; 
        hitMeanLL = zeros(size(hitVals,1), length(timBorders)); 
        missMeanLL = zeros(size(hitVals,1), length(timBorders)); 
        for tt = 1:length(timBorders)-1
            ti = tim>=timBorders(tt) & tim<=timBorders(tt+1); 
            hitMeanLL(:,tt) = arrayfun(@(x) wmean(LLtim, squeeze(mean(hitVals(x,:,ti),3)) - min(squeeze(mean(hitVals(x,:,ti),3))) ), 1:size(hitVals,1));
            missMeanLL(:,tt) = arrayfun(@(x) wmean(LLtim, squeeze(mean(missVals(x,:,ti),3)) - min(squeeze(mean(missVals(x,:,ti),3))) ), 1:size(missVals,1));
        end

        LLdat.hitTim_sub = hitMeanTim; 
        LLdat.missTim_sub = missMeanTim; 
        LLdat.hitLL_sub = hitMeanLL; 
        LLdat.missLL_sub = missMeanLL; 




        %% retrieval latency

        tim = LLdat.retTim; 
        LLtim = -150:150; 
        hitVals = permute(squeeze(regRes2(1,:,:,:)), [3,1,2]);
        missVals = permute(squeeze(regRes2(2,:,:,:)), [3,1,2]);
        hitMeanTim = arrayfun(@(x) wmean(tim', squeeze(mean(hitVals(x,:,:),2)) - min(squeeze(mean(hitVals(x,:,:),2))) ), 1:size(hitVals,1));
        
        missMeanTim = arrayfun(@(x) wmean(tim', squeeze(mean(missVals(x,:,:),2)) - min(squeeze(mean(missVals(x,:,:),2))) ), 1:size(missVals,1));
        

        timBorders = [-500, 0, 500, 1000, 1500, 2000]; 
        hitMeanLL = zeros(size(hitVals,1), length(timBorders)); 
        missMeanLL = zeros(size(hitVals,1), length(timBorders)); 
        for tt = 1:length(timBorders)-1
            ti = tim>=timBorders(tt) & tim<=timBorders(tt+1); 
            hitMeanLL(:,tt) = arrayfun(@(x) wmean(LLtim, squeeze(mean(hitVals(x,:,ti),3)) - min(squeeze(mean(hitVals(x,:,ti),3))) ), 1:size(hitVals,1));
            missMeanLL(:,tt) = arrayfun(@(x) wmean(LLtim, squeeze(mean(missVals(x,:,ti),3)) - min(squeeze(mean(missVals(x,:,ti),3))) ), 1:size(missVals,1));
        end

        LLdat.hitTim_ret = hitMeanTim; 
        LLdat.missTim_ret = missMeanTim; 
        LLdat.hitLL_ret = hitMeanLL; 
        LLdat.missLL_ret = missMeanLL; 
        



        save([statFiles(idx).folder '/' statFiles(idx).name], 'LLdat', '-v7.3');

end