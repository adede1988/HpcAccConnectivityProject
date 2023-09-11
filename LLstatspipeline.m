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





        %% retrieval! 
        

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
        
       
        save([statFiles(idx).folder '/' statFiles(idx).name], 'LLdat', '-v7.3'); 







end