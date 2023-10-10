function [] = HFBstatspipeline(statFiles, idx, figPath)

        HFBdat = load([statFiles(idx).folder '/' statFiles(idx).name]).HFBdat; 

        
        regRes = HFBdat.regRes; 
        regRes2 =HFBdat.regRes2; 
        regSubs = HFBdat.regSubs; 
        regSubIDs = HFBdat.regSubIDs; 
        aggTargs = HFBdat.aggTargs; 
        reg1 = HFBdat.reg1; 
        regd = HFBdat.d; 
        regacc = HFBdat.acc; 
        chani = HFBdat.chani; 
        realID = HFBdat.realID; 

        
        %% encoding mean difference

        hitVals = permute(squeeze(regRes(1,:,:)), [2,1]);
        missVals = permute(squeeze(regRes(2,:,:)), [2,1]);
        
        hmSort = [ones(size(hitVals,1),1); zeros(size(hitVals,1),1)];
        hmSort = hmSort==1; 
        tVals = zeros(size(hitVals,2),1);
    
        for ti = 1:size(hitVals,2)
    
            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, hmSort, [realID; realID], ...
                'VariableNames', {'HFB', 'hitMiss', 'sub'}); 
            lme = fitlme(modDat, 'HFB ~  hitMiss +  (1|sub)'); 
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
        disp('calculated encoding t-vals ')
       
        perms = 1000; 
        tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            slice = nullTs(:,ii); 
            for ti = 1:size(hitVals,2)
                curdat = [hitVals(:,ti); missVals(:,ti)]; 
        
                curdat = randsample(curdat, length(curdat), false); 
    
                modDat = table(curdat, hmSort, [realID; realID], ...
                'VariableNames', {'HFB', 'hitMiss', 'sub'}); 
                lme = fitlme(modDat, 'HFB ~  hitMiss +  (1|sub)'); 
                slice(ti) = lme.Coefficients(2,4);
            end
            nullTs(:,ii) = slice; 
        end
        toc
        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        disp('calculated encoding permutation ')


       HFBdat.tVals_sub = tVals; 
       HFBdat.hitVals_sub = mean(hitVals); 
       HFBdat.missVals_sub = mean(missVals);
       HFBdat.p_sub = p; 


       %% retrieval mean difference


        hitVals = permute(squeeze(regRes2(1,:,:)), [2,1]);
        missVals = permute(squeeze(regRes2(2,:,:)), [2,1]);
        
        hmSort = [ones(size(hitVals,1),1); zeros(size(hitVals,1),1)];
        hmSort = hmSort==1; 
        tVals = zeros(size(hitVals,2),1);
    
        for ti = 1:size(hitVals,2)
    
            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, hmSort, [realID; realID], ...
                'VariableNames', {'HFB', 'hitMiss', 'sub'}); 
            lme = fitlme(modDat, 'HFB ~  hitMiss +  (1|sub)'); 
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
    
       
        disp('calculated retrieval t-vals ')
        perms = 1000; 
        tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            slice = nullTs(:,ii); 
            for ti = 1:size(hitVals,2)
                curdat = [hitVals(:,ti); missVals(:,ti)]; 
        
                curdat = randsample(curdat, length(curdat), false); 
    
                modDat = table(curdat, hmSort, [realID; realID], ...
                'VariableNames', {'HFB', 'hitMiss', 'sub'}); 
                lme = fitlme(modDat, 'HFB ~  hitMiss +  (1|sub)'); 
                slice(ti) = lme.Coefficients(2,4);
            end
            nullTs(:,ii) = slice; 
        end
        toc
        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 

        disp('calculated retrieval permutation')

       HFBdat.tVals_ret = tVals; 
       HFBdat.hitVals_ret = mean(hitVals); 
       HFBdat.missVals_ret = mean(missVals);
       HFBdat.p_ret = p;

        %% encoding latency

        tim = HFBdat.encTim; 
        hitVals = permute(squeeze(regRes(1,:,:)), [2,1]);
        missVals = permute(squeeze(regRes(2,:,:)), [2,1]);
        
        hitMeanTim = arrayfun(@(x) wmean(tim, hitVals(x,:) - min(hitVals(x,:)) ), 1:size(hitVals,1));
        missMeanTim = arrayfun(@(x) wmean(tim, missVals(x,:) - min(missVals(x,:)) ), 1:size(missVals,1));
        

        HFBdat.hitTim_sub = hitMeanTim; 
        HFBdat.missTim_sub = missMeanTim; 

        %% retrieval latency

        tim = HFBdat.retTim; 
        hitVals = permute(squeeze(regRes2(1,:,:)), [2,1]);
        missVals = permute(squeeze(regRes2(2,:,:)), [2,1]);
        
        hitMeanTim = arrayfun(@(x) wmean(tim, hitVals(x,:) - min(hitVals(x,:)) ), 1:size(hitVals,1));
        missMeanTim = arrayfun(@(x) wmean(tim, missVals(x,:) - min(missVals(x,:)) ), 1:size(missVals,1));
        

        HFBdat.hitTim_ret = hitMeanTim; 
        HFBdat.missTim_ret = missMeanTim; 



        save([statFiles(idx).folder '/' statFiles(idx).name], 'HFBdat', '-v7.3');



end