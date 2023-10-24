        %matrices of channels X timepoints for hits and misses separately
        hitVals = permute(squeeze(regRes(1,:,:)), [2,1]); 
        missVals = permute(squeeze(regRes(2,:,:)), [2,1]);
        
        %vector that is 2*numchanels long and has 1s and 0s for hits and
        %misses
        hmSort = [ones(size(hitVals,1),1); zeros(size(hitVals,1),1)];
        hmSort = hmSort==1; 
        
        %empty vector that is 1Xtimepoints for t values
        tVals = zeros(size(hitVals,2),1);
    
        %loop on timepoints
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