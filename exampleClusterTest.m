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
    
            %it is implicitly assumed here that hitVals and missVals are in
            %the same channel order. If this is not true, then the subject
            %IDs won't line up correctly. 
            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, hmSort, [realID; realID], ...
                'VariableNames', {'HFB', 'hitMiss', 'sub'}); 
            lme = fitlme(modDat, 'HFB ~  hitMiss +  (1|sub)'); 
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
        disp('calculated encoding t-vals ')
       
        perms = 1000; 
       
        nullTs = squeeze(zeros([size(tVals), perms])); 

        %using the "slice" variable makes the code ready for parfor usage
        %in outer loop if speedup through parallelization is desired 
        for ii = 1:perms
            slice = nullTs(:,ii); 
            for ti = 1:size(hitVals,2)

                curdat = [hitVals(:,ti); missVals(:,ti)]; 
                %shuffle within each channel independently. This avoids
                %having both observations of a channel be hits or both
                %misses. Instead, it's going to do a coin flip for each
                %channel to determine if the first or second observation is
                %the "hit" (in truth, it's always the first)
                hmShuff = zeros(length(hmSort), 1); 
                for chan = 1:length(realID)
                    chani = [chan, chan+length(realID)]; 
                    if rand() > .5 %coin flip
                        hmShuff(chani(1)) = 1; %first obs of chan is HIT
                    else
                        hmShuff(chani(2)) = 1; %second obs of chan is HIT
                    end
                end
                hmShuff = hmShuff ==1; %convert to bool
    
          
                modDat = table(curdat, hmShuff, [realID; realID], ...
                'VariableNames', {'HFB', 'hitMiss', 'sub'}); 
                lme = fitlme(modDat, 'HFB ~  hitMiss +  (1|sub)'); 
                slice(ti) = lme.Coefficients(2,4);
            end
            nullTs(:,ii) = slice; 
        end
        
        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 