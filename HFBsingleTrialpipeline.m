function [] = HFBsingleTrialpipeline(statFiles, fileIdx, statType)

        HFBdat = load([statFiles(fileIdx).folder '/' statFiles(fileIdx).name]).statInfo; 
        perms = 2000; 
    
        %% encoding mean difference

        %matrices of channels X timepoints for hits and misses separately
        hitVals = HFBdat.hits'; 
        missVals = HFBdat.misses';
        tim = HFBdat.tim; 
      
        hitVals(:, tim<-450 | tim>3000) = []; 
        missVals(:, tim<-450 | tim>3000) = []; 
        tim(tim<-450 | tim>3000) = []; 
        


        
        %vector that is 2*numchanels long and has 1s and 0s for hits and
        %misses
        hmSort = [ones(size(hitVals,1),1); zeros(size(missVals,1),1)];
        hmSort = hmSort==1; 
        
        %random effects variables
        realID = cellfun(@(x, y, z) [x,'_', num2str(y),'_', num2str(z)], HFBdat.hitSub, ...
            num2cell(HFBdat.hitTriali), num2cell(HFBdat.hitChi), 'UniformOutput',false);

        realID2 = cellfun(@(x, y, z) [x,'_', num2str(y),'_', num2str(z)], HFBdat.missSub, ...
            num2cell(HFBdat.missTriali), num2cell(HFBdat.missChi), 'UniformOutput',false);

        realID = [realID; realID2]; 

        chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.hitSub, ...
             num2cell(HFBdat.hitChi), 'UniformOutput',false);

        chanID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.missSub, ...
             num2cell(HFBdat.missChi), 'UniformOutput',false);

        trialID = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.hitSub, ...
             num2cell(HFBdat.hitTriali), 'UniformOutput',false);

        trialID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.missSub, ...
             num2cell(HFBdat.missTriali), 'UniformOutput',false);


        chanID = [chanID; chanID2]; 
        trialID = [trialID; trialID2]; 
        subID = [HFBdat.hitSub; HFBdat.missSub]; 
        sub_uni = unique(subID); 


        if statType == 1
        %empty vector that is 1Xtimepoints for t values
        tVals = zeros(size(hitVals,2),1);
    
        %loop on timepoints
        for ti = 1:size(hitVals,2)
    

            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, categorical(hmSort), chanID, trialID, subID, realID, ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'trial', 'sub', 'realID'}); 
           
            lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan:sub)  + (1|trial)');
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
        disp('calculated encoding t-vals ')
       

      
        
        tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            if(mod(ii, 100)) ==0 
                disp(['..........................' num2str(ii) ...
                    ' time:' num2str(round(toc/60,1))])
            end

            hmShuff = zeros(length(hmSort),1); 
            for sub = 1:length(sub_uni)
                curi = find(ismember(subID, sub_uni{sub}));
                curHM = hmSort(curi);
                %index into the curi set of current trials indicating
                %unique trial identifiers 
                curiUni = unique(cellfun(@(x) find(...
                    strcmp(x, trialID(curi)),1), ...
                    trialID(curi))); 
                
                %get unique trials
                trials_uni = trialID(curi(curiUni)); 
                curHits = sum(curHM(curiUni)); 
                curMisses = length(curiUni) - curHits;  

                %do the coin flipping!
                for trial =1:length(trials_uni)
                    trialidx = find(ismember(trialID,...
                                            trials_uni{trial}));
                    if rand() > curHits/(curHits+curMisses) 
                        hmShuff(trialidx) = 0; 
                        curMisses = curMisses - 1; 
                    else
                        hmShuff(trialidx) = 1; 
                        curHits = curHits - 1; 
                    end

                end

            end


            
            slice = nullTs(:,ii); 
            for ti = 1:size(hitVals,2)
                curdat = [hitVals(:,ti); missVals(:,ti)]; 
                
                %shuffle such that the number of hits and misses per channel is preserved
                %but which individual trials are hits and which are misses
                %gets shuffled
                
                %to do this, I need to loop through channel IDs (sub + chi)
                %and for each channel id, get its indices across hits and
                %misses, then shuffle hit miss values across those indices
%                 hmShuff = zeros(length(hmSort),1); 
%                 for chan = 1:length(chan_uni)
%                     curi = find(ismember(chanID, chan_uni{chan}));
%                     curHits = find(diff(curi)>2); 
%                     curMisses = length(curi) - curHits; 
%                     %do the coin flipping!
%                     for trial =1:length(curi)
%                         if rand() > curHits/(curHits+curMisses)
%                             hmShuff(curi(trial)) = 0; 
%                             curMisses = curMisses - 1; 
%                         else
%                             hmShuff(curi(trial)) = 1; 
%                             curHits = curHits - 1; 
%                         end
% 
%                     end
% 
%                 end

                %alternate shuffle: shuffle such that all observations of
                %each trial are randomly changed to either hit or miss
                %while maintaining the total number of hits and misses
                %within an individual person
                


                modDat = table(curdat, categorical(hmShuff), chanID, trialID, subID, realID, ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'trial', 'sub', 'realID'}); 
           
                lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan:sub)  + (1|trial)');
                slice(ti) = lme.Coefficients(2,4);
            end
            nullTs(:,ii) = slice; 
        end
        toc
        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        disp('calculated encoding permutation ')


       HFBdat.tVals_image = tVals; 
       HFBdat.hitVals_image = mean(hitVals); 
       HFBdat.missVals_image = mean(missVals);
       HFBdat.p_image = p; 

       statInfo = HFBdat; 
       save([statFiles(fileIdx).folder '/' 'stat' num2str(statType) statFiles(fileIdx).name], 'statInfo', '-v7.3');
        
        else
        
        %% HFB peak aligned analysis
        hitVals = HFBdat.hits'; 
        missVals = HFBdat.misses';
        tim = HFBdat.tim; 

        w = 20; 

        hitLatIdx = arrayfun(@(x) find(x==tim), HFBdat.hitLat); 
        hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 

        peakHit = arrayfun(@(x) hitVals(x,hitLatIdx(x)-w:...
                                         hitLatIdx(x)+w), [1:length(hitLatIdx)], ...
                                         'uniformoutput', false);
        peakHit = cell2mat(peakHit.');

        missLatIdx = arrayfun(@(x) find(x==tim), HFBdat.missLat); 
        missLatIdx(missLatIdx>length(tim)-w) = length(tim) - w; 

        peakMiss = arrayfun(@(x) missVals(x,missLatIdx(x)-w:...
                                         missLatIdx(x)+w), [1:length(missLatIdx)], ...
                                         'uniformoutput', false);
        peakMiss = cell2mat(peakMiss.');



        tVals = zeros(size(peakHit,2),1);
    
        %loop on timepoints
        for ti = 1:size(peakHit,2)
    

            curdat = [peakHit(:,ti); peakMiss(:,ti)];
    
            modDat = table(curdat, categorical(hmSort), chanID, trialID, subID, realID, ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'trial', 'sub', 'realID'}); 
           
            lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan:sub)  + (1|trial)');
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
        disp('calculated encoding t-vals ')

        tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            if(mod(ii, 100)) ==0 
                disp(['..........................' num2str(ii) ...
                    ' time:' num2str(round(toc/60,1))])
            end
            
            hmShuff = zeros(length(hmSort),1); 
            for sub = 1:length(sub_uni)
                curi = find(ismember(subID, sub_uni{sub}));
                curHM = hmSort(curi);
                %index into the curi set of current trials indicating
                %unique trial identifiers 
                curiUni = unique(cellfun(@(x) find(...
                    strcmp(x, trialID(curi)),1), ...
                    trialID(curi))); 
                
                %get unique trials
                trials_uni = trialID(curi(curiUni)); 
                curHits = sum(curHM(curiUni)); 
                curMisses = length(curiUni) - curHits;  

                %do the coin flipping!
                for trial =1:length(trials_uni)
                    trialidx = find(ismember(trialID,...
                                            trials_uni{trial}));
                    if rand() > curHits/(curHits+curMisses) 
                        hmShuff(trialidx) = 0; 
                        curMisses = curMisses - 1; 
                    else
                        hmShuff(trialidx) = 1; 
                        curHits = curHits - 1; 
                    end

                end

            end

            hmShuff = hmShuff ==1; 
            
            slice = nullTs(:,ii); 
            for ti = 1:size(peakHit,2)
                curdat = [peakHit(:,ti); peakMiss(:,ti)];
                
           
                


                 modDat = table(curdat, categorical(hmShuff), chanID, ...
                     trialID, subID, realID, ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', ...
                'trial', 'sub', 'realID'}); 
           
            lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan:sub)  + (1|trial)');
             slice(ti) = lme.Coefficients(2,4);
            end
            nullTs(:,ii) = slice; 
            
        end
        toc
        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        disp('calculated encoding permutation ')


       HFBdat.tVals_HFB = tVals; 
       HFBdat.hitVals_HFB = mean(hitVals); 
       HFBdat.missVals_HFB = mean(missVals);
       HFBdat.p_HFB = p; 

      
       statInfo = HFBdat; 
       save([statFiles(fileIdx).folder '/' 'stat' num2str(statType) statFiles(fileIdx).name], 'statInfo', '-v7.3');

        end

end