function [] = HFBsingleTrialpipeline(statFiles, fileIdx, statType)

        HFBdat = load([statFiles(fileIdx).folder '/' statFiles(fileIdx).name]).statInfo; 
        perms = 200; 
    
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
%         hmSort = [ones(size(hitVals,1),1); zeros(size(missVals,1),1)];
%         hmSort = hmSort==1; 
        
        %random effects variables
%         realID = cellfun(@(x, y, z) [x,'_', num2str(y),'_', num2str(z)], HFBdat.hitSub, ...
%             num2cell(HFBdat.hitTriali), num2cell(HFBdat.hitChi), 'UniformOutput',false);
% 
%         realID2 = cellfun(@(x, y, z) [x,'_', num2str(y),'_', num2str(z)], HFBdat.missSub, ...
%             num2cell(HFBdat.missTriali), num2cell(HFBdat.missChi), 'UniformOutput',false);
% 
%         realID = [realID; realID2]; 

        chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.hitSub, ...
             num2cell(HFBdat.hitChi), 'UniformOutput',false);

        chanID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.missSub, ...
             num2cell(HFBdat.missChi), 'UniformOutput',false);
% 
%         trialID = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.hitSub, ...
%              num2cell(HFBdat.hitTriali), 'UniformOutput',false);
% 
%         trialID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.missSub, ...
%              num2cell(HFBdat.missTriali), 'UniformOutput',false);


%         chanID = [chanID; chanID2]; 
%         trialID = [trialID; trialID2]; 
        subID = [HFBdat.hitSub; HFBdat.missSub]; 
        sub_uni = unique(subID); 
        chanUni = unique(chanID); 

        



        if statType == 1


        %average across trials within channels
        hitTrials = hitVals; 
        missTrials = missVals; 
        hitVals = zeros(length(chanUni), length(tim)); 
        missVals = zeros(length(chanUni), length(tim)); 
        sub_out = cell(size(chanUni)); 

        eliminate = []; 

        for chan = 1:length(chanUni)

            hitidx = ismember(chanID, chanUni(chan));
            missidx = ismember(chanID2, chanUni(chan)); 
            if sum(hitidx)>0 && sum(missidx)>0
                hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
                missVals(chan, :) = mean(missTrials(missidx, :), 1); 
                tmp = subID(hitidx); 
                sub_out{chan} = tmp{1}; 
            else
                eliminate = [eliminate, chan]; 
            end
        end
       
        hitVals(eliminate,:) = []; 
        missVals(eliminate, :) = []; 
        sub_out(eliminate) = []; 
        chanUni(eliminate) = []; 
        hmSort = [ones(length(chanUni),1); zeros(length(chanUni),1)]; 




        %empty vector that is 1Xtimepoints for t values
        tVals = zeros(size(hitVals,2),1);
    
        %loop on timepoints
        for ti = 1:size(hitVals,2)
    

            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, categorical(hmSort), ...
                [chanUni; chanUni], [sub_out; sub_out], ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
           
            lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan)');
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
            for chan = 1:length(chanUni)
                curi = find(ismember([chanUni; chanUni], chanUni{chan}));
                curHM = hmSort(curi);

                if rand() > .5
                    hmShuff(curi(1)) = 1; 
                else
                    hmShuff(curi(2)) = 1; 
                end

            end


            
            slice = nullTs(:,ii); 
            for ti = 1:size(hitVals,2)
                curdat = [hitVals(:,ti); missVals(:,ti)]; 
         
                curdat = [hitVals(:,ti); missVals(:,ti)];
    
                modDat = table(curdat, categorical(hmShuff), ...
                    [chanUni; chanUni], [sub_out; sub_out], ...
                    'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
               
                lme = fitlme(modDat, ...
                    'HFB ~ hitMiss + (1|chan)');
                slice(ti) = lme.Coefficients(2,4); 
                     
            end
            nullTs(:,ii) = slice; 
        end
        
        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        disp('calculated encoding permutation ')


       HFBdat.tVals_image = tVals; 
       HFBdat.hitVals_image = mean(hitVals); 
       HFBdat.missVals_image = mean(missVals);
       HFBdat.p_image = p; 

       statInfo = HFBdat; 
       save([statFiles(fileIdx).folder '/out/' 'stat' num2str(statType) '_' statFiles(fileIdx).name], 'statInfo', '-v7.3');
        
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





        %average across trials within channels
        hitTrials = peakHit; 
        missTrials = peakMiss; 
        hitVals = zeros(length(chanUni), 41); 
        missVals = zeros(length(chanUni), 41); 
        sub_out = cell(size(chanUni)); 


        eliminate = []; 

        for chan = 1:length(chanUni)

            hitidx = ismember(chanID, chanUni(chan));
            missidx = ismember(chanID2, chanUni(chan)); 
            if sum(hitidx)>0 && sum(missidx)>0
                hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
                missVals(chan, :) = mean(missTrials(missidx, :), 1); 
                tmp = subID(hitidx); 
                sub_out{chan} = tmp{1}; 
            else
                eliminate = [eliminate, chan]; 
            end
        end
       
        hitVals(eliminate,:) = []; 
        missVals(eliminate, :) = []; 
        sub_out(eliminate) = []; 
        chanUni(eliminate) = []; 
        hmSort = [ones(length(chanUni),1); zeros(length(chanUni),1)]; 

% 
%         for chan = 1:length(chanUni)
% 
%             hitidx = ismember(chanID, chanUni(chan));
%             missidx = ismember(chanID2, chanUni(chan)); 
%             hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
%             missVals(chan, :) = mean(missTrials(missidx, :), 1); 
%             tmp = subID(hitidx); 
%             sub_out{chan} = tmp{1}; 
% 
%         end
       




        %empty vector that is 1Xtimepoints for t values
        tVals = zeros(size(hitVals,2),1);
    
        %loop on timepoints
        for ti = 1:size(hitVals,2)
    

            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, categorical(hmSort), ...
                [chanUni; chanUni], [sub_out; sub_out], ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
           
            lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan)');
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
        disp('calculated encoding t-vals ')





        %% fix below here


           tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            if(mod(ii, 100)) ==0 
                disp(['..........................' num2str(ii) ...
                    ' time:' num2str(round(toc/60,1))])
            end

            hmShuff = zeros(length(hmSort),1); 
            for chan = 1:length(chanUni)
                curi = find(ismember([chanUni; chanUni], chanUni{chan}));
              

                if rand() > .5
                    hmShuff(curi(1)) = 1; 
                else
                    hmShuff(curi(2)) = 1; 
                end



            end


            
            slice = nullTs(:,ii); 
            for ti = 1:size(hitVals,2)
             
                curdat = [hitVals(:,ti); missVals(:,ti)];
    
                modDat = table(curdat, categorical(hmShuff), ...
                    [chanUni; chanUni], [sub_out; sub_out], ...
                    'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
               
                lme = fitlme(modDat, ...
                    'HFB ~ hitMiss + (1|chan)');
                slice(ti) = lme.Coefficients(2,4); 
                     
            end
            nullTs(:,ii) = slice; 
        end
        


        %% 



        
        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        disp('calculated encoding permutation ')


       HFBdat.tVals_HFB = tVals; 
       HFBdat.hitVals_HFB = mean(hitVals); 
       HFBdat.missVals_HFB = mean(missVals);
       HFBdat.p_HFB = p; 

      
       statInfo = HFBdat; 
       save([statFiles(fileIdx).folder '/out/' 'stat' num2str(statType) '_' statFiles(fileIdx).name], 'statInfo', '-v7.3');

        end

end