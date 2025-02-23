function [] = HFBsingleTrialpipeline(statFiles, fileIdx, statType, permi)
        rng(permi)
        HFBdat = load([statFiles(fileIdx).folder '/' statFiles(fileIdx).name]).statInfo; 
        perms = 50; 
    
        %% encoding mean difference

        %matrices of channels X timepoints for hits and misses separately
        hitVals = HFBdat.hits'; 
        missVals = HFBdat.misses';
        tim = HFBdat.tim; 
      
        hitVals(:, tim<-450 | tim>3000) = []; 
        missVals(:, tim<-450 | tim>3000) = []; 
        tim(tim<-450 | tim>3000) = []; 
        

      

        chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.hitSub, ...
             num2cell(HFBdat.hitChi), 'UniformOutput',false);

        chanID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.missSub, ...
             num2cell(HFBdat.missChi), 'UniformOutput',false);

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
            
                hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
                missVals(chan, :) = mean(missTrials(missidx, :), 1); 
                tmp = subID(hitidx); 
                sub_out{chan} = tmp{1}; 
            if sum(hitidx)>0 && sum(missidx)>0
            else
                eliminate = [eliminate, chan]; 
            end
        end
       

        hmSort = [ones(length(chanUni),1); zeros(length(chanUni),1)]; 




        %empty vector that is 1Xtimepoints for t values
        tVals = zeros(size(hitVals,2),1);
    
        %loop on timepoints
        for ti = 1:size(hitVals,2)
            try

            curdat = [hitVals(:,ti); missVals(:,ti)];
    
            modDat = table(curdat, categorical(hmSort), ...
                [chanUni; chanUni], [sub_out; sub_out], ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
           
            lme = fitlme(modDat, ...
                'HFB ~ hitMiss + (1|chan:sub) ');
            tVals(ti) = lme.Coefficients(2,4); 
            catch 
            end
        end 
        disp('calculated encoding t-vals ')
       

       
        
        tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            ii
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
                try
                curdat = [hitVals(:,ti); missVals(:,ti)]; 
         
    
                modDat = table(curdat, categorical(hmShuff), ...
                    [chanUni; chanUni], [sub_out; sub_out], ...
                    'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
               
                lme = fitlme(modDat, ...
                    'HFB ~ hitMiss + (1|chan:sub)');
                slice(ti) = lme.Coefficients(2,4); 
                catch
                end
                     
            end
            nullTs(:,ii) = slice; 
        end
        


        outDat = struct;

       outDat.tVals = tVals; 
       outDat.hitVals = hitVals; 
       outDat.missVals = missVals;
       outDat.eliminate = eliminate; 
       outDat.nulls = nullTs; 
 
       save([statFiles(fileIdx).folder '/out/'...
           'stat' num2str(statType) '_' num2str(permi) '_' statFiles(fileIdx).name], ...
           'outDat', '-v7.3');

  
        else
        
        %% HFB peak aligned analysis
        hitVals = HFBdat.hits'; 
        missVals = HFBdat.misses';
        tim = HFBdat.tim; 

        w = 20; 

        hitLatIdx = arrayfun(@(x) find(x==tim), HFBdat.hitLat); 
        hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
        hitLatIdx(hitLatIdx<=w) =  w+1; 

        peakHit = arrayfun(@(x) hitVals(x,hitLatIdx(x)-w:...
                                         hitLatIdx(x)+w), [1:length(hitLatIdx)], ...
                                         'uniformoutput', false);
        peakHit = cell2mat(peakHit.');

        missLatIdx = arrayfun(@(x) find(x==tim), HFBdat.missLat); 
        missLatIdx(missLatIdx>length(tim)-w) = length(tim) - w; 
        missLatIdx(missLatIdx<=w) = w+1; 

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
            
                hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
                missVals(chan, :) = mean(missTrials(missidx, :), 1); 
                tmp = subID(hitidx); 
                sub_out{chan} = tmp{1}; 
            if sum(hitidx)>0 && sum(missidx)>0
            else
                eliminate = [eliminate, chan]; 
            end
        end
        
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
                    'HFB ~ hitMiss + (1|chan:sub)');
            tVals(ti) = lme.Coefficients(2,4); 
    
        end 
        disp('calculated encoding t-vals ')







           tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            ii
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
                    'HFB ~ hitMiss + (1|chan:sub)');
                slice(ti) = lme.Coefficients(2,4); 
                     
            end
            nullTs(:,ii) = slice; 
        end
        

       outDat = struct;

       outDat.tVals = tVals; 
       outDat.hitVals = hitVals; 
       outDat.missVals = missVals;
       outDat.eliminate = eliminate; 
       outDat.nulls = nullTs; 
 
       save([statFiles(fileIdx).folder '/out/'...
           'stat' num2str(statType) '_' num2str(permi) '_' statFiles(fileIdx).name], ...
           'outDat', '-v7.3');

        end

end