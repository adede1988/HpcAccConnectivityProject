function [] = TFsingleTrialpipeline(statFiles, fileIdx, statType, permi)

        HFBdat = load([statFiles(fileIdx).folder '/' statFiles(fileIdx).name]).statInfo; 
        perms = 100; 
    
        %% encoding mean difference

        %matrices of channels X timepoints for hits and misses separately
        hitVals = HFBdat.hits; 
        missVals = HFBdat.misses;
        tim = HFBdat.tim; 
      
        hitVals(tim<-450 | tim>3000, :, :) = []; 
        missVals(tim<-450 | tim>3000, :, :) = []; 
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
        %chan X time X frequency
        hitVals = zeros([length(chanUni), ...
            size(hitTrials,1), size(hitTrials,3)]); 
        missVals = zeros([length(chanUni), ...
            size(hitTrials,1), size(hitTrials,3)]); 
        sub_out = cell(size(chanUni)); 


        eliminate = []; 

        for chan = 1:length(chanUni)

            hitidx = ismember(chanID, chanUni(chan));
            missidx = ismember(chanID2, chanUni(chan)); 
            
                hitVals(chan, :,:) = mean(hitTrials(:, hitidx,:), 2); 
                missVals(chan, :,:) = mean(missTrials(:,missidx, :), 2); 
                tmp = subID(hitidx); 
                sub_out{chan} = tmp{1}; 
            if sum(hitidx)>0 && sum(missidx)>0
            else
                eliminate = [eliminate, chan]; 
            end
        end



        hmSort = [ones(length(chanUni),1); zeros(length(chanUni),1)]; 




         %empty vector that is timepoints X frequency
        tVals = zeros(size(hitVals,[2,3]));
    
        %loop on timepoints
        for ti = 1:size(hitVals,2)
            slice = tVals(ti,:); 
            for fi = 1:100

            curdat = [squeeze(hitVals(:,ti,fi));...
                      squeeze(missVals(:, ti, fi))];

            modDat = table(curdat, categorical(hmSort), ...
                [chanUni; chanUni], [sub_out; sub_out], ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
           
            lme = fitlme(modDat, ...
                    'HFB ~ hitMiss + (1|chan)');
            slice(fi) = lme.Coefficients(2,4); 
            end
            tVals(ti,:) = slice; 
        end 
        disp('calculated encoding t-vals ')
       

       
        
      


           tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            ii
            if(mod(ii, 10)) ==0 
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


            
            curSet = nullTs(:,:,ii); 
                    %loop on timepoints
            for ti = 1:size(hitVals,2)
                slice = curSet(ti,:); 
                for fi = 1:100
    
                curdat = [squeeze(hitVals(:,ti,fi));...
                          squeeze(missVals(:, ti, fi))];
    
                modDat = table(curdat, categorical(hmShuff), ...
                    [chanUni; chanUni], [sub_out; sub_out], ...
                    'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
               
                lme = fitlme(modDat, ...
                        'HFB ~ hitMiss + (1|chan)');
                slice(fi) = lme.Coefficients(2,4); 
                end
                curSet(ti,:) = slice; 
            end 
            nullTs(:,:,ii) = curSet;

        end


        
        
%         [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
%         disp('calculated encoding permutation ')
       outDat = struct;

       outDat.tVals = tVals; 
       outDat.hitVals = mean(hitVals); 
       outDat.missVals = mean(missVals);
%        HFBdat.p_image = p; 
       outDat.eliminate = eliminate; 
       outDat.nulls = nullTs; 
 
       save([statFiles(fileIdx).folder '/out/'...
           'stat' num2str(statType) '_' num2str(permi) '_' statFiles(fileIdx).name], ...
           'outDat', '-v7.3');
        
        else
        
        %% HFB peak aligned analysis
        hitVals = HFBdat.hits; 
        missVals = HFBdat.misses;
        tim = HFBdat.tim; 
      
%         hitVals(tim<-450 | tim>3000, :, :) = []; 
%         missVals(tim<-450 | tim>3000, :, :) = []; 
%         tim(tim<-450 | tim>3000) = []; 

        w = 20; 

        hitLatIdx = arrayfun(@(x) find(x==tim), HFBdat.hitLat); 
        hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
        hitLatIdx(hitLatIdx<=w) =  w+1; 

        peakHit = arrayfun(@(x) squeeze(hitVals(hitLatIdx(x)-w:...
                                         hitLatIdx(x)+w, x, :)),  [1:length(hitLatIdx)], ...
                                         'uniformoutput', false);
        peakHit = cat(3, peakHit{:});

        missLatIdx = arrayfun(@(x) find(x==tim), HFBdat.missLat); 
        missLatIdx(missLatIdx>length(tim)-w) = length(tim) - w; 
        missLatIdx(missLatIdx<=w) = w+1; 

        peakMiss = arrayfun(@(x) squeeze(missVals(missLatIdx(x)-w:...
                                         missLatIdx(x)+w, x, :)), [1:length(missLatIdx)], ...
                                         'uniformoutput', false);
        peakMiss = cat(3, peakMiss{:});





        %average across trials within channels
        hitTrials = peakHit; 
        missTrials = peakMiss; 
        %chan X time X frequency
        hitVals = zeros([length(chanUni), 41, size(hitTrials,2)]); 
        missVals = zeros([length(chanUni), 41, size(hitTrials,2)]); 
        sub_out = cell(size(chanUni)); 


        eliminate = []; 

        for chan = 1:length(chanUni)

            hitidx = ismember(chanID, chanUni(chan));
            missidx = ismember(chanID2, chanUni(chan)); 
            
                hitVals(chan, :,:) = mean(hitTrials(:, :,hitidx), 3); 
                missVals(chan, :,:) = mean(missTrials(:, :,missidx), 3); 
                tmp = subID(hitidx); 
                sub_out{chan} = tmp{1}; 
            if sum(hitidx)>0 && sum(missidx)>0
            else
                eliminate = [eliminate, chan]; 
            end
        end

        hmSort = [ones(length(chanUni),1); zeros(length(chanUni),1)]; 





        %empty vector that is 1Xtimepoints for t values
        tVals = zeros(size(hitVals,[2,3]));
    
        %loop on timepoints
        for ti = 1:size(hitVals,2)
            slice = tVals(ti,:); 
            for fi = 1:100

            curdat = [squeeze(hitVals(:,ti,fi));...
                      squeeze(missVals(:, ti, fi))];

            modDat = table(curdat, categorical(hmSort), ...
                [chanUni; chanUni], [sub_out; sub_out], ...
                'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
           
            lme = fitlme(modDat, ...
                    'HFB ~ hitMiss + (1|chan)');
            slice(fi) = lme.Coefficients(2,4); 
            end
            tVals(ti,:) = slice; 
        end 
        disp('calculated encoding t-vals ')







           tic
        nullTs = squeeze(zeros([size(tVals), perms])); 
        for ii = 1:perms
            ii
            if(mod(ii, 10)) ==0 
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


            
            curSet = nullTs(:,:,ii); 
                    %loop on timepoints
            for ti = 1:size(hitVals,2)
                slice = curSet(ti,:); 
                for fi = 1:100
    
                curdat = [squeeze(hitVals(:,ti,fi));...
                          squeeze(missVals(:, ti, fi))];
    
                modDat = table(curdat, categorical(hmShuff), ...
                    [chanUni; chanUni], [sub_out; sub_out], ...
                    'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
               
                lme = fitlme(modDat, ...
                        'HFB ~ hitMiss + (1|chan)');
                slice(fi) = lme.Coefficients(2,4); 
                end
                curSet(ti,:) = slice; 
            end 
            nullTs(:,:,ii) = curSet;

        end


        
%         [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
        disp('calculated encoding permutation ')


       outDat = struct;

       outDat.tVals = tVals; 
       outDat.hitVals = mean(hitVals); 
       outDat.missVals = mean(missVals);
%        HFBdat.p_image = p; 
       outDat.eliminate = eliminate; 
       outDat.nulls = nullTs; 
 
       save([statFiles(fileIdx).folder '/out/'...
           'stat' num2str(statType) '_' permi '_' statFiles(fileIdx).name], ...
           'outDat', '-v7.3');
        end

end