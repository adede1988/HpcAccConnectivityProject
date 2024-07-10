function [] = connectivitypipeline(statFiles, reg1, reg2, encRet,...
    statType, permi, datPre)
    rng(permi)
    perms = 50; 
    frex = logspace(log10(2), log10(25), 20);
    regions = {'acc', 'dlPFC', 'hip', 'lTemp', 'iTemp', 'mtl', ...
    'pcc', 'pPFC', 'vis'}; 
    phase = {'enc' , 'ret'}; 
    %% get data from single channel files
    %initialize 
    hitVals = cell(length(statFiles), 1); 
    missVals = cell(length(statFiles),1); 

    chiVals = zeros(length(statFiles),2);
    subVals = cell(length(statFiles),1); 
    %loop pair files 
    for ii = 1: length(statFiles)
        ii
        pairPaths = load([statFiles(ii).folder '/' ...
            statFiles(ii).name]).outVar;
        %check it's not the same channel to itself! 
        if ~strcmp(pairPaths.setName, pairPaths.pairName)
            %load the chanDat
            setChan = load([datPre 'CHANDAT/finished/' ...
                pairPaths.setName]).chanDat; 
%                 pairChan = load([pairPaths.pairFolder '/' ...
%                     pairPaths.pairName]).chanDat; 
            %extract channel index values
            set_chi = split(pairPaths.setName, '_');
            subVals{ii} = set_chi{2}; 
            set_chi = split(set_chi{3}, '.');
            set_chi = str2num(set_chi{1});
            pair_chi = split(pairPaths.pairName, '_');
            pair_chi = split(pair_chi{3}, '.');
            pair_chi = str2num(pair_chi{1});

            chiVals(ii, 1) = set_chi; 
            chiVals(ii, 2) = pair_chi; 
            %get the target connectivity data 
            %encRet: 1 = subsequent, 2 = retrieval
            %statType: 1 = image lock, 0 = HFB lock
            if encRet == 1 && statType == 1
                tim = setChan.HFB.encMulTim; 
                missVals{ii} = squeeze(...
                    setChan.ISPC.subMiss(pair_chi, :, :, 2));
                hitVals{ii} = squeeze(...
                    setChan.ISPC.subHit(pair_chi, :, :, 2));
            elseif encRet == 1 && statType == 0 
                tim = [-500:25:500];
                missVals{ii} = squeeze(...
                    setChan.ISPC.subMiss(pair_chi, 1:41, :, 4));
                hitVals{ii} = squeeze(...
                    setChan.ISPC.subHit(pair_chi, 1:41, :, 4));
            elseif encRet == 2 && statType == 1 
                tim = setChan.HFB.onMulTim; 
                missVals{ii} = squeeze(...
                    setChan.ISPC.miss_on(pair_chi, :, :, 2));
                hitVals{ii} = squeeze(...
                    setChan.ISPC.hit_on(pair_chi, :, :, 2));
            elseif encRet == 2 && statType == 0 
                tim = [-500:25:500];
                missVals{ii} = squeeze(...
                    setChan.ISPC.miss_on(pair_chi, 1:41, :, 4));
                hitVals{ii} = squeeze(...
                    setChan.ISPC.hit_on(pair_chi, 1:41, :, 4));
            end


        end


    end
    %eliminate empties caused by self pairs
    hitVals(cellfun(@(x) isempty(x), hitVals(:,1)), :) = []; 
    missVals(cellfun(@(x) isempty(x), missVals(:,1)), :) = [];
    subVals(cellfun(@(x) isempty(x), subVals(:,1)), :) = [];
    chiVals(chiVals(:,1)==0, :) = []; 

    %concatenate data 
    hitVals = cat(3, hitVals{:}); 
    missVals = cat(3, missVals{:}); 

    %if this is a within region connectivity measure AND it's image locked,
    %then eliminate double counts
    pairIDs = cell(length(subVals), 1); 
    eliminate = []; 
    for ii = 1:length(subVals)
        if ii > 1
            allPairSoFar = unique(pairIDs(1:ii-1)); 
        else 
            allPairSoFar = {'aksj'}; 
        end
        if chiVals(ii, 1) < chiVals(ii, 2)
            pairIDs{ii} = [subVals{ii} '_' num2str(chiVals(ii,1)) '_' ...
                            num2str(chiVals(ii,2))];
        elseif chiVals(ii,1) > chiVals(ii, 2)
            pairIDs{ii} = [subVals{ii} '_' num2str(chiVals(ii,2)) '_' ...
                            num2str(chiVals(ii,1))];
        else
            disp('self connection made it through!')
            alskdj %error here
        end


        if sum(cellfun(@(x) strcmp(pairIDs{ii}, x), allPairSoFar))>0
            %this pair is already done! 
            eliminate = [eliminate, ii]; 
        end




    end

    if ~isempty(eliminate) %get rid of repeated pairs! 
        hitVals(:,:,eliminate) = []; 
        missVals(:,:,eliminate) = []; 
        subVals(eliminate) = []; 
        chiVals(eliminate, :) = []; 
        pairIDs(eliminate) = []; 

    end

    hmSort = [ones(length(pairIDs),1); zeros(length(pairIDs),1)]; 

    clear setChan
    
    %% do stats! 

    disp('ready to calcualte observed values')


    %empty matrix that is timepoints X frequency
    tVals = zeros(size(hitVals,[1,2]));
    
    %loop on timepoints
    for ti = 1:size(hitVals,1)
        disp(['observed time: ' num2str(ti)])
        slice = tVals(ti,:); 
        for fi = 1:size(hitVals,2)

        curdat = [squeeze(hitVals(ti,fi,:));...
                  squeeze(missVals(ti, fi, :))];

        modDat = table(curdat, categorical(hmSort), ...
            [pairIDs; pairIDs], [subVals; subVals], ...
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
            for chan = 1:length(pairIDs)
                curi = find(ismember([pairIDs; pairIDs], pairIDs{chan}));
              

                if rand() > .5
                    hmShuff(curi(1)) = 1; 
                else
                    hmShuff(curi(2)) = 1; 
                end



            end


            
            curSet = nullTs(:,:,ii); 
            %loop on timepoints
            for ti = 1:size(hitVals,1)
                slice = curSet(ti,:); 
                for fi = 1:size(hitVals,2)
        
                curdat = [squeeze(hitVals(ti,fi,:));...
                          squeeze(missVals(ti, fi, :))];
        
                modDat = table(curdat, categorical(hmShuff), ...
                    [pairIDs; pairIDs], [subVals; subVals], ...
                    'VariableNames', {'HFB', 'hitMiss', 'chan', 'sub'}); 
               
                lme = fitlme(modDat, ...
                        'HFB ~ hitMiss + (1|chan)');
                slice(fi) = lme.Coefficients(2,4); 
                end
                curSet(ti,:) = slice; 
            end 
            nullTs(:,:,ii) = curSet;

        end
        toc

        
        
       outDat = struct;

       outDat.tVals = tVals; 
       %store out all the hit and miss vals by pair once
       if permi > 1
       outDat.hitVals = mean(hitVals,3); 
       outDat.missVals = mean(missVals,3);
       else
       outDat.hitVals = hitVals; 
       outDat.missVals = missVals;
       end
       outDat.pairIDs = pairIDs;  
       outDat.subVals = subVals; 
       outDat.chiVals = chiVals; 
       outDat.reg1 = reg1; 
       outDat.reg2 = reg2; 
       outDat.statType = statType; 
       outDat.tim = tim; 
       outDat.frex = frex; 
       outDat.nulls = nullTs; 
 
       save([statFiles(1).folder '/out/'...
           regions{reg1} '_' regions{reg2} '_' ...
           phase{encRet} '_' 'stat' num2str(statType) '_' ...
           num2str(permi) '.mat'], ...
           'outDat', '-v7.3');
        
       

end