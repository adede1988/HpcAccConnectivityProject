function [] = lowFreqConPipeline_twoVar(cndFiles, idx, outprefix)


connectionDat = load([cndFiles(idx).folder '/' cndFiles(idx).name]).connectionDat;
disp(cndFiles(idx).name)  

%using raw PPC data
lowDat = connectionDat.lowBand(:, :,2); %raw ppc
highDat = connectionDat.highBand(:, :,2); %raw ppc

% lowDat = lowDat - min(lowDat, [], 2); 
% highDat = highDat - min(highDat, [], 2); 
% 
% lowDat = lowDat ./ max(lowDat, [], 2); 
% highDat = highDat ./ max(highDat, [], 2); 
% 
% %log transform to achieve normality 
% lowDat = log10((lowDat - min(lowDat,[],'all')) + .1);
% highDat = log10((highDat - min(highDat,[], 'all')) + .1); 

%rank transform the d' data
%rank transform the d' data
if ~isfield(connectionDat, 'dOrig')
    connectionDat.dOrig = connectionDat.d; 
    connectionDat.subDOrig = connectionDat.subD; 
else
    connectionDat.d = connectionDat.dOrig; 
    connectionDat.subD = connectionDat.subDOrig; 
end
[dVals, rankd] = sort(connectionDat.subD); 
dVals = dVals(rankd); 
% rankd = norminv((rankd - .5)/length(rankd) ); %rankit transform
% if ~isempty(find(diff(dVals)==0))
%     connectionDat.repeatedDvals = 'yes'; 
% else
%     connectionDat.repeatedDvals = 'no'; 
% end
for sub = 1:length(rankd)
    connectionDat.d(connectionDat.dOrig==dVals(sub)) = rankd(sub); 
    connectionDat.subD(connectionDat.subDOrig==dVals(sub)) = rankd(sub); 
end
connectionDat.subD = norminv((connectionDat.subD -.5) / length(connectionDat.subD)); 
connectionDat.d = norminv((connectionDat.d - .5) / length(connectionDat.subD)); 



% double up the allSubs and d variables
connectionDat.d = [connectionDat.d, connectionDat.d]; 
connectionDat.allSubs = [connectionDat.allSubs, connectionDat.allSubs]; 


%t vals for hitMiss, memory, and interaction
%order of effects: 
%1: hit/miss
%2: memory (rank d')
%3: interaction
lowtVals = zeros(size(lowDat,2),3); 
hightVals = lowtVals; 

%do the initial model fit: 
%connectivity ~ memory + (1|subject) 
% where n is electrode pair nested in subject
for ti = 1:size(lowDat, 2)
    [conVals, rankCon] = sort(lowDat(:,ti)); 
    ranks = zeros(length(rankCon),1); 
    for ri = 1:length(rankCon)
        ranks(ri) = find(rankCon==ri); 
    end


    modDat = table(norminv((ranks - .5) / length(ranks)), connectionDat.hmSort, connectionDat.d', connectionDat.allSubs', ...
        'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory * hitMiss +  (1|sub)'); 
    lowtVals(ti,:) = lme.Coefficients(2:4,4); 
  
    [conVals, rankCon] = sort(highDat(:,ti)); 
    ranks = zeros(length(rankCon),1); 
    for ri = 1:length(rankCon)
        ranks(ri) = find(rankCon==ri); 
    end
   
    modDat = table(norminv((ranks - .5) / length(ranks)), connectionDat.hmSort,connectionDat.d', connectionDat.allSubs', ...
        'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
    hightVals(ti,:) = lme.Coefficients(2:4,4);

end

connectionDat.lowtVals = lowtVals; 
connectionDat.hightVals = hightVals; 

disp('raw t values calculated')


%do permutation test, refitting the model 1000 times shuffling the memory
%performance values randomly each time keeping the data such that each
%subject has one value for memory performance over all their electrodes
perms = 1000; 
lownullTs = zeros([size(hightVals), perms]); 
highnullTs = lownullTs; 

subidx = cellfun(@(y) cellfun(@(x) strcmp(x,y), connectionDat.allSubs), ...
                 connectionDat.uniqueSubs, 'UniformOutput', false); 

for ii = 1:perms
%     if mod(ii, 100)==0
%         disp(['...........................' num2str(ii) 'permutations complete'])
%     end
    shuffd = zeros(size(connectionDat.d));
    shuffHM = zeros(size(connectionDat.hmSort)); 
    shuffVals = randsample(connectionDat.subD, length(connectionDat.subD), false); 
    for sub = 1:length(subidx)
        shuffd(subidx{sub}) = shuffVals(sub);
        tmpHM = zeros(sum(subidx{sub}),1); 
        tmpHM(1:sum(subidx{sub})/2) = 1; 
        tmpHM = tmpHM==1; 
        shuffHM(subidx{sub}) = randsample(tmpHM, length(tmpHM), false); 
    end
    sliceT = lownullTs(:,:,ii); 
    sliceT2 = highnullTs(:,:,ii); 
    for ti = 1:size(lowDat,2)
        [conVals, rankCon] = sort(lowDat(:,ti)); 
        ranks = zeros(length(rankCon),1); 
        for ri = 1:length(rankCon)
            ranks(ri) = find(rankCon==ri); 
        end
    
    
        modDat = table(norminv((ranks - .5) / length(ranks)), shuffHM, shuffd', connectionDat.allSubs',...
            'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
        sliceT(ti,:) = lme.Coefficients(2:4,4); 

        [conVals, rankCon] = sort(highDat(:,ti)); 
        ranks = zeros(length(rankCon),1); 
        for ri = 1:length(rankCon)
            ranks(ri) = find(rankCon==ri); 
        end
       
        modDat = table(norminv((ranks - .5) / length(ranks)), shuffHM, shuffd', connectionDat.allSubs', ...
            'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'});  
        lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
        sliceT2(ti,:) = lme.Coefficients(2:4,4); 
    end
    lownullTs(:,:,ii) = sliceT; 
    highnullTs(:,:,ii) = sliceT2; 
end

connectionDat.low975 = prctile(lownullTs, 97.5, 3); 
connectionDat.low025 = prctile(lownullTs, 2.5, 3); 
connectionDat.high975 = prctile(highnullTs, 97.5, 3); 
connectionDat.high025 = prctile(highnullTs, 2.5, 3);
lowp = zeros(size(connectionDat.lowtVals)); 
highp = lowp; 
[~,lowp(:,1), ~] = cluster_test(lowtVals(:,1), lownullTs(:,1,:)); 
[~,lowp(:,2), ~] = cluster_test(lowtVals(:,1), lownullTs(:,2,:));
[~,lowp(:,3), ~] = cluster_test(lowtVals(:,1), lownullTs(:,3,:));
[~,highp(:,1), ~] = cluster_test(hightVals(:,1), highnullTs(:,1,:)); 
[~,highp(:,2), ~] = cluster_test(hightVals(:,1), highnullTs(:,2,:)); 
[~,highp(:,3), ~] = cluster_test(hightVals(:,1), highnullTs(:,3,:)); 

connectionDat.lowp = lowp; 
connectionDat.highp = highp;
save([cndFiles(idx).folder '/'  outprefix '_' cndFiles(idx).name], 'connectionDat')
disp('permutation on full set complete')

%cross validate using leave one out 
if (sum(connectionDat.lowp<.05, 'all') > 0) || (sum(connectionDat.highp<.05,'all') > 0)
    %low / high X subjectLeftOut X time X hit/miss; memory; interaction
    crossP = ones([2, connectionDat.subN, length(connectionDat.tim), 3]); 
    if ~isfield(connectionDat, 'lowCrossT')
        connectionDat.lowCrossT = zeros(connectionDat.subN, length(connectionDat.tim), 3); 
        connectionDat.highCrossT = zeros(connectionDat.subN, length(connectionDat.tim), 3); 
        connectionDat.lowCrossP = squeeze(crossP(1,:,:,:)); 
        connectionDat.highCrossP = squeeze(crossP(2,:,:,:)); 
    end
    %try leaving out each subject one at a time to confirm if results
    %are the same
    for subOut = 1:length(connectionDat.uniqueSubs)
        %check if this subject is already done
        if sum(connectionDat.lowCrossT(subOut,:,1)==0) == length(connectionDat.tim) %all zeros! need to calculate
        disp(['doing stats without: ' num2str(subOut) ' of ' num2str(length(connectionDat.uniqueSubs))])
        subSelect = cellfun(@(x) ~strcmp(x, connectionDat.uniqueSubs{subOut}), connectionDat.allSubs); 
        cur = lowDat(subSelect,:);
        cur2 = highDat(subSelect,:); 
        curd = connectionDat.d(subSelect); 
        curhm = connectionDat.hmSort(subSelect); 
        cursub = connectionDat.allSubs(subSelect); 
        tmpT = zeros(size(lowDat,2),3); %low frequency
        tmpT2 = tmpT; %high frequency
        for ti = 1:size(lowDat,2)
            [conVals, rankCon] = sort(cur(:,ti)); 
            ranks = zeros(length(rankCon),1); 
            for ri = 1:length(rankCon)
                ranks(ri) = find(rankCon==ri); 
            end
        
        
            modDat = table(norminv((ranks - .5) / length(ranks)), curhm, curd', cursub', ...
                'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
            lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
            tmpT(ti,:) = lme.Coefficients(2:4,4); %t-value associated with memory!

            [conVals, rankCon] = sort(cur2(:,ti)); 
            ranks = zeros(length(rankCon),1); 
            for ri = 1:length(rankCon)
                ranks(ri) = find(rankCon==ri); 
            end
        
        
            modDat = table(norminv((ranks - .5) / length(ranks)), curhm, curd', cursub', ...
                'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
            lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
            tmpT2(ti,:) = lme.Coefficients(2:4,4); %t-value associated with memory!
        end

        connectionDat.lowCrossT(subOut,:,:) = tmpT; 
        connectionDat.highCrossT(subOut,:,:) = tmpT2; 
        
        lownullTs = zeros([size(tmpT), perms]); 
        highnullTs = lownullTs;  
        curUSubs = connectionDat.uniqueSubs;
        curUSubs(subOut) = [];
        cursubD = connectionDat.subD; 
        cursubD(subOut) = []; 
        subidx = cellfun(@(y) cellfun(@(x) strcmp(x,y), cursub), ...
                 curUSubs, 'UniformOutput', false); 
        
        for ii = 1:perms
%             if mod(ii, 100)==0
%                 disp(['...........................' num2str(ii) 'permutations complete'])
%             end
            shuffd = zeros(length(cursub),1);
            shuffHM = shuffd; 
            shuffVals = randsample(cursubD, length(cursubD), false); 
            for sub = 1:length(subidx)
                shuffd(subidx{sub}) = shuffVals(sub); 
                tmpHM = zeros(sum(subidx{sub}),1); 
                tmpHM(1:sum(subidx{sub})/2) = 1; 
                tmpHM = tmpHM==1; 
                shuffHM(subidx{sub}) = randsample(tmpHM, length(tmpHM), false); 
            end
            slice1 = lownullTs(:,:,ii); 
            slice2 = highnullTs(:,:,ii); 
            for ti = 1:size(lowDat,2)
                [conVals, rankCon] = sort(cur(:,ti)); 
                ranks = zeros(length(rankCon),1); 
                for ri = 1:length(rankCon)
                    ranks(ri) = find(rankCon==ri); 
                end
            
            
                modDat = table(norminv((ranks - .5) / length(ranks)), shuffHM, shuffd, cursub',...
                    'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
                lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
                slice1(ti,:) = lme.Coefficients(2:4,4); %t-value associated with memory!
                
                [conVals, rankCon] = sort(cur2(:,ti)); 
                ranks = zeros(length(rankCon),1); 
                for ri = 1:length(rankCon)
                    ranks(ri) = find(rankCon==ri); 
                end
            
            
                modDat = table(norminv((ranks - .5) / length(ranks)), shuffHM, shuffd, cursub',...
                    'VariableNames', {'connectivity', 'hitMiss', 'memory', 'sub'}); 
                lme = fitlme(modDat, 'connectivity ~ memory * hitMiss + (1|sub)'); 
                slice2(ti,:) = lme.Coefficients(2:4,4); %t-value associated with memory!
            end
            lownullTs(:,:, ii) = slice1; 
            highnullTs(:,:, ii) = slice2; 

        end
      
        %low frequency hitMiss
        [~, connectionDat.lowCrossP(subOut,:,1), ~] = cluster_test(tmpT(:,1), lownullTs(:,1,:)); 
        %high frequency hitMiss
        [~, connectionDat.highCrossP(subOut,:,1), ~] = cluster_test(tmpT2(:,1), highnullTs(:,1,:)); 
        %low frequency memory
        [~, connectionDat.lowCrossP(subOut,:,2), ~] = cluster_test(tmpT(:,2), lownullTs(:,2,:)); 
        %high frequency memory
        [~, connectionDat.highCrossP(subOut,:,2), ~] = cluster_test(tmpT2(:,2), highnullTs(:,2,:));
        %low frequency interaction
        [~, connectionDat.lowCrossP(subOut,:,3), ~] = cluster_test(tmpT(:,3), lownullTs(:,3,:)); 
        %high frequency interaction
        [~, connectionDat.highCrossP(subOut,:,3), ~] = cluster_test(tmpT2(:,3), highnullTs(:,3,:));
        save([cndFiles(idx).folder '/' outprefix '_' cndFiles(idx).name], 'connectionDat')
        end
    end
else
    connectionDat.lowCrossP = ones([connectionDat.subN, length(connectionDat.tim),3]); 
    connectionDat.highCrossP = ones([connectionDat.subN, length(connectionDat.tim),3]); 
end

save([cndFiles(idx).folder '/' outprefix '_' cndFiles(idx).name], 'connectionDat')

















end