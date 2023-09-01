function [] = lowFreqConPipeline(cndFiles, idx)


connectionDat = load([cndFiles(idx).folder '/' cndFiles(idx).name]).connectionDat;
disp(cndFiles(idx).name)  

%using z-scored PPC data
lowDat = connectionDat.lowBand(connectionDat.hmSort, :,4); %raw ppc
highDat = connectionDat.highBand(connectionDat.hmSort, :,4); %raw ppc

%log transform to achieve normality 
lowDat = log10((lowDat - min(lowDat,[],'all')) + .1);
highDat = log10((highDat - min(highDat,[], 'all')) + .1); 

%rank transform the d' data
connectionDat.dOrig = connectionDat.d; 
connectionDat.subDOrig = connectionDat.subD; 
[dVals, rankd] = sort(connectionDat.subD); 
rankd = norminv((rankd - .5)/length(rankd) ); %rankit transform
if ~isempty(find(diff(dVals)==0))
    connectionDat.repeatedDvals = 'yes'; 
else
    connectionDat.repeatedDvals = 'no'; 
end
for sub = 1:length(rankd)
    connectionDat.d(connectionDat.d==connectionDat.subD(sub)) = rankd(sub); 
end
connectionDat.subD = rankd; 



lowtVals = zeros(size(lowDat,2),1);
hightVals = lowtVals; 

%do the initial model fit: 
%connectivity ~ memory + (1|subject) 
% where n is electrode pair nested in subject
parfor ti = 1:size(lowDat, 2)
   
    modDat = table(lowDat(:,ti), connectionDat.d', connectionDat.allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
    lowtVals(ti) = lme.Coefficients(2,4); %t-value associated with memory! 

   
    modDat = table(highDat(:,ti), connectionDat.d', connectionDat.allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
    hightVals(ti) = lme.Coefficients(2,4); %t-value associated with memory! 

end

connectionDat.lowtVals = lowtVals; 
connectionDat.hightVals = hightVals; 

disp('raw t values calculated')


%do permutation test, refitting the model 1000 times shuffling the memory
%performance values randomly each time keeping the data such that each
%subject has one value for memory performance over all their electrodes
perms = 1000; 
lownullTs = zeros([length(hightVals), perms]); 
highnullTs = lownullTs; 

subidx = cellfun(@(y) cellfun(@(x) strcmp(x,y), connectionDat.allSubs), ...
                 connectionDat.uniqueSubs, 'UniformOutput', false); 

parfor ii = 1:perms
    shuffd = zeros(size(connectionDat.d));
    shuffVals = randsample(connectionDat.subD, length(connectionDat.subD), false); 
    for sub = 1:length(subidx)
        shuffd(subidx{sub}) = shuffVals(sub); 
    end
    sliceT = lownullTs(:,ii); 
    sliceT2 = highnullTs(:,ii); 
    for ti = 1:size(lowDat,2)
        modDat = table(lowDat(:,ti), shuffd', connectionDat.allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
        sliceT(ti) = lme.Coefficients(2,4); %t-value associated with memory! 
        modDat = table(highDat(:,ti), shuffd', connectionDat.allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
        sliceT2(ti) = lme.Coefficients(2,4); %t-value associated with memory! 
    end
    lownullTs(:,ii) = sliceT; 
    highnullTs(:,ii) = sliceT2; 
end

connectionDat.low975 = prctile(lownullTs, 97.5, 2); 
connectionDat.low025 = prctile(lownullTs, 2.5, 2); 
connectionDat.high975 = prctile(highnullTs, 97.5, 2); 
connectionDat.high025 = prctile(highnullTs, 2.5, 2);
[~,connectionDat.lowp, ~] = cluster_test(lowtVals, lownullTs); 
[~,connectionDat.highp, ~] = cluster_test(hightVals, highnullTs); 

disp('permutation on full set complete')

if (sum(connectionDat.lowp<.05) > 0) || (sum(connectionDat.highp<.05) > 0)
    %low / high X subjectLeftOut X time
    permP = ones([2, connectionDat.subN, length(connectionDat.tim)]); 
    connectionDat.lowPermT = zeros(connectionDat.subN, length(connectionDat.tim)); 
    connectionDat.highPermT = zeros(connectionDat.subN, length(connectionDat.tim)); 
    connectionDat.lowPermP = squeeze(permP(1,:,:)); 
    connectionDat.highPermP = squeeze(permP(2,:,:)); 
    %try leaving out each subject one at a time to confirm if results
    %are the same
    for subOut = 1:length(connectionDat.uniqueSubs)
        disp(['doing permutation without: ' num2str(subOut) ' of ' num2str(length(connectionDat.uniqueSubs))])
        subSelect = cellfun(@(x) ~strcmp(x, connectionDat.uniqueSubs{subOut}), connectionDat.allSubs); 
        cur = lowDat(subSelect,:);
        cur2 = highDat(subSelect,:); 
        curd = connectionDat.d(subSelect); 
        cursub = connectionDat.allSubs(subSelect); 
        tmpT = zeros(size(lowDat,2),1); %low frequency
        tmpT2 = tmpT; %high frequency
        for ti = 1:size(lowDat,2)
            modDat = table(cur(:,ti), curd', cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
            lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
            tmpT(ti) = lme.Coefficients(2,4); %t-value associated with memory!

            modDat = table(cur2(:,ti), curd', cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
            lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
            tmpT2(ti) = lme.Coefficients(2,4); %t-value associated with memory!
        end

        connectionDat.lowPermT(subOut,:) = tmpT; 
        connectionDat.highPermT(subOut,:) = tmpT2; 
        
        lownullTs = zeros([length(hightVals), perms]); 
        highnullTs = lownullTs;  
        curUSubs = connectionDat.uniqueSubs;
        curUSubs(subOut) = [];
        cursubD = connectionDat.subD; 
        cursubD(subOut) = []; 
        subidx = cellfun(@(y) cellfun(@(x) strcmp(x,y), cursub), ...
                 curUSubs, 'UniformOutput', false); 
        
        for ii = 1:perms
            shuffd = zeros(length(cursub),1);
            shuffVals = randsample(cursubD, length(cursubD), false); 
            for sub = 1:length(subidx)
                shuffd(subidx{sub}) = shuffVals(sub); 
            end
            slice1 = lownullTs(:,ii); 
            slice2 = highnullTs(:,ii); 
            for ti = 1:size(lowDat,2)
                modDat = table(cur(:,ti), shuffd, cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
                lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
                slice1(ti) = lme.Coefficients(2,4); %t-value associated with memory!
    
                modDat = table(cur2(:,ti), shuffd, cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
                lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
                slice2(ti) = lme.Coefficients(2,4); %t-value associated with memory!
            end
            lownullTs(:,ii) = slice1; 
            highnullTs(:,ii) = slice2; 

        end
      
        %low frequency
        [~, connectionDat.lowPermP(subOut,:), ~] = cluster_test(tmpT, lownullTs); 
        %high frequency
        [~, connectionDat.highPermP(subOut,:), ~] = cluster_test(tmpT2, highnullTs); 
        save([cndFiles(idx).folder '/' cndFiles(idx).name], 'connectionDat')
    end
else
    connectionDat.lowPermP = ones([connectionDat.subN, length(connectionDat.tim)]); 
    connectionDat.highPermP = ones([connectionDat.subN, length(connectionDat.tim)]); 
end

save([cndFiles(idx).folder '/' cndFiles(idx).name], 'connectionDat')

















end