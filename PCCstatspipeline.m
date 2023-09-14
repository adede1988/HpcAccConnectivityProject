function [] = PCCstatspipeline(cndFiles, idx, outprefix)


connectionDat = load([cndFiles(idx).folder '/' cndFiles(idx).name]).connectionDat;
disp(cndFiles(idx).name)  


%% ENCODING! 

%using raw PPC data HITS ONLY
lowDat = connectionDat.lowBand(connectionDat.hmSort, :,2); %raw ppc
highDat = connectionDat.highBand(connectionDat.hmSort, :,2); %raw ppc



%t vals for memory
lowtVals = zeros(size(lowDat,2),1); 
hightVals = lowtVals; 

%do the initial model fit: 
%connectivity ~ memory + (1|subject) 
% where n is electrode pair nested in subject
for ti = 1:size(lowDat, 2)
   
    modDat = table(lowDat(:,ti), connectionDat.d, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
    lowtVals(ti) = lme.Coefficients(2,4); 
  
    modDat = table(highDat(:,ti), connectionDat.d, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
    lowtVals(ti) = lme.Coefficients(2,4); 

end

connectionDat.lowtVals = lowtVals; 
connectionDat.hightVals = hightVals; 

disp('raw t values calculated')


%do permutation test, refitting the model 1000 times shuffling the memory
%performance values randomly each time keeping the data such that each
%subject has one value for memory performance over all their electrodes
perms = 10; 
lownullTs = zeros([size(hightVals), perms]); 
highnullTs = lownullTs; 

subidx = cellfun(@(y) cellfun(@(x) strcmp(x,y), connectionDat.allSubs), ...
                 connectionDat.uniqueSubs, 'UniformOutput', false); 

for ii = 1:perms
%     if mod(ii, 100)==0
%         disp(['...........................' num2str(ii) 'permutations complete'])
%     end
    shuffd = zeros(size(connectionDat.d));
   
    shuffVals = randsample(connectionDat.subD, length(connectionDat.subD), false); 
    for sub = 1:length(subidx)
        shuffd(subidx{sub}) = shuffVals(sub);
    end
    sliceT = lownullTs(:,ii); 
    sliceT2 = highnullTs(:,ii); 
    for ti = 1:size(lowDat,2)
         modDat = table(lowDat(:,ti), shuffd, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
        sliceT(ti) = lme.Coefficients(2,4); 

         modDat = table(highDat(:,ti), shuffd, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
        sliceT2(ti) = lme.Coefficients(2,4);  
    end
    lownullTs(:,ii) = sliceT; 
    highnullTs(:,ii) = sliceT2; 
end

connectionDat.low975 = prctile(lownullTs, 97.5, 2); 
connectionDat.low025 = prctile(lownullTs, 2.5, 2); 
connectionDat.high975 = prctile(highnullTs, 97.5, 2); 
connectionDat.high025 = prctile(highnullTs, 2.5, 2);

[~,connectionDat.lowp_sub, ~] = cluster_test(lowtVals, lownullTs); 

[~,connectionDat.highp_sub, ~] = cluster_test(hightVals, highnullTs); 


save([cndFiles(idx).folder '/'  outprefix '_' cndFiles(idx).name], 'connectionDat')
disp('permutation on encoding complete')


%% retrieval 

%using raw PPC data HITS ONLY
lowDat = connectionDat.lowBand2(connectionDat.hmSort, :,2); %raw ppc
highDat = connectionDat.highBand2(connectionDat.hmSort, :,2); %raw ppc



%t vals for memory
lowtVals = zeros(size(lowDat,2),1); 
hightVals = lowtVals; 

%do the initial model fit: 
%connectivity ~ memory + (1|subject) 
% where n is electrode pair nested in subject
for ti = 1:size(lowDat, 2)
   
    modDat = table(lowDat(:,ti), connectionDat.d, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
    lowtVals(ti) = lme.Coefficients(2,4); 
  
    modDat = table(highDat(:,ti), connectionDat.d, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
    lowtVals(ti) = lme.Coefficients(2,4); 

end

connectionDat.lowtVals_ret = lowtVals; 
connectionDat.hightVals_ret = hightVals; 

disp('raw t values calculated')


%do permutation test, refitting the model 1000 times shuffling the memory
%performance values randomly each time keeping the data such that each
%subject has one value for memory performance over all their electrodes
perms = 10; 
lownullTs = zeros([size(hightVals), perms]); 
highnullTs = lownullTs; 

subidx = cellfun(@(y) cellfun(@(x) strcmp(x,y), connectionDat.allSubs), ...
                 connectionDat.uniqueSubs, 'UniformOutput', false); 

for ii = 1:perms
%     if mod(ii, 100)==0
%         disp(['...........................' num2str(ii) 'permutations complete'])
%     end
    shuffd = zeros(size(connectionDat.d));
   
    shuffVals = randsample(connectionDat.subD, length(connectionDat.subD), false); 
    for sub = 1:length(subidx)
        shuffd(subidx{sub}) = shuffVals(sub);
    end
    sliceT = lownullTs(:,ii); 
    sliceT2 = highnullTs(:,ii); 
    for ti = 1:size(lowDat,2)
         modDat = table(lowDat(:,ti), shuffd, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
        sliceT(ti) = lme.Coefficients(2,4); 

         modDat = table(highDat(:,ti), shuffd, connectionDat.allSubs, ...
        'VariableNames', {'connectivity',  'memory', 'sub'}); 
        lme = fitlme(modDat, 'connectivity ~ memory  +  (1|sub)'); 
        sliceT2(ti) = lme.Coefficients(2,4);  
    end
    lownullTs(:,ii) = sliceT; 
    highnullTs(:,ii) = sliceT2; 
end

connectionDat.low975_ret = prctile(lownullTs, 97.5, 2); 
connectionDat.low025_ret = prctile(lownullTs, 2.5, 2); 
connectionDat.high975_ret = prctile(highnullTs, 97.5, 2); 
connectionDat.high025_ret = prctile(highnullTs, 2.5, 2);

[~,connectionDat.lowp_ret, ~] = cluster_test(lowtVals, lownullTs); 

[~,connectionDat.highp_ret, ~] = cluster_test(hightVals, highnullTs); 


save([cndFiles(idx).folder '/'  outprefix '_' cndFiles(idx).name], 'connectionDat')
disp('permutation on retrieval complete')


















end