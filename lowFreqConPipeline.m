function [] = lowFreqConPipeline(cndFiles, idx)


connectionDat = load([cndFiles(idx).folder '/' cndFiles(idx).name]).connectionDat;
disp(cndFiles(idx).name)  
%% implementation of linear mixed effects model! 

lowDat = connectionDat.lowBand(connectionDat.hmSort, :,4); %z-score hit data only
highDat = connectionDat.highBand(connectionDat.hmSort, :,4); %z-score hit data only

%log transform 
lowDat = log10((lowDat - min(lowDat,[],'all')) + .1);
highDat = log10((highDat - min(highDat,[], 'all')) + .1); 

lowtVals = zeros(size(lowDat,1),1);
hightVals = lowtVals; 


for ti = 1:size(lowDat, 2)
   
    modDat = table(lowDat(:,ti), connectionDat.d', connectionDat.allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
    lowtVals(ti) = lme.Coefficients(2,4); %t-value associated with memory! 

   
    modDat = table(highDat(:,ti), connectionDat.d', connectionDat.allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
    lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
    hightVals(ti) = lme.Coefficients(2,4); %t-value associated with memory! 

end

connectionDat.lottVals = lowtVals; 
connectionDat.hightVals = hightVals; 

disp('raw t values calculated')
%         %potentially use this to get a normally distributed rank value for
%         %d'
% %         rankit = INVNORMAL ((Rank of X - 0.5)/ n)
%         
% 

perms = 1000; 
lownullTs = zeros([length(hightVals), perms]); 
highnullTs = lownullTs; 
for ii = 1:perms
    
    shuffd = randsample(connectionDat.d, length(connectionDat.d), false); 
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


[~,connectionDat.lowp, ~] = cluster_test(lowtVals, lownullTs); 
[~,connectionDat.highp, ~] = cluster_test(hightVals, highnullTs); 

disp('permutation on full set complete')

if (sum(connectionDat.lowp<.05) > 0) || (sum(connectionDat.highp<.05) > 0)
    %low / high X subjectLeftOut X time
    permP = ones([2, connectionDat.subN, length(connectionDat.tim)]); 
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
        tmpT = zeros(size(lowDat,1),1); %low frequency
        tmpT2 = tmpT; %high frequency
        for ti = 1:size(lowDat,2)
            modDat = table(cur(:,ti), curd', cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
            lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
            tmpT(ti) = lme.Coefficients(2,4); %t-value associated with memory!

            modDat = table(cur2(:,ti), curd', cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
            lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
            tmpT2(ti) = lme.Coefficients(2,4); %t-value associated with memory!
        end

       
        lownullTs = zeros([length(hightVals), perms]); 
        highnullTs = lownullTs;  
        for ii = 1:perms
            shuffd = randsample(curd, length(curd), false); 
            slice1 = lownullTs(:,ii); 
            slice2 = highnullTs(:,ii); 
            for ti = 1:size(lowDat,2)
                modDat = table(cur(:,ti), shuffd', cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
                lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
                slice1(ti) = lme.Coefficients(2,4); %t-value associated with memory!
    
                modDat = table(cur2(:,ti), shuffd', cursub', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
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


toc














end