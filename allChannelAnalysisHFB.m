% channel HFB assessment


%github access token: 
%  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
%  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


%local paths: 

codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath(genpath([codePre 'mni2atlas']))
% addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\fieldtrip-20230118')
% ft_defaults
datFolder = [datPre 'CHANDAT/finished']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 


%% loop over channels and split saved files by output variables

targVariable = 'TF';
noneFlag = true; 
errorChans = cell(37,1); 
allDat = cell(37,1); 
parfor sub = 1:37
    sub
    if isempty(allDat{sub})
        [errorChans{sub}, allDat{sub}] = getAllChanDat(chanFiles, sub, targVariable, noneFlag); %includes age and memory filter
    end
end


%% loop over channels and get single channel HFB responses

HFBdat = load('R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\hip.mat').HFBdat; 
regions = {HFBdat.aggTargs.lab}; 
regions(2) = []; 
% highfrex = linspace(70, 150, 81); 
% highnumfrex = length(highfrex);
%get the key information out for all channels, but leave empty if it's not
%a region of interest! 
%1: HFB trials X time hit, 
%2: HFB trials X time miss, 
%3: trial RTs Hit, 
%4: trial RTs Miss,
%5: region, 
%6: time
%7: subID
%8: chani
%9: TF trials X time LOW THETA HIT
%10: TF trials X time LOW THETA MISS
%11: TF trials X time HIGH theta HIT
%12: TF trials X time HIGH theta MISS
allResENC = cell(length(chanFiles), 12); 
allResRET = allResENC; 
parfor ii = 1:length(chanFiles)
ii
chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
lab = chanDat.labels{chanDat.chi,3}; 
if sum(cellfun(@(x) strcmp(x, lab), regions)) == 1  && sum(chanDat.reactiveRes==1)>0 
    %the label of this channel matches one of the ROIs and is reactive
    
    % ENCODING
    out = cell(12,1); 
    out{1} = chanDat.HFB.subHit;  
    out{2} = chanDat.HFB.subMiss; 
    out{3} = chanDat.encInfo(chanDat.use & chanDat.hits, 4); 
    out{4} = chanDat.encInfo(chanDat.use & chanDat.misses, 4); 
    out{5} = lab; 
    out{6} = chanDat.HFB.encMulTim; 
    out{7} = chanDat.subID; 
    out{8} = chanDat.chi; 
    out{9} = mean(chanDat.TF2.subHit(:,:,1:23), 3); 
    out{10}= mean(chanDat.TF2.subMiss(:,:,1:23), 3); 
    out{11}= mean(chanDat.TF2.subHit(:,:,24:38), 3); 
    out{12}= mean(chanDat.TF2.subMiss(:,:,24:38), 3); 


    allResENC(ii,:) = out; 

    % RETRIEVAL
    out = cell(6,1); 
    out{1} = chanDat.HFB.hit_on;  %pow(:, chanDat.use & chanDat.hits); 
    out{2} = chanDat.HFB.miss_on; %pow(:, chanDat.use & chanDat.misses);
    out{3} = chanDat.retInfo(chanDat.retInfo(:,1)==1, 3); 
    out{4} = chanDat.retInfo(chanDat.retInfo(:,1)==2, 3); 
    out{5} = lab; 
    out{6} = chanDat.HFB.onMulTim; 
    out{7} = chanDat.subID; 
    out{8} = chanDat.chi; 
    out{9} = mean(chanDat.TF2.hit_on(:,:,1:23), 3); 
    out{10}= mean(chanDat.TF2.miss_on(:,:,1:23), 3); 
    out{11}= mean(chanDat.TF2.hit_on(:,:,24:38), 3); 
    out{12}= mean(chanDat.TF2.miss_on(:,:,24:38), 3); 

    allResRET(ii,:) = out; 

end


end

%ditch the empties
test = allResENC; 
allResENC = test; 
allResENC(cellfun(@(x) isempty(x), allResENC(:,1)), :) = []; 
allResRET(cellfun(@(x) isempty(x), allResRET(:,1)), :) = [];
%get the latency information and put it into the 13th and 14th columns
tim = allResENC{1,6}; 
allResENC(:,13) = cellfun(@(x,y) gausLat(x, tim, y), allResENC(:,1), allResENC(:,3), ...
    'uniformoutput', false );
allResENC(:,14) = cellfun(@(x,y) gausLat(x, tim, y), allResENC(:,2), allResENC(:,4), ...
    'uniformoutput', false );
tim = allResRET{1,6}; 
allResRET(:,13) = cellfun(@(x,y) gausLat(x, tim, y), allResRET(:,1), allResRET(:,3), ...
    'uniformoutput', false );
allResRET(:,14) = cellfun(@(x,y) gausLat(x, tim, y), allResRET(:,2), allResRET(:,4), ...
    'uniformoutput', false );

%pulling out .csv file for stats in R
%chan X stats
%col 1: subID
%col 2: chi
%col 3: peak latency
%col 4: encode v. retrieve
%col 5: hit v. miss
%col 6: region
%col 7: RT
%col 8: peak latency / RT
aovDat = table;
aovDat.subID = repmat("askj", 75000,1); 
aovDat.chi = zeros(75000,1); 
aovDat.peakLat = zeros(75000,1); 
aovDat.encRet = repmat("askj", 75000,1); 
aovDat.hitMiss = repmat("askj", 75000,1); 
aovDat.reg = repmat("askj", 75000,1); 
aovDat.RT = zeros(75000,1); 
aovDat.adjTime = zeros(75000,1); 

aovi = 1; 

for ii = 1:length(allResENC)
    ii
    %hits
    L = length(allResENC{ii,9});
    aovDat.subID(aovi:aovi+L-1) = allResENC{ii,7};
    aovDat.chi(aovi:aovi+L-1) = allResENC{ii,8};
    aovDat.peakLat(aovi:aovi+L-1) = allResENC{ii,9};
    aovDat.encRet(aovi:aovi+L-1) = 'Enc';
    aovDat.hitMiss(aovi:aovi+L-1) = 'Hit';
    aovDat.reg(aovi:aovi+L-1) = allResENC{ii,5};
    aovDat.RT(aovi:aovi+L-1) = allResENC{ii,3};
    aovDat.adjTime(aovi:aovi+L-1) = aovDat.peakLat(aovi:aovi+L-1) ./ aovDat.RT(aovi:aovi+L-1);
    %misses
    aovi = aovi + L; 
    L = length(allResENC{ii,10});
    aovDat.subID(aovi:aovi+L-1) = allResENC{ii,7};
    aovDat.chi(aovi:aovi+L-1) = allResENC{ii,8};
    aovDat.peakLat(aovi:aovi+L-1) = allResENC{ii,10};
    aovDat.encRet(aovi:aovi+L-1) = 'Enc';
    aovDat.hitMiss(aovi:aovi+L-1) = 'Miss';
    aovDat.reg(aovi:aovi+L-1) = allResENC{ii,5};
    aovDat.RT(aovi:aovi+L-1) = allResENC{ii,4};
    aovDat.adjTime(aovi:aovi+L-1) = aovDat.peakLat(aovi:aovi+L-1) ./ aovDat.RT(aovi:aovi+L-1);

    aovi = aovi + L; 

end


for ii = 1:length(allResRET)
    ii
    %hits retrieval
    L = length(allResRET{ii,9});
    aovDat.subID(aovi:aovi+L-1) = allResRET{ii,7};
    aovDat.chi(aovi:aovi+L-1) = allResRET{ii,8};
    aovDat.peakLat(aovi:aovi+L-1) = allResRET{ii,9};
    aovDat.encRet(aovi:aovi+L-1) = 'Ret';
    aovDat.hitMiss(aovi:aovi+L-1) = 'Hit';
    aovDat.reg(aovi:aovi+L-1) = allResRET{ii,5};
    aovDat.RT(aovi:aovi+L-1) = allResRET{ii,3};
    aovDat.adjTime(aovi:aovi+L-1) = aovDat.peakLat(aovi:aovi+L-1) ./ aovDat.RT(aovi:aovi+L-1);
    %misses retrieval
    aovi = aovi + L; 
    L = length(allResRET{ii,10});
    aovDat.subID(aovi:aovi+L-1) = allResRET{ii,7};
    aovDat.chi(aovi:aovi+L-1) = allResRET{ii,8};
    aovDat.peakLat(aovi:aovi+L-1) = allResRET{ii,10};
    aovDat.encRet(aovi:aovi+L-1) = 'Ret';
    aovDat.hitMiss(aovi:aovi+L-1) = 'Miss';
    aovDat.reg(aovi:aovi+L-1) = allResRET{ii,5};
    aovDat.RT(aovi:aovi+L-1) = allResRET{ii,4};
    aovDat.adjTime(aovi:aovi+L-1) = aovDat.peakLat(aovi:aovi+L-1) ./ aovDat.RT(aovi:aovi+L-1);

     aovi = aovi + L; 

end

writetable(aovDat, 'R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject\trialLatDat.csv')


%% cross correlogram concept didn't really work very well. 

%cross correlograms? I think the cross correlograms don't really work
%these are not spikes, and areas have different time courses with HFB, so
%the plots just aren't very compelling. 
figure
subIDs = unique(allResENC(:,7));
targReg = [1,2,5,7,9]; 
for ri = 1:5
    for ri2 = 1:5
        subplot(5,5,(ri-1)*5+ri2)
        reg = targReg(ri); 
        reg2 = targReg(ri2); 
        %need to loop on subjects and build up my +- 500 crosscor plots
        hitCross = []; 
        missCross = []; 
        for sub = 1:length(subIDs)
            cur = allResRET(cellfun(@(x) strcmp(x, subIDs{sub}), allResRET(:,7)),: );
            tim = allResRET{1,6}; 
            cur1 = find(cellfun(@(x) strcmp(x, regions{reg}), cur(:,5))); 
            cur2 = find(cellfun(@(x) strcmp(x, regions{reg2}), cur(:,5))); 
            for cc1 = 1:length(cur1)
                for cc2 = 1:length(cur2)
                    if cur1(cc1) ~= cur2(cc2) 
                    %get snippets of region 2 data based on region 1
                    %latency
                    hitLats = cur{cur1(cc1),9};
                    hits = cur{cur2(cc2),1}; 
                    hits(:, hitLats==-1) = []; 
                    hitLats(hitLats==-1) = []; 

                    missLats = cur{cur1(cc1), 10}; 
                    misses = cur{cur2(cc2),2}; 
                    misses(:,missLats==-1)= []; 
                    missLats(missLats==-1) = [];

                    hitIdx = arrayfun(@(x) find(x==tim), hitLats); 
                    peakHit = arrayfun(@(x) hits(hitIdx(x)-20:...
                                                         hitIdx(x)+20,x), [1:length(hitIdx)], ...
                                                         'uniformoutput', false);
                    peakHit = cell2mat(peakHit);
                    
                    missIdx = arrayfun(@(x) find(x==tim), missLats); 
                    peakMiss = arrayfun(@(x) misses(missIdx(x)-20:...
                                                     missIdx(x)+20,x), [1:length(missIdx)], ...
                                                     'uniformoutput', false);
                    peakMiss = cell2mat(peakMiss);
                    
                    hitCross = [hitCross, peakHit]; 
                    missCross = [missCross, peakMiss]; 
                    end
                end
            end
        end
        
        plot(mean(hitCross,2))
        hold on 
        plot(mean(missCross,2))
        title([regions{reg} ' ' regions{reg2}])

    end
end

%% key visualizations 

%encoding visualizations 
% subset by region
regionalLatencies = cell(9,2); %regions X hit/miss
for reg = 1:9

%     regionalLatencies(reg,:) = visualizeHFBsingleTrialDat(allResENC, reg, regions, 'sub');
%     visualizeHFBsingleTrialDat_RTsort(allResENC, reg, regions, 'sub');
    visualizeTFsingleTrialDat(allResENC, reg, regions, 'sub', 1); 
    visualizeTFsingleTrialDat(allResENC, reg, regions, 'sub', 2);
    visualizeTFsingleTrialDat(allResRET, reg, regions, 'ret', 1); 
    visualizeTFsingleTrialDat(allResRET, reg, regions, 'ret', 2);  

end

% make a cool figure with the cumulative distributions of the peak
% latencies
figure
subplot 121
hold on 
targRegs = [1,2,5,6,7, 9];
colorVals = {'magenta','yellow' ,  'green',  'black','red', 'blue'};
for reg = 1:length(targRegs)
    lats = regionalLatencies{targRegs(reg), 1}; 
  
    latsPrc = arrayfun(@(x) prctile(lats, x), [0:.5:100]);
    plot(flip(latsPrc), flip([0:.5:100]), 'color', colorVals{reg}, 'linewidth', 2)

end
legend(regions(targRegs))
title('subsequent hits')

subplot 122
hold on 
for reg = 1:length(targRegs)
    lats = regionalLatencies{targRegs(reg), 2}; 

    latsPrc = arrayfun(@(x) prctile(lats, x), [0:.5:100]);
    plot(flip(latsPrc), flip([0:.5:100]), 'color', colorVals{reg}, 'linewidth', 2)

end
legend(regions(targRegs))
title('subsequent misses')



%visualize the retrieval data
regionalLatencies2 = cell(9,2); %regions X hit/miss
for reg = 1:9

    regionalLatencies2(reg,:) = visualizeHFBsingleTrialDat(allResRET, reg, regions, 'ret');
    visualizeHFBsingleTrialDat_RTsort(allResRET, reg, regions, 'ret');
end

% make a cool figure with the cumulative distributions of the peak
% latencies
figure
subplot 121
hold on 
targRegs = [1,2,5,6,7, 9];
colorVals = {'magenta','yellow' ,  'green',  'black','red', 'blue'};
for reg = 1:length(targRegs)
    lats = regionalLatencies2{targRegs(reg),1}; 
    latsPrc = arrayfun(@(x) prctile(lats, x), [0:.5:100]);
    plot(flip(latsPrc), flip([0:.5:100]), 'color', colorVals{reg}, 'linewidth', 2)

end
legend(regions(targRegs))
title('retrieval hits')

subplot 122
hold on 
for reg = 1:length(targRegs)
    lats = regionalLatencies2{targRegs(reg),2}; 
    latsPrc = arrayfun(@(x) prctile(lats, x), [0:.5:100]);
    plot(flip(latsPrc), flip([0:.5:100]), 'color', colorVals{reg}, 'linewidth', 2)

end
legend(regions(targRegs))
title('retrieval misses')




















%% PPC vs PLV testing 
% 
% %chans X  3 (trialNum, PPC, PLV)
% allPPCPLV = zeros(5000, 3); 
% allPPCPLV2 = zeros(5000, 3); 
% pi = 1; 
% for sub = 1:64
%     sub
%     if ~isempty(allDat{sub})
%         nChan = size(allDat{sub}.ISPC.subMiss,1); 
%         nsubMiss = sum(allDat{sub}.misses & allDat{sub}.use); 
%         nsubHit = sum(allDat{sub}.hits & allDat{sub}.use); 
%         ti = 60; %choose one time point
%         fi = 3; %choose one frequency
%         for ii = 1:nChan
%             for jj = 1:nChan
%                 if ii>jj
%                 allPPCPLV(pi, 1) = nsubMiss;
%                 allPPCPLV(pi, 2) = allDat{sub}.ISPC.subMiss(ii,jj, ti, fi, 2); %get PPC
%                 allPPCPLV(pi, 3) = allDat{sub}.ISPC.subMiss(ii,jj, ti, fi, 1); %get PLV
%                 allPPCPLV2(pi, 1) = nsubHit;
%                 allPPCPLV2(pi, 2) = allDat{sub}.ISPC.subHit(ii, jj, ti, fi, 2); %get PPC
%                 allPPCPLV2(pi, 3) = allDat{sub}.ISPC.subHit(ii, jj, ti, fi, 1); %get PLV
%                 pi = pi + 1; 
%                 end
%             end
%         end
% 
%     end
% end
% 
% 
% tCounts = unique(allPPCPLV(:,1));
% tCounts2 = unique(allPPCPLV2(:,1)); 
% sumPPCPLV = zeros(length(tCounts), 3); 
% sumPPCPLV2 = zeros(length(tCounts2), 3); 
% for ii = 1:length(tCounts)
%     tmp = allPPCPLV(allPPCPLV(:,1)==tCounts(ii), :); 
%     sumPPCPLV(ii,1) = tCounts(ii); 
%     sumPPCPLV(ii,2) = mean(tmp(:,2)); 
%     sumPPCPLV(ii,3) = mean(tmp(:,3)); 
% 
% end
% 
% for ii = 1:length(tCounts2)
%     tmp = allPPCPLV2(allPPCPLV2(:,1)==tCounts2(ii), :); 
%     sumPPCPLV2(ii,1) = tCounts2(ii); 
%     sumPPCPLV2(ii,2) = mean(tmp(:,2)); 
%     sumPPCPLV2(ii,3) = mean(tmp(:,3)); 
% 
% end
% 
% figure
% hold off
% scatter(sumPPCPLV(:,1), sumPPCPLV(:,2), 'filled')
% hold on 
% scatter(sumPPCPLV(:,1), sumPPCPLV(:,3), 'filled')
% ylabel("connectivity")
% xlabel("number of trials")
% title('PLV vs. PPC as a function of trial count misses')
% 
% figure
% hold off
% scatter(sumPPCPLV2(:,1), sumPPCPLV2(:,2), 'filled')
% hold on 
% scatter(sumPPCPLV2(:,1), sumPPCPLV2(:,3), 'filled')
% ylabel("connectivity")
% xlabel("number of trials")
% title('PLV vs. PPC as a function of trial count hits')
% 
% 
% figure 
% hold off
% scatter(allPPCPLV2(allPPCPLV2(:,1)>40, 2), allPPCPLV2(allPPCPLV2(:,1)>40, 3).^2, 'filled')
% hold on 
% 
% scatter(allPPCPLV2(allPPCPLV2(:,1)<30, 2), allPPCPLV2(allPPCPLV2(:,1)<30, 3).^2, 'filled')


%% get age, sex, memory performance, subject list
age = []; 
sex = cell(50,1); 
acc = []; 
d = []; 
subID = cell(50,1); 
si = 1; 
for sub = 1:length(allDat)
    if ~isempty(allDat{sub})
        age(si) = allDat{sub}.age;
        sex{si} = allDat{sub}.sex{1}; 
        T = sum(allDat{sub}.retInfo(:,1)==1 | allDat{sub}.retInfo(:,1)==2); 
        Hr = sum(allDat{sub}.retInfo(:,1)==1) / T; 
        T = sum(allDat{sub}.retInfo(:,1)==3 | allDat{sub}.retInfo(:,1)==4); 
        F = sum(allDat{sub}.retInfo(:,1)==4);
        if F == 0
            acc(si) = Hr - F/T; 
            F = 1; 
        else
            acc(si) = Hr - F/T; 
        end
        Fr = F / T; 
       
        d(si) = norminv(Hr) - norminv(Fr); 

        subID{si} = allDat{sub}.subID; 
        if si==25
            disp(sub)
        end
        si = si+1; 


    end

end


figure
subplot 121
scatter(age, d, 50, 'filled')
xlabel('age')
ylabel("memory(d')")
title('memory as a function of age')

subplot 122

scatter(age, acc, 50, 'filled')
xlabel('age')
ylabel("memory (HR - FR)")
title('memory as a function of age')

% subplot 223
% histogram(d)
% xlabel("memory (d')")
% ylabel("count")
% title('distribution of memory performance')
% 
% subplot 224
% scatter (d, acc, 50, 'filled')
% xlabel("memory (d')")
% ylabel("memory (acc)")
% title("comparison of memory measures")



%% define regions of interest


allBrod = getAllLabs(allDat);
% 
% targBrod = allBrod; %([allBrod.subN]>5); 
% targBrod(cellfun(@(x) strcmp('ERROR', x), {targBrod.lab})) = []; 
% 
% aggTargs = struct; 
% aggTargs(1).lab = {'BA8', 'BA9'};
% aggTargs(1).ROI = 'dlPFC';
% aggTargs(2).lab = {'BA44','BA45', 'BA46'}; 
% aggTargs(2).ROI = 'mlPFC'; 
% aggTargs(3).lab = {'BA47','BA10', 'BA11'}; 
% aggTargs(3).ROI = 'piPFC'; 
% aggTargs(4).lab = {'BA24', 'BA32', 'BA33', 'BA25'}; 
% aggTargs(4).ROI = 'ACC'; 
% aggTargs(5).lab = {'BA21', 'BA22', 'Fusiform (37)'}; 
% aggTargs(5).ROI = 'lTemp'; 
% aggTargs(6).lab = {'BA7', 'BA40', 'BA39'}; 
% aggTargs(6).ROI = 'Par'; 
% aggTargs(7).lab = {'BA19', 'VisualAssoc (18)', 'PrimVisual (17)'};
% aggTargs(7).ROI = 'Vis'; 
% aggTargs(8).lab = {'BA20'}; 
% aggTargs(8).ROI = 'iTemp'; 
% aggTargs(9).lab = {'BA6', 'PrimMotor (4)'}; 
% aggTargs(9).ROI = 'motor'; 
% aggTargs(10).lab = {'Parahip (36)', 'Hippocampus (54)'};
% aggTargs(10).ROI = 'MTL'; 
% aggTargs(11).lab = {'BA23', 'BA31'};
% aggTargs(11).ROI = 'PCC'; 
% 
% for reg = 1:11
%     aggTargs(reg).count = 0; 
%     aggTargs(reg).subIDs = cell(1,1); 
%     aggTargs(reg).countR = 0; 
%     aggTargs(reg).subIDsR = cell(1,1); 
% 
% end
% reactCount = zeros(length(chanFiles), 11); 
% counter = zeros(length(chanFiles), 11); 
% subIDs = cell(length(chanFiles),11); 
% subIDsreact = subIDs; 
% parfor ii = 1:length(chanFiles)
%     ii
%     chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
%     if isfield(chanDat, 'HFB') && isfield(chanDat, 'leadLag4') && isfield(chanDat, 'ISPC')%check for complete processing
%     T = sum(chanDat.retInfo(:,1)==1 | chanDat.retInfo(:,1)==2); 
%     Hr = sum(chanDat.retInfo(:,1)==1) / T; 
%     T = sum(chanDat.retInfo(:,1)==3 | chanDat.retInfo(:,1)==4); 
%     F = sum(chanDat.retInfo(:,1)==4);
%     if F == 0
%         acc =  Hr - F/T; 
%         F = 1; 
%     else
%         acc =  Hr - F/T; 
%     end
%     Fr = F / T; 
%    
%     d = norminv(Hr) - norminv(Fr); 
% 
%     if acc>0 && chanDat.age > 13 %memory
%     if ~strcmp(chanDat.subID, "SLCH010") && ~strcmp(chanDat.subID, "OH30") %subjects not ready for full integration yet
%         %find which ROI (if any) it's in
%         for reg = 1:11
%             for b = 1:length(aggTargs(reg).lab)
%                 
%                 if strcmp(chanDat.brodmann, "ERROR")
%                     zachLabs = readtable('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject/brodmann_by_subj.csv');
%                     zachLabs = zachLabs(cell2mat(cellfun(@(x) strcmp(x, chanDat.subID), {zachLabs.subj}, 'uniformoutput', false )), :);
%                     chanDat.brodmann = zachLabs.brodmann{chanDat.chi};
%                     parsave([chanFiles(ii).folder '/' chanFiles(ii).name], chanDat)
%                 end
%            
%                 reactive = reactiveTest_100(chanDat.HFB);
%                 if strcmp(aggTargs(reg).lab{b}, chanDat.brodmann)
%                     counter(ii, reg) = 1; 
%                     subIDs{ii, reg} = chanDat.subID; 
%                     if sum(reactive==1)>0
%                         reactCount(ii,reg) = 1; 
%                         subIDsreact{ii, reg} = chanDat.subID; 
%                     end
%                 end
%             end
% 
%         end
% 
% 
%     end
%     end
%     end
%     
% end
% 
% for reg = 1:11
%     aggTargs(reg).count = sum(counter(:,reg));
%     aggTargs(reg).countR = sum(reactCount(:,reg)); 
%     curSubs = subIDs(:,reg);
%     curSubsR = subIDsreact(:,reg); 
%     curSubs(cellfun(@(x) isempty(x), curSubs)) = [];
%     curSubsR(cellfun(@(x) isempty(x), curSubsR)) = []; 
%     aggTargs(reg).subIDs = length(unique(curSubs)); 
%     aggTargs(reg).subIDsR = length(unique(curSubsR)); 
% 
% end
% 


%% extract stats on leadLag HFB analysis for sending to cluster

%data files are too big, so I can't work with them all loaded in at the
%same time anymore


%  getSigConnection(allBrod, allDat); %encoding data
%  getSigConnection2(allBrod, allDat): %retrieval data


%% make time mask to match other conditions to the time of the leadlag HFB data

LL = load('R:\MSS\Johnson_Lab\dtf8829\SUMDAT\UCD_DA8_leadLag.mat').metaDat; 
HFB = load('R:\MSS\Johnson_Lab\dtf8829\SUMDAT\UCD_DA8_HFB.mat').metaDat; 
tim = LL.leadLag.encTim;
encTim = tim; 
timHFB = HFB.HFB.encMulTim; 

for ii = 1:length(timHFB)
    if ~ismember(timHFB(ii), tim)
        timHFB(ii) = 99999; 
    end
end
timMask = zeros(size(timHFB)); 
timMask(timHFB<99999) = 1; 

tim = LL.leadLag.retTim;
retTim = tim; 
timHFB = HFB.HFB.onMulTim; 

for ii = 1:length(timHFB)
    if ~ismember(timHFB(ii), tim)
        timHFB(ii) = 99999; 
    end
end
timMask2 = zeros(size(timHFB)); 
timMask2(timHFB<99999) = 1; 

%% 

[sigHFBSub, sigHFBRet] = getSigHFB(allBrod, allDat, timMask, timMask2, encTim, retTim); 
[sigTFSub, sigTFRet] = getSigTF(allBrod, allDat, timMask, timMask2, encTim, retTim);

% aggTargs = getSigISPC2(aggTargs, allDat, timMask); 

[conN, conID] = getSigISPC(allBrod, allDat, timMask, timMask2, encTim, retTim);

%% get all region - region PPC to demonstrate choice of bands

%hit/miss X pair X time X frequency
allConnections = zeros(2, 20000, sum(timMask), 20);
pi = 1; 
for sub = 1:64
    sub
    if ~isempty(allDat{sub})
        nChan = size(allDat{sub}.ISPC.subMiss,1); 
        for ii = 1:nChan
            for jj = 1:nChan
                allConnections(1,pi, :, :) = squeeze(allDat{sub}.ISPC.subHit(ii,jj,find(timMask),:,2)); 
                allConnections(2,pi, :, :) = squeeze(allDat{sub}.ISPC.subMiss(ii,jj,find(timMask),:,2));
                pi = pi+1; 
            end
        end
        

    end
end
 frex = logspace(log10(2), log10(25), 20); 
figure
imagesc(allDat{3}.leadLag.encTim, [], squeeze(mean(allConnections(1,:,:,:), 2))')
yticks([5:5:20])
yticklabels(round(frex([5:5:20])) )
set(gca, 'ydir', 'normal')
ylabel('frequency')
xlabel('time (ms)')
title('overall mean TF PPC plot subsequent hit')
caxis([0.005, .02])

figure
imagesc(allDat{3}.leadLag.encTim, [], squeeze(mean(allConnections(2,:,:,:), 2))')
yticks([5:5:20])
yticklabels(round(frex([5:5:20])) )
set(gca, 'ydir', 'normal')
ylabel('frequency')
xlabel('time (ms)')
title('overall mean TF PPC plot subsequent miss')
caxis([0.005, .02])
%% patch to get center of mass latency information and put it into the HFB summary data files
[sigHFBSub, sigHFBRet] = getSigHFB_latencyPatch(aggTargs, allDat, timMask, timMask2);

%% scratch d' calculation needs to be integrated into get all later
% dprime = []; 
% acc = []; 
% for sub = 1:length(allDat)
%     if ~isempty(allDat{sub}) && allDat{sub}.age > 16
%         T = sum(allDat{sub}.retInfo(:,1)==1 | allDat{sub}.retInfo(:,1)==2); 
%         Hr = sum(allDat{sub}.retInfo(:,1)==1) / T; 
%         T = sum(allDat{sub}.retInfo(:,1)==3 | allDat{sub}.retInfo(:,1)==4); 
%         F = sum(allDat{sub}.retInfo(:,1)==4);
%         if F == 0
%             acc = [acc, Hr - F/T]; 
%             F = 1; 
%         else
%             acc = [acc, Hr - F/T]; 
%         end
%         Fr = F / T; 
%        
%         allDat{sub}.d = norminv(Hr) - norminv(Fr); 
%         dprime = [dprime, allDat{sub}.d];
% 
%     end
% end




% [allCon, allConRet, allConN] = getAverageConnection(aggTargs, allDat);

 %; ./ sum(allConN([13,24], :),'all'); 
% MTL_LTL_hit = squeeze(sum(allCon(10, 5, 1,:,:), [1,2]) )  ./ allConN(10, 5); 
% MTL_LTL_miss = squeeze(sum(allCon(10, 5, 2,:,:), [1,2]) ) ./ allConN(10, 5);
% imagesc(MTL_LTL_hit)



% for ii = 1:11
%     for jj = 1:11
%         plotMat = (squeeze(allCon(ii,jj,1,:,:)) - squeeze(allCon(ii,jj,2,:,:))) ./ allConN(ii,jj); 
%         figure
%         imagesc(tim, [-150:50:150], plotMat)
%         yline(0)
%         xline(0)
%         caxis([-.05, .05])
%         colorbar
%         title([aggTargs(ii).ROI ' to ' aggTargs(jj).ROI])
%     end
% end


hitOverMiss = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'hitOverMiss_sig.avi'],''));
open(hitOverMiss); 
f = figure;
f.Position = [100 100 1000 600];
for tt = 1:100
    makeTimePointPlot(sigConSub, aggTargs,sigHFBSub, tt, tim)
    frame = getframe(gcf);
    writeVideo(hitOverMiss, frame); 
end

close(hitOverMiss)





% 
% hitOverMiss = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'hitOverMissRet.avi'],''));
% open(hitOverMiss); 
% f = figure;
% f.Position = [100 100 600 600];
% for tt = 1:100
%     makeTimePointPlot(allConRet, allConN, aggTargs, tt, tim)
%     frame = getframe(gcf);
%     writeVideo(hitOverMiss, frame); 
% end
% 
% close(hitOverMiss)





MTL_all_hit = squeeze(sum(allCon([20,13],[1:12,14:19,21:33],1,:,:), [1,2])) ./ sum(allConN([20,13],[1:12,14:19,21:33]), 'all'); 
MTL_all_miss = squeeze(sum(allCon([20,13],[1:12,14:19,21:33],2,:,:), [1,2])) ./ sum(allConN([20,13],[1:12,14:19,21:33]), 'all'); 

imagesc(all_all_hit)

imagesc(allDat{1}.leadLag.encTim, -150:150, all_all')


figure %MTL connectivity
%where do we have paired electrodes recorded with the MTL? 
plot(sum(allConN([13,24], :),1))
xticks([1:size(allConN,1)])
xticklabels({targBrod.lab})
ylabel('number of pair electrodes')
xlabel('regions of pair electrodes')
title('regions recorded simultaneously with MTL')


% hits > misses MTL connectivity
tim = allDat{1}.leadLag.encTim;
MTL_all = squeeze(sum(allCon([13,20], :, 2,:,:), [1,2]) );
p_bin = leadLagFigure((MTL_all_hit - MTL_all_miss)', tim, 'MTL', 'MTL to ALL subsequent hit > subsequent miss'); 
[vals, locs] = max(MTL_all); 
[valAll, loc2] = max(vals); 
offset = loc2; 
timePnt = locs(loc2); 
LLtim = [-150:150];
targTim = tim(timePnt); 
targOff = LLtim(offset);


%misses > hits MTL connectivity
MTL_all = squeeze(sum(allCon([13,20], :, :,:,2), [1,2]) );
p_bin = leadLagFigure(MTL_all, tim, 'MTL', 'MTL to ALL subsequent hit < subsequent miss'); 

[vals, locs] = max(MTL_all); 
[valAll, loc2] = max(vals); 
offset = loc2; 
timePnt = locs(loc2); 
targTim2 = tim(timePnt); 
targOff2 = LLtim(offset);

out = getExampleConnection(targBrod, allDat, [13,20], targTim, targOff,targTim2, targOff2); 

%hits > misses ALL connectivity
all_all = squeeze(sum(allCon(:, :, :,:,1), [1,2]) );
p_bin = leadLagFigure((all_all_hit - all_all_miss)', tim, 'ALL', 'ALL to ALL subsequent hit > subsequent miss'); 

%misses > hits ALL connectivity
all_all = squeeze(sum(allCon(:, :, :,:,2), [1,2]) );
p_bin = leadLagFigure(all_all, tim, 'ALL', 'ALL to ALL subsequent hit < subsequent miss'); 



figure
MTL_all = squeeze(sum(allCon([13,20], :, :,:,2), [1,2]) );
imagesc(allDat{1}.leadLag.encTim, -150:150,  MTL_all')
title('MTL to ALL subsequent hit < subsequent miss')
ylabel('MTL lags                    MTL leads')


figure %all to all connectivity
all_all = squeeze(sum(allCon(:, :, :,:,1), [1,2]) ); %; ./ sum(allConN,'all'); 
imagesc(allDat{1}.leadLag.encTim, -150:150, all_all')
ylabel('lagging connections                    leading connections')
title('ALL to ALL subsequent hit > subsequent miss')
figure %all to all connectivity
all_all = squeeze(sum(allCon(:, :, :,:,2), [1,2]) ); %; ./ sum(allConN,'all'); 
imagesc(allDat{1}.leadLag.encTim, -150:150, all_all')
ylabel('lagging connections                    leading connections')
title('ALL to ALL subsequent hit < subsequent miss')


%% get the leadLag data for all over 16 year olds

sumDat = dir([datPre 'SUMDAT']); 
LLidx = find(cellfun(@(x) length(x)>0, strfind({sumDat.name}, '_LL')));
HFBidx = find(cellfun(@(x) length(x)>0, strfind({sumDat.name}, '_HFB'))); 

LL = struct;
LLi = 1; 

for sub = 1:length(HFBidx)
    sub
    HFB = load([sumDat(HFBidx(sub)).folder '\' sumDat(HFBidx(sub)).name]).HFBout;
    if HFB.age >=16 %age check
        cur = load([sumDat(LLidx(sub)).folder '\' sumDat(LLidx(sub)).name]).LLout;
        if LLi == 1
            LL = cur; 
            LLi = LLi+1; 
        else
            LL(LLi) = cur; 
            LLi = LLi + 1; 
        end

        


    end


end

%% Permutation testing of lead lag data within pairwise connections to find significant relationships within condition




LLstats = rmfield(LL, 'LL');
LLstats(1).encTim = LL(1).LL.encTim;
LLstats(1).subMiss = 1; 
LLstats(1).subHit = 1; 
LLstats(1).contrastEnc = 1; 
LLstats(1).retTim = LL(1).LL.retTim; 
LLstats(1).missRet = 1; 
LLstats(1).hitRet = 1; 
LLstats(1).contrastRet = 1; 


parfor sub = 1:length(LL)

    %encoding conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tim = LL(sub).LL.encTim;
    LLstats(sub).encTim = tim; 
    %subMiss
    allObs = LL(sub).LL.subMiss; 
    LLStats(sub).subMiss = conditionCluTest(allObs, tim); 

    %subHit
    allObs = LL(sub).LL.subHit; 
    LLStats(sub).subHit = conditionCluTest(allObs, tim); 

    %contrastEnc
    allObs = LL(sub).LL.subHit - LL(sub).LL.subMiss;
    LLStats(sub).contrastEnc = conditionCluTest(allObs, tim); 
    

    %retrieval conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tim = LL(sub).LL.retTim;
    LLstats(sub).retTim = tim; 
    %subMiss
    allObs = LL(sub).LL.missRet; 
    LLStats(sub).missRet = conditionCluTest(allObs, tim); 

    %subHit
    allObs = LL(sub).LL.hitRet; 
    LLStats(sub).hitRet = conditionCluTest(allObs, tim); 

    %contrastEnc
    allObs = LL(sub).LL.hitRet - LL(sub).LL.missRet;
    LLStats(sub).contrastRet = conditionCluTest(allObs, tim); 






end



%% scratch
sumMat = zeros(size(obs));
count = 0; 
for ii = 1:size(cluStats,1)
    for jj = 1:size(cluStats,1)
        cur = squeeze(cluStats(ii,jj,2,:,:)); 
        ccidx = find(~isnan(cur(:,1)));
        for cc = 1:length(ccidx)
            count = count+1; 
            sumMat(tim>=cur(ccidx(cc),3) & tim<=cur(ccidx(cc),4), LLtim>=cur(ccidx(cc),6) & LLtim<=cur(ccidx(cc),7) ) = ...
                sumMat(tim>=cur(ccidx(cc),3) & tim<=cur(ccidx(cc),4), LLtim>=cur(ccidx(cc),6) & LLtim<=cur(ccidx(cc),7) )+1; 
        end

    end
end




%%


%get rid of the not currently working subjects
allChanDat([allChanDat.goodSub]==0) = []; 
%lose the kids
allChanDat([allChanDat.age]<16) = []; 

%% how many available electrodes in each region? 
regions = unique({allChanDat.brodmann}); 
LRcombo = cell(length(regions),1); 
for reg = 1:length(regions)
   LRcombo{reg} = split(regions{reg}, '-');
   LRcombo{reg} = LRcombo{reg}{end}; 
end
LRcombo = unique(LRcombo); 

regCount = struct; 
for reg = 1:length(LRcombo)
    regCount(reg).name = LRcombo{reg}; 
    regIDX = find(cellfun(@(x) ~isempty(x), cellfun(@(x) strfind(x, LRcombo{reg}), {allChanDat.brodmann}, 'uniformoutput', false)) );
    regCount(reg).n = length(unique({allChanDat(regIDX).subID}));
    regCount(reg).count = length(regIDX);
    regCount(reg).subs = unique({allChanDat(regIDX).subID});



end


%% I want to know about inter-regional connections, but not all regions are simulrecorded in all patients (FAR FROM IT!)
%I will create a regions X regions X time X lag matrix, each cell will get
%a value that is a mean of all the available connections
areaCount = sum([regCount.n]>=5);
regCount([regCount.n]<5) = []; 

subIDs = unique({allChanDat.subID});

%region 1 (ii) X region 2 (jj) X time X leadLag
subMiss = zeros([areaCount, areaCount, size(allChanDat(1).leadLag.subMiss, [3,4])]);
subHit = subMiss; 
%countDat: region 1 (ii) X region 2 (jj)
connectionCount = zeros(areaCount, areaCount); 
for ii = 1:areaCount
    ii
    for jj = 1:areaCount
        if ii ~= jj
            %need to ID all subs who have both ii and jj regions
            idx_ii = find(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(x, regCount(ii).name), {allChanDat.brodmann}, 'uniformoutput', false)) );
            sub_ii = regCount(ii).subs;
%             allSubii = {allChanEncDat(idx_ii).subID};
            idx_jj = find(cellfun(@(y) ~isempty(y), cellfun(@(x) strfind(x, regCount(jj).name), {allChanDat.brodmann}, 'uniformoutput', false)) );
            sub_jj = regCount(jj).subs;
%             allSubjj = {allChanEncDat(idx_jj).subID};
            cc=0;
            for ss = 1:length(sub_ii)
                
                if sum(cellfun(@(x) strcmp(x, sub_ii{ss}), sub_jj )) > 0 %is a subject implanted in both regions? 
                    %get their leadLag info
                    subIDX = find(cellfun(@(x) strcmp(sub_ii{ss}, x), {allChanDat.subID}));
                    leadLag = allChanDat(subIDX(1)).leadLag; 
                    subReg = {allChanDat(subIDX).brodmann};
                    LRcombo = cell(length(subReg),1); 
                    for reg = 1:length(subReg)
                       LRcombo{reg} = split(subReg{reg}, '-');
                       LRcombo{reg} = LRcombo{reg}{end}; 
                    end
                    
                    subReg_ii = find(cellfun(@(x) strcmp(x, regCount(ii).name), LRcombo)); 
                    subReg_jj = find(cellfun(@(x) strcmp(x, regCount(jj).name), LRcombo)); 
                    
                    for iii = 1:length(subReg_ii)
                        for jjj = 1:length(subReg_jj)
                            subMiss(ii,jj,:,:) = subMiss(ii,jj,:,:) + leadLag.subMiss(subReg_ii(iii), subReg_jj(jjj), :, :); 
                            subHit(ii,jj,:,:) = subHit(ii,jj,:,:) + leadLag.subHit(subReg_ii(iii), subReg_jj(jjj), :, :);
                            cc = cc+1; 
                        end
                    end
                    



                end
            end
            
            subMiss(ii,jj,:,:) = subMiss(ii,jj,:,:) ./ cc; 
            subHit(ii, jj,:,:) = subHit(ii,jj,:,:) ./ cc; 

            connectionCount(ii,jj) = cc; 

        end
    end
end

hitLead = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'hitLead.avi'],''));
open(hitLead); 
hitLag = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'hitLag.avi'],''));
open(hitLag);
missLead = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'missLead.avi'],''));
open(missLead);
missLag = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'missLag.avi'],''));
open(missLag); 

hitLeadLagDiff = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'hitLeadLagDiff.avi'],''));
open(hitLeadLagDiff); 

f = figure;
f.Position = [100 100 600 600];
for ii = 1:length(allChanDat(1).leadLag.encTim)
ii
    test = squeeze(mean(subHit(:,:,ii,151:170),4)); 
    test(isnan(test)) = 0; 
    imagesc(test)
    axis square
    xticks([1:areaCount])
    xticklabels({regCount.name})
    yticks([1:areaCount])
    yticklabels({regCount.name})
    caxis([-.25, .25])
    title(['subHit lead map at: ' num2str(allChanEncDat(1).leadLag.encTim(ii))])
    frame = getframe(gcf);
    writeVideo(hitLead, frame); 
    writeVideo(hitLead, frame); 
    writeVideo(hitLead, frame); 
    writeVideo(hitLead, frame); 


    test2 = squeeze(mean(subHit(:,:,ii,130:150),4)); 
    test2(isnan(test2)) = 0; 
    imagesc(test2)
    axis square
    xticks([1:areaCount])
    xticklabels({regCount.name})
    yticks([1:areaCount])
    yticklabels({regCount.name})
    caxis([-.25, .25])
    title(['subHit follow map at: ' num2str(allChanEncDat(1).leadLag.encTim(ii))])
    frame = getframe(gcf);
    writeVideo(hitLag, frame); 
    writeVideo(hitLag, frame); 
    writeVideo(hitLag, frame); 
    writeVideo(hitLag, frame); 

    
    imagesc(test - test2)
    axis square
    xticks([1:areaCount])
    xticklabels({regCount.name})
    yticks([1:areaCount])
    yticklabels({regCount.name})
    caxis([.05, .15])
    title(['subHit lead over Lag: ' num2str(allChanDat(1).leadLag.encTim(ii))])
    frame = getframe(gcf);
    writeVideo(hitLeadLagDiff, frame); 
    writeVideo(hitLeadLagDiff, frame); 
    writeVideo(hitLeadLagDiff, frame); 
    writeVideo(hitLeadLagDiff, frame); 



    test = squeeze(mean(subMiss(:,:,ii,151:170),4)); 
    test(isnan(test)) = 0; 
    imagesc(test)
    axis square
    xticks([1:areaCount])
    xticklabels({regCount.name})
    yticks([1:areaCount])
    yticklabels({regCount.name})
    caxis([-.25, .25])
    title(['subMiss lead map at: ' num2str(allChanEncDat(1).leadLag.encTim(ii))])
    frame = getframe(gcf);
    writeVideo(missLead, frame); 
    writeVideo(missLead, frame); 
    writeVideo(missLead, frame); 
    writeVideo(missLead, frame); 


    test = squeeze(mean(subMiss(:,:,ii,130:150),4)); 
    test(isnan(test)) = 0; 
    imagesc(test)
    axis square
    xticks([1:areaCount])
    xticklabels({regCount.name})
    yticks([1:areaCount])
    yticklabels({regCount.name})
    caxis([-.25, .25])
    title(['subMiss follow map at: ' num2str(allChanEncDat(1).leadLag.encTim(ii))])
    frame = getframe(gcf);
    writeVideo(missLag, frame); 
    writeVideo(missLag, frame); 
    writeVideo(missLag, frame); 
    writeVideo(missLag, frame); 
% % 
% 
% export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'leadLagAll' num2str(ii) '.jpg'],''), '-r300')



end

close(hitLead); 
close(hitLag); 
close(missLead); 
close(missLag); 
close(hitLeadLagDiff); 
% 
% 
% for ii = 1:length(allChanDat(1).leadLag.encTim)
%     figure('visible', false)
% plot(squeeze(mean(subMiss(13,:,ii,151:170),4)), 'linewidth', 3, 'color', 'red')
% xticks([1:areaCount])
% xticklabels({regCount.name})
% hold on 
% plot(squeeze(mean(subMiss(13,:,ii,130:150),4)), 'linewidth', 3, 'color', 'magenta')
% plot(squeeze(mean(subHit(13,:,ii,151:170),4)), 'linewidth', 3, 'color', 'green')
% plot(squeeze(mean(subHit(13,:,ii,130:150),4)), 'linewidth', 3, 'color', 'blue')
% legend({'lead miss', 'lag miss', 'lead hit', 'lag hit'})
% title(['hippocampus leadlag at: ' num2str(allChanDat(1).leadLag.encTim(ii))])
% ylim([-.25,.25])
% export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'leadLagHip' num2str(ii) '.jpg'],''), '-r300')
% end
% 
% 


%% can I make a 3d video of all the HFB activity? 
surface_template_l = load('R:\MSS\Johnson_Lab\smg1656\Recons\ft_templates/surface_pial_left.mat'); % from ft
surface_template_r = load('R:\MSS\Johnson_Lab\smg1656\Recons\ft_templates/surface_pial_right.mat');
% HFBvideo = VideoWriter(['G:\My Drive\Johnson\MTL_PFC_networkFigs\' 'HFBvid.avi']);
% open(HFBvideo); 

times = [21, 31, 41, 45, 49,53,57,61,71,81];


for ii = 1:length(times) 
    f = figure('visible', false, 'Position', [0 0 1200 1200]);
   
    ii
    subplot 221
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    for chan = 1:length(allChanDat)
        chani = allChanDat(chan).chi; 
        val = mean(allChanDat(chan).HFB.subMiss(times(ii),:));
        alpha = 0.5;
        h = scatter3(allChanDat(chan).elecpos(chani,1),allChanDat(chan).elecpos(chani,2),allChanDat(chan).elecpos(chani,3), ...
            val^2, 'filled');
        if val<0
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
        else
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        end

    end
    view([0,90])
    title(['subMiss ' num2str(allChanDat(1).leadLag.encTim(times(ii)))])
    hold off
     subplot 222
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    for chan = 1:length(allChanDat)
        chani = allChanDat(chan).chi; 
        val = mean(allChanDat(chan).HFB.subMiss(times(ii),:));
        alpha = 0.5;
        h = scatter3(allChanDat(chan).elecpos(chani,1),allChanDat(chan).elecpos(chani,2),allChanDat(chan).elecpos(chani,3), ...
            val^2, 'filled');
        if val<0
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
        else
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        end

    end
    view([90,0])
    title(['subMiss ' num2str(allChanDat(1).leadLag.encTim(times(ii)))])
    hold off

    subplot 223
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    for chan = 1:length(allChanDat)
        chani = allChanDat(chan).chi; 
        val = mean(allChanDat(chan).HFB.subHit(times(ii),:));
        alpha = 0.5;
        h = scatter3(allChanDat(chan).elecpos(chani,1),allChanDat(chan).elecpos(chani,2),allChanDat(chan).elecpos(chani,3), ...
            val^2, 'filled');
        if val<0
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
        else
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        end

    end
    view([0, 90])
    title(['subHit ' num2str(allChanDat(1).leadLag.encTim(times(ii)))])
    hold off

    subplot 224
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    for chan = 1:length(allChanDat)
        chani = allChanDat(chan).chi; 
        val = mean(allChanDat(chan).HFB.subHit(times(ii),:));
        alpha = 0.5;
        h = scatter3(allChanDat(chan).elecpos(chani,1),allChanDat(chan).elecpos(chani,2),allChanDat(chan).elecpos(chani,3), ...
            val^2, 'filled');
        if val<0
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
        else
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        end

    end
    view([90, 0])
    title(['subHit ' num2str(allChanDat(1).leadLag.encTim(times(ii)))])
    hold off
    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\' 'HFB_brainFigs' num2str(ii) '.jpg'], '-r300')
%     frame = getframe(gcf);
%     writeVideo(HFBvideo, frame); 
end


close(HFBvideo);


%% making difference score figures

times = [21, 31, 41, 45, 49,53,57,61,71,81, 91, 101, 111, 121, 131, 141, 151, 161, 171];


for ii = 1:length(times) 
    f = figure('visible', false, 'Position', [0 0 2000 1000]);
   
    ii
    subplot 121
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    for chan = 1:length(allChanDat)
        chani = allChanDat(chan).chi; 
        [h,~,~,stats] = ttest2(allChanDat(chan).HFB.subHit(times(ii),:), allChanDat(chan).HFB.subMiss(times(ii),:));
        val = stats.tstat; 
        if h==1
        alpha = 0.5;
        h = scatter3(allChanDat(chan).elecpos(chani,1),allChanDat(chan).elecpos(chani,2),allChanDat(chan).elecpos(chani,3), ...
            (val*5)^2, 'filled');
        if val<0
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
        else
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        end
        end
    end
    view([0,90])
    title(['sub Hit - sub Miss ' num2str(allChanDat(1).leadLag.encTim(times(ii)))])
    hold off
     subplot 122
    hold on 
    ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.1);
    for chan = 1:length(allChanDat)
        chani = allChanDat(chan).chi; 
        [h,~,~,stats] = ttest2(allChanDat(chan).HFB.subHit(times(ii),:), allChanDat(chan).HFB.subMiss(times(ii),:));
        val = stats.tstat; 
        if h==1
        alpha = 0.5;
        h = scatter3(allChanDat(chan).elecpos(chani,1),allChanDat(chan).elecpos(chani,2),allChanDat(chan).elecpos(chani,3), ...
            (val*5)^2, 'filled');
        if val<0
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
        else
            set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        end
        end
    end
    view([90,0])
    title(['sub Hit - sub Miss ' num2str(allChanDat(1).leadLag.encTim(times(ii)))])
    hold off

    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\' 'HFB_brainFigs_diff_t' num2str(ii) '.jpg'], '-r300')
%     frame = getframe(gcf);
%     writeVideo(HFBvideo, frame); 
end



%% go over the reactive channels and note their ROI membership

%*********************************************************************************************************%
%regions of interest need to be fixed up a bit. If it's based on mni, then
%there should be no hippocampal electrodes because it's likely a surface
%only case, so for participants with no notes, any hippocampal electrodes
%should be reclassified as parahippocampal gyrus. 
%IMPLEMENTED!
%*********************************************************************************************************%

for chan = 1:length(allChanDat)
    try
    if sum(sum(allChanDat(chan).roiNote)) == 0 
        roi = allChanDat(chan).roimni(allChanDat(chan).chi, :);
        if roi(2) == 1
            roi(2) = 0; %ECoG channels cannot be in the hippocampus
            roi(3) = 1; %any ECoG channels flagged as H were probably in the PHG
        end
    else
        roi = allChanDat(chan).roiNote(allChanDat(chan).chi, :);
    end

    allChanDat(chan).dlPFC = roi(1); 
    allChanDat(chan).hip = roi(2); 
    allChanDat(chan).phg = roi(3); 
    allChanDat(chan).acc = roi(4); 

    catch

    allChanDat(chan).dlPFC = 0; 
    allChanDat(chan).hip = 0; 
    allChanDat(chan).phg = 0; 
    allChanDat(chan).acc = 0;     

    end






end



%% trim off the children data

allChanDat([allChanDat.age]<16) = []; 
save('R:\MSS\Johnson_Lab\dtf8829\allChanEncDat.mat', 'allChanDat', '-v7.3')

%% start here for processed data! 

allChanDat = load('R:\MSS\Johnson_Lab\dtf8829\allChanEncDat.mat');
allChanDat = allChanDat.allChanEncDat; 

%% pull out the ROI specific data



% dlPFCEnc = allChanEncDat([allChanEncDat.dlPFC] ==1); 
% hipEnc = allChanEncDat([allChanEncDat.hip] == 1);
% phgEnc = allChanEncDat([allChanEncDat.phg] == 1);
% accEnc = allChanEncDat([allChanEncDat.acc] == 1);
% 
% ROIDat = {dlPFCEnc, hipEnc, phgEnc, accEnc}; 
RoiNames = {'dlPFC', 'hip', 'phg', 'acc'}; 

%% check the reactivity, looking for only upward going reactivity! 


HFBConditions = fieldnames(allChanDat(1).HFB); 

for chan = 1:length(allChanDat)
    
    %enc onset
    conditions = [1 2]; 
    comboHFB = []; 
    for con = 1:length(conditions)
        comboHFB = [comboHFB allChanDat(chan).HFB.(HFBConditions{conditions(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, allChanDat(chan).HFB.encMulTim, [-50,2000] ); 
    allChanDat(chan).HFBenc = test;

    %enc RT
    conditions = [5 6]; 
    comboHFB = []; 
    for con = 1:length(conditions)
        comboHFB = [comboHFB allChanDat(chan).HFB.(HFBConditions{conditions(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, allChanDat(chan).HFB.encMulTim, [-1500,500] ); 
    allChanDat(chan).HFBencRT = test;

    %ret on
    conditions = [9 10 11 12]; 
    comboHFB = []; 
    for con = 1:length(conditions)
        comboHFB = [comboHFB allChanDat(chan).HFB.(HFBConditions{conditions(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, allChanDat(chan).HFB.encMulTim, [-50,2000] ); 
    allChanDat(chan).HFBretOn = test;

    %ret rt
    conditions = [15 16 17 18]; 
    comboHFB = []; 
    for con = 1:length(conditions)
        comboHFB = [comboHFB allChanDat(chan).HFB.(HFBConditions{conditions(con)})];
    end

    test = mean(comboHFB,2); 
    test = checkForThreshold(test, allChanDat(chan).HFB.encMulTim, [-1500,500] ); 
    allChanDat(chan).HFBretRT = test;



end

%% add column on whether each channel is in ANY ROI

for chan = 1:length(allChanDat)
    allChanDat(chan).anyROI = sum([allChanDat(chan).dlPFC, allChanDat(chan).hip, allChanDat(chan).phg, allChanDat(chan).acc]);
    allChanDat(chan).anyReact = any([allChanDat(chan).HFBenc>0, allChanDat(chan).HFBretOn>0,allChanDat(chan).HFBretRT>0,allChanDat(chan).HFBencRT>0]);



end

%how many ROI channels are there: 
sum([allChanDat.anyROI]>0)
%how many ROI AND reactive channels are there: 
sum([allChanDat.anyROI]>0 & [allChanDat.anyReact]>0)


%% get demographics 
test = allChanDat([allChanDat.anyROI]>0 & [allChanDat.anyReact]>0); 
IDs = unique({test.subID}); 

sex = zeros(length(IDs),1); 
age = zeros(length(IDs),1); 
for sub = 1:length(IDs)
    temp = test(cellfun(@(x) strcmp(IDs{sub}, x), {test.subID} ) );
    if strcmp(temp(1).sex, 'M')
        sex(sub) = 1; %male
    else 
        sex(sub) = 2; %female
    end

    age(sub) = temp(1).age; 
end
%lose DA8 because still haven't fixed up the sample rate issue
% sex(3) = []; 
% age(3) = []; 


%% heatmaps of hits vs. misses and activity w/i regions across all trials
%loop over ROIs, combine all trials within a region, plot sorted by RT and
%split by hit v. miss

%chan X stats
%col 1: subID
%col 2: chi
%col 3: center of mass
%col 4: peak latency
%col 5: peak value
%col 6: encode v. retrieve
%col 7: condition
%col 8: region
%col 9: mean RT
%col 10: mean ( centerOfMass / RT)
aovDat = table;
aovDat.subID = repmat("askj", 1000,1); 
aovDat.chi = zeros(1000,1); 
aovDat.centerOfMass = zeros(1000,1); 
aovDat.peakLat = zeros(1000,1); 
aovDat.peakVal = zeros(1000,1); 
aovDat.encRet = repmat("askj", 1000,1); 
aovDat.cond = repmat("askj", 1000,1); 
aovDat.reg = repmat("askj", 1000,1); 
aovDat.RT = zeros(1000,1); 
aovDat.adjTime = zeros(1000,1); 

aovi = 1; 


for rr = 1:4
    curDat = allChanDat([allChanDat.(RoiNames{rr})] == 1); 
    
    hitRates = arrayfun(@(x) sum(curDat(x).retInfo(:,1)==1) / sum(curDat(x).retInfo(:,1)<5), 1:length(curDat));
    faRates = arrayfun(@(x) sum(curDat(x).retInfo(:,1)==4) / sum(curDat(x).retInfo(:,1)<5), 1:length(curDat));
    dPrime = norminv(hitRates) - norminv(faRates);
    
%     curDat(dPrime<1.5) = []; 


%     curDat(cellfun(@(x) strcmp('DA8', x), {curDat.subID})) = []; 
    curDatSum = HFBSummary(curDat, 1); 

    %need to add in the subsequent memory RT locked responses RESOLVED
    
    
    


    [aovDat, aovi] = plotConditionCompare(curDatSum, [1,2], RoiNames{rr}, linspace(-50,2500, 30), aovDat, aovi);
    [aovDat, aovi] = plotConditionCompare(curDatSum, [3,4], RoiNames{rr}, linspace(-1500,500, 30), aovDat, aovi);
    [aovDat, aovi] = plotConditionCompare(curDatSum, [7,5], RoiNames{rr}, linspace(-50,2000, 30), aovDat, aovi);
    [aovDat, aovi] = plotConditionCompare(curDatSum, [11,9], RoiNames{rr}, linspace(-1500,500, 30), aovDat, aovi);
    



end

writetable(aovDat, join(['R:\MSS\Johnson_Lab\dtf8829\' 'latencyAOVdat.csv'],''))
  





%% doing leadlag for subset of subjects with contacts in two ROIs
%goal here is to save out the subject level files of ROI active channels. 
%Then use the HPC with leadLagWrapper.m and leadLagPipeline.m to do the
%analysis


%down select to ROI channels
allChanEncROI = allChanDat([allChanDat.dlPFC]==1 | [allChanDat.hip]==1 | [allChanDat.phg]==1 | [allChanDat.acc]==1);
 allChanEncROI(cellfun(@(x) strcmp('DA8', x), {allChanEncROI.subID})) = []; 
%get subject IDs
IDs = unique({allChanEncROI.subID}); 

for sub = 1:length(IDs)
    sub
    subMask = cellfun(@(x) strcmp(IDs{sub}, x), {allChanEncROI.subID}); 
    curSub = allChanEncROI(subMask);  
    roiMat = [[curSub.dlPFC]; [curSub.hip]; [curSub.phg]; [curSub.acc]]';
    
    %eliminat double counted electrodes (eCog Hip/Phg)
    doubles = find(sum(roiMat,2) == 2); 
    roiMat(doubles,2) = 0; 



    curSub(1).activeRoiMat = roiMat; 
  

    save(['R:\MSS\Johnson_Lab\dtf8829\HFBCONDAT\' 'HFBCONDAT_' IDs{sub} '.mat'], 'curSub', '-v7.3')
end


%% scratch connectivity code
% 
   

     
%         for rr = 1:4
%             if rr ~= curRoi && sum(roiMat(:,rr)) ~= 0
%                 %if both of these are passed, then it's possible there were
%                 %simultaneous channels! But, they need to be in the
%                 %visually responsive set
%                 
%                 %get the numbers of the candidates
%                 partnerIDX = find(roiMat(:,rr));
%                 subID = allChanEncDat(chan).subID; 
% 
% 
%                 
% 
%             end
% 
% 
% 
% 
% 
%         end
% 
% 
% 
% 
% 
% 
% 
%   
% 
% 
% 
% 
% 
% 
% 
% end












%  if ~isstring(allChanEncDat(chan).leadLag) %check if there's any leadLag to be found
%         curChani = allChanEncDat(chan).chi; 
%         %get the Roi matrix
%         if sum(sum(allChanEncDat(chan).roiNote)) == 0
%             roiMat = allChanEncDat(chan).roimni; 
%         else
%             roiMat = allChanEncDat(chan).roiNote; 
%         end










%% Leadlag aggregating 

Ei = 1; 
errorChans = []; 
for chan = 247:length(chanFiles)
    chan
    chanDat = load([chanFiles(chan).folder '/' chanFiles(chan).name]).chanDat; 
    if isfield(chanDat, 'leadLag')
        try
    if ~isstring(chanDat.leadLag)
        if Ei == 1
            allChanDat = chanDat; 
            Ei = Ei+1; 
        else
            allChanDat(Ei) = chanDat; 
            Ei = Ei + 1; 
        end
    end
        catch
            errorChans = [errorChans chan]; 
        end
    else 
        errorChans = [errorChans chan]; 
    end
end

%% get ROI memberships

for chan = 1:length(allChanDat)
    try
    if sum(sum(allChanDat(chan).roiNote)) == 0 
        roi = allChanDat(chan).roimni(allChanDat(chan).chi, :);
    else
        roi = allChanDat(chan).roiNote(allChanDat(chan).chi, :);
    end

    allChanDat(chan).dlPFC = roi(1); 
    allChanDat(chan).hip = roi(2); 
    allChanDat(chan).phg = roi(3); 
    allChanDat(chan).acc = roi(4); 

    catch

    end



end


%% pull out structs for each region separately and plot its mean leadLag
% RoiNames = {'dlPFC', 'hip', 'phg', 'acc'}; 
% for rr = 1:4
% 
%     curDat = allChanEncDat([allChanEncDat.(RoiNames{rr})] == 1); 
% 
%     %count the channels
%     allLeadLagRoi = []; 
%     for chan = 1:length(curDat)
%         allLeadLagRoi = [allLeadLagRoi; curDat(chan).leadLagRoi]; 
%         
%     end
% 
%     %preallocate
%     allLeadLagHit = zeros([size(allLeadLagRoi,1), size(curDat(1).leadLag.encHit, [2,3]) ] ) ; 
%     allLeadLagMiss = zeros([size(allLeadLagRoi,1), size(curDat(1).leadLag.encHit, [2,3]) ] ) ; 
%     lli = 1; 
%     for chan = 1:length(curDat)
%         allLeadLagHit(lli:size(curDat(chan).leadLag.encHit,1)+lli-1, :, :) = curDat(chan).leadLag.encHit; 
%         allLeadLagMiss(lli:size(curDat(chan).leadLag.encMiss,1)+lli-1, :, :) = curDat(chan).leadLag.encMiss; 
%         lli = lli + size(curDat(chan).leadLag.encMiss,1);
%     end
% 
%     for rr2 = 1:4 %loop on the region whose connection with respect to the current is being tested
%        tmpHit = allLeadLagHit(allLeadLagRoi(:,rr2)==1, :, :);  
%        tmpMiss = allLeadLagMiss(allLeadLagRoi(:,rr2)==1, :, :);
%     
%        figure 
%        subplot 211
%        imagesc(squeeze(mean(tmpHit,1)))
%        xticks([101:100:size(curDat(1).HFB.subMiss,1)])
%        xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
%        xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
%        yticks([0:100:400])
%        ylim([-.05, 400])
%        yticklabels([-200:100:200])
%        title([RoiNames{rr} ' hit correlation with ' RoiNames{rr2}])
%        colorbar
%        caxis([0,.15])
%        ylabel('lags               leads')
%        xlabel('time relative to stim onset (ms)')
% 
%        subplot 212
%        imagesc(squeeze(mean(tmpMiss,1)))
%        xticks([101:100:size(curDat(1).HFB.subMiss,1)])
%        xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
%        xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
%        yticks([0:100:400])
%        ylim([-.05, 400])
%        yticklabels([-200:100:200])
%        title([RoiNames{rr} ' miss correlation with ' RoiNames{rr2}])
%        colorbar
%        caxis([0,.15])
%        ylabel('lags               leads') 
%        xlabel('time relative to stim onset (ms)')
% 
%     end
% 
% end
% 
% % 
% % dlPFCEnc = allChanEncDat([allChanEncDat.dlPFC] ==1); 
% % hipEnc = allChanEncDat([allChanEncDat.hip] == 1);
% % phgEnc = allChanEncDat([allChanEncDat.phg] == 1);
% % accEnc = allChanEncDat([allChanEncDat.acc] == 1);
% % 

% for rr = 1:4
%     curDat = allChanEncDat; %ROIDat{rr}; 
% 
% 
%     for chan = 1:length(curDat)
%     
%         
%         chan
% 
% 
%         figure('visible', false)
%         subplot 311
%         imagesc(curDat(chan).HFB.subMiss')
%         caxis([-10,10])
%         xticks([101:100:size(curDat(chan).HFB.subMiss,1)])
%         xticklabels(curDat(chan).HFB.encMulTim([101:100:size(curDat(chan).HFB.subMiss,1)]))
%         xline(find(curDat(chan).HFB.encMulTim>=0,1), '--', 'linewidth', 2, 'color', 'red')
%         xline(find(curDat(chan).HFB.encMulTim>=3000,1), '--', 'linewidth', 2, 'color', 'red')
%         ylabel("miss trials")
%         title(curDat(chan).subID)
%     
%         subplot 312
%         imagesc(curDat(chan).HFB.subHit')
%         caxis([-10,10])
%         xticks([101:100:size(curDat(chan).HFB.subMiss,1)])
%         xticklabels(curDat(chan).HFB.encMulTim([101:100:size(curDat(chan).HFB.subMiss,1)]))
%         xline(find(curDat(chan).HFB.encMulTim>=0,1), '--', 'linewidth', 2, 'color', 'red')
%         xline(find(curDat(chan).HFB.encMulTim>=3000,1), '--', 'linewidth', 2, 'color', 'red')
%         if curDat(chan).dlPFC == 1
%             titReg = 'dlpfc'; 
%         elseif curDat(chan).hip == 1
%             titReg = 'hip'; 
%         elseif curDat(chan).phg == 1
%             titReg = 'phg'; 
%         elseif curDat(chan).acc == 1
%             titReg = 'acc'; 
%         else
%             titReg = 'other';
%         end
%         
%        
%         title([titReg ' chanNum: ' num2str(curDat(chan).chi)])
%         ylabel("hit trials")
%     
% 
%         %t test at each time point
%         
% 
%         subplot 313
%         plot(mean(curDat(chan).HFB.subHit,2), 'color', 'blue', 'linewidth', 3)
%         hold on 
%         plot(mean(curDat(chan).HFB.subMiss,2), 'color', 'red', 'linewidth', 3)
% 
%         test = arrayfun(@(x) ttest2(curDat(chan).HFB.subMiss(x,:), curDat(chan).HFB.subHit(x,:)), 1:size(curDat(chan).HFB.subHit,1) );
%         test(test==1) = max(max(mean(curDat(chan).HFB.subHit,2)), max(mean(curDat(chan).HFB.subMiss,2)))+1; 
%         xVals = [1:length(test)];
%         xVals(test==0) = []; 
%         test(test==0) = []; 
%         breakPoints = [1 find(diff(xVals)>1) length(xVals)];
%         for bb=1:length(breakPoints)-1
%             if breakPoints(bb+1) - breakPoints(bb) > 5 %only plot significance if it's sustained for 50ms
%                 scatter(xVals(breakPoints(bb)+1:breakPoints(bb+1)), test(breakPoints(bb)+1:breakPoints(bb+1)), 'k', 'filled')
%             end
%         end
%         
% 
% 
%         xlim([0, length(mean(curDat(chan).HFB.subMiss,2))])
%         xticks([101:100:size(curDat(chan).HFB.subMiss,1)])
%         xticklabels(curDat(chan).HFB.encMulTim([101:100:size(curDat(chan).HFB.subMiss,1)]))
%         xline(find(curDat(chan).HFB.encMulTim>=0,1), '--', 'linewidth', 2, 'color', 'red')
%         xline(find(curDat(chan).HFB.encMulTim>=3000,1), '--', 'linewidth', 2, 'color', 'red')
% %         legend({"hits", "misses"}, 'location', 'best', 'autoupdate', 'off')
%         yline(-1.96, '--')
%         yline(1.96, '--')
%         ylabel('mean z-score')
%         xlabel('time (ms)')
%     
%          export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\singleChanHFB\' titReg '_' curDat(chan).subID '_' num2str(curDat(chan).chi) '.jpg'],''), '-r300')
%     
%     end
% 

% end

%% what about just looking at all data






% 
% 
% %% average line plotting locked to stim onset and RT
% for rr = 1:4
%     curDat = allChanEncDat([allChanEncDat.(RoiNames{rr})] == 1); 
%     allHit = []; 
%     allHitRT = []; 
%     allMiss = []; 
%     allMissRT = []; 
%     for ii = 1:length(curDat)
%         comboHFB = [curDat(ii).HFB.subHit curDat(ii).HFB.subMiss]; 
%         test = mean(comboHFB,2); 
%       
% 
%         test = checkForThreshold(test, curDat(ii).HFB.encMulTim, [-450 2500]);
%         if test
% 
%         allHit = [allHit curDat(ii).HFB.subHit];
%         allHitRT = [allHitRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).hits, 4)]; 
%         allMiss = [allMiss curDat(ii).HFB.subMiss]; 
%         allMissRT = [allMissRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).misses, 4)];
% 
%         end
%     end
% 
%     % aligned to stim onset: *********************************************
%     meanHit = mean(allHit,2);
%     hitStd = std(allHit,[],2) ./ sqrt(size(allHit,2)); 
%     hitLow = meanHit - hitStd; %prctile(allHit, 2.5, 2); 
%     hitHigh = meanHit + hitStd; %prctile(allHit, 97.5, 2); 
% 
%     meanMiss = mean(allMiss,2); 
%     missStd = std(allMiss,[],2) ./ sqrt(size(allMiss,2)); 
%     missLow = meanMiss - missStd; %prctile(allHit, 2.5, 2); 
%     missHigh = meanMiss + missStd; %prctile(allHit, 97.5, 2); 
% 
% 
%     colors = {[75, 122, 71]./255, [236, 146, 72]./255};
% 
% 
%     figure
%     subplot 211
%     plot(meanHit, 'color', colors{1}, 'linewidth', 2)
%     hold on 
%     xticks([101:100:size(curDat(1).HFB.subMiss,1)])
%     xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
%     xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
%     x = [1:length(meanHit)]; 
%     x = [x flip(x)]; 
%     y = [hitLow' flip(hitHigh')];
%     h = fill(x,y,colors{1},'LineStyle','none'); 
%     set(h, 'facealpha', .5)
%     xlim([5,length(meanHit)-5])
% 
%     plot(meanMiss, 'color', colors{2}, 'linewidth', 2)
%     hold on 
%    
%     x = [1:length(meanMiss)]; 
%     x = [x flip(x)]; 
%     y = [missLow' flip(missHigh')];
%     h = fill(x,y,colors{2},'LineStyle','none'); 
%     set(h, 'facealpha', .5)
% 
%     title(["mean hit and miss for " RoiNames{rr}])
%     xlabel('time relative to image onset (ms)')
%     ylabel('z-scored HFB')
% 
%     % aligned to RT: ******************************************************
%     
%     
% 
%     allHit_align = nan(4000, size(allHit,2)); 
%     allHitRT_down = 2500 - (round(allHitRT/5) + 200); 
%     
%     allMiss_align = nan(4000, size(allMiss,2)); 
%     allMissRT_down = 2500 - (round(allMissRT/5) + 200); 
% 
%     for ii = 1:size(allHit,2)
%         if(allHitRT_down(ii)>0)
%         allHit_align(allHitRT_down(ii):allHitRT_down(ii)+length(meanHit)-1, ii) = allHit(:,ii); 
%         end
% 
%     end
% 
%     for ii = 1:size(allMiss,2)
% 
%         if(allMissRT_down(ii)>0)
% 
%             allMiss_align(allMissRT_down(ii):allMissRT_down(ii)+length(meanMiss)-1,ii) = allMiss(:,ii); 
%         end
% 
%     end
% 
%     %trim down Hits
%     test = sum(isnan(allHit_align),2);
%     RT_tim = [-12495:5:7500];
%     test = find(test > size(allHit_align,2)/2 );
%     RT_tim(test) = []; 
%     allHit_align(test,:) = []; 
%     
%     %use the same time to trim for misses
%     allMiss_align(test,:) = []; 
% 
%     [~, order] = sort(allMissRT); 
% 
%     meanHit = mean(allHit_align,2, 'omitnan');
%     hitStd = std(allHit_align,[],2, 'omitnan') ./ sqrt(size(allHit_align,2)); 
%     hitLow = meanHit - hitStd; %prctile(allHit, 2.5, 2); 
%     hitHigh = meanHit + hitStd; %prctile(allHit, 97.5, 2); 
% 
% 
%     meanMiss = mean(allMiss_align,2, 'omitnan'); 
%     missStd = std(allMiss_align,[],2, 'omitnan') ./ sqrt(size(allMiss_align,2)); 
%     missLow = meanMiss - missStd; %prctile(allHit, 2.5, 2); 
%     missHigh = meanMiss + missStd; %prctile(allHit, 97.5, 2); 
% 
%  
%     subplot 212
%     plot(meanHit, 'color', colors{1}, 'linewidth', 2)
%     hold on 
%     xticks([100:100:length(RT_tim)])
%     xticklabels(RT_tim([100:100:length(RT_tim)]))
%     xline(find(RT_tim>=0,1), '--', 'linewidth', 4, 'color', 'red')
%     x = [1:length(meanHit)]; 
%     x = [x flip(x)]; 
%     y = [hitLow' flip(hitHigh')];
%     h = fill(x,y,colors{1},'LineStyle','none'); 
%     set(h, 'facealpha', .5)
%     xlim([5,length(RT_tim)-5])
% 
%     plot(meanMiss, 'color', colors{2}, 'linewidth', 2)
%     hold on 
%    
%     x = [1:length(meanMiss)]; 
%     x = [x flip(x)]; 
%     y = [missLow' flip(missHigh')];
%     h = fill(x,y,colors{2},'LineStyle','none'); 
%     set(h, 'facealpha', .5)
% 
%     xlabel('time relative to response time (ms)')
%     ylabel('z-scored HFB')
% 
% 
% 
% end
% 
% %% pull out the trial by trial peak time of HFB activity
% 
% for rr = 1:4
%     curDat = allChanEncDat([allChanEncDat.(RoiNames{rr})] == 1); 
%     allHit = []; 
%     allHitRT = []; 
%     allMiss = []; 
%     allMissRT = []; 
%     for ii = 1:length(curDat)
%         comboHFB = [curDat(ii).HFB.subHit curDat(ii).HFB.subMiss]; 
%         test = mean(comboHFB,2); 
%       
% 
%         test = checkForThreshold(test, curDat(ii).HFB.encMulTim, [-450 2500]);
%         if test
% 
%         allHit = [allHit curDat(ii).HFB.subHit];
%         allHitRT = [allHitRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).hits, 4)]; 
%         allMiss = [allMiss curDat(ii).HFB.subMiss]; 
%         allMissRT = [allMissRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).misses, 4)];
% 
%         end
%     end
% 
% 
%     %Latency data now on a trial by trial basis 
% 
% 
% 
% 
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 









%    b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)']]/255;
%     b2w2r(129:end, 1) = 1; 
%     b2w2r(1:128, 3) = 1; 
% 
%     [sortedRT, order] = sort(allHitRT); 
%     figure
%     subplot 121
%     imagesc(allHit(:,order)')
%     caxis([-7,7])
%     hold on 
%     plot(sortedRT/5 + find(curDat(1).HFB.encMulTim>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
%     xticks([101:100:size(curDat(1).HFB.subMiss,1)])
%     xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
%     xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
%     title([RoiNames{rr} ' all hit trials'])
% %     colormap(jet)
% 
% 
%     subplot 122
%     [sortedRT, order] = sort(allMissRT); 
%     imagesc(allMiss(:,order)')
%     caxis([-7,7])
%     hold on 
%     plot(sortedRT/5 + find(curDat(1).HFB.encMulTim>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
%     xticks([101:100:size(curDat(1).HFB.subMiss,1)])
%     xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
%     xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
%     title([RoiNames{rr} ' all miss trials'])
% %     colormap(jet)
















