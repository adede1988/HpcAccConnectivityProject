% channel HFB assessment


%github access token: 
%  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
%  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


%local paths: 

% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'R:\MSS\Johnson_Lab\dtf8829\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])


datFolder = [datPre 'CHANDAT']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 


%% loop over channels and grab all the summary info for responsive channels
Ei = 1; 
Oi = 1; 
Ri = 1; 
errorChans = []; 
curRawDir = 'someStuff'; 
for chan = 1:length(chanFiles)
    chan
    chanDat = load([chanFiles(chan).folder '/' chanFiles(chan).name]).chanDat; 
    try
        chanDat = rmfield(chanDat, {'leadLag', 'leadLagRoi', 'chansRoi'});
    catch
        disp(['chan: ' num2str(chan) ' did not have leadlag'])
    end

    %make sure that the encoding info is available in the data
    if ~isfield(chanDat, 'encInfo')
        if strcmp(curRawDir, chanDat.dataDir)
            chanDat.encInfo = rawDat.trialinfo;
        else
            rawDat = load([chanDat.dataDir '/' chanDat.encDatFn]).data; 
            curRawDir = chanDat.dataDir; 
            chanDat.encInfo = rawDat.trialinfo; 
            save([chanFiles(chan).folder '/' chanFiles(chan).name], 'chanDat')
        end
    end

    try
    %check for HFB encoding reactivity
    if chanDat.HFBenc == 1  || chanDat.HFBretOn == 1 || chanDat.HFBretRT == 1
        if Ei == 1
            allChanEncDat = chanDat; 
            Ei = Ei+1; 
        else
            allChanEncDat(Ei) = chanDat; 
            Ei = Ei + 1; 
        end
    end
    catch
        errorChans = [errorChans chan]; 
    end

    


end



%% go over the reactive channels and note their ROI membership

%*********************************************************************************************************%
%regions of interest need to be fixed up a bit. If it's based on mni, then
%there should be no hippocampal electrodes because it's likely a surface
%only case, so for participants with no notes, any hippocampal electrodes
%should be reclassified as parahippocampal gyrus. 
%IMPLEMENTED!
%*********************************************************************************************************%

for chan = 1:length(allChanEncDat)
    try
    if sum(sum(allChanEncDat(chan).roiNote)) == 0 
        roi = allChanEncDat(chan).roimni(allChanEncDat(chan).chi, :);
        if roi(2) == 1
            roi(2) = 0; %ECoG channels cannot be in the hippocampus
            roi(3) = 1; %any ECoG channels flagged as H were probably in the PHG
        end
    else
        roi = allChanEncDat(chan).roiNote(allChanEncDat(chan).chi, :);
    end

    allChanEncDat(chan).dlPFC = roi(1); 
    allChanEncDat(chan).hip = roi(2); 
    allChanEncDat(chan).phg = roi(3); 
    allChanEncDat(chan).acc = roi(4); 

    catch

    allChanEncDat(chan).dlPFC = 0; 
    allChanEncDat(chan).hip = 0; 
    allChanEncDat(chan).phg = 0; 
    allChanEncDat(chan).acc = 0;     

    end






end


%% trim off the children data

allChanEncDat([allChanEncDat.age]<16) = []; 
save('R:\MSS\Johnson_Lab\dtf8829\allChanEncDat.mat', 'allChanEncDat', '-v7.3')

%% pull out the ROI specific data


% 
% dlPFCEnc = allChanEncDat([allChanEncDat.dlPFC] ==1); 
% hipEnc = allChanEncDat([allChanEncDat.hip] == 1);
% phgEnc = allChanEncDat([allChanEncDat.phg] == 1);
% accEnc = allChanEncDat([allChanEncDat.acc] == 1);
% 
% ROIDat = {dlPFCEnc, hipEnc, phgEnc, accEnc}; 
RoiNames = {'dlPFC', 'hip', 'phg', 'acc'}; 

%% heatmaps of hits vs. misses and activity w/i regions across all trials
%loop over ROIs, combine all trials within a region, plot sorted by RT and
%split by hit v. miss

for rr = 1:4
    curDat = allChanEncDat([allChanEncDat.(RoiNames{rr})] == 1); 
    allHit = []; 
    allHitRT = []; 
    allMiss = []; 
    allMissRT = []; 
    for ii = 1:length(curDat)
        comboHFB = [curDat(ii).HFB.subHit curDat(ii).HFB.subMiss]; 
        test = mean(comboHFB,2); 
      

        test = checkForThreshold(test, curDat(ii).HFB.encMulTim, [-450 2500]); %this is getting rid of negatives! 
        if test

        allHit = [allHit curDat(ii).HFB.subHit];
        allHitRT = [allHitRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).hits, 4)]; 
        allMiss = [allMiss curDat(ii).HFB.subMiss]; 
        allMissRT = [allMissRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).misses, 4)];

        end
    end

    b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)'; linspace(255,0,128)']]/255;
    b2w2r(129:end, 1) = 1; 
    b2w2r(1:128, 3) = 1; 

    [sortedRT, order] = sort(allHitRT); 
    figure
    subplot 121
    imagesc(allHit(:,order)')
    caxis([-7,7])
    hold on 
    plot(sortedRT/5 + find(curDat(1).HFB.encMulTim>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
    xticks([101:100:size(curDat(1).HFB.subMiss,1)])
    xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
    xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
    title([RoiNames{rr} ' all hit trials'])
%     colormap(jet)


    subplot 122
    [sortedRT, order] = sort(allMissRT); 
    imagesc(allMiss(:,order)')
    caxis([-7,7])
    hold on 
    plot(sortedRT/5 + find(curDat(1).HFB.encMulTim>=0,1), [1:length(sortedRT)], 'color', 'red', 'linewidth', 2)
    xticks([101:100:size(curDat(1).HFB.subMiss,1)])
    xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
    xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
    title([RoiNames{rr} ' all miss trials'])
%     colormap(jet)
end


%% average line plotting locked to stim onset and RT
for rr = 1:4
    curDat = allChanEncDat([allChanEncDat.(RoiNames{rr})] == 1); 
    allHit = []; 
    allHitRT = []; 
    allMiss = []; 
    allMissRT = []; 
    for ii = 1:length(curDat)
        comboHFB = [curDat(ii).HFB.subHit curDat(ii).HFB.subMiss]; 
        test = mean(comboHFB,2); 
      

        test = checkForThreshold(test, curDat(ii).HFB.encMulTim, [-450 2500]);
        if test

        allHit = [allHit curDat(ii).HFB.subHit];
        allHitRT = [allHitRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).hits, 4)]; 
        allMiss = [allMiss curDat(ii).HFB.subMiss]; 
        allMissRT = [allMissRT; curDat(ii).encInfo(curDat(ii).use & curDat(ii).misses, 4)];

        end
    end

    % aligned to stim onset: *********************************************
    meanHit = mean(allHit,2);
    hitStd = std(allHit,[],2) ./ sqrt(size(allHit,2)); 
    hitLow = meanHit - hitStd; %prctile(allHit, 2.5, 2); 
    hitHigh = meanHit + hitStd; %prctile(allHit, 97.5, 2); 

    meanMiss = mean(allMiss,2); 
    missStd = std(allMiss,[],2) ./ sqrt(size(allMiss,2)); 
    missLow = meanMiss - missStd; %prctile(allHit, 2.5, 2); 
    missHigh = meanMiss + missStd; %prctile(allHit, 97.5, 2); 


    colors = {[75, 122, 71]./255, [236, 146, 72]./255};


    figure
    subplot 211
    plot(meanHit, 'color', colors{1}, 'linewidth', 2)
    hold on 
    xticks([101:100:size(curDat(1).HFB.subMiss,1)])
    xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
    xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
    x = [1:length(meanHit)]; 
    x = [x flip(x)]; 
    y = [hitLow' flip(hitHigh')];
    h = fill(x,y,colors{1},'LineStyle','none'); 
    set(h, 'facealpha', .5)
    xlim([5,length(meanHit)-5])

    plot(meanMiss, 'color', colors{2}, 'linewidth', 2)
    hold on 
   
    x = [1:length(meanMiss)]; 
    x = [x flip(x)]; 
    y = [missLow' flip(missHigh')];
    h = fill(x,y,colors{2},'LineStyle','none'); 
    set(h, 'facealpha', .5)

    title(["mean hit and miss for " RoiNames{rr}])
    xlabel('time relative to image onset (ms)')
    ylabel('z-scored HFB')

    % aligned to RT: ******************************************************
    
    

    allHit_align = nan(4000, size(allHit,2)); 
    allHitRT_down = 2500 - (round(allHitRT/5) + 200); 
    
    allMiss_align = nan(4000, size(allMiss,2)); 
    allMissRT_down = 2500 - (round(allMissRT/5) + 200); 

    for ii = 1:size(allHit,2)
        if(allHitRT_down(ii)>0)
        allHit_align(allHitRT_down(ii):allHitRT_down(ii)+length(meanHit)-1, ii) = allHit(:,ii); 
        end

    end

    for ii = 1:size(allMiss,2)

        if(allMissRT_down(ii)>0)

            allMiss_align(allMissRT_down(ii):allMissRT_down(ii)+length(meanMiss)-1,ii) = allMiss(:,ii); 
        end

    end

    %trim down Hits
    test = sum(isnan(allHit_align),2);
    RT_tim = [-12495:5:7500];
    test = find(test > size(allHit_align,2)/2 );
    RT_tim(test) = []; 
    allHit_align(test,:) = []; 
    
    %use the same time to trim for misses
    allMiss_align(test,:) = []; 

    [~, order] = sort(allMissRT); 

    meanHit = mean(allHit_align,2, 'omitnan');
    hitStd = std(allHit_align,[],2, 'omitnan') ./ sqrt(size(allHit_align,2)); 
    hitLow = meanHit - hitStd; %prctile(allHit, 2.5, 2); 
    hitHigh = meanHit + hitStd; %prctile(allHit, 97.5, 2); 


    meanMiss = mean(allMiss_align,2, 'omitnan'); 
    missStd = std(allMiss_align,[],2, 'omitnan') ./ sqrt(size(allMiss_align,2)); 
    missLow = meanMiss - missStd; %prctile(allHit, 2.5, 2); 
    missHigh = meanMiss + missStd; %prctile(allHit, 97.5, 2); 

 
    subplot 212
    plot(meanHit, 'color', colors{1}, 'linewidth', 2)
    hold on 
    xticks([100:100:length(RT_tim)])
    xticklabels(RT_tim([100:100:length(RT_tim)]))
    xline(find(RT_tim>=0,1), '--', 'linewidth', 4, 'color', 'red')
    x = [1:length(meanHit)]; 
    x = [x flip(x)]; 
    y = [hitLow' flip(hitHigh')];
    h = fill(x,y,colors{1},'LineStyle','none'); 
    set(h, 'facealpha', .5)
    xlim([5,length(RT_tim)-5])

    plot(meanMiss, 'color', colors{2}, 'linewidth', 2)
    hold on 
   
    x = [1:length(meanMiss)]; 
    x = [x flip(x)]; 
    y = [missLow' flip(missHigh')];
    h = fill(x,y,colors{2},'LineStyle','none'); 
    set(h, 'facealpha', .5)

    xlabel('time relative to response time (ms)')
    ylabel('z-scored HFB')



end

%% pull out the trial by trial peak time of HFB activity

for rr = 1:4






end



%% doing leadlag for subset of subjects with contacts in two ROIs
%goal here is to save out the subject level files of ROI active channels. 
%Then use the HPC with leadLagWrapper.m and leadLagPipeline.m to do the
%analysis


%down select to ROI channels
allChanEncROI = allChanEncDat([allChanEncDat.dlPFC]==1 | [allChanEncDat.hip]==1 | [allChanEncDat.phg]==1 | [allChanEncDat.acc]==1);
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
%     %initialize the output for all channels to all channels 
%     %chan X chan X trial X time X offset
%     leadLagHit = zeros(length(curSub), length(curSub), size(curSub(1).HFB.subHit,2), size(curSub(1).HFB.subHit,1), 41);
%     leadLagMiss = zeros(length(curSub), length(curSub), size(curSub(1).HFB.subHit,2), size(curSub(1).HFB.subHit,1), 41);
%     
% 
%     for chan = 1:length(curSub)
%         for chan2 = 1:length(curSub)
%             for offSet = -20:20 %a time step is 5ms, so -20:20 is -100ms:100ms
%                 %Hit trials 
%                 HFB1 = curSub(chan).HFB.subHit; 
%                 HFB2 = curSub(chan2).HFB.subHit; 
%                 if offSet<0
%                     HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
%                 elseif offSet>0
%                     HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
%                 end
%                 
%                 %hard code 30 time step (150ms +- around time point)
%                 test = reshape(cell2mat(arrayfun(@(x) ... %x is loop on trials
%                     arrayfun(@(y) ... %y is loop on time
%                         corr(HFB1(y-30:y+30,x), HFB2(y-30:y+30,x)), ...
%                     31:size(HFB1,1)-30), ...
%                 1:size(HFB1,2), 'uniformoutput', false)),[], size(HFB2,2));
%                 
%                 leadLagHit(chan, chan2, :, 31:size(test,1)+30, offSet+21) = test'; 
% 
% 
%             end
%         end
%     end
% 
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
            allChanEncDat = chanDat; 
            Ei = Ei+1; 
        else
            allChanEncDat(Ei) = chanDat; 
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

for chan = 1:length(allChanEncDat)
    try
    if sum(sum(allChanEncDat(chan).roiNote)) == 0 
        roi = allChanEncDat(chan).roimni(allChanEncDat(chan).chi, :);
    else
        roi = allChanEncDat(chan).roiNote(allChanEncDat(chan).chi, :);
    end

    allChanEncDat(chan).dlPFC = roi(1); 
    allChanEncDat(chan).hip = roi(2); 
    allChanEncDat(chan).phg = roi(3); 
    allChanEncDat(chan).acc = roi(4); 

    catch

    end



end


%% pull out structs for each region separately and plot its mean leadLag
RoiNames = {'dlPFC', 'hip', 'phg', 'acc'}; 
for rr = 1:4

    curDat = allChanEncDat([allChanEncDat.(RoiNames{rr})] == 1); 

    %count the channels
    allLeadLagRoi = []; 
    for chan = 1:length(curDat)
        allLeadLagRoi = [allLeadLagRoi; curDat(chan).leadLagRoi]; 
        
    end

    %preallocate
    allLeadLagHit = zeros([size(allLeadLagRoi,1), size(curDat(1).leadLag.encHit, [2,3]) ] ) ; 
    allLeadLagMiss = zeros([size(allLeadLagRoi,1), size(curDat(1).leadLag.encHit, [2,3]) ] ) ; 
    lli = 1; 
    for chan = 1:length(curDat)
        allLeadLagHit(lli:size(curDat(chan).leadLag.encHit,1)+lli-1, :, :) = curDat(chan).leadLag.encHit; 
        allLeadLagMiss(lli:size(curDat(chan).leadLag.encMiss,1)+lli-1, :, :) = curDat(chan).leadLag.encMiss; 
        lli = lli + size(curDat(chan).leadLag.encMiss,1);
    end

    for rr2 = 1:4 %loop on the region whose connection with respect to the current is being tested
       tmpHit = allLeadLagHit(allLeadLagRoi(:,rr2)==1, :, :);  
       tmpMiss = allLeadLagMiss(allLeadLagRoi(:,rr2)==1, :, :);
    
       figure 
       subplot 211
       imagesc(squeeze(mean(tmpHit,1)))
       xticks([101:100:size(curDat(1).HFB.subMiss,1)])
       xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
       xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
       yticks([0:100:400])
       ylim([-.05, 400])
       yticklabels([-200:100:200])
       title([RoiNames{rr} ' hit correlation with ' RoiNames{rr2}])
       colorbar
       caxis([0,.15])
       ylabel('lags               leads')
       xlabel('time relative to stim onset (ms)')

       subplot 212
       imagesc(squeeze(mean(tmpMiss,1)))
       xticks([101:100:size(curDat(1).HFB.subMiss,1)])
       xticklabels(curDat(1).HFB.encMulTim([101:100:size(curDat(1).HFB.subMiss,1)]))
       xline(find(curDat(1).HFB.encMulTim>=0,1), '--', 'linewidth', 4, 'color', 'green')
       yticks([0:100:400])
       ylim([-.05, 400])
       yticklabels([-200:100:200])
       title([RoiNames{rr} ' miss correlation with ' RoiNames{rr2}])
       colorbar
       caxis([0,.15])
       ylabel('lags               leads') 
       xlabel('time relative to stim onset (ms)')

    end

end


dlPFCEnc = allChanEncDat([allChanEncDat.dlPFC] ==1); 
hipEnc = allChanEncDat([allChanEncDat.hip] == 1);
phgEnc = allChanEncDat([allChanEncDat.phg] == 1);
accEnc = allChanEncDat([allChanEncDat.acc] == 1);


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
         export_fig(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\singleChanHFB\' titReg '_' curDat(chan).subID '_' num2str(curDat(chan).chi) '.jpg'],''), '-r300')
%     
%     end
% 

% end

%% what about just looking at all data















































