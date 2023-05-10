% channel HFB assessment


%github access token: 
%  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
%  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


%local paths: 

% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'C:\Users\dtf8829\Documents\QuestConnect\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])


datFolder = [datPre 'CHANDAT']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 


%% loop over channels and grab all the summary info for responsive channels
Ei = 1; 
Oi = 1; 
Ri = 1; 
errorChans = []; 
for chan = 1:length(chanFiles)
    chan
    chanDat = load([chanFiles(chan).folder '/' chanFiles(chan).name]).chanDat; 
    try
    %check for HFB encoding reactivity
    if chanDat.HFBenc == 1
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


%% pull out the ROI specific data

dlPFCEnc = allChanEncDat([allChanEncDat.dlPFC] ==1); 
hipEnc = allChanEncDat([allChanEncDat.hip] == 1);
phgEnc = allChanEncDat([allChanEncDat.phg] == 1);
accEnc = allChanEncDat([allChanEncDat.acc] == 1);

ROIDat = {dlPFCEnc, hipEnc, phgEnc, accEnc}; 
RoiNames = {'dlpfc', 'hip', 'phg', 'acc'}; 

for rr = 1:4
    curDat = ROIDat{rr}; 


    for chan = 1:length(curDat)
    
        figure
        subplot 311
        imagesc(movmean(curDat(chan).HFB.subMiss,50)')
%         caxis([-2,2])
        xticks([500:500:4500])
        xticklabels([-500:500:4000])
        xline(1000, '--', 'linewidth', 2, 'color', 'red')
        ylabel("miss trials")
        title(curDat(chan).subID)
    
        subplot 312
        imagesc(movmean(curDat(chan).HFB.subHit,50)')
%         caxis([-2,2])
        xticks([500:500:4500])
        xticklabels([-500:500:4000])
        xline(1000, '--', 'linewidth', 2, 'color', 'red')
        title([RoiNames{rr} '; chanNum: ' num2str(curDat(chan).chi)])
        ylabel("hit trials")
    
        subplot 313
        plot(mean(curDat(chan).HFB.subHit,2), 'color', 'blue', 'linewidth', 3)
        hold on 
        plot(mean(curDat(chan).HFB.subMiss,2), 'color', 'red', 'linewidth', 3)
        xlim([0, length(mean(curDat(chan).HFB.subMiss,2))])
        xticks([500:500:4500])
        xticklabels([-500:500:4000])
        xline(1000, '--', 'linewidth', 2, 'color', 'red')
        legend({"hits", "misses"})
    
    
    
    end


end

%% what about just looking at all data







