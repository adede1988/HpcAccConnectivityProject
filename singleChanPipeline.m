function [] = singleChanPipeline(chanFiles, idx, codePre)

%% set frequency parameters
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);

highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = linspace(5, 7, highnumfrex); 

%% load the data

try %try loading the processed file
    chanDat = load([chanFiles(idx).folder '/' chanFiles(idx).name]).chanDat; 
catch
    chanDat = load([chanFiles(idx).folder '/CHANRAW/' chanFiles(idx).name]).chanDat; % go raw if it's not working!
end

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])

%% get channel labels from Zach's CSV

zachLabs = readtable([codePre 'HpcAccConnectivityProject/brodmann_by_subj.csv']);
zachLabs = zachLabs(cell2mat(cellfun(@(x) strcmp(x, chanDat.subID), {zachLabs.subj}, 'uniformoutput', false )), :);
zachLabs(isnan(zachLabs.x), :) = []; 
%check that order is the same 
try

    if sum(abs(zachLabs.x - chanDat.elecpos(:,1))) < .00001
        zachLabs = zachLabs(chanDat.chi, :); 
        chanDat.brodmann = zachLabs.brodmann{1}; 
        chanDat.ogChan = zachLabs.Channel{1}; 
        if isfield(chanDat, 'aal_lab')
            chanDat = rmfield(chanDat, 'aal_lab');
        end
    else
        disp(['MISMATCH WITH ZACH COORDINATES!!!'])
        chanDat.brodmann = 'ERROR'; 
        chanDat.aal_lab = 'ERROR'; 
    end
catch
        for xx = 1:length(chanDat.elecpos(:,1))
        labIdx = find(arrayfun(@(x) abs(chanDat.elecpos(xx,1)-x), zachLabs.x)<.00001);
        if isempty(labIdx)
            disp(xx)
        end

        end


%     end

    disp(['MISMATCH WITH ZACH COORDINATES!!!'])
    chanDat.brodmann = 'ERROR'; 
    chanDat.aal_lab = 'ERROR';

end

%% check for encoding info
% if ~isfield(chanDat, 'encInfo')
%     dataDirPath = split(chanDat.dataDir, 'Johnson_Lab');
%     dat = load(fullfile(['/projects/p31578' dataDirPath{2} '/' chanDat.encDatFn])).data; 
%     chanDat.encInfo = dat.trialinfo; 
%     clear dat
%     save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat')
% end

%% time points! hard code

chanDat.enctim = [-1000:3500];
chanDat.enctimRT = [-2000:500];
chanDat.retOtim = [-1000:3000];
chanDat.retRtim = [-2000:500];

%% ROI membership

    if sum(sum(chanDat.roiNote)) == 0 
        roi = chanDat.roimni(chanDat.chi, :);
        if roi(2) == 1
            roi(2) = 0; %ECoG channels cannot be in the hippocampus
            roi(3) = 1; %any ECoG channels flagged as H were probably in the PHG
        end
    else
        roi = chanDat.roiNote(chanDat.chi, :);
    end

    chanDat.dlPFC = roi(1); 
    chanDat.hip = roi(2); 
    chanDat.phg = roi(3); 
    chanDat.acc = roi(4); 

    



%% High frequency Broadband 

%note hardcoded baselines: encoding: -450 : -50   ms
%                          retOn   : -450 : -50   ms
%                          retRT   : -2000: -1600 ms

if ~isfield(chanDat, 'HFB')
    HFB = getHFB(chanDat, highfrex); 

    chanDat.HFB = HFB; 
    
    clear HFB 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])

else
    disp('HFB already done')
end


%% Lead lag analysis




if ~isfield(chanDat, 'leadLag4')

        leadLag = struct; 
        reactive = reactiveTest(chanDat.HFB);
        includedChans = []; 
        start = 1; 


          
        leadLagEncTim = chanDat.enctim(501:25:end-500);
        leadLagRetTim = chanDat.retOtim(501:25:end-500);

        %correlation matrices
        %partner chan, hit/miss, offset, time
        outCluStats = nan(length(chanFiles), 2, 301, length(leadLagEncTim)); %subsequent memory
        outCluStats2 = nan(length(chanFiles), 2, 301, length(leadLagRetTim)); %retrieval




    if sum(reactive==1)>0

    
    %trial index values
    subMiss = find(chanDat.use & chanDat.misses); 
    subHit = find(chanDat.use & chanDat.hits); 
    miss_on = find(chanDat.retInfo(:,1)==2);
    hit_on = find(chanDat.retInfo(:,1)==1); 
   






    if start<length(chanFiles) %is there even any work left to be done? 

    for chan = start:length(chanFiles)
        disp(['leadLag analysis with channel: ' num2str(chan) ' of ' num2str(length(chanFiles))])
        tic

        chanDat2 = load([chanFiles(chan).folder '/CHANRAW/' chanFiles(chan).name]).chanDat; 
        chanDat2.HFB = getHFB(chanDat2, highfrex);
        reactive2 = reactiveTest(chanDat2.HFB);
        if sum(reactive2==1)>0
            includedChans = [includedChans chan]; 
        %need to grab the HFB data at higher resolution, so calculate from
        %scratch 
    
        %ENCODING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHAN 1
        [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.enc, highfrex, chanDat.fsample, chanDat.enctim, 1);  
        pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
        pow = cell2mat(pow); %organize
        highnumfrex = length(mulFrex); 
        pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
        pow = squeeze(mean(pow, 3)); %take the mean over frequencies
    
        HFB1 = pow; 
        clear pow
    
        %CHAN 2
        [pow2, mulTim, mulFrex] = getChanMultiTF(chanDat2.enc, highfrex, chanDat.fsample, chanDat.enctim, 1);  
        pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow2,3), 'UniformOutput',false ); %z-score
        pow2 = cell2mat(pow2); %organize
        highnumfrex = length(mulFrex); 
        pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/highnumfrex, []); %organize
        pow2 = squeeze(mean(pow2, 3)); %take the mean over frequencies
        
        disp('encoding: ')
        outCluStats(chan,:,:,:) = getLL(HFB1, pow2, subMiss, subHit, chanDat.enctim, leadLagEncTim);


      


    
        %RERTRIEVAL DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHAN 1
        [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retOn, highfrex, chanDat.fsample, chanDat.retOtim, 1);  
        pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
        pow = cell2mat(pow); %organize
        highnumfrex = length(mulFrex); 
        pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
        pow = squeeze(mean(pow, 3)); %take the mean over frequencies
    
        HFB1 = pow; 
        clear pow
    
        %CHAN 2
        [pow2, mulTim, mulFrex] = getChanMultiTF(chanDat2.retOn, highfrex, chanDat.fsample, chanDat.retOtim, 1);  
        pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow2,3), 'UniformOutput',false ); %z-score
        pow2 = cell2mat(pow2); %organize
        highnumfrex = length(mulFrex); 
        pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/highnumfrex, []); %organize
        pow2 = squeeze(mean(pow2, 3)); %take the mean over frequencies
        

        disp('retrieval: ')
        outCluStats2(chan,:,:,:) = getLL(HFB1, pow2, miss_on, hit_on, chanDat.retOtim, leadLagRetTim);
    
        x = num2str(toc/60); 
        disp(['.......................................................' x])
        leadLag.inclChan = includedChans; 
        leadLag.encTim = leadLagEncTim; 
        leadLag.retTim = leadLagRetTim; 
        leadLag.subMem = outCluStats; 
        leadLag.retMem = outCluStats2; 
        chanDat.leadLag4 = leadLag; 
        disp('interim save')
        save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat');
        end
    end

    leadLag.inclChan = includedChans; 
    leadLag.encTim = leadLagEncTim; 
    leadLag.retTim = leadLagRetTim; 
    leadLag.subMem = outCluStats; 
    leadLag.retMem = outCluStats2; 
    chanDat.leadLag4 = leadLag; 
    
    clear HFB1 HFB2 pow2 leadLag subMiss subHit miss_on hit_on 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
    
    else %there was no work to be done
        disp('channel leadLag already complete!')
    end
    else
        chanDat.leadLag4 = 1; 
        disp('non-reactive channel save')
        save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
        disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
    end
else
    disp('already done with leadLag')

end


%% Lead lag analysis




if ~isfield(chanDat, 'leadLag3')
    if isfield(chanDat, 'leadLag')
        chanDat = rmfield(chanDat, 'leadLag');
    end
  
% 
%     if isfield(chanDat, 'leadLag3') %check for previous work! 
%         chanDat = rmfield(chanDat, 'leadLag3');
%         leadLag = struct; 
%         reactive = reactiveTest(chanDat.HFB);
%         includedChans = []; 
%         start = 1; 
%         outCluStats = nan(length(chanFiles), 2, 30, 15); %subsequent memory
%         outCluStats2 = nan(length(chanFiles), 2, 30, 15); %retrieval
% %         leadLag = chanDat.leadLag3; 
%         
% %         if isstruct(leadLag) %is this a reactive channel with previous leadLag calculation? 
% %             start = max(leadLag.inclChan); %start value for looping below
% %             includedChans = leadLag.inclChan; %previous work done on these channels
% %             outCluStats = leadLag.subMem; 
% %             outCluStats2 = leadLag.retMem; 
% %             reactive = [1,1,1,1]; 
% %         else
% %             reactive = [0,0,0,0]; %if it's not a struct, then this channel itself is not reactive, so skip it
% %         end
% 
%     else %if no work has been done, then start from scratch here

        leadLag = struct; 
        reactive = reactiveTest(chanDat.HFB);
        includedChans = []; 
        start = 1; 
        outCluStats = nan(length(chanFiles), 2, 200, 15); %subsequent memory
        outCluStats2 = nan(length(chanFiles), 2, 200, 15); %retrieval
%     end



    if sum(reactive==1)>0
    %chan X time X offSet
    leadLagEncTim = chanDat.enctim(501:25:end-500);
%     subMiss = zeros([length(chanFiles), length(leadLagEncTim), length([-150:150])]);
%     subHit = subMiss; 
    
    leadLagRetTim = chanDat.retOtim(501:25:end-500); 
%     miss_on = zeros([length(chanFiles), length(leadLagRetTim), length([-150:150])]);
%     hit_on = miss_on; 

    
    %trial index values
    subMiss = find(chanDat.use & chanDat.misses); 
    subHit = find(chanDat.use & chanDat.hits); 
    miss_on = find(chanDat.retInfo(:,1)==2);
    hit_on = find(chanDat.retInfo(:,1)==1); 
    %out stats will be chan X pos/neg X clust X stat: 
    %stat 1: num points
    %stat 2: mean time
    %stat 3: min time
    %stat 4: max time
    %stat 5: mean LL
    %stat 6: min LL
    %stat 7: max LL
    %stat 8: mean correlation Hit
    %stat 9: min correlation Hit
    %stat 10: max correlation Hit
    %stat 11: median correlation Hit
    %stat 12: mean correlation Miss
    %stat 13: min correlation Miss
    %stat 14: max correlation Miss
    %stat 15: median correlation Miss






    if start<length(chanFiles) %is there even any work left to be done? 

    for chan = start:length(chanFiles)
        disp(['leadLag analysis with channel: ' num2str(chan) ' of ' num2str(length(chanFiles))])
        tic

        chanDat2 = load([chanFiles(chan).folder '/CHANRAW/' chanFiles(chan).name]).chanDat; 
        chanDat2.HFB = getHFB(chanDat2, highfrex);
        reactive2 = reactiveTest(chanDat2.HFB);
        if sum(reactive2==1)>0
            includedChans = [includedChans chan]; 
        %need to grab the HFB data at higher resolution, so calculate from
        %scratch 
    
        %ENCODING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHAN 1
        [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.enc, highfrex, chanDat.fsample, chanDat.enctim, 1);  
        pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
        pow = cell2mat(pow); %organize
        highnumfrex = length(mulFrex); 
        pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
        pow = squeeze(mean(pow, 3)); %take the mean over frequencies
    
        HFB1 = pow; 
        clear pow
    
        %CHAN 2
        [pow2, mulTim, mulFrex] = getChanMultiTF(chanDat2.enc, highfrex, chanDat.fsample, chanDat.enctim, 1);  
        pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow2,3), 'UniformOutput',false ); %z-score
        pow2 = cell2mat(pow2); %organize
        highnumfrex = length(mulFrex); 
        pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/highnumfrex, []); %organize
        pow2 = squeeze(mean(pow2, 3)); %take the mean over frequencies
        
        disp('encoding: ')
        outCluStats = conditionCluTest(HFB1, pow2, subMiss, subHit, outCluStats, chanDat.enctim, leadLagEncTim, chan);


      


    
        %RERTRIEVAL DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHAN 1
        [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retOn, highfrex, chanDat.fsample, chanDat.retOtim, 1);  
        pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
        pow = cell2mat(pow); %organize
        highnumfrex = length(mulFrex); 
        pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
        pow = squeeze(mean(pow, 3)); %take the mean over frequencies
    
        HFB1 = pow; 
        clear pow
    
        %CHAN 2
        [pow2, mulTim, mulFrex] = getChanMultiTF(chanDat2.retOn, highfrex, chanDat.fsample, chanDat.retOtim, 1);  
        pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow2,3), 'UniformOutput',false ); %z-score
        pow2 = cell2mat(pow2); %organize
        highnumfrex = length(mulFrex); 
        pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/highnumfrex, []); %organize
        pow2 = squeeze(mean(pow2, 3)); %take the mean over frequencies
        

        disp('retrieval: ')
        outCluStats2 = conditionCluTest(HFB1, pow2, miss_on, hit_on, outCluStats2, chanDat.retOtim, leadLagRetTim, chan);
    
        x = num2str(toc/60); 
        disp(['.......................................................' x])
        leadLag.inclChan = includedChans; 
        leadLag.encTim = leadLagEncTim; 
        leadLag.retTim = leadLagRetTim; 
        leadLag.subMem = outCluStats; 
        leadLag.retMem = outCluStats2; 
        chanDat.leadLag3 = leadLag; 
        disp('interim save')
        save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat');
        end
    end

    leadLag.inclChan = includedChans; 
    leadLag.encTim = leadLagEncTim; 
    leadLag.retTim = leadLagRetTim; 
    leadLag.subMem = outCluStats; 
    leadLag.retMem = outCluStats2; 
    chanDat.leadLag3 = leadLag; 
    
    clear HFB1 HFB2 pow2 leadLag subMiss subHit miss_on hit_on 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
    
    else %there was no work to be done
        disp('channel leadLag already complete!')
    end
    else
        chanDat.leadLag3 = 1; 
        disp('non-reactive channel save')
        save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
        disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
    end
else
    disp('already done with leadLag')

end




%% time frequency decomposition, extract TF summaries for target trial types: 

% subsequent hit / subsequent miss (encoding data)
% hit / miss / CR / FA (retrieval locked to onset data)
% hit / miss / CR / FA (retrieval locked to response data)

if ~isfield(chanDat, 'TF')
    disp('working on TF')
    %to keep size down, don't put large variables into the chanDat struct!
    TFout = struct;
    %ENCODING DATA: ***********************************************************
    pow = log10(abs(getChanTrialTF(chanDat.enc, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctim>=-450,1), find(chanDat.enctim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean misses: 
    TFout.subMiss = squeeze(mean(pow(:,chanDat.use & chanDat.misses, :), 2)); 
    %get mean hits: 
    TFout.subHit = squeeze(mean(pow(:,chanDat.use & chanDat.hits, :), 2));
    %clean up
    clear pow
    disp('encoding done')

    %ENCODING DATA RESPONSE: ***********************************************************
    pow = log10(abs(getChanTrialTF(chanDat.encRT, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctimRT>=-2000,1), find(chanDat.enctimRT>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean misses: 
    TFout.subMissRT = squeeze(mean(pow(:,chanDat.use & chanDat.misses, :), 2)); 
    %get mean hits: 
    TFout.subHitRT = squeeze(mean(pow(:,chanDat.use & chanDat.hits, :), 2));
    %clean up
    clear pow
    disp('encoding RT done')

    %RETRIEVAL STIM ONSET: ****************************************************
    pow = log10(abs(getChanTrialTF(chanDat.retOn, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retOtim>=-450,1), find(chanDat.retOtim>=-50,1)]), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean hit: 
    TFout.hit_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==1, :), 2)); 
    %get mean CRs: 
    TFout.cr_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==3, :), 2));
    %get mean miss: 
    TFout.miss_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==2, :), 2));
    %get mean FA: 
    TFout.fa_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==4, :), 2));
    %clean up
    clear pow
    disp('retrieval 1 done')
    
    %RETRIEVAL RESPONSE LOCKED: ***********************************************
    pow = log10(abs(getChanTrialTF(chanDat.retRT, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retRtim>=-2000,1), find(chanDat.retRtim>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean hit: 
    TFout.hit_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==1, :), 2)); 
    %get mean CRs: 
    TFout.cr_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==3, :), 2));
    %get mean miss: 
    TFout.miss_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==2, :), 2));
    %get mean FA: 
    TFout.fa_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==4, :), 2));
    %clean up
    clear pow
    disp('retrieval 2 done')
    

    %reduce the size! 
    multim = chanDat.HFB.encMulTim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.enctim,1), multim);
    TFout.subMiss = TFout.subMiss(encIDX, :); 
    TFout.subHit = TFout.subHit(encIDX, :);

    multim = chanDat.HFB.encRT_tim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.enctimRT,1), multim);
    TFout.subMissRT = TFout.subMissRT(encIDX, :); 
    TFout.subHitRT = TFout.subHitRT(encIDX, :);

    multim = chanDat.HFB.onMulTim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.retOtim,1), multim);
    TFout.hit_on = TFout.hit_on(encIDX, :); 
    TFout.miss_on = TFout.miss_on(encIDX, :);
    TFout.fa_on = TFout.fa_on(encIDX, :); 
    TFout.cr_on = TFout.cr_on(encIDX, :);

    multim = chanDat.HFB.rtMulTim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.retRtim,1), multim);
    TFout.hit_rt = TFout.hit_rt(encIDX, :); 
    TFout.miss_rt = TFout.miss_rt(encIDX, :);
    TFout.fa_rt = TFout.fa_rt(encIDX, :); 
    TFout.cr_rt = TFout.cr_rt(encIDX, :);





    chanDat.TF = TFout; 
    
    clear TFout 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
else
    disp('TF already done, skipping')
end

%% get ISPC and PPC values 

if ~isfield(chanDat, 'ISPC')
    disp('working on ISPC')
    ISPCout = struct; 
    %store the downsample index (di) 
    multim = chanDat.HFB.encMulTim; 
    ISPCout.encdi = arrayfun(@(x) find(x<=chanDat.enctim,1), multim);
    
    multim = chanDat.HFB.encRT_tim; 
    ISPCout.encRdi = arrayfun(@(x) find(x<=chanDat.enctimRT,1), multim);

    multim = chanDat.HFB.onMulTim; 
    ISPCout.ondi = arrayfun(@(x) find(x<=chanDat.retOtim,1), multim);
    
    multim = chanDat.HFB.rtMulTim; 
    ISPCout.rtdi = arrayfun(@(x) find(x<=chanDat.retRtim,1), multim);
    
    %preallocate: 
    %channels X time X frequencies X ISPC/PPC
    frex = logspace(log10(2), log10(25), 20); 
    numfrex = length(frex); 
    ISPCout.subMiss = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 4); 
    ISPCout.subHit = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 4); 

    ISPCout.subMissR = zeros(length(chanFiles), length(ISPCout.encRdi), length(frex), 4); 
    ISPCout.subHitR = zeros(length(chanFiles), length(ISPCout.encRdi), length(frex), 4); 
    
    ISPCout.hit_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.cr_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.miss_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.fa_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);


    ISPCout.hit_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.cr_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.miss_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.fa_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);


    %will need to loop channels
    %NOTE: all trial types must have at least two trials! 
    for chan = 1:length(chanFiles)
        tic
%         chan
        if chan > idx %don't do repeat work! 
        chanDat2 = load([chanFiles(chan).folder '/CHANRAW/' chanFiles(chan).name]).chanDat; 

        %ENCODING DATA: ***********************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMiss(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.misses);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHit(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.hits);
        end

        %ENCODING DATA RT: ***********************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMissR(chan,:,:,:) = getChanISPC(chanDat.encRT, chanDat2.encRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encRdi, chanDat.use & chanDat.misses);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHitR(chan,:,:,:) = getChanISPC(chanDat.encRT, chanDat2.encRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encRdi, chanDat.use & chanDat.hits);
        end

        %RETRIEVAL STIM ONSET: ****************************************************
        if sum(chanDat.retInfo(:,1)==1) > 1
        ISPCout.hit_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==1);
        end
        if sum(chanDat.retInfo(:,1)==3) > 1
        ISPCout.cr_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==3);
        end
        if sum(chanDat.retInfo(:,1)==2) > 1
        ISPCout.miss_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==2);
        end
        if sum(chanDat.retInfo(:,1)==4) > 1
        ISPCout.fa_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==4);
        end
        
        %RETRIEVAL RESPONSE LOCKED: ****************************************************
        if sum(chanDat.retInfo(:,1)==1) > 1
        ISPCout.hit_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==1);
        end
        if sum(chanDat.retInfo(:,1)==3) > 1
        ISPCout.cr_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==3);
        end
        if sum(chanDat.retInfo(:,1)==2) > 1
        ISPCout.miss_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==2);
        end
        if sum(chanDat.retInfo(:,1)==4) > 1
        ISPCout.fa_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==4);
        end
        

        end
        disp(['channel: ' num2str(chan) ' took ' num2str(round(toc/60,1)) ' minutes'])
    end
    chanDat.ISPC = ISPCout; 
%     chanDat = rmfield(chanDat, 'sizeReduce'); 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])



else
    disp('connectivity already done, skipping')
end












end








%% files are getting too large. Need to shrink them by grabbing epoch means for connectivity
%shrinking is accomplished by averaging across 150ms temporal epochs

% if ~isfield(chanDat, 'sizeReduce')
%     disp('shrinking connectivity data')
%     temp = chanDat; 
%     %choose some epochs: 
%     chanDat.encepoch = -450:150:3001;
%     chanDat.onepoch = -450:150:2001; 
%     chanDat.rtepoch = -1950:150:500; 
%     
%     %ENCODING DATA: ***********************************************************
%     %shrink connectivity dat: 
%     tim = chanDat.enctim; 
%     di = chanDat.ISPCout.encdi; 
%     tim = tim(di); 
%     epoch = chanDat.encepoch; 
%     chanDat.ISPCout.subMiss = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.subMiss(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.subHit = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.subHit(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     %shrink power dat: 
%     tim = chanDat.enctim; 
%     chanDat.TFout.subMiss = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.subMiss(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.subHit = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.subHit(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
% 
% 
% 
%     %RETRIEVAL STIM ONSET: ***********************************************************
%     %shrink connectivity dat: 
%     tim = chanDat.retOtim; 
%     di = chanDat.ISPCout.ondi; 
%     tim = tim(di); 
%     epoch = chanDat.onepoch; 
%     chanDat.ISPCout.hit_on = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.hit_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.cr_on = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.cr_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.miss_on = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.miss_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.fa_on = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.fa_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     %shrink power dat: 
%     tim = chanDat.retOtim; 
%     chanDat.TFout.hit_on = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.hit_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.cr_on = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.cr_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.miss_on = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.miss_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.fa_on = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.fa_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%         
% 
%     %RETRIEVAL RESPONSE LOCKED: ****************************************************
%     %shrink connectivity dat: 
%     tim = chanDat.retRtim; 
%     di = chanDat.ISPCout.rtdi; 
%     tim = tim(di); 
%     epoch = chanDat.rtepoch; 
%     chanDat.ISPCout.hit_rt = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.hit_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.cr_rt = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.cr_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.miss_rt = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.miss_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     chanDat.ISPCout.fa_rt = cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.ISPCout.fa_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
%                 1:length(epoch)-1 , 'uniformoutput', false));
%     %shrink power dat: 
%     tim = chanDat.retRtim; 
%     chanDat.TFout.hit_rt = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.hit_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.cr_rt = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.cr_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.miss_rt = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.miss_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
%     chanDat.TFout.fa_rt = reshape(cell2mat(arrayfun(@(x) ...
%                 mean(chanDat.TFout.fa_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
%                 1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
% 
%     chanDat.sizeReduce = true; 
%     disp('attempting saving')
%     save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
%     disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
% 
% 
% 
% else
%     disp('size reduction already done')
% end













