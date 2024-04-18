function [] = singleChanPipeline(chanFiles, idx, subFiles, codePre)

%% set frequency parameters
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = logspace(log10(3),log10(10),numfrex)./(2*pi*frex);

highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = logspace(log10(10),log10(20),highnumfrex)./(2*pi*highfrex);

%% load the data
folderName = split(chanFiles(idx).folder, 'CHANRAW'); 
folderName = folderName{1}; 
try %try loading the processed file
    
    chanDat = load([folderName 'finished/' chanFiles(idx).name]).chanDat; 
catch
    chanDat = load([chanFiles(idx).folder '/' chanFiles(idx).name]).chanDat; % go raw if it's not working!
end

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])


%% time points! hard code

chanDat.enctim = [-1000:3500];
chanDat.enctimRT = [-2000:500];
chanDat.retOtim = [-1000:3000];
chanDat.retRtim = [-2000:500];




%% High frequency Broadband 

%note hardcoded baselines: encoding: -450 : -50   ms
%                          retOn   : -450 : -50   ms
%                          retRT   : -2000: -1600 ms

if ~isfield(chanDat, 'reactiveRes')
    HFB = getHFB(chanDat, highfrex); 

    chanDat.HFB = HFB; 
    chanDat.reactiveRes = reactiveTest_100(chanDat.HFB);
    % downsample for size
    HFB = chanDat.HFB;
    HFB_names = fieldnames(HFB); 
    
    datidx = {[1,2], [5,6], [9,10,11,12], [15,16,17,18]}; 
    timidx = [3,7,13,19]; 

    for dati = 1:4
        curtim = HFB.(HFB_names{timidx(dati)});
        HFB.(HFB_names{timidx(dati)}) = curtim(1:5:end);
        curDi = datidx{dati}; 
        for fi = 1:length(curDi)
            cur = HFB.(HFB_names{curDi(fi)});

            HFB.(HFB_names{curDi(fi)}) = cur(1:5:end,:);
        end
    end
    chanDat.HFB = HFB; 
    clear HFB 
    % done downsample 
    disp('attempting saving')
    save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name])

else
    disp('HFB already done')
end

%% get the HFB latencies and save them in order to reference other variables to them

if ~isfield(chanDat, 'gjfds')
    HFB_lat = struct; 
   

    %encoding
    HFB_lat.subHit = gausLat(chanDat.HFB.subHit, ...
        chanDat.HFB.encMulTim, ...
        chanDat.encInfo(chanDat.use & chanDat.hits, 4), ...
        1);

    HFB_lat.subMiss = gausLat(chanDat.HFB.subMiss, ...
        chanDat.HFB.encMulTim, ...
        chanDat.encInfo(chanDat.use & chanDat.misses, 4), ...
        1);

    %retrieval
    triali = chanDat.retInfo(:,1)==1;
    HFB_lat.retHit = gausLat(chanDat.HFB.hit_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);

    triali = chanDat.retInfo(:,1)==2;
    HFB_lat.retMiss = gausLat(chanDat.HFB.miss_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);
    
    triali = chanDat.retInfo(:,1)==3;
    HFB_lat.retCR = gausLat(chanDat.HFB.cr_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);

    triali = chanDat.retInfo(:,1)==4;
    HFB_lat.retFA = gausLat(chanDat.HFB.fa_on, ...
        chanDat.HFB.onMulTim, ...
        chanDat.retInfo(triali, 3), ...
        1);

    chanDat.HFB_lat = HFB_lat; 



end






%% single trial time frequency

if ~isfield(chanDat, 'TF2')
    disp('working on TF')
    %to keep size down, don't put large variables into the chanDat struct!
    TFout = struct;
    %ENCODING DATA: ***********************************************************
    pow = getChanTrialTF(chanDat.enc, frex, numfrex, stds, chanDat.fsample); %get power time series for all trials/frequencies
    phase = angle(pow); 
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctim>=-450,1), find(chanDat.enctim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean misses: 
    TFout.subMiss = squeeze(pow(:,chanDat.use & chanDat.misses, :)); 
    TFout.subMiss_p = squeeze(phase(:,chanDat.use & chanDat.misses, :)); 
    %get mean hits: 
    TFout.subHit = squeeze(pow(:,chanDat.use & chanDat.hits, :));
    TFout.subHit_p = squeeze(phase(:,chanDat.use & chanDat.hits, :));
    %clean up
    clear pow
    disp('encoding done')

    %ENCODING DATA RESPONSE: *********************************************************** 
    pow = getChanTrialTF(chanDat.encRT, frex, numfrex, stds, chanDat.fsample); %get power time series for all trials/frequencies
    phase = angle(pow);
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctimRT>=-2000,1), find(chanDat.enctimRT>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean misses: 
    TFout.subMissRT = squeeze(pow(:,chanDat.use & chanDat.misses, :)); 
    TFout.subMissRT_p = squeeze(phase(:,chanDat.use & chanDat.misses, :)); 
    %get mean hits: 
    TFout.subHitRT = squeeze(pow(:,chanDat.use & chanDat.hits, :));
    TFout.subHitRT_p = squeeze(phase(:,chanDat.use & chanDat.hits, :));
    %clean up
    clear pow
    disp('encoding RT done')

    %RETRIEVAL STIM ONSET: ****************************************************
    pow = getChanTrialTF(chanDat.retOn, frex, numfrex, stds, chanDat.fsample); %get power time series for all trials/frequencies
    phase = angle(pow);
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retOtim>=-450,1), find(chanDat.retOtim>=-50,1)]), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean hit: 
    TFout.hit_on = squeeze(pow(:,chanDat.retInfo(:,1)==1, :)); 
    TFout.hit_on_p = squeeze(phase(:,chanDat.retInfo(:,1)==1, :)); 
    %get mean CRs: 
    TFout.cr_on = squeeze(pow(:,chanDat.retInfo(:,1)==3, :));
    TFout.cr_on_p = squeeze(phase(:,chanDat.retInfo(:,1)==3, :));
    %get mean miss: 
    TFout.miss_on = squeeze(pow(:,chanDat.retInfo(:,1)==2 , :));
    TFout.miss_on_p = squeeze(phase(:,chanDat.retInfo(:,1)==2 , :));
    %get mean FA: 
    TFout.fa_on = squeeze(pow(:,chanDat.retInfo(:,1)==4 , :));
    TFout.fa_on_p = squeeze(phase(:,chanDat.retInfo(:,1)==4 , :));
    %clean up
    clear pow
    disp('retrieval 1 done')
    
    %RETRIEVAL RESPONSE LOCKED: ***********************************************
    pow = getChanTrialTF(chanDat.retRT, frex, numfrex, stds, chanDat.fsample); %get power time series for all trials/frequencies
    phase = angle(pow);
    pow = abs(pow).^2;
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retRtim>=-2000,1), find(chanDat.retRtim>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean hit: 
    TFout.hit_rt = squeeze(pow(:,chanDat.retInfo(:,1)==1 , :)); 
    TFout.hit_rt_p = squeeze(phase(:,chanDat.retInfo(:,1)==1 , :));
    %get mean CRs: 
    TFout.cr_rt = squeeze(pow(:,chanDat.retInfo(:,1)==3 , :));
    TFout.cr_rt_p = squeeze(phase(:,chanDat.retInfo(:,1)==3 , :));
    %get mean miss: 
    TFout.miss_rt = squeeze(pow(:,chanDat.retInfo(:,1)==2 , :));
    TFout.miss_rt_p = squeeze(phase(:,chanDat.retInfo(:,1)==2 , :));
    %get mean FA: 
    TFout.fa_rt = squeeze(pow(:,chanDat.retInfo(:,1)==4 , :));
    TFout.fa_rt_p = squeeze(phase(:,chanDat.retInfo(:,1)==4 , :));
    %clean up
    clear pow
    disp('retrieval 2 done')
    

    %reduce the size! 
    multim = chanDat.HFB.encMulTim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.enctim,1), multim);
    TFout.subMiss = TFout.subMiss(encIDX, :, :); 
    TFout.subHit = TFout.subHit(encIDX, :, :);
    TFout.subMiss_p = TFout.subMiss_p(encIDX, :, :); 
    TFout.subHit_p = TFout.subHit_p(encIDX, :, :);

    multim = chanDat.HFB.encRT_tim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.enctimRT,1), multim);
    TFout.subMissRT = TFout.subMissRT(encIDX, :, :); 
    TFout.subHitRT = TFout.subHitRT(encIDX, :, : );
    TFout.subMissRT_p = TFout.subMissRT_p(encIDX, :, :); 
    TFout.subHitRT_p = TFout.subHitRT_p(encIDX, :, : );

    multim = chanDat.HFB.onMulTim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.retOtim,1), multim);
    TFout.hit_on = TFout.hit_on(encIDX, :,:); 
    TFout.miss_on = TFout.miss_on(encIDX, :,:);
    TFout.fa_on = TFout.fa_on(encIDX, :,:); 
    TFout.cr_on = TFout.cr_on(encIDX, :,:);
    TFout.hit_on_p = TFout.hit_on_p(encIDX, :,:); 
    TFout.miss_on_p = TFout.miss_on_p(encIDX, :,:);
    TFout.fa_on_p = TFout.fa_on_p(encIDX, :,:); 
    TFout.cr_on_p = TFout.cr_on_p(encIDX, :,:);

    multim = chanDat.HFB.rtMulTim; 
    encIDX = arrayfun(@(x) find(x<=chanDat.retRtim,1), multim);
    TFout.hit_rt = TFout.hit_rt(encIDX, :,:); 
    TFout.miss_rt = TFout.miss_rt(encIDX, :,:);
    TFout.fa_rt = TFout.fa_rt(encIDX, :,:); 
    TFout.cr_rt = TFout.cr_rt(encIDX, :,:);
    TFout.hit_rt_p = TFout.hit_rt_p(encIDX, :,:); 
    TFout.miss_rt_p = TFout.miss_rt_p(encIDX, :,:);
    TFout.fa_rt_p = TFout.fa_rt_p(encIDX, :,:); 
    TFout.cr_rt_p = TFout.cr_rt_p(encIDX, :,:);



    chanDat.TF2 = TFout; 
    
    clear TFout 
     disp('attempting saving')
    save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name])
else
    disp('TF already done, skipping')
end






%% get ISPC and PPC values 

if ~isfield(chanDat, 'hjk')
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
    stds = logspace(log10(3),log10(5),numfrex)./(2*pi*frex);

    ISPCout.subMiss = zeros(length(subFiles), length(ISPCout.encdi), length(frex), 4); 
    ISPCout.subHit = zeros(length(subFiles), length(ISPCout.encdi), length(frex), 4); 
    ISPCout.subMissRT = zeros(length(subFiles), length(ISPCout.encRdi), length(frex), 4);
    ISPCout.subHitRT = zeros(length(subFiles), length(ISPCout.encRdi), length(frex), 4);
    
    
    ISPCout.hit_on = zeros(length(subFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.miss_on = zeros(length(subFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.cr_on = zeros(length(subFiles), length(ISPCout.ondi), length(frex), 4);
    ISPCout.miss_onRT = zeros(length(subFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.hit_onRT = zeros(length(subFiles), length(ISPCout.rtdi), length(frex), 4);
    ISPCout.cr_onRT = zeros(length(subFiles), length(ISPCout.rtdi), length(frex), 4);




    %will need to loop channels
    %NOTE: all trial types must have at least two trials! 
    for chan = 1:length(subFiles)
        tic
       
        chanDat2 = load([subFiles(chan).folder '/' subFiles(chan).name]).chanDat; 

        %ENCODING DATA: ***********************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMiss(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, ...
            chanDat.use & chanDat.misses, chanDat.HFB_lat.subMiss, ...
            chanDat.HFB.encMulTim, true);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHit(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, ...
            chanDat.use & chanDat.hits, chanDat.HFB_lat.subHit, ...
            chanDat.HFB.encMulTim, true);
        end

     
        %ENCODING DATA BEHAVIOR RESPONSE: *************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMissRT(chan,:,:,:) = getChanISPC(chanDat.encRT, chanDat2.encRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encRdi, ...
            chanDat.use & chanDat.misses, chanDat.HFB_lat.subMiss, ...
            chanDat.HFB.encRT_tim, false);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHitRT(chan,:,:,:) = getChanISPC(chanDat.encRT, chanDat2.encRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encRdi, ...
            chanDat.use & chanDat.hits, chanDat.HFB_lat.subHit, ...
            chanDat.HFB.encRT_tim, false);
        end



        %RETRIEVAL STIM ONSET: ****************************************************
        if sum(chanDat.retInfo(:,1)==1) > 1
        ISPCout.hit_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, ...
            chanDat.retInfo(:,1)==1, chanDat.HFB_lat.retHit, ...
            chanDat.HFB.onMulTim, true);
        end
        if sum(chanDat.retInfo(:,1)==2) > 1
        ISPCout.miss_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, ...
            chanDat.retInfo(:,1)==2, chanDat.HFB_lat.retMiss, ...
            chanDat.HFB.onMulTim, true);
        end
        if sum(chanDat.retInfo(:,1)==3) > 1
        ISPCout.cr_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, ...
            chanDat.retInfo(:,1)==3, chanDat.HFB_lat.retCR, ...
            chanDat.HFB.onMulTim, true);
        end
     
        
        %RETRIEVAL BEHAVIOR RESPONSE: ****************************************************
        if sum(chanDat.retInfo(:,1)==1) > 1
        ISPCout.hit_onRT(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, ...
            chanDat.retInfo(:,1)==1, chanDat.HFB_lat.retHit, ...
            chanDat.HFB.rtMulTim, false);
        end
        if sum(chanDat.retInfo(:,1)==2) > 1
        ISPCout.miss_onRT(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, ...
            chanDat.retInfo(:,1)==2, chanDat.HFB_lat.retMiss, ...
            chanDat.HFB.rtMulTim, false);
        end
        if sum(chanDat.retInfo(:,1)==3) > 1
        ISPCout.cr_onRT(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, ...
            chanDat.retInfo(:,1)==3, chanDat.HFB_lat.retCR, ...
            chanDat.HFB.rtMulTim, false);
        end
     
        

        
        disp(['channel: ' num2str(chan) ' took ' num2str(round(toc/60,1)) ' minutes'])
    end
    chanDat.ISPC = ISPCout; 
    disp('attempting saving')
    save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name])



else
    disp('connectivity already done, skipping')
end




%% final save out

save([folderName 'finished/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' folderName 'finished/' chanFiles(idx).name]) 
% delete([chanFiles(idx).folder '/' chanFiles(idx).name]) 







end















