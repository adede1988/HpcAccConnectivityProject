function [] = singleChanPipeline(chanFiles, idx)

%% set frequency parameters
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);

%% load the data
chanDat = load([chanFiles(idx).folder '/' chanFiles(idx).name]).chanDat; 

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])

%% needs time points! hard code
chanDat.enctim = [-1000:3500];
chanDat.retOtim = [-1000:3000];
chanDat.retRtim = [-2000:500];


%% time frequency decomposition, extract TF summaries for target trial types: 
% subsequent hit / subsequent miss (encoding data)
% hit / miss / CR / FA (retrieval locked to onset data)
% hit / miss / CR / FA (retrieval locked to response data)

if ~isfield(chanDat, 'TFout')
    disp('working on TF')
    %to keep size down, don't put large variables into the chanDat struct!
    TFout = struct;
    %ENCODING DATA: ***********************************************************
    pow = abs(getChanTrialTF(chanDat.enc, frex, numfrex, stds, chanDat.fsample)).^2; %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
    %get mean misses: 
    TFout.subMiss = squeeze(mean(pow(:,chanDat.use & chanDat.misses, :), 2)); 
    %get mean hits: 
    TFout.subHit = squeeze(mean(pow(:,chanDat.use & chanDat.hits, :), 2));
    %clean up
    clear pow
    disp('encoding done')

    %RETRIEVAL STIM ONSET: ****************************************************
    pow = abs(getChanTrialTF(chanDat.retOn, frex, numfrex, stds, chanDat.fsample)).^2; %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
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
    pow = abs(getChanTrialTF(chanDat.retRT, frex, numfrex, stds, chanDat.fsample)).^2; %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
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
    
    chanDat.TFout = TFout; 
    
    clear TFout 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
else
    disp('TF already done, skipping')
end


if ~isfield(chanDat, 'ISPCout')
    disp('working on ISPC')
    ISPCout = struct; 
    %store the downsample index (di) 
    ISPCout.encdi = 1:20:length(chanDat.enctim);
    ISPCout.ondi = 1:20:length(chanDat.retOtim); 
    ISPCout.rtdi = 1:20:length(chanDat.retRtim); 
    %channels X time X frequencies X ISPC/PPC
    ISPCout.subMiss = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 2); 
    ISPCout.subHit = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 2); 
    
    ISPCout.hit_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 2);
    ISPCout.cr_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 2);
    ISPCout.miss_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 2);
    ISPCout.fa_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 2);


    ISPCout.hit_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 2);
    ISPCout.cr_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 2);
    ISPCout.miss_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 2);
    ISPCout.fa_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 2);


    %will need to loop channels
    for chan = 1:length(chanFiles)
        chan
        if chan ~= idx %skip self connection
        chanDat2 = load([chanFiles(chan).folder '/' chanFiles(chan).name]).chanDat; 

        %ENCODING DATA: ***********************************************************
        ISPCout.subMiss(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.misses);
        ISPCout.subHit(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.hits);

        %RETRIEVAL STIM ONSET: ****************************************************
        ISPCout.hit_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==1);
        ISPCout.cr_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==3);
        ISPCout.miss_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==2);
        ISPCout.fa_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==4);
        
        %RETRIEVAL RESPONSE LOCKED: ****************************************************
        ISPCout.hit_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==1);
        ISPCout.cr_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==3);
        ISPCout.miss_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==2);
        ISPCout.fa_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==4);
        

        end
    end
    
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])



else
    disp('connectivity already done, skipping')
end














end