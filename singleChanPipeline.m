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

if ~isfield(chanDat, 'ISPCboot')
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

%% get ISPC and PPC values 

if ~isfield(chanDat, 'ISPCboot')
    disp('working on ISPC')
    ISPCout = struct; 
    %store the downsample index (di) 
    ISPCout.encdi = 1:20:length(chanDat.enctim);
    ISPCout.ondi = 1:20:length(chanDat.retOtim); 
    ISPCout.rtdi = 1:20:length(chanDat.retRtim); 
    
    %preallocate: 
    %channels X time X frequencies X ISPC/PPC
    ISPCout.subMiss = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 4); 
    ISPCout.subHit = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 4); 
    
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
        chanDat2 = load([chanFiles(chan).folder '/' chanFiles(chan).name]).chanDat; 

        %ENCODING DATA: ***********************************************************
        if sum(chanDat.use & chanDat.misses)>1
        ISPCout.subMiss(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.misses);
        end
        if sum(chanDat.use & chanDat.hits)>1
        ISPCout.subHit(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
            frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.hits);
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
    chanDat.ISPCout = ISPCout; 
    chanDat.ISPCboot = true; 
    chanDat = rmfield(chanDat, 'sizeReduce'); 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])



else
    disp('connectivity already done, skipping')
end


%% files are getting too large. Need to shrink them by grabbing epoch means for connectivity
%shrinking is accomplished by averaging across 150ms temporal epochs

if ~isfield(chanDat, 'sizeReduce')
    disp('shrinking connectivity data')
    temp = chanDat; 
    %choose some epochs: 
    chanDat.encepoch = -450:150:3001;
    chanDat.onepoch = -450:150:2001; 
    chanDat.rtepoch = -2000:150:500; 
    
    %ENCODING DATA: ***********************************************************
    %shrink connectivity dat: 
    tim = chanDat.enctim; 
    di = chanDat.ISPCout.encdi; 
    tim = tim(di); 
    epoch = chanDat.encepoch; 
    chanDat.ISPCout.subMiss = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.subMiss(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.subHit = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.subHit(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    %shrink power dat: 
    tim = chanDat.enctim; 
    chanDat.TFout.subMiss = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.subMiss(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.subHit = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.subHit(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';



    %RETRIEVAL STIM ONSET: ***********************************************************
    %shrink connectivity dat: 
    tim = chanDat.retOtim; 
    di = chanDat.ISPCout.ondi; 
    tim = tim(di); 
    epoch = chanDat.onepoch; 
    chanDat.ISPCout.hit_on = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.hit_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.cr_on = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.cr_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.miss_on = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.miss_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.fa_on = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.fa_on(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    %shrink power dat: 
    tim = chanDat.retOtim; 
    chanDat.TFout.hit_on = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.hit_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.cr_on = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.cr_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.miss_on = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.miss_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.fa_on = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.fa_on(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
        

    %RETRIEVAL RESPONSE LOCKED: ****************************************************
    %shrink connectivity dat: 
    tim = chanDat.retRtim; 
    di = chanDat.ISPCout.rtdi; 
    tim = tim(di); 
    epoch = chanDat.rtepoch; 
    chanDat.ISPCout.hit_rt = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.hit_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.cr_rt = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.cr_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.miss_rt = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.miss_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    chanDat.ISPCout.fa_rt = cell2mat(arrayfun(@(x) ...
                mean(chanDat.ISPCout.fa_rt(:,tim>=epoch(x) & tim < epoch(x+1), :,:), 2),...
                1:length(epoch)-1 , 'uniformoutput', false));
    %shrink power dat: 
    tim = chanDat.retRtim; 
    chanDat.TFout.hit_rt = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.hit_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.cr_rt = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.cr_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.miss_rt = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.miss_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';
    chanDat.TFout.fa_rt = reshape(cell2mat(arrayfun(@(x) ...
                mean(chanDat.TFout.fa_rt(tim>=epoch(x) & tim < epoch(x+1), :), 1),...
                1:length(epoch)-1, 'uniformoutput', false )), [length(frex), length(epoch)-1 ])';

    chanDat.sizeReduce = true; 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])



else
    disp('size reduction already done')
end























end