function [] = singleChanPipeline(chanFiles, idx, datPre)

%% set frequency parameters
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);

highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = linspace(5, 7, highnumfrex); 

%% load the data

% try %try loading the processed file
%     chanDat = load([chanFiles(idx).folder '/' chanFiles(idx).name]).chanDat; 
% catch
    chanDat = load([chanFiles(idx).folder '/CHANRAW/' chanFiles(idx).name]).chanDat; % go raw if it's not working!
% end

disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi)])


%% check for encoding info
% if ~isfield(chanDat, 'encInfo')
%     dataDirPath = split(chanDat.dataDir, 'Johnson_Lab');
%     dat = load(fullfile(['/projects/p31578' dataDirPath{2} '/' chanDat.encDatFn])).data; 
%     chanDat.encInfo = dat.trialinfo; 
%     clear dat
%     save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat')
% end

%% needs time points! hard code
if chanDat.fsample == 1000
    chanDat.enctim = [-1000:3500];
    chanDat.retOtim = [-1000:3000];
    chanDat.retRtim = [-2000:500];
else
    step = 1000 / chanDat.fsample; 
    chanDat.enctim = [-1000:step:3500];
    chanDat.retOtim = [-1000:step:3000];
    chanDat.retRtim = [-2000:step:500];
end



%% High frequency Broadband 

%note hardcoded baselines: encoding: -450 : -50   ms
%                          retOn   : -450 : -50   ms
%                          retRT   : -2000: -1600 ms

% if ~isfield(chanDat, 'HFBenc')

    HFB = struct; 
    %ENCODING DATA: ***********************************************************
    chanDat.HFBenc = 0; %note if it's a reactive channel, assume not
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.enc, highfrex, chanDat.fsample, chanDat.enctim); 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(mulFrex); 
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3)); 

%     %can I make an RT aligned set? 
%     pow_align = nan(4000, size(pow,2)); 
%     RTs = chanDat.encInfo(:, 4); 
%     RT_down = 2500 - (round(RTs/5) + 200); 
% 
%     for tt = 1:size(pow,2)
%         pow_align(RT_down(tt):RT_down(tt)+size(pow,1)-1, tt) = pow(:,tt); 
%     end
% 
%     test = sum(isnan(pow_align),2);
%     RT_tim = [-12495:5:7500];
%     test = find(RT_tim<-2000 | RT_tim>2000);
%     RT_tim(test) = []; 
%     pow_align(test,:) = []; 

    
    test = mean(pow(:,chanDat.use), 2);
    test = test>1.96;
    testidx = find(test([find(mulTim>=-450,1):find(mulTim>=2500,1)]));
    if ~isempty(testidx)
        breakPoints = [1 find(diff(testidx)>1)']; 
        if length(breakPoints)==1 && length(testidx)>10
            chanDat.HFBenc = 1; 
        else
            breakPoints = [breakPoints, length(testidx)]; 
            for bb = 1:length(breakPoints)-1
                if breakPoints(bb+1) - breakPoints(bb) > 10
                    chanDat.HFBenc = 1; 
                end
            end

        end
    end
    %get mean misses: 
    HFB.subMiss = pow(:,chanDat.use & chanDat.misses); 
    %get mean hits: 
    HFB.subHit = pow(:,chanDat.use & chanDat.hits);
    %onset locked time
    HFB.encMulTim = mulTim; 
    HFB.enconFrex = mulFrex; 

    %ENCODING RT aligned DATA: ***********************************************************
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.encRT, highfrex, chanDat.fsample, chanDat.retRtim); 
    %hack for long RT trials: 
        nanIdx = find(isnan(pow(30,:,1)));
        pow(:,nanIdx,:) = mean(pow, 2, 'omitnan'); 
    %
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(mulFrex); 
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3)); 


    %get misses RT locked
    HFB.subMissRT = pow(:,chanDat.use & chanDat.misses);
    %get hits RT locked
    HFB.subHitRT = pow(:,chanDat.use & chanDat.hits);


    %RT locked time
    HFB.encRT_tim = mulTim; 
    HFB.encRTFrex = mulFrex; 
    %clean up
    clear pow
    disp('encoding done')


    %RETRIEVAL STIM ONSET: ****************************************************
    chanDat.HFBretOn = 0; 
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retOn, highfrex, chanDat.fsample, chanDat.retOtim); 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    highnumfrex = length(mulFrex); 
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize

    test = mean(pow(:,chanDat.retInfo(:,1)>0 & chanDat.retInfo(:,1)<5,:), [2,3]);
    test = abs(test)>1.96;
    testidx = find(test([find(mulTim>=-450,1): find(mulTim>=2000,1)]));
    if ~isempty(testidx)
        breakPoints = [1 find(diff(testidx)>1)']; 
        if length(breakPoints)==1 && length(testidx)>10
            chanDat.HFBretOn = 1; 
        else
            breakPoints = [breakPoints, length(testidx)]; 
            for bb = 1:length(breakPoints)-1
                if breakPoints(bb+1) - breakPoints(bb) > 10
                    chanDat.HFBretOn = 1; 
                end
            end

        end
    end
    %get mean hit: 
    HFB.hit_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==1, :), [3])); 
    %get mean CRs: 
    HFB.cr_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==3, :), [3]));
    %get mean miss: 
    HFB.miss_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==2, :), [3]));
    %get mean FA: 
    HFB.fa_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==4, :), [3]));
    %multi taper stats    
    HFB.onMulTim = mulTim; 
    HFB.onFrex = mulFrex; 
    %clean up
    clear pow
    disp('retrieval 1 done')
    
    %RETRIEVAL RESPONSE LOCKED: ***********************************************
    chanDat.HFBretRT = 0; 
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retRT, highfrex, chanDat.fsample, chanDat.retRtim); 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-2000,1), find(mulTim>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    highnumfrex = length(mulFrex); 
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
    test = mean(pow(:,chanDat.retInfo(:,1)>0 & chanDat.retInfo(:,1)<5,:), [2,3]);
    test = abs(test)>1.96;
    testidx = find(test([find(mulTim>=-1600,1): find(mulTim>=300,1)]));

    if ~isempty(testidx)
        breakPoints = [1 find(diff(testidx)>1)']; 
        if length(breakPoints)==1 && length(testidx)>10
            chanDat.HFBretRT = 1; 
        else
            breakPoints = [breakPoints, length(testidx)]; 
            for bb = 1:length(breakPoints)-1
                if breakPoints(bb+1) - breakPoints(bb) > 10
                    chanDat.HFBretRT = 1; 
                end
            end

        end
    end


    %get mean hit: 
    HFB.hit_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==1, :), [3])); 
    %get mean CRs: 
    HFB.cr_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==3, :), [3]));
    %get mean miss: 
    HFB.miss_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==2, :), [3]));
    %get mean FA: 
    HFB.fa_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==4, :), [3]));
     %multi taper stats    
    HFB.rtMulTim = mulTim; 
    HFB.rtFrex = mulFrex; 
    %clean up
    clear pow
    disp('retrieval 2 done')
% 
% 
    chanDat.HFB = HFB; 
    
    clear HFB 
    disp('attempting saving')
    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
    disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])

% else
%     disp('HFB already done')
% end


%% Lead lag analysis
% x = 5
% if sum(sum(chanDat.roiNote)) == 0 
%     roiMat = chanDat.roimni; 
% else
%     roiMat = chanDat.roiNote; 
% end
% 
% chanDat.chansRoi = sum(roiMat,2);
% roiIDX = find(chanDat.chansRoi==1); 
% chanDat.leadLagRoi = roiMat(roiIDX,:);
% if ismember(chanDat.chi, roiIDX)
% 
% %look at lead v. lag across + - 200 ms
% %chan X lead/lag X time
% leadLag = struct; 
% encHit = zeros([sum(chanDat.chansRoi==1), 401, size(chanDat.HFB.subHit,1)]);
% encMiss = zeros([sum(chanDat.chansRoi==1), 401, size(chanDat.HFB.subHit,1)]);
% for chan = 1:length(roiIDX)
%     chan
%     chanDat2 = load([chanFiles(roiIDX(chan)).folder '/' chanFiles(roiIDX(chan)).name]).chanDat; 
% 
%     % ENCODING DATA! ******************************************************
%     %HITS: 
%     for offSet = -200:200
%         HFB1 = chanDat.HFB.subHit;
%         HFB2 = chanDat2.HFB.subHit; 
%         tim = chanDat.HFB.encMulTim; 
%         if offSet<0
%             HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
%         elseif offSet>0
%             HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
%         end
%         
%         encHit(chan, offSet+201, :) = arrayfun(@(x) corr(HFB1(x,:)', HFB2(x,:)'), [1:size(HFB1,1)]);
%         
%     end
%     leadLag.encHit = encHit; 
% 
%     %MISSES: 
%     for offSet = -200:200
%         HFB1 = chanDat.HFB.subMiss;
%         HFB2 = chanDat2.HFB.subMiss; 
%         tim = chanDat.HFB.encMulTim; 
%         if offSet<0
%             HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
%         elseif offSet>0
%             HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
%         end
%         
%         encMiss(chan, offSet+201, :) = arrayfun(@(x) corr(HFB1(x,:)', HFB2(x,:)'), [1:size(HFB1,1)]);
%         
%     end
%     leadLag.encMiss = encMiss; 
% 
% 
% 
%     
% 
% 
% 
% end
% disp('saving leadlag')
% chanDat.leadLag = leadLag; 
% save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
% 
% 
% else
% 
% chanDat.leadLag = "no ROI electrodes"; 
% save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
% 
% %THIS CHANNEL IS NOT IN A ROI, so SKIP LEAD LAG! 
% 
% 
% end






%% time frequency decomposition, extract TF summaries for target trial types: 

% subsequent hit / subsequent miss (encoding data)
% hit / miss / CR / FA (retrieval locked to onset data)
% hit / miss / CR / FA (retrieval locked to response data)

% if ~isfield(chanDat, 'ISPCboot')
%     disp('working on TF')
%     %to keep size down, don't put large variables into the chanDat struct!
%     TFout = struct;
%     %ENCODING DATA: ***********************************************************
%     pow = abs(getChanTrialTF(chanDat.enc, frex, numfrex, stds, chanDat.fsample)).^2; %get power time series for all trials/frequencies
%     pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
%     pow = cell2mat(pow); %organize
%     pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
%     %get mean misses: 
%     TFout.subMiss = squeeze(mean(pow(:,chanDat.use & chanDat.misses, :), 2)); 
%     %get mean hits: 
%     TFout.subHit = squeeze(mean(pow(:,chanDat.use & chanDat.hits, :), 2));
%     %clean up
%     clear pow
%     disp('encoding done')
% 
%     %RETRIEVAL STIM ONSET: ****************************************************
%     pow = abs(getChanTrialTF(chanDat.retOn, frex, numfrex, stds, chanDat.fsample)).^2; %get power time series for all trials/frequencies
%     pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
%     pow = cell2mat(pow); %organize
%     pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
%     %get mean hit: 
%     TFout.hit_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==1, :), 2)); 
%     %get mean CRs: 
%     TFout.cr_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==3, :), 2));
%     %get mean miss: 
%     TFout.miss_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==2, :), 2));
%     %get mean FA: 
%     TFout.fa_on = squeeze(mean(pow(:,chanDat.retInfo(:,1)==4, :), 2));
%     %clean up
%     clear pow
%     disp('retrieval 1 done')
%     
%     %RETRIEVAL RESPONSE LOCKED: ***********************************************
%     pow = abs(getChanTrialTF(chanDat.retRT, frex, numfrex, stds, chanDat.fsample)).^2; %get power time series for all trials/frequencies
%     pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
%     pow = cell2mat(pow); %organize
%     pow = reshape(pow, size(pow,1), size(pow,2)/100, []); %organize
%     %get mean hit: 
%     TFout.hit_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==1, :), 2)); 
%     %get mean CRs: 
%     TFout.cr_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==3, :), 2));
%     %get mean miss: 
%     TFout.miss_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==2, :), 2));
%     %get mean FA: 
%     TFout.fa_rt = squeeze(mean(pow(:,chanDat.retInfo(:,1)==4, :), 2));
%     %clean up
%     clear pow
%     disp('retrieval 2 done')
%     
%     chanDat.TFout = TFout; 
%     
%     clear TFout 
%     disp('attempting saving')
%     save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
%     disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
% else
%     disp('TF already done, skipping')
% end

%% get ISPC and PPC values 

% if ~isfield(chanDat, 'ISPCboot')
%     disp('working on ISPC')
%     ISPCout = struct; 
%     %store the downsample index (di) 
%     ISPCout.encdi = 1:20:length(chanDat.enctim);
%     ISPCout.ondi = 1:20:length(chanDat.retOtim); 
%     ISPCout.rtdi = 1:20:length(chanDat.retRtim); 
%     
%     %preallocate: 
%     %channels X time X frequencies X ISPC/PPC
%     ISPCout.subMiss = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 4); 
%     ISPCout.subHit = zeros(length(chanFiles), length(ISPCout.encdi), length(frex), 4); 
%     
%     ISPCout.hit_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
%     ISPCout.cr_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
%     ISPCout.miss_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
%     ISPCout.fa_on = zeros(length(chanFiles), length(ISPCout.ondi), length(frex), 4);
% 
% 
%     ISPCout.hit_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
%     ISPCout.cr_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
%     ISPCout.miss_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
%     ISPCout.fa_rt = zeros(length(chanFiles), length(ISPCout.rtdi), length(frex), 4);
% 
% 
%     %will need to loop channels
%     %NOTE: all trial types must have at least two trials! 
%     for chan = 1:length(chanFiles)
%         tic
% %         chan
%         if chan > idx %don't do repeat work! 
%         chanDat2 = load([chanFiles(chan).folder '/CHANRAW/' chanFiles(chan).name]).chanDat; 
% 
%         %ENCODING DATA: ***********************************************************
%         if sum(chanDat.use & chanDat.misses)>1
%         ISPCout.subMiss(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.misses);
%         end
%         if sum(chanDat.use & chanDat.hits)>1
%         ISPCout.subHit(chan,:,:,:) = getChanISPC(chanDat.enc, chanDat2.enc, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.encdi, chanDat.use & chanDat.hits);
%         end
% 
%         %RETRIEVAL STIM ONSET: ****************************************************
%         if sum(chanDat.retInfo(:,1)==1) > 1
%         ISPCout.hit_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==1);
%         end
%         if sum(chanDat.retInfo(:,1)==3) > 1
%         ISPCout.cr_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==3);
%         end
%         if sum(chanDat.retInfo(:,1)==2) > 1
%         ISPCout.miss_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==2);
%         end
%         if sum(chanDat.retInfo(:,1)==4) > 1
%         ISPCout.fa_on(chan,:,:,:) = getChanISPC(chanDat.retOn, chanDat2.retOn, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.ondi, chanDat.retInfo(:,1)==4);
%         end
%         
%         %RETRIEVAL RESPONSE LOCKED: ****************************************************
%         if sum(chanDat.retInfo(:,1)==1) > 1
%         ISPCout.hit_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==1);
%         end
%         if sum(chanDat.retInfo(:,1)==3) > 1
%         ISPCout.cr_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==3);
%         end
%         if sum(chanDat.retInfo(:,1)==2) > 1
%         ISPCout.miss_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==2);
%         end
%         if sum(chanDat.retInfo(:,1)==4) > 1
%         ISPCout.fa_rt(chan,:,:,:) = getChanISPC(chanDat.retRT, chanDat2.retRT, ...
%             frex, numfrex, stds, chanDat.fsample, ISPCout.rtdi, chanDat.retInfo(:,1)==4);
%         end
%         
% 
%         end
%         disp(['channel: ' num2str(chan) ' took ' num2str(round(toc/60,1)) ' minutes'])
%     end
%     chanDat.ISPCout = ISPCout; 
%     chanDat.ISPCboot = true; 
% %     chanDat = rmfield(chanDat, 'sizeReduce'); 
%     disp('attempting saving')
%     save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 
%     disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
% 
% 
% 
% else
%     disp('connectivity already done, skipping')
% end


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























end