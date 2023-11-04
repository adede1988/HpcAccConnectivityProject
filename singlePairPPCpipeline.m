
function [] = singlePairPPCpipeline(chanFiles, savFolder)

tic
chanDat2 = load([chanFiles.folder '/' chanFiles.name2]).chanDat; 
chanDat = load([chanFiles.folder '/' chanFiles.name]).chanDat; 
pairDat = chanFiles; 
pairDat.chi1 = chanDat.chi;
clear chanDat
pairDat.chi2 = chanDat2.chi;
fn = [chanDat2.subID '_' num2str(pairDat.chi1) '_' num2str(pairDat.chi2) '_'...
    pairDat.chan1Reg{1} '_' pairDat.chan2Reg{1} '.mat'];

%% set frequency parameters
frex = linspace(2,10,10);
numfrex = length(frex); 
stds = linspace(2,3,numfrex);

highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = linspace(5, 7, highnumfrex); 

%% load the data, but try seeing if this file has already been worked on

if ~exist([savFolder fn], 'file')

pairDat.frex = frex; 
pairDat.stds = stds; 
pairDat.chan1Reg = pairDat.chan1Reg{1}; 
pairDat.chan2Reg = pairDat.chan2Reg{1};

%get the HFB latencies and reactivity from channel 2

HFB2 = getHFB(chanDat2, highfrex); %get HFB
pairDat.react2 = sum(reactiveTest_100(HFB2)==1)>0; %reactivity
HFB2 = HFBdown(HFB2); %down sample
%latency calculations for encoding
tim = HFB2.encMulTim;
RT = chanDat2.encInfo(chanDat2.use & chanDat2.misses, 4); 
pairDat.subMissLat2 = gausLat(HFB2.subMiss, tim, RT);
RT = chanDat2. encInfo(chanDat2.use & chanDat2.hits, 4); 
pairDat.subHitLat2 = gausLat(HFB2.subHit, tim, RT);
%latency calculations for retrieval
tim = HFB2.onMulTim;
RT = chanDat2.retInfo(chanDat2.retInfo(:,1)==1, 3); 
pairDat.retHitLat2 = gausLat(HFB2.hit_on, tim, RT);
RT = chanDat2.retInfo(chanDat2.retInfo(:,1)==2, 3); 
pairDat.retMissLat2 = gausLat(HFB2.miss_on, tim, RT);

enc2 = chanDat2.enc;
ret2 = chanDat2.retOn; 

clear chanDat2 HFB2 RT tim

chanDat = load([chanFiles.folder '/' chanFiles.name]).chanDat; 


%% get time info set
chanDat.enctim = [-1000:3500];
chanDat.enctimRT = [-2000:500];
chanDat.retOtim = [-1000:3000];
chanDat.retRtim = [-2000:500];

%% get the HFB latency and reaction info for channel 1

HFB2 = getHFB(chanDat, highfrex); %get HFB
pairDat.react1 = sum(reactiveTest_100(HFB2)==1)>0; %reactivity
HFB2 = HFBdown(HFB2); %down sample
%latency calculations for encoding
tim = HFB2.encMulTim;
RT = chanDat.encInfo(chanDat.use & chanDat.misses, 4); 
pairDat.subMissLat1 = gausLat(HFB2.subMiss, tim, RT);
RT = chanDat. encInfo(chanDat.use & chanDat.hits, 4); 
pairDat.subHitLat1 = gausLat(HFB2.subHit, tim, RT);
%latency calculations for retrieval
tim = HFB2.onMulTim;
RT = chanDat.retInfo(chanDat.retInfo(:,1)==1, 3); 
pairDat.retHitLat1 = gausLat(HFB2.hit_on, tim, RT);
RT = chanDat.retInfo(chanDat.retInfo(:,1)==2, 3); 
pairDat.retMissLat1 = gausLat(HFB2.miss_on, tim, RT);

 

%get the downsamples for PPC calculation 
pairDat.encdi = arrayfun(@(x) find(x<=chanDat.enctim,1), HFB2.encMulTim);
pairDat.ondi = arrayfun(@(x) find(x<=chanDat.retOtim,1), HFB2.onMulTim);
pairDat.encTim = HFB2.encMulTim; 
pairDat.retTim = HFB2.onMulTim; 

clear HFB2 RT tim 

pairDat.encInfo = chanDat.encInfo; 
pairDat.retInfo = chanDat.retInfo; 
pairDat.misses = chanDat.misses; 
pairDat.use = chanDat.use; 
pairDat.hits = chanDat.hits; 


disp(['data loaded: ' chanDat.subID ' ' num2str(chanDat.chi) ' ' num2str(pairDat.chi2)])



 

%% PPC over time at single trial level


disp('working on PPC_time')


%NOTE: all trial types must have at least two trials! 

if ~isfield(pairDat, 'subMissPPC')
pairDat.toc1 = toc; 
%ENCODING DATA: ***********************************************************
if sum(chanDat.use & chanDat.misses)>1 
pairDat.subMissPPC = getChanPPC_time(chanDat.enc, enc2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.encdi, chanDat.use & chanDat.misses);
end
pairDat.toc2 = toc; 

save([savFolder '/' fn], 'pairDat'); 
end

if ~isfield(pairDat, 'subHitPPC')
if sum(chanDat.use & chanDat.hits)>1
pairDat.subHitPPC= getChanPPC_time(chanDat.enc, enc2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.encdi, chanDat.use & chanDat.hits);
end
pairDat.toc2 = toc; 
save([savFolder '/' fn], 'pairDat'); 
end
% TF
pow = log10(abs(getChanTrialTF(chanDat.enc, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctim>=-450,1), find(chanDat.enctim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
pow = cell2mat(pow); %organize
pow = reshape(pow, size(pow,1), size(pow,2)/10, []); %organize
%get mean misses: 
pairDat.subMissTF1 = pow(pairDat.encdi,chanDat.use & chanDat.misses, :); 
%get mean hits: 
pairDat.subHitTF1 = pow(pairDat.encdi,chanDat.use & chanDat.hits, :);

pow = log10(abs(getChanTrialTF(enc2, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctim>=-450,1), find(chanDat.enctim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
pow = cell2mat(pow); %organize
pow = reshape(pow, size(pow,1), size(pow,2)/10, []); %organize
%get mean misses: 
pairDat.subMissTF2 = pow(pairDat.encdi,chanDat.use & chanDat.misses, :); 
%get mean hits: 
pairDat.subHitTF2 = pow(pairDat.encdi,chanDat.use & chanDat.hits, :);

pairDat.toc3 = toc; 


%RETRIEVAL STIM ONSET: ****************************************************
if ~isfield(pairDat, 'hit_onPPC')
if sum(chanDat.retInfo(:,1)==1) > 1
pairDat.hit_onPPC = getChanPPC_time(chanDat.retOn, ret2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.ondi, chanDat.retInfo(:,1)==1);
end
pairDat.toc4 = toc; 
save([savFolder '/' fn], 'pairDat'); 
end

if ~isfield(pairDat, 'miss_onPPC')
if sum(chanDat.retInfo(:,1)==2) > 1
pairDat.miss_onPPC = getChanPPC_time(chanDat.retOn, ret2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.ondi, chanDat.retInfo(:,1)==2);
end
pairDat.toc5 = toc; 
save([savFolder '/' fn], 'pairDat'); 
end
% TF
pow = log10(abs(getChanTrialTF(chanDat.retOn, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retOtim>=-450,1), find(chanDat.retOtim>=-50,1)]), 1:size(pow,3), 'UniformOutput',false ); %z-score
pow = cell2mat(pow); %organize
pow = reshape(pow, size(pow,1), size(pow,2)/10, []); %organize
%get mean hit: 
pairDat.hit_onTF1 = pow(pairDat.ondi,chanDat.retInfo(:,1)==1, :); 

%get mean miss: 
pairDat.miss_onTF1 = pow(pairDat.ondi,chanDat.retInfo(:,1)==2 , :);


pow = log10(abs(getChanTrialTF(ret2, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retOtim>=-450,1), find(chanDat.retOtim>=-50,1)]), 1:size(pow,3), 'UniformOutput',false ); %z-score
pow = cell2mat(pow); %organize
pow = reshape(pow, size(pow,1), size(pow,2)/10, []); %organize
%get mean hit: 
pairDat.hit_onTF2 = pow(pairDat.ondi,chanDat.retInfo(:,1)==1, :); 

%get mean miss: 
pairDat.miss_onTF2 = pow(pairDat.ondi,chanDat.retInfo(:,1)==2 , :);

pairDat.toc6 = toc; 


else %% this file has already been created and run to here, just adding on the leadlag analysis
    tic
    chanDat2 = load([chanFiles.folder '/' chanFiles.name2]).chanDat; 
    chanDat = load([chanFiles.folder '/' chanFiles.name]).chanDat; 
    %% get time info set
chanDat.enctim = [-1000:3500];
chanDat.enctimRT = [-2000:500];
chanDat.retOtim = [-1000:3000];
chanDat.retRtim = [-2000:500];
    pairDat = load([savFolder fn]).pairDat; 
    disp('doing encoding lead lag')
    % encoding lead lag
    pow = log10(abs(getChanTrialTF(chanDat.enc, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctim>=-450,1), ...
        find(chanDat.enctim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/10, []); %organize

    missidx = find(chanDat.use & chanDat.misses); 
    hitidx = find(chanDat.use & chanDat.hits); 

    pow2 = log10(abs(getChanTrialTF(chanDat2.enc, frex, numfrex, stds, chanDat2.fsample)).^2); %get power time series for all trials/frequencies
    pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), [find(chanDat.enctim>=-450,1),...
        find(chanDat.enctim>=-50,1)] ), 1:size(pow2,3), 'UniformOutput',false ); %z-score
    pow2 = cell2mat(pow2); %organize
    pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/10, []); %organize

    leadLagEncTim = chanDat.enctim(501:25:end-500);
    leadLagRetTim = chanDat.retOtim(501:25:end-500);
    LLsubMiss = zeros(length(frex), length(missidx), 301, length(leadLagEncTim)); 
    LLsubHit = zeros(length(frex), length(hitidx), 301, length(leadLagEncTim)); 
    for fi = 1:length(frex)
        disp(num2str(fi))
        [LLsubHit(fi,:,:,:), LLsubMiss(fi,:,:,:)] = getLL_trials(squeeze(pow(:,:,1)), squeeze(pow2(:,:,2)), missidx, hitidx, ...
            chanDat.enctim, leadLagEncTim); 
    end
    pairDat.LLsubMiss = LLsubMiss; 
    pairDat.LLsubHit = LLsubHit; 

    % retrieval lead lag
    disp('doing retrieval leadlag')
    pow = log10(abs(getChanTrialTF(chanDat.retOn, frex, numfrex, stds, chanDat.fsample)).^2); %get power time series for all trials/frequencies
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retOtim>=-450,1), ...
        find(chanDat.retOtim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/10, []); %organize

    miss_on = find(chanDat.retInfo(:,1)==2);
    hit_on = find(chanDat.retInfo(:,1)==1);  

    pow2 = log10(abs(getChanTrialTF(chanDat2.retOn, frex, numfrex, stds, chanDat2.fsample)).^2); %get power time series for all trials/frequencies
    pow2 = arrayfun(@(x) myChanZscore(pow2(:,:,x), [find(chanDat.retOtim>=-450,1),...
        find(chanDat.retOtim>=-50,1)] ), 1:size(pow2,3), 'UniformOutput',false ); %z-score
    pow2 = cell2mat(pow2); %organize
    pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/10, []); %organize

    LLretMiss = zeros(length(frex), length(miss_on), 301, length(leadLagRetTim)); 
    LLretHit = zeros(length(frex), length(hit_on), 301, length(leadLagRetTim)); 
    for fi = 1:length(frex)
        disp(num2str(fi))
        [LLretHit(fi,:,:,:), LLretMiss(fi,:,:,:)] = getLL_trials(squeeze(pow(:,:,fi)), ...
            squeeze(pow2(:,:,fi)), miss_on, hit_on, ...
            chanDat.retOtim, leadLagRetTim); 
    end
    pairDat.LLretMiss = LLretMiss; 
    pairDat.LLretHit = LLretHit; 



    save([savFolder '/' fn], 'pairDat'); 


% 
%     tmp = pairDat.retHitLat1;
%     tmp2 = LLretHit(:,tmp>-1,:,:); 
%     tmp(tmp==-1) = []; 
%     tmp(tmp>=2000) = 1975;
%     test = arrayfun(@(x) squeeze(tmp2(:, x,:, ...
%         find(leadLagRetTim>=tmp(x),1)-20 : ...
%         find(leadLagRetTim>=tmp(x),1)+20)), ...
%         1:length(tmp), 'uniformoutput', false);
%     test2 = reshape(cell2mat(test), [10, 301, length(test), 41]); 




%     test = arrayfun(@(x) squeeze(hitTemp(x, :, ...
%         find(pairDat.encTim>=pairDat.subHitLat1(x),1)-20 : ...
%         find(pairDat.encTim>=pairDat.subHitLat1(x),1)+20)), ...
%         1:36, 'uniformoutput', false);
%     test2 = reshape(cell2mat(test), [301, 41, 36]); 


    
    
toc



end





%     chanDat = rmfield(chanDat, 'sizeReduce'); 
disp('attempting saving')
save([savFolder '/' fn], 'pairDat'); 
disp(['save success: ' savFolder '/' chanFiles.name])





end
