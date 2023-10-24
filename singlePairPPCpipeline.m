
function [] = singlePairPPCpipeline(chanFiles, savFolder)

%% set frequency parameters
frex = linspace(2,10,10);
numfrex = length(frex); 
stds = linspace(2,3,numfrex);

highfrex = linspace(70, 150, 81); 
highnumfrex = length(highfrex); 
highstds = linspace(5, 7, highnumfrex); 

%% load the data
pairDat = chanFiles; 
pairDat.frex = frex; 
pairDat.stds = stds; 
pairDat.chan1Reg = pairDat.chan1Reg{1}; 
pairDat.chan2Reg = pairDat.chan2Reg{1};

%get the HFB latencies and reactivity from channel 2
chanDat2 = load([chanFiles.folder '/' chanFiles.name2]).chanDat; 
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
pairDat.chi2 = chanDat2.chi;

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

pairDat.chi1 = chanDat.chi; 

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

fn = [chanDat.subID '_' num2str(pairDat.chi1) '_' num2str(pairDat.chi2) '_' pairDat.chan1Reg '_' pairDat.chan2Reg '.mat'];


 

%% PPC over time at single trial level


disp('working on PPC_time')


%NOTE: all trial types must have at least two trials! 



%ENCODING DATA: ***********************************************************
if sum(chanDat.use & chanDat.misses)>1 
pairDat.subMissPPC = getChanPPC_time(chanDat.enc, enc2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.encdi, chanDat.use & chanDat.misses);
end
save([savFolder '/' fn], 'pairDat'); 

if sum(chanDat.use & chanDat.hits)>1
pairDat.subHitPPC= getChanPPC_time(chanDat.enc, enc2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.encdi, chanDat.use & chanDat.hits);
end
save([savFolder '/' fn], 'pairDat'); 

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



%RETRIEVAL STIM ONSET: ****************************************************
if sum(chanDat.retInfo(:,1)==1) > 1
pairDat.hit_onPPC = getChanPPC_time(chanDat.retOn, ret2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.ondi, chanDat.retInfo(:,1)==1);
end
save([savFolder '/' fn], 'pairDat'); 

if sum(chanDat.retInfo(:,1)==2) > 1
pairDat.miss_onPPC = getChanPPC_time(chanDat.retOn, ret2, ...
    frex, numfrex, stds, chanDat.fsample, pairDat.ondi, chanDat.retInfo(:,1)==2);
end
save([savFolder '/' fn], 'pairDat'); 

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





%     chanDat = rmfield(chanDat, 'sizeReduce'); 
disp('attempting saving')
save([savFolder '/' fn], 'pairDat'); 
disp(['save success: ' savFolder '/' chanFiles.name])





end
