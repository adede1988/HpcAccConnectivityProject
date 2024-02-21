
function HFB = getHFB(chanDat, highfrex)

    chanDat.enctim = [-1000:3500];
    chanDat.enctimRT = [-2000:500];
    chanDat.retOtim = [-1000:3000];
    chanDat.retRtim = [-2000:500];
    HFB = struct; 
    %ENCODING DATA: ***********************************************************
    chanDat.HFBenc = 0; %note if it's a reactive channel, assume not
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.enc, highfrex, chanDat.fsample, chanDat.enctim, 5); 
   
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(mulFrex); 
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3, 'omitnan')); 
    
   
    %get mean misses: 
    HFB.subMiss = pow(:,chanDat.use & chanDat.misses); 
    %get mean hits: 
    HFB.subHit = pow(:,chanDat.use & chanDat.hits);
    %onset locked time
    HFB.encMulTim = mulTim; 
    HFB.enconFrex = mulFrex; 

    %ENCODING RT aligned DATA: ***********************************************************
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.encRT, highfrex, chanDat.fsample, chanDat.retRtim, 5); 
    %hack for long RT trials: 
        nanIdx = find(isnan(pow(100,:,1)));
        if ~isempty(nanIdx)
            for nani = 1:length(nanIdx)
                pow(:,nanIdx(nani),:) = mean(pow, 2, 'omitnan'); 
            end
        end
    %
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-2000,1), find(mulTim>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(mulFrex); 
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3, 'omitnan')); 


    %get misses RT locked
    HFB.subMissRT = pow(:,chanDat.use & chanDat.misses);
    %get hits RT locked
    HFB.subHitRT = pow(:,chanDat.use & chanDat.hits);


    %RT locked time
    HFB.encRT_tim = mulTim; 
    HFB.encRTFrex = mulFrex; 
    %clean up
    clear pow
%     disp('encoding done')


    %RETRIEVAL STIM ONSET: ****************************************************
    chanDat.HFBretOn = 0; 
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retOn, highfrex, chanDat.fsample, chanDat.retOtim, 5); 
 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    highnumfrex = length(mulFrex); 
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
    pow = squeeze(mean(pow, 3, 'omitnan')); 
    
    %get mean hit: 
    HFB.hit_on = pow(:,chanDat.retInfo(:,1)==1, :); 
    %get mean CRs: 
    HFB.cr_on = pow(:,chanDat.retInfo(:,1)==3, :);
    %get mean miss: 
    HFB.miss_on = pow(:,chanDat.retInfo(:,1)==2, :);
    %get mean FA: 
    HFB.fa_on = pow(:,chanDat.retInfo(:,1)==4, :);
    %multi taper stats    
    HFB.onMulTim = mulTim; 
    HFB.onFrex = mulFrex; 
    %clean up
    clear pow
%     disp('retrieval 1 done')
    
    %RETRIEVAL RESPONSE LOCKED: ***********************************************
    chanDat.HFBretRT = 0; 
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retRT, highfrex, chanDat.fsample, chanDat.retRtim, 5); 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-2000,1), find(mulTim>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    highnumfrex = length(mulFrex); 
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
    %take the mean across trials and frequencies to get a test timeseries
    pow = squeeze(mean(pow, 3, 'omitnan')); 


    %get mean hit: 
    HFB.hit_rt = pow(:,chanDat.retInfo(:,1)==1, :); 
    %get mean CRs: 
    HFB.cr_rt = pow(:,chanDat.retInfo(:,1)==3, :);
    %get mean miss: 
    HFB.miss_rt = pow(:,chanDat.retInfo(:,1)==2, :);
    %get mean FA: 
    HFB.fa_rt = pow(:,chanDat.retInfo(:,1)==4, :);
     %multi taper stats    
    HFB.rtMulTim = mulTim; 
    HFB.rtFrex = mulFrex; 
    %clean up
    clear pow
%     disp('retrieval 2 done')


    end