
function HFB = getHFB(chanDat, highfrex, stdVals)

   
    HFB = struct; 
    %ENCODING DATA: ***********************************************************
   

    pow = getChanTrialTF(chanDat.enc, highfrex, length(highfrex), ...
        stdVals, chanDat.fsample); %get power time series for all trials/frequencies
   
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctim>=-450,1), find(chanDat.enctim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/length(highfrex), []); %organize


    %take the mean across frequencies
    pow = squeeze(mean(pow, 3, 'omitnan')); 
    %take every 5th time point
    pow = pow(1:5:end, :); 
    curTim = chanDat.enctim(1:5:end); 
   
    %get misses: 
    HFB.subMiss = pow(:,chanDat.use & chanDat.misses); 
    %get hits: 
    HFB.subHit = pow(:,chanDat.use & chanDat.hits);
    %onset locked time
    HFB.encMulTim = curTim; 
    HFB.enconFrex = highfrex; 

    
    %ENCODING RT aligned DATA: ***********************************************************
    pow = getChanTrialTF(chanDat.encRT, highfrex, length(highfrex), ...
        stdVals, chanDat.fsample); %get power time series for all trials/frequencies
   
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.enctimRT>=-2000,1),...
                                                  find(chanDat.enctimRT>=-1600,1)] ), ...
                                                  1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/length(highfrex), []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3, 'omitnan')); 
    %take every 5th time point
    pow = pow(1:5:end, :); 
    curTim = chanDat.enctimRT(1:5:end); 

    %get misses RT locked
    HFB.subMissRT = pow(:,chanDat.use & chanDat.misses);
    %get hits RT locked
    HFB.subHitRT = pow(:,chanDat.use & chanDat.hits);


    %RT locked time
    HFB.encRT_tim = curTim; 
    HFB.encRTFrex = highfrex; 
    %clean up
    clear pow
%     disp('encoding done')


    %RETRIEVAL STIM ONSET: ****************************************************
    
    pow = getChanTrialTF(chanDat.retOn, highfrex, length(highfrex), ...
        stdVals, chanDat.fsample); %get power time series for all trials/frequencies
   
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retOtim>=-450,1),...
                                                  find(chanDat.retOtim>=-50,1)] ), ...
                                                  1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/length(highfrex), []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3, 'omitnan')); 
    %take every 5th time point
    pow = pow(1:5:end, :); 
    curTim = chanDat.retOtim(1:5:end); 
    
    %get mean hit: 
    HFB.hit_on = pow(:,chanDat.retInfo(:,1)==1, :); 
    %get mean CRs: 
    HFB.cr_on = pow(:,chanDat.retInfo(:,1)==3, :);
    %get mean miss: 
    HFB.miss_on = pow(:,chanDat.retInfo(:,1)==2, :);
    %get mean FA: 
    HFB.fa_on = pow(:,chanDat.retInfo(:,1)==4, :);
    %multi taper stats    
    HFB.onMulTim = curTim; 
    HFB.onFrex = highfrex; 
    %clean up
    clear pow
    
    %RETRIEVAL RESPONSE LOCKED: ***********************************************
  
    pow = getChanTrialTF(chanDat.retRT, highfrex, length(highfrex), ...
        stdVals, chanDat.fsample); %get power time series for all trials/frequencies
   
    pow = abs(pow).^2; 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(chanDat.retRtim>=-2000,1),...
                                                  find(chanDat.retRtim>=-1600,1)] ), ...
                                                  1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/length(highfrex), []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3, 'omitnan')); 
    %take every 5th time point
    pow = pow(1:5:end, :); 
    curTim = chanDat.retRtim(1:5:end); 


    %get mean hit: 
    HFB.hit_rt = pow(:,chanDat.retInfo(:,1)==1, :); 
    %get mean CRs: 
    HFB.cr_rt = pow(:,chanDat.retInfo(:,1)==3, :);
    %get mean miss: 
    HFB.miss_rt = pow(:,chanDat.retInfo(:,1)==2, :);
    %get mean FA: 
    HFB.fa_rt = pow(:,chanDat.retInfo(:,1)==4, :);
     %multi taper stats    
    HFB.rtMulTim = curTim; 
    HFB.rtFrex = highfrex; 
    %clean up
    clear pow


    end