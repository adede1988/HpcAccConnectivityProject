
function HFB = getHFB(chanDat, highfrex)

    chanDat.enctim = [-1000:3500];
    chanDat.enctimRT = [-2000:500];
    chanDat.retOtim = [-1000:3000];
    chanDat.retRtim = [-2000:500];
    HFB = struct; 
    %ENCODING DATA: ***********************************************************
    chanDat.HFBenc = 0; %note if it's a reactive channel, assume not
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.enc, highfrex, chanDat.fsample, chanDat.enctim, 25); 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-450,1), find(mulTim>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    highnumfrex = length(mulFrex); 
    pow = cell2mat(pow); %organize
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize

    %take the mean across frequencies
    pow = squeeze(mean(pow, 3)); 
    
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
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.encRT, highfrex, chanDat.fsample, chanDat.retRtim, 25); 
    %hack for long RT trials: 
        nanIdx = find(isnan(pow(30,:,1)));
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
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retOn, highfrex, chanDat.fsample, chanDat.retOtim, 25); 
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
    [pow, mulTim, mulFrex] = getChanMultiTF(chanDat.retRT, highfrex, chanDat.fsample, chanDat.retRtim, 25); 
    pow = arrayfun(@(x) myChanZscore(pow(:,:,x), [find(mulTim>=-2000,1), find(mulTim>=-1600,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
    pow = cell2mat(pow); %organize
    highnumfrex = length(mulFrex); 
    pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
    %take the mean across trials and frequencies to get a test timeseries
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


    end