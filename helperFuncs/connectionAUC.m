function [aovDat, ai] = connectionAUC(aovDat, ai, curDat, ...
                regions, phase, stat, datPre)

reg1 = find(cellfun(@(x) strcmp(x, curDat.reg1), regions)); 
reg2 = find(cellfun(@(x) strcmp(x, curDat.reg2), regions)); 

if isfield(curDat.([phase '_' stat '_clust']), 'pos_clusters')
for cc = 1:length(curDat.([phase '_' stat '_clust']).pos_clusters)
if curDat.([phase '_' stat '_clust']).pos_clusters(cc).p < .05 %sig check
    for pp = 1:length(curDat.subVals) 
        
        %loop pairs to calculate AUC for each channel pair
        if curDat.chiVals(pp,1) < 10
            chi = ['00' num2str(curDat.chiVals(pp,1))]; 
        elseif curDat.chiVals(pp,1) < 100
            chi = ['0' num2str(curDat.chiVals(pp,1))]; 
        else
            chi = num2str(curDat.chiVals(pp,1)); 
        end
        ch1 = load([datPre 'CHANDAT/finished/chanDat_' ...
            curDat.subVals{pp} '_' chi '.mat']).chanDat; 
        if curDat.chiVals(pp,2) < 10
            chi = ['00' num2str(curDat.chiVals(pp,2))]; 
        elseif curDat.chiVals(pp,2) < 100
            chi = ['0' num2str(curDat.chiVals(pp,2))]; 
        else
            chi = num2str(curDat.chiVals(pp,2)); 
        end
        ch2 = load([datPre 'CHANDAT/finished/chanDat_' ...
            curDat.subVals{pp} '_' chi '.mat']).chanDat;

        %some cases are flipped because in
        %connectionIntegrate.m line 25 I had to flip regions
        %sometimes for the image locked data since timing
        %doesn't matter for these, so I avoided repeat
        %calculation by flipping for half of the pairs, but I
        %then wrote down the flipped chi values for these and
        %it does matter here where I'm reconstructing file
        %names for HFB locked connectivity
        if ~strcmp(ch1.labels{curDat.chiVals(pp,1),3}, regions{reg1})
            tmp = ch1; 
            ch1 = ch2; 
            ch2 = tmp; 
            tmp = curDat.chiVals(pp,1); 
            curDat.chiVals(pp,1) = curDat.chiVals(pp,2); 
            curDat.chiVals(pp,2) = tmp; 
        end

        %what's the mean frequency of this connection? 
        confrex = logspace(log10(2), log10(25), 20); 
        weight = sum(curDat.([phase '_' stat '_clust']).pos_clusters(cc).inds);
        meanFrex = sum((confrex .* weight)) / sum(weight);
        [~, freqi] = min(abs(confrex - meanFrex));
       
        
        %what's the mean time index for this connection? 
        weight = sum(curDat.([phase '_' stat '_clust']).pos_clusters(cc).inds'); 
        ti = round(sum([1:length(weight)] .* weight) / sum(weight)); 

        stds = logspace(log10(3),log10(5),20)./(2*pi*confrex);

        %set up the inputs for getting angle difs and PPC
        if strcmp(stat, 'HFB')
            HFBflag = true; 
        else
            HFBflag = false; 
        end

        if strcmp(phase, 'enc')
            dat1 = ch1.enc; 
            dat2 = ch2.enc; 
            di = ch1.ISPC.encdi; 
            hitTrials = ch1.use & ch1.hits;
            missTrials = ch1.use & ch1.misses;
            HFBLatHit = ch1.HFB_lat.subHit;
            HFBLatMiss = ch1.HFB_lat.subMiss;
            outTim = ch1.HFB.encMulTim;
        else
            dat1 = ch1.retOn; 
            dat2 = ch2.retOn; 
            di = ch1.ISPC.ondi; 
            hitTrials = ch1.retInfo(:,1)==1;
            missTrials = ch1.retInfo(:,1)==2;
            HFBLatHit = ch1.HFB_lat.retHit;
            HFBLatMiss = ch1.HFB_lat.retMiss;
            outTim = ch1.HFB.onMulTim;

        end
        
        %get angle difs and PPC values 
        [hitDifs, ppcHit] = getClustAngDat(dat1, dat2, confrex, ...
            20, stds,...
            1000, di, hitTrials, ...
            HFBLatHit, ...
            outTim, HFBflag, freqi, ti);

        [missDifs, ppcMiss] = getClustAngDat(dat1, dat2, confrex, ...
            20, stds,...
            1000, di, missTrials, ...
            HFBLatMiss, ...
            outTim, HFBflag, freqi, ti);

        %get the mean hit angle dif, then convert all angle
        %difs into distance from the mean in order to make
        %non-circular for ROC analysis
        %the assumption here is that the mean angle difference for the hits
        %is the good angle difference for memory performance and deviations
        %from this mean difference should predict worse memory
        mean_cos = mean(cos(hitDifs));
        mean_sin = mean(sin(hitDifs));
        mean_angle = atan2(mean_sin, mean_cos);
    
        hitDists = abs(angdiff(hitDifs, repmat(mean_angle, size(hitDifs))));
        missDists = abs(angdiff(missDifs, repmat(mean_angle, size(missDifs))));
        %criteria go from 0 to pi because on a circle, you can never be
        %more than pi away from anywhere (if both directions are
        %equivalent). A difference of 0 indicates that the observed value
        %for a particular trial is exactly on the mean difference, which
        %should predict optimal memory
        crit = [0:pi/20:pi]; 

        HR = arrayfun(@(x) sum(hitDists<x) / length(hitDists), crit);
        FA = arrayfun(@(x) sum(missDists<x) / length(missDists), crit);
        
        AUC = 0; 
        for ci = 1:length(crit)
            if(ci>1)
              baseDist = FA(ci) - FA(ci-1);
              %top triangle area:
              tria = ((HR(ci) - HR(ci-1)) * baseDist) / 2;
              %base rectangle area:
              reca = HR(ci-1) * baseDist;
              AUC = AUC + tria + reca;
             else 
              baseDist = FA(ci);
              %top triangle area:
              tria = ((HR(ci)) * baseDist) / 2;
              AUC = AUC + tria;
            end
        end
        tim = outTim; 
        aovDat.subID(ai) = curDat.subVals{pp}; 
        aovDat.chi(ai) = curDat.chiVals(pp,1); 
        aovDat.chi2(ai) = curDat.chiVals(pp,2); %only applicable to connectivity
        aovDat.HFB_Image(ai) = stat; %HFB / Image timing
        aovDat.statType(ai) = 'Con'; %HFB, TF power, ITPC, connectivity
        aovDat.encRet(ai) = phase; %enc / ret
        aovDat.reg1(ai) = regions{reg1}; 
        aovDat.reg2(ai) = regions{reg2}; %only applicable for connectivity 
        aovDat.AUC(ai) = AUC; %area under the curve to differentiate hit/miss
        aovDat.hitVal(ai) = ppcHit; %raw stat for hits
        aovDat.missVal(ai) = ppcMiss; %raw stat for misses
        aovDat.tVal(ai) = curDat.([phase '_' stat '_tVal'])(ti,freqi);
        aovDat.time(ai) = tim(ti); % time in ms 
        aovDat.freq(ai) = meanFrex; % mean frequency (HFB uses 120)
        ai = ai+1; 
    end
end %is the current cluster p<.05? 
end %loop on the clusters
end %is there a positive clust field? 



end