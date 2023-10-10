function [sigTFSub, sigTFRet] = getSigTF(aggTargs, allDat, timMask, timMask2, encTim, retTim)
frex = logspace(log10(2),log10(80),100);
%reg  X tVals / HFB diff / p values  X time X frequency
sigTFSub = zeros([length(aggTargs),  3,  sum(timMask==1), 100] );
sigTFRet = 1; 

for reg1 = 1:length(aggTargs)
    reg1
    regRes = nan([2, sum(timMask==1), 100, 100]);
    regRes2 = nan([2, sum(timMask2==1), 100, 100]); 
    regSubs = nan([100,1]);
    regSubIDs = cell(100,1); 
    regd = regSubs; 
    regacc = regSubs; 
    chani = regd; 
    submissRT = regd; 
    subhitRT = regd; 
    retmissRT = regd; 
    rethitRT = regd; 
    ri = 1; 
    for sub = 1:length(allDat)
        if ~isempty(allDat{sub})
            T = sum(allDat{sub}.retInfo(:,1)==1 | allDat{sub}.retInfo(:,1)==2); 
            Hr = sum(allDat{sub}.retInfo(:,1)==1) / T; 
            T = sum(allDat{sub}.retInfo(:,1)==3 | allDat{sub}.retInfo(:,1)==4); 
            F = sum(allDat{sub}.retInfo(:,1)==4);
            if F == 0
                acc = Hr - F/T; 
                F = 1; 
            else
                acc = Hr - F/T; 
            end
            Fr = F / T; 
           
            d = norminv(Hr) - norminv(Fr); 

            smRT = mean(allDat{sub}.encInfo(allDat{sub}.misses & allDat{sub}.use, 4));
            shRT = mean(allDat{sub}.encInfo(allDat{sub}.hits & allDat{sub}.use, 4));
            rmRT = mean(allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==2, 3));
            rhRT = mean(allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==1, 3));

            sm = allDat{sub}.TF.subMiss(:,timMask==1,:); 
            sh = allDat{sub}.TF.subHit(:,timMask==1,:); 
            rm = allDat{sub}.TF.hit_on(:,timMask2==1,:); 
            rh = allDat{sub}.TF.miss_on(:,timMask2==1,:);
             b = allDat{sub}.meetLabs(:,3); 
            ID = allDat{sub}.subID;
            reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
            if ~isempty(reg1i)
                for i1 = 1:length(reg1i)
                
                    regRes(1,:,:,ri) = sh(reg1i(i1),:,:); %mean over trials hits
                    regRes(2,:,:,ri) = sm(reg1i(i1),:,:); %mean over trials misses
                    regRes2(1,:,:,ri) = rh(reg1i(i1),:,:); %mean over trials hits
                    regRes2(2,:,:,ri) = rm(reg1i(i1),:,:); %mean over trials misses
                    regSubs(ri) = sub; 
                    regSubIDs{ri} = ID; 
                    regd(ri) = d; 
                    regacc(ri) = acc; 
                    chani(ri) = i1; 
                    submissRT(ri) = smRT; 
                    subhitRT(ri) = shRT; 
                    retmissRT(ri) = rmRT; 
                    rethitRT(ri) = rhRT;
                    ri = ri+1; 
                end
    
            end
        end
    end

    if ri <= 100
        regRes(:,:,:,ri:end) = []; 
        regRes2(:,:,:,ri:end) = []; 
        regSubs(ri:end) = []; 
        regSubIDs(ri:end) = [];
        regd(ri:end) = []; 
        regacc(ri:end) = []; 
        chani(ri:end) = [];
        submissRT(ri:end) = []; 
        subhitRT(ri:end) = []; 
        retmissRT(ri:end) = []; 
        rethitRT(ri:end) = [];
    end

    realID = cell(length(regd),1); 
    for ii = 1:length(regd)
        realID{ii} = [regSubIDs{ii} num2str(chani(ii))];

    end

    TFdat = struct; 
    TFdat.regRes = regRes; 
    TFdat.regRes2 = regRes2; 
    TFdat.regSubs = regSubs; 
    TFdat.regSubIDs = regSubIDs; 
    TFdat.aggTargs = aggTargs; 
    TFdat.reg1 = reg1; 
    TFdat.d = regd; 
    TFdat.acc = regacc; 
    TFdat.chani = chani; 
    TFdat.realID = realID; 
    TFdat.n_sub = length(unique(regSubs));
    TFdat.n_pair = length(regSubs); 
    TFdat.encTim = encTim;  
    TFdat.retTim = retTim;  
    TFdat.submissRT = submissRT; 
    TFdat.subhitRT = subhitRT; 
    TFdat.retmissRT = retmissRT; 
    TFdat.rethitRT = rethitRT;

    %% encoding mean difference

    hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
    
    tVals = myArrayT(hitVals, missVals,1);
    perms = 1000; 
    
    nullTs = squeeze(zeros([size(tVals), perms])); 
    parfor ii = 1:perms
    
        nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
    end
    [h, p, clusterinfo] = cluster_test(tVals, nullTs); 


    TFdat.tVals_sub = tVals; 
    TFdat.hitVals_sub = mean(hitVals); 
    TFdat.missVals_sub = mean(missVals);
    TFdat.p_sub = p;

    figure('visible', false, 'Position', [0 0 1000 800])
    subplot 311
    imagesc(TFdat.encTim, [], squeeze(mean(hitVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([aggTargs(reg1).lab ': subsequent hit z-score power'])
    colorbar


    subplot 312
    imagesc(TFdat.encTim, [], squeeze(mean(missVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([aggTargs(reg1).lab ': subsequent miss z-score power'])
    colorbar

    subplot 313
    imagesc( tVals')
    test = addRedOutline(p, .05, 'red');
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    xticks([1:20:141])
    xticklabels(TFdat.encTim([1:20:141]))
    title([aggTargs(reg1).lab ': sub hit v. miss t-values'])
    colorbar
    

    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\' 'sub_' aggTargs(reg1).lab '.jpg'])



    %% retrieval mean difference
    hitVals = permute(squeeze(regRes2(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(regRes2(2,:,:,:)), [3,1,2]);
    
    tVals = myArrayT(hitVals, missVals,1);
    perms = 1000; 
    
    nullTs = squeeze(zeros([size(tVals), perms])); 
    parfor ii = 1:perms
    
        nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
    end
    [h, p, clusterinfo] = cluster_test(tVals, nullTs); 


    TFdat.tVals_ret = tVals; 
    TFdat.hitVals_ret = mean(hitVals); 
    TFdat.missVals_ret = mean(missVals);
    TFdat.p_ret = p;

    figure('visible', false, 'Position', [0 0 1000 800])
    subplot 311
    imagesc(TFdat.encTim, [], squeeze(mean(hitVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([aggTargs(reg1).lab ': retrieval hit z-score power'])
    colorbar


    subplot 312
    imagesc(TFdat.encTim, [], squeeze(mean(missVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([aggTargs(reg1).lab ': retrieval miss z-score power'])
    colorbar

    subplot 313
    imagesc( tVals')
    test = addRedOutline(p, .05, 'red');
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    xticks([1:20:141])
    xticklabels(TFdat.encTim([1:20:141]))
    title([aggTargs(reg1).lab ': ret hit v. miss t-values'])
    colorbar
    

    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\' 'ret_' aggTargs(reg1).lab '.jpg'])

    save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\TF_KEY_STATS\' aggTargs(reg1).lab  '.mat'], 'TFdat', '-v7.3')

    
% 
%     sigTFSub(reg1, 1, :, :) = tVals; 
%     sigTFSub(reg1, 2, :, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
%     sigTFSub(reg1, 3, :, :) = p; 
    disp(['........................' num2str(round(toc/60, 1))])



end
end

 