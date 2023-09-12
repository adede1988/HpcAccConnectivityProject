function [sigHFBSub, sigHFBRet] = getSigHFB(aggTargs, allDat, timMask, timMask2)

%reg  X tVals / HFB diff / p values  X time
sigHFBSub = zeros([length(aggTargs),  3,  sum(timMask==1)] );
sigHFBRet = zeros([length(aggTargs),  3,  sum(timMask2==1)] );

for reg1 = 1:length(aggTargs)
    reg1
    regRes = nan([2, sum(timMask==1), 100]);
    regRes2 = nan([2, sum(timMask2==1), 100]); 
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


            sm = allDat{sub}.HFB.subMiss(:,timMask==1,:); 
            sh = allDat{sub}.HFB.subHit(:,timMask==1,:); 
            rm = allDat{sub}.HFB.hit_on(:,timMask2==1,:); 
            rh = allDat{sub}.HFB.miss_on(:,timMask2==1,:); 
            b = allDat{sub}.brodmann; 
            ID = allDat{sub}.subID;
            reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
            if ~isempty(reg1i)
                for i1 = 1:length(reg1i)
                
                    regRes(1,:,ri) = mean(sh(reg1i(i1),:,:),3); %mean over trials hits subsequent
                    regRes(2,:,ri) = mean(sm(reg1i(i1),:,:),3); %mean over trials misses subsequent
                    regRes2(1,:,ri) = mean(rh(reg1i(i1),:,:),3); %mean over trials hits
                    regRes2(2,:,ri) = mean(rm(reg1i(i1),:,:),3); %mean over trials misses
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
        regRes(:,:,ri:end) = []; 
        regRes2(:,:,ri:end) = []; 
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

    HFBdat = struct; 
    HFBdat.regRes = regRes; 
    HFBdat.regRes2 = regRes2; 
    HFBdat.regSubs = regSubs; 
    HFBdat.regSubIDs = regSubIDs; 
    HFBdat.aggTargs = aggTargs; 
    HFBdat.reg1 = reg1; 
    HFBdat.d = regd; 
    HFBdat.acc = regacc; 
    HFBdat.chani = chani; 
    HFBdat.realID = realID; 
    HFBdat.n_sub = length(unique(regSubs));
    HFBdat.n_pair = length(regSubs); 
    HFBdat.encTim = allDat{3}.leadLag.encTim;  
    HFBdat.retTim = allDat{3}.leadLag.retTim;  
    HFBdat.submissRT = submissRT; 
    HFBdat.subhitRT = subhitRT; 
    HFBdat.retmissRT = retmissRT; 
    HFBdat.rethitRT = rethitRT;



    save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\' aggTargs(reg1).ROI  '.mat'], 'HFBdat', '-v7.3')



%    
% 
%     figure('visible', false, 'Position', [0 0 1000 800])
%     subplot 211
%     plot(allDat{1}.leadLag.encTim, mean(hitVals,1))
%     hold on 
%     plot(allDat{1}.leadLag.encTim,mean(missVals,1))
%     title(aggTargs(reg1).ROI)
% 
%     subplot 212
%     plot(allDat{1}.leadLag.encTim, mean(hitVals,1) - mean(missVals,1), 'color', 'green', 'linewidth', 3)
%     hold on 
%     plot(allDat{1}.leadLag.encTim, prctile(nullTs,97.5, 2), 'color', 'k')
%     plot(allDat{1}.leadLag.encTim, prctile(nullTs,2.5, 2), 'color', 'k')
%     P_plot = p; 
%     P_plot(p<.05) = max(prctile(nullTs,97.5, 2))*1.1; 
%     P_plot(p>.05) = nan;
%     plot(allDat{1}.leadLag.encTim, P_plot, 'color', 'red', 'linewidth', 3)
%     title('hit - miss difference and 95% confidence interval of shuffle')
% 
%     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\' aggTargs(reg1).ROI '_HFB' '.jpg'])
% 
%     
% 
%     sigHFBSub(reg1, 1, :) = tVals; 
%     sigHFBSub(reg1, 2, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
%     sigHFBSub(reg1, 3, :) = p; 
%     disp(['........................' num2str(round(toc/60, 1))])



end
end

 
%         for sub = 1:length(allDat)
%            
%             if ~isempty(allDat{sub})
%                 c = allDat{sub}.leadLag; 
%                 b = allDat{sub}.brodmann; 
%                 
%                 reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
%                 reg2i = find(cellfun(@(x) sum(strcmp(aggTargs(reg2).lab, x)), b));
% 
%                 if ~isempty(reg1i) && ~isempty(reg2i) 
%                     for i1 = 1:length(reg1i)
%                         for i2 = 1:length(reg2i)
%                             if reg1i(i1) ~= reg2i(i2)
%                             regRes(:, :, :, ri) = c.subMem(reg1i(i1), reg2i(i2), :, :, :); 
%                             ri = ri + 1; 
%                             end
%                         end
%                     end
% 
% 
%                 end
% 
%                
% 
% 
%             end
%         end
%        
% 
%         hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
%         missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
%         tVals = myArrayT(hitVals, missVals,1);
%         perms = 100; 
% 
%         nullTs = zeros([size(tVals), perms]); 
%         parfor ii = 1:perms
% 
%             nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
%         end
% 
%         [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
%             
%        
%     end
% end
% 
% 
% 










% end