function [sigHFBSub, sigHFBRet] = getSigHFB_latencyPatch(aggTargs, allDat, timMask, timMask2)



%reg  X tVals / HFB diff / p values  X time
sigHFBSub = zeros([length(aggTargs),  3,  sum(timMask==1)] );
sigHFBRet = zeros([length(aggTargs),  3,  sum(timMask2==1)] );

for reg1 = 1:length(aggTargs)
    reg1
    HFBdat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\' aggTargs(reg1).ROI  '.mat']).HFBdat;
%     HFBdat.hitTim_sub_old = HFBdat.hitTim_sub; 
%     HFBdat.missTim_sub_old = HFBdat.missTim_sub;
%     HFBdat.hitTim_ret_old = HFBdat.hitTim_ret; 
%     HFBdat.missTim_ret_old = HFBdat.missTim_ret;
    ri = 1; 
    for sub = 1:length(allDat)
        if ~isempty(allDat{sub})

            
            sm = allDat{sub}.HFB.subMiss(:,timMask==1,:); 
            sh = allDat{sub}.HFB.subHit(:,timMask==1,:); 
            rh = allDat{sub}.HFB.hit_on(:,timMask2==1,:); 
            rm = allDat{sub}.HFB.miss_on(:,timMask2==1,:); 
            b = allDat{sub}.brodmann; 
            ID = allDat{sub}.subID;
            reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
            if ~isempty(reg1i)
                for i1 = 1:length(reg1i)
                    sm_cm = []; 
                    sh_cm = []; 
                    rm_cm = []; 
                    rh_cm = []; 
                    RT = allDat{sub}.encInfo(allDat{sub}.use & allDat{sub}.misses,4);
                    for tt = 1:size(sm,3) %trial loop subMiss
                        threshTestWin = [-50 2500]; %window to look in 
                        curTrial = sm(reg1i(i1), :, tt); %single sub miss trial
                        test = curTrial( HFBdat.encTim>=threshTestWin(1) & HFBdat.encTim<=RT(tt));
                        testTim = HFBdat.encTim(HFBdat.encTim>=threshTestWin(1) & HFBdat.encTim<=RT(tt));
                        test = (test - min(test)); 
                        test = test./ max(test);
                        test(test<.6) = 0; 
                        testMean = wmean([1:length(testTim)], test); 
                        sm_cm = [sm_cm testTim(round(testMean))];
                    end
                    RT = allDat{sub}.encInfo(allDat{sub}.use & allDat{sub}.hits,4);
                    for tt = 1:size(sh,3)%trial loop subHit
                        threshTestWin = [-50 2500]; %window to look in 
                        curTrial = sh(reg1i(i1), :, tt); %single sub miss trial
                        test = curTrial( HFBdat.encTim>=threshTestWin(1) & HFBdat.encTim<=RT(tt));
                        testTim = HFBdat.encTim(HFBdat.encTim>=threshTestWin(1) & HFBdat.encTim<=RT(tt));
                        test = (test - min(test)); 
                        test = test./ max(test);
                        test(test<.6) = 0; 
                        testMean = wmean([1:length(testTim)], test); 
                        sh_cm = [sh_cm testTim(round(testMean))];
                    end
                    RT = allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==2,3);
                    for tt = 1:size(rm,3) %trial loop retrieval miss
                        threshTestWin = [-50 2500]; %window to look in 
                        curTrial = rm(reg1i(i1), :, tt); %single sub miss trial
                        test = curTrial( HFBdat.retTim>=threshTestWin(1) & HFBdat.retTim<=RT(tt));
                        testTim = HFBdat.retTim(HFBdat.retTim>=threshTestWin(1) & HFBdat.retTim<=RT(tt));
                        test = (test - min(test)); 
                        test = test./ max(test);
                        test(test<.6) = 0; 
                        testMean = wmean([1:length(testTim)], test); 
                        rm_cm = [rm_cm testTim(round(testMean))];
                    end
                    RT = allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==1,3);
                    for tt = 1:size(rh,3)%trial loop retrieval hit
                        threshTestWin = [-50 2500]; %window to look in 
                        curTrial = rh(reg1i(i1), :, tt); %single sub miss trial
                        test = curTrial( HFBdat.retTim>=threshTestWin(1) & HFBdat.retTim<=RT(tt));
                        testTim = HFBdat.retTim(HFBdat.retTim>=threshTestWin(1) & HFBdat.retTim<=RT(tt));
                        test = (test - min(test)); 
                        test = test./ max(test);
                        test(test<.6) = 0; 
                        testMean = wmean([1:length(testTim)], test); 
                        rh_cm = [rh_cm testTim(round(testMean))];
                    end
                    
                    HFBdat.hitTim_sub(ri) = mean(sh_cm); 
                    HFBdat.missTim_sub(ri) = mean(sm_cm); 
                    HFBdat.hitTim_ret(ri) = mean(rh_cm); 
                    HFBdat.missTim_ret(ri) = mean(rm_cm);
                 

                    ri = ri+1; 
                end
    
            end
        end
    end

   




    save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\' aggTargs(reg1).ROI  '.mat'], 'HFBdat', '-v7.3')






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