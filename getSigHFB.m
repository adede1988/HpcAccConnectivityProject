function [sigHFBSub, sigHFBRet] = getSigHFB(aggTargs, allDat, timMask)

%reg  X tVals / HFB diff / p values  X time
sigHFBSub = zeros([length(aggTargs),  3,  sum(timMask==1)] );
sigHFBRet = 1; 

for reg1 = 1:length(aggTargs)
    reg1
    regRes = nan([2, sum(timMask==1), 100]);
    ri = 1; 
    for sub = 1:length(allDat)
        if ~isempty(allDat{sub})
            sm = allDat{sub}.HFB.subMiss(:,timMask==1,:); 
            sh = allDat{sub}.HFB.subHit(:,timMask==1,:); 
            b = allDat{sub}.brodmann; 
            reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
            if ~isempty(reg1i)
                for i1 = 1:length(reg1i)
                
                    regRes(1,:,ri) = mean(sh(reg1i(i1),:,:),3); %mean over trials hits
                    regRes(2,:,ri) = mean(sm(reg1i(i1),:,:),3); %mean over trials misses
                    ri = ri+1; 
                end
    
            end
        end
    end

    if ri <= 100
        regRes(:,:,ri:end) = []; 
    end


    hitVals = permute(squeeze(regRes(1,:,:)), [2,1]);
    missVals = permute(squeeze(regRes(2,:,:)), [2,1]);
    
    tVals = myArrayT(hitVals, missVals,1);
    perms = 100; 
    
    nullTs = squeeze(zeros([size(tVals), perms])); 
    parfor ii = 1:perms
    
        nullTs(:,ii) = myArrayT(hitVals, missVals, 2);
    end
    [h, p, clusterinfo] = cluster_test(tVals', nullTs); 

    figure('visible', false, 'Position', [0 0 1000 800])
    subplot 211
    plot(allDat{1}.leadLag.encTim, mean(hitVals,1))
    hold on 
    plot(allDat{1}.leadLag.encTim,mean(missVals,1))
    title(aggTargs(reg1).ROI)

    subplot 212
    plot(allDat{1}.leadLag.encTim, mean(hitVals,1) - mean(missVals,1), 'color', 'green', 'linewidth', 3)
    hold on 
    plot(allDat{1}.leadLag.encTim, prctile(nullTs,97.5, 2), 'color', 'k')
    plot(allDat{1}.leadLag.encTim, prctile(nullTs,2.5, 2), 'color', 'k')
    P_plot = p; 
    P_plot(p<.05) = max(prctile(nullTs,97.5, 2))*1.1; 
    P_plot(p>.05) = nan;
    plot(allDat{1}.leadLag.encTim, P_plot, 'color', 'red', 'linewidth', 3)
    title('hit - miss difference and 95% confidence interval of shuffle')

    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\' aggTargs(reg1).ROI '_HFB' '.jpg'])

    

    sigHFBSub(reg1, 1, :) = tVals; 
    sigHFBSub(reg1, 2, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
    sigHFBSub(reg1, 3, :) = p; 
    disp(['........................' num2str(round(toc/60, 1))])



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