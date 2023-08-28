function [sigTFSub, sigTFRet] = getSigTF(aggTargs, allDat, timMask)
frex = logspace(log10(2),log10(80),100);
%reg  X tVals / HFB diff / p values  X time X frequency
sigTFSub = zeros([length(aggTargs),  3,  sum(timMask==1), 100] );
sigTFRet = 1; 

for reg1 = 1:length(aggTargs)
    reg1
    regRes = nan([2, sum(timMask==1), 100, 100]);
    ri = 1; 
    for sub = 1:length(allDat)
        if ~isempty(allDat{sub})
            sm = allDat{sub}.TF.subMiss(:,timMask==1,:); 
            sh = allDat{sub}.TF.subHit(:,timMask==1,:); 
            b = allDat{sub}.brodmann; 
            reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
            if ~isempty(reg1i)
                for i1 = 1:length(reg1i)
                
                    regRes(1,:,:,ri) = sh(reg1i(i1),:,:); %mean over trials hits
                    regRes(2,:,:,ri) = sm(reg1i(i1),:,:); %mean over trials misses
                    ri = ri+1; 
                end
    
            end
        end
    end

    if ri <= 100
        regRes(:,:,:,ri:end) = []; 
    end


    hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
    missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
    
    tVals = myArrayT(hitVals, missVals,1);
    perms = 100; 
    
    nullTs = squeeze(zeros([size(tVals), perms])); 
    parfor ii = 1:perms
    
        nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
    end
    [h, p, clusterinfo] = cluster_test(tVals, nullTs); 

    figure('visible', false, 'Position', [0 0 1000 800])
    subplot 311
    imagesc(allDat{1}.leadLag.encTim, [], squeeze(mean(hitVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([aggTargs(reg1).ROI ': subsequent hit z-score power'])
    colorbar


    subplot 312
    imagesc(allDat{1}.leadLag.encTim, [], squeeze(mean(missVals,1))')
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    title([aggTargs(reg1).ROI ': subsequent miss z-score power'])
    colorbar

    subplot 313
    imagesc( tVals')
    test = addRedOutline(p);
    yticks([10:20:100])
    yticklabels(round(frex([10:20:100])))
    set(gca, 'YDir','normal')
    xticks([1:20:141])
    xticklabels(allDat{1}.leadLag.encTim([1:20:141]))
    title([aggTargs(reg1).ROI ': sub hit v. miss t-values'])
    colorbar
    

    export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\' aggTargs(reg1).ROI '.jpg'])

    

    

    sigTFSub(reg1, 1, :, :) = tVals; 
    sigTFSub(reg1, 2, :, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
    sigTFSub(reg1, 3, :, :) = p; 
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