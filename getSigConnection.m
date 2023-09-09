function [sigConSub, sigConRet, aggTargs] = getSigConnection(aggTargs, allDat)
dimVals = size(allDat{3}.leadLag.subMem, [4,5]);
dimVals2 = size(allDat{3}.leadLag.retMem, [4,5]); 
%reg X reg X tVals / corr diff / p values X offset X time
sigConSub = zeros([length(aggTargs), length(aggTargs), 3,  dimVals] );
sigConRet = zeros([length(aggTargs), length(aggTargs), 3,  dimVals2] ); 

for reg1 = 1:length(aggTargs)
    reg1
    for reg2 = 1:length(aggTargs)
        tic
        disp(['working on: ' num2str(reg1) ' X ' num2str(reg2)])
        regRes = nan([2, dimVals, 100]);
        regRes2 = nan([2, dimVals2, 100]); %retrieval 
        ri = 1; 
        for sub = 1:length(allDat)
           
            if ~isempty(allDat{sub})
                c = allDat{sub}.leadLag; 
                b = allDat{sub}.brodmann; 
                
                reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
                reg2i = find(cellfun(@(x) sum(strcmp(aggTargs(reg2).lab, x)), b));

                if ~isempty(reg1i) && ~isempty(reg2i) 
                    for i1 = 1:length(reg1i)
                        for i2 = 1:length(reg2i)
                            if reg1i(i1) ~= reg2i(i2)
                            regRes(:, :, :, ri) = c.subMem(reg1i(i1), reg2i(i2), :, :, :); 
                            regRes2(:,:, :, ri) = c.retMem(reg1i(i1), reg2i(i2), :, :, :); 
                            ri = ri + 1; 
                            end
                        end
                    end


                end

               


            end
        end
        if ri < 100
            regRes(:,:,:,:,:,ri:end) = []; 
            regRes2(:,:,:,:,:,ri:end) = [];
        end

        hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
        missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
        tVals = myArrayT(hitVals, missVals,1);
        perms = 1000; 

        nullTs = zeros([size(tVals), perms]); 
        parfor ii = 1:perms

            nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
        end

        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
            
        sigConSub(reg1, reg2, 1, :, :) = tVals; 
        sigConSub(reg1, reg2, 2, :, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
        sigConSub(reg1, reg2, 3, :, :) = p; 

        hitVals = permute(squeeze(regRes2(1,:,:,:)), [3,1,2]);
        missVals = permute(squeeze(regRes2(2,:,:,:)), [3,1,2]);
        tVals = myArrayT(hitVals, missVals,1);
        perms = 1000; 

        nullTs = zeros([size(tVals), perms]); 
        parfor ii = 1:perms

            nullTs(:,:,ii) = myArrayT(hitVals, missVals, 2);
        end

        [h, p, clusterinfo] = cluster_test(tVals, nullTs); 
            
        sigConRet(reg1, reg2, 1, :, :) = tVals; 
        sigConRet(reg1, reg2, 2, :, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
        sigConRet(reg1, reg2, 3, :, :) = p; 





        disp(['........................' num2str(round(toc/60, 1))])
    end
end













end