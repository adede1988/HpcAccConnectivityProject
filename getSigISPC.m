function [conN, conID] = getSigISPC(aggTargs, allDat, timMask, timMask2, encTim, retTim)
frex = logspace(log10(2), log10(25), 20); 

%ISPC stats are ISPC, PPC, ISPC z-scored, PPC z-scored
conN = zeros(length(aggTargs)); 
conID = cell(length(aggTargs)); 
for reg1 = 1:length(aggTargs)
%    reg1 = 10;                        % DEBUG HARDCODE OF REGION!
    for reg2 = 1:length(aggTargs)
%        reg2 = 1;                    % DEBUG HARDCODE OF REGION!!

        tic
        disp(['working on: ' num2str(reg1) ' X ' num2str(reg2)])
  
        %hit/miss X frequency X time  X pair (get PCC only)
        regRes = nan([2,  sum(timMask==1), length(frex),  100]);
        regRes2= nan([2, sum(timMask2==1), length(frex),  100]); 
        regSubs = nan([100,1]);
        regSubIDs = cell(100,1); 
        chani = cell(100,1); 
        regd = regSubs; 
        regacc = regSubs; 
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
               
               

                c = allDat{sub}.ISPC; 
                c = rmfield(c, {'subMissR', 'subHitR', 'cr_on', 'fa_on', 'hit_rt', 'miss_rt', 'cr_rt', 'fa_rt'});
                c.subMiss = squeeze(c.subMiss(:,:,:,:,2)); 
                c.subHit = squeeze(c.subHit(:,:,:,:,2)); 
                c.hit_on = squeeze(c.hit_on(:,:,:,:,2)); 
                c.miss_on = squeeze(c.miss_on(:,:,:,:,2)); 


                smRT = mean(allDat{sub}.encInfo(allDat{sub}.misses & allDat{sub}.use, 4));
                shRT = mean(allDat{sub}.encInfo(allDat{sub}.hits & allDat{sub}.use, 4));
                rmRT = mean(allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==2, 3));
                rhRT = mean(allDat{sub}.retInfo(allDat{sub}.retInfo(:,1)==1, 3));
                b = allDat{sub}.meetLabs(:,3); 
                ID = allDat{sub}.subID;
                reg1i = find(cellfun(@(x) sum(strcmp(aggTargs(reg1).lab, x)), b));
                reg2i = find(cellfun(@(x) sum(strcmp(aggTargs(reg2).lab, x)), b));

                if ~isempty(reg1i) && ~isempty(reg2i) 
                    for i1 = 1:length(reg1i)
                        for i2 = 1:length(reg2i)
                            if reg1i(i1) ~= reg2i(i2)
                                sh = c.subHit(reg1i(i1), reg2i(i2), timMask==1, :);
                                sm = c.subMiss(reg1i(i1), reg2i(i2), timMask==1, :);
                                rh = c.hit_on(reg1i(i1), reg2i(i2), timMask2==1, :);
                                rm = c.miss_on(reg1i(i1), reg2i(i2), timMask2==1, :);
                                
                                regRes(1, :, :,  ri) = sh; 
                                regRes(2, :, :,  ri) = sm; 
                                regRes2(1, :, :,  ri) = rh; 
                                regRes2(2, :, :,  ri) = rm; 
                                clear sh sm rh rm
                                regd(ri) = d; 
                                conN(reg1, reg2) = conN(reg1, reg2) + 1; 
                                regSubs(ri) = sub; 
                                regSubIDs{ri} = ID; 
                                regacc(ri) = acc; 
                                submissRT(ri) = smRT; 
                                subhitRT(ri) = shRT; 
                                retmissRT(ri) = rmRT; 
                                rethitRT(ri) = rhRT; 
                                chani{ri} = [num2str(reg1i(i1)) '_' num2str(reg2i(i2))];


                                
                                ri = ri + 1; 
                            


                            end
                        end
                    end
                clear c

                end

               

          
            end
        end
        try
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
        
        subIDs = unique(regSubIDs); 
        if length(subIDs) >= 5  %subject threshold!!! 

            %encoding
            hitVals = permute(squeeze(regRes(1,:,:,:)), [3,1,2]);
            missVals = permute(squeeze(regRes(2,:,:,:)), [3,1,2]);
            combo = [hitVals; missVals]; 

            %rewriting to go after the specific bands! 
            %low band, frex indices 1-3
            %high band, frex indices 11-12
            lowBand = squeeze(mean(combo(:,:,1:3), 3)); 
            highBand = squeeze(mean(combo(:,:,11:12), 3)); 

            sub_d = regd(cellfun(@(y) find(cellfun(@(x) strcmp(y,x), regSubIDs) ,1) , subIDs));
    
            hmSort = [ones(length(regd),1); zeros(length(regd),1)]; 
            hmSort = hmSort == 1; 

              
            connectionDat = struct; 
            connectionDat.reg1 = aggTargs(reg1).lab;
            connectionDat.reg2 = aggTargs(reg2).lab;
            connectionDat.reg1i = reg1; 
            connectionDat.reg2i = reg2; 
            connectionDat.N = length(regSubIDs);
            connectionDat.subN = length(subIDs); 
            connectionDat.uniqueSubs = subIDs; 
            connectionDat.subD = sub_d; 
            connectionDat.d = regd; 
            connectionDat.allSubs = regSubIDs; 
            connectionDat.lowBand = lowBand; 
            connectionDat.highBand = highBand; 
            connectionDat.tim = encTim; 
            connectionDat.hmSort = hmSort; 
            connectionDat.regacc= regacc; 
            connectionDat.chani = chani;
            connectionDat.submissRT = submissRT; 
            connectionDat.subhitRT = subhitRT; 
            connectionDat.retmissRT = retmissRT; 
            connectionDat.rethitRT = rethitRT;


            %retrieval
            hitVals = permute(squeeze(regRes2(1,:,:,:)), [3,1,2]);
            missVals = permute(squeeze(regRes2(2,:,:,:)), [3,1,2]);
            combo = [hitVals; missVals]; 

            %rewriting to go after the specific bands! 
            %low band, frex indices 1-3
            %high band, frex indices 11-12
            lowBand = squeeze(mean(combo(:,:,1:3), 3)); 
            highBand = squeeze(mean(combo(:,:,11:12), 3)); 

              
            connectionDat.lowBand2 = lowBand; 
            connectionDat.highBand2 = highBand; 
            connectionDat.tim2 = retTim; 

            save(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\PCC_KEY_STATS\' connectionDat.reg1 '_' connectionDat.reg2 '.mat'], 'connectionDat')

        end
        catch
            'what happened?'
        end
    end
end
end
        

% allReg = nan([60000, sum(timMask==1), length(frex)]); 
% allRegi = 1; 
% 
% 
% for sub = 1:length(allDat)
%    
%     if ~isempty(allDat{sub})
% 
%         if allDat{sub}.age>16  % AGE FILTER!!!!!!
%             c = allDat{sub}.ISPC;
%             for nulli = 1:size(c.subMiss,1)
%                 allReg(allRegi:allRegi+size(c.subMiss,1)-1 ,:,:) = squeeze(c.subMiss(nulli, :, timMask==1, :, 2));
%                 allRegi = allRegi + size(c.subMiss,1); 
%                 allReg(allRegi:allRegi+size(c.subMiss,1)-1 ,:,:)  = squeeze(c.subHit(nulli, :, timMask==1, :, 2)); 
%                 allRegi = allRegi + size(c.subHit,1); 
%             end
% 
%         end
%         
%     end
% end
% clear c curSamp 
% 
% 
% if allRegi<=60000
%     allReg(allRegi:end, :, :) = []; 
% end
% 
% 
% 'wait'






%         %% implementation of linear mixed effects model! 
%         
%         lowDat = lowBand(hmSort, :); %log transformed z-score hit data only
%         highDat = highBand(hmSort, :); %log transofmred z-score hit data only
%         tVals = zeros(size(lowDat,2),1);
% 
%         
% 
%         tic
%         parfor ti = 1:size(lowDat, 2)
%             ti
%             modDat = table(lowDat(:,ti), alld', allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
%             lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
%             tVals(ti) = lme.Coefficients(2,4); %t-value associated with memory! 
% 
%         end
%         toc
%         %potentially use this to get a normally distributed rank value for
%         %d'
% %         rankit = INVNORMAL ((Rank of X - 0.5)/ n)
%         
% 
%         perms = 20; 
%         nullTs = zeros([length(tVals), perms]); 
% 
%         parfor ii = 1:perms
% 
%             shuffd = randsample(alld, length(alld), false); 
%             sliceT = nullTs(:,ii); 
%             for ti = 1:size(lowDat,2)
%                 modDat = table(lowDat(:,ti), shuffd', allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
%                 lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
%                 sliceT(ti) = lme.Coefficients(2,4); %t-value associated with memory! 
%             end
%             nullTs(:,ii) = sliceT; 
%         end
% 
% 
%         [~,p, ~] = cluster_test(tVals, nullTs); 
% 
% 
%         if sum(p<.05) >0
%             permP = zeros([length(subIDs), length(p)]); 
%             'wait'
%             %try leaving out each subject one at a time to confirm if results
%             %are the same
%             for subOut = 1:length(subIDs)
%                 subOut
%                 allSubs = conID{reg1, reg2};
%                 subSelect = cellfun(@(x) ~strcmp(x, subIDs{subOut}), allSubs); 
%                 cur = lowDat(subSelect,:);
%                 curd = alld(subSelect); 
%                 cursub = allSubs(subSelect); 
%                 tVals = 1;
%                 for ti = 1:size(lowDat,2)
%                     modDat = table(cur(:,ti), shuffd', allSubs', 'VariableNames', {'connectivity', 'memory', 'sub'}); 
%                     lme = fitlme(modDat, 'connectivity ~ memory + (1|sub)'); 
%                     sliceT(ti) = lme.Coefficients(2,4); %t-value associated with memory!
%                 end
% 
%                 %get the overall significance
%                 Fvals = myArrayAnova(cur, curHM,  curd, 1);
%               
%                 perms = 1000; 
%         
%                 nullTs = zeros([size(Fvals), perms]); 
%                 parfor ii = 1:perms
%         
%                     nullTs(:,:,:, ii) = myArrayAnova(cur, curHM,  curd, 2);
%         
%                 end
%                  %main effect hit/miss
%                 [~, permP(subOut,1,:,:), ~] = cluster_test(squeeze(Fvals(1,:,:)), squeeze(nullTs(1,:,:,:)) ); 
%                 %main effect hit/miss
%                 [~, permP(subOut,2,:,:), ~] = cluster_test(squeeze(Fvals(2,:,:)), squeeze(nullTs(2,:,:,:)) ); 
%                 %interaction
%                 [~, permP(subOut,3,:,:), ~] = cluster_test(squeeze(Fvals(3,:,:)), squeeze(nullTs(3,:,:,:)) ); 
%             end
%         else
%             permP = ones([length(subIDs), size(allP)]); 
%         end
% 
%         %main effect of hit miss
% %         [h, p, clusterinfo] = cluster_test(squeeze(Fvals(1,:,:)), squeeze(nullTs(1,:,:,:)) ); 
%             
% % 
% %     dSplit_half = alld > med_d;
%     figure('visible', true, 'Position', [0 0 1000 2000])
% 
%     %hits for good memory performers
%     subplot 521
%     p1 = makeScaledHeat(combo(dSplit' & hmSort, :, :), allDat{1}.leadLag.encTim, frex, adjFact);
%     allMax = max(squeeze(mean(p1 ,1)), [], 'all'); 
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' sub Hit connectivity'])
% 
%     subplot 522
%     p1 = makeScaledHeat(combo(~dSplit' & hmSort, :, :), allDat{1}.leadLag.encTim, frex, adjFact);
%     allMax = max([allMax, max(squeeze(mean(p1 ,1)), [], 'all')]);
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' sub Hit connectivity'])
%  
%     subplot 523
%     p1 = makeScaledHeat(combo(dSplit' & ~hmSort, :, :), allDat{1}.leadLag.encTim, frex, adjFact);
%     allMax = max([allMax, max(squeeze(mean(p1 ,1)), [], 'all')]);
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' sub Miss connectivity'])
% 
% 
%     subplot 524
%     p1 = makeScaledHeat(combo(~dSplit' & ~hmSort, :, :), allDat{1}.leadLag.encTim, frex, adjFact);
%     allMax = max([allMax, max(squeeze(mean(p1 ,1)), [], 'all')]);
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' sub Miss connectivity'])
% 
%     caxis([0,allMax*.75])
%     subplot 521
%     caxis([0, allMax*.75])
%     subplot 522
%     caxis([0, allMax*.75])
%     subplot 523
%     caxis([0, allMax*.75])
% 
% 
% 
% 
%     subplot(5,2,5)
%     %main effect of hit miss
% %     [h, p, clusterinfo] = cluster_test(squeeze(Fvals(1,:,:)), squeeze(nullTs(1,:,:,:)) ); 
%     p = squeeze(permP(:,1,:,:));
%     p = squeeze(sum(p<.05,1)); 
%     imagesc( p')
%     caxis([0,length(subIDs)])
%     pOut = ones(size(p)); 
%     pOut(p == length(subIDs)) = .01; 
%     maxAgree = max(p, [], 'all'); %how many different subjects can be left out and hold significance
%     mask = zeros(size(p));
%     mask(p==maxAgree) = 1; 
%     curF = squeeze(Fvals(1,:,:));
%     curF(mask==0) = 0; 
%     [mfv, mfi] = max(curF'); 
%     [mtv, mti] = max(mfv); 
%     mfi = mfi(mti); 
%     test = addRedOutline(pOut, .05);
%     hold on 
%     scatter(mti, mfi, 50, 'red', 'filled')
%     yticks([1:2:20])
%     yticklabels(round(frex(1:2:20)))
%     set(gca, 'YDir','normal')
%     xticks([1:20:141])
%     xticklabels(allDat{1}.leadLag.encTim([1:20:141]))
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' hit - miss leave one out sig count'])
%     colorbar
% 
%     subplot(5,2,6)
%     hold off
%     plotVals = squeeze(combo(:,mti,mfi)); 
%     plotVals = 10.^plotVals; 
%     plotVals = plotVals + adjFact - .01; 
%     b = boxchart(hmSort*1, plotVals);
%     b.MarkerStyle = 'none'; 
%     xticks([0,1])
%     xticklabels({'subMiss', 'subHit'})
%     title(['hit v. miss at: ' num2str(round(frex(mfi))) ' Hz  X ' ...
%         num2str(allDat{1}.leadLag.encTim(mti)) ' ms; p=' num2str(round(p(mti,mfi),2)) ])
%     
%     PLH = plotVals(hmSort); 
%     PLM = plotVals(~hmSort); 
%     randVals = (rand(length(PLH),1)-.5)*.5;
%     hold on 
%     scatter(randVals, PLM, 10,  'blue')
%     scatter(randVals+1, PLH, 10, 'blue')
%     
%     for pi = 1:length(PLH)
%         plot([0+randVals(pi),1+randVals(pi)], [PLM(pi),PLH(pi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
%         
% 
%     end
%     yline(0, 'linewidth', 3, 'linestyle', '--', 'color', [169,169,169]./256)
%   
%     
% 
% 
%     subplot(5,2,7)
%     p = squeeze(permP(:,2,:,:));
%     p = squeeze(sum(p<.05,1)); 
%     imagesc( p')
%     caxis([0,length(subIDs)])
%     pOut = ones(size(p)); 
%     pOut(p == length(subIDs)) = .01; 
%     maxAgree = max(p, [], 'all'); %how many different subjects can be left out and hold significance
%     mask = zeros(size(p));
%     mask(p==maxAgree) = 1; 
%     curF = squeeze(Fvals(1,:,:));
%     curF(mask==0) = 0; 
%     [mfv, mfi] = max(curF'); 
%     [mtv, mti] = max(mfv); 
%     mfi = mfi(mti); 
%     test = addRedOutline(pOut, .05);
%     hold on 
%     scatter(mti, mfi, 50, 'red', 'filled')
%     yticks([1:2:20])
%     yticklabels(round(frex(1:2:20)))
%     set(gca, 'YDir','normal')
%     xticks([1:20:141])
%     xticklabels(allDat{1}.leadLag.encTim([1:20:141]))
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' good/bad mem leave one out sig count'])
%     colorbar
% 
%     subplot(5,2,8)
%     hold off
%     plotVals = squeeze(combo(:, mti, mfi));
%     plotVals = 10.^plotVals; 
%     plotVals = plotVals + adjFact - .01; 
%     scatter([alld,alld], plotVals)
%     subMeans = zeros(length(subIDs),1); 
%    
%     curIDs = conID{reg1, reg2};
%     for si = 1:length(subMeans)
%         idx = cellfun(@(x) strcmp(x, subIDs{si}), curIDs); 
%         subMeans(si) = mean(plotVals(idx)); 
%     end
%    hold on 
%    scatter(sub_d, subMeans, 100, 'red', 'filled')
%    yline(0, 'linewidth', 3, 'linestyle', '--', 'color', [169,169,169]./256)
%     title(['good/bad mem at: ' num2str(round(frex(mfi))) ' Hz  X ' ...
%         num2str(allDat{1}.leadLag.encTim(mti)) ' ms; p=' num2str(round(p(mti,mfi),2)) ])
%     xlim([min(sub_d)-.2, max(sub_d)+.2])
% 
% 
%     subplot(5,2,9)
%     p = squeeze(permP(:,3,:,:));
%     p = squeeze(sum(p<.05,1)); 
%     imagesc( p')
%     caxis([0,length(subIDs)])
%     pOut = ones(size(p)); 
%     pOut(p == length(subIDs)) = .01; 
%     maxAgree = max(p, [], 'all'); %how many different subjects can be left out and hold significance
%     mask = zeros(size(p));
%     mask(p==maxAgree) = 1; 
%     curF = squeeze(Fvals(1,:,:));
%     curF(mask==0) = 0; 
%     [mfv, mfi] = max(curF'); 
%     [mtv, mti] = max(mfv); 
%     mfi = mfi(mti);
%     test = addRedOutline(pOut, .05);
%     hold on 
%     scatter(mti, mfi, 50, 'red', 'filled')
%     yticks([1:2:20])
%     yticklabels(round(frex(1:2:20)))
%     set(gca, 'YDir','normal')
%     xticks([1:20:141])
%     xticklabels(allDat{1}.leadLag.encTim([1:20:141]))
%     title([aggTargs(reg1).ROI ' to ' aggTargs(reg2).ROI ' interaction leave one out sig count'])
%     colorbar
% 
%     subplot(5,2,10)
%     hold off
%     curd = [alld,alld]; 
%     curd(hmSort==0) = []; 
%     plotVals = squeeze(combo(:, mti, mfi));
%     plotVals = 10.^plotVals; 
%     plotVals = plotVals + adjFact - .01; 
% 
% 
%     scatter(curd+.02, plotVals(hmSort), 'filled')
%     hold on 
%     curd = [alld,alld]; 
%     curd(hmSort==1) = []; 
%     scatter(curd-.02, plotVals(~hmSort), 'filled')
%     xlim([min(sub_d)-.2, max(sub_d)+.2])
%     yline(0, 'linewidth', 3, 'linestyle', '--', 'color', [169,169,169]./256)
%    
%     title(['interaction at: ' num2str(round(frex(mfi))) ' Hz  X ' ...
%         num2str(allDat{1}.leadLag.encTim(mti)) ' ms; p=' num2str(round(p(mti,mfi),2)) ])
% 
%     export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\ISPC_regional\' aggTargs(reg1).ROI '_' aggTargs(reg2).ROI '.jpg'])
%         end
%         catch 
%             'whoops'
%         end
% 
% 
% 
% %         sigConSub(reg1, reg2, 1, :, :) = tVals; 
% %         sigConSub(reg1, reg2, 2, :, :) = squeeze(mean(hitVals,1)) - squeeze(mean(missVals,1)); 
% %         sigConSub(reg1, reg2, 3, :, :) = p; 
% %         disp(['........................' num2str(round(toc/60, 1))])
% end %DEBUG 
%     end
% end
% 
% 
% % aovDat(aovDat.clu==0, :) = []; 
% % 
% % writetable(aovDat, 'R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject\PPC_regConnectionDat.csv')
% 
% 
% 
% 
% 
% 
% 
% end
% 
% 






%% look at cross correlation between TF points across all connections
% allCors = zeros(prod(size(allReg,[2,3]))); 
% 
% for ii = 1:size(allReg, 2)
%     ii
%    
%     for jj = 1:size(allReg, 3)
%         cur = allReg(:, ii, jj); 
%         cur = reshape(repmat(cur, size(allReg, [2,3])), size(allReg));
%         cur = squeeze(myArrayCorr(cur, allReg)); 
%         allCors((ii-1)*size(allReg,3)+jj, :) = cur(:); 
%     end
%   
% end
% 
% 
% seedSpread = zeros(size(allReg, [2,3])); 
% 
% for ii = 1:size(allReg,2)
%     for jj = 1:size(allReg,3)
%         cur = reshape(allCors((ii-1)*size(allReg,3)+jj, :), size(allReg, [2,3]));  
%    
%         seedSpread(ii,jj) = sum(cur>.3, 'all');
% 
% 
%     end
% end
% yScale = frex; 
% xScale = allDat{1}.leadLag.encTim; 
% figure
% imagesc(squeeze(mean(allReg,1))')
% yticks([1:2:length(yScale)])
% yticklabels(round(yScale(1:2:length(yScale))))
% xticks([1:20:length(xScale)])
% xticklabels(xScale([1:20:length(xScale)]))
% set(gca, 'YDir','normal')
% colorbar
% figure 
% imagesc(squeeze(var(allReg,1))')
% yticks([1:2:length(yScale)])
% yticklabels(round(yScale(1:2:length(yScale))))
% xticks([1:20:length(xScale)])
% xticklabels(xScale([1:20:length(xScale)]))
% set(gca, 'YDir','normal')
% colorbar

% goal here is to group the clusters of high threshold connections
% together! 

% %select the high 
% test = squeeze(mean(allReg,1)); 
% test = test(:); 
% p = squeeze(mean(allReg,1)) > prctile(test, 95); 
% 
% pidxtt = repmat([1:size(allReg,2)], [size(allReg,3),1])';
% pidxfi = repmat([1:size(allReg,3)], [size(allReg,2),1]);
% pidxtt = pidxtt(:); 
% pidxfi = pidxfi(:); 
% p = p(:);
% test = find(p);
% pidxtt = pidxtt(test);
% pidxfi = pidxfi(test); 
% clusSim = zeros(length(test)); 
% for ii = 1:length(test)
%     for jj = 1:length(test)
%         cur = reshape(allCors((pidxtt(ii)-1)*size(allReg,3)+pidxfi(ii),:), size(allReg, [2,3])); 
%         clusSim(ii,jj) = cur(pidxtt(jj), pidxfi(jj));
%     end
% end
% 
% cluVals = DBscanDynamicEpi(clusSim, 20, 3, 1, 1);
% 
% 
% p = reshape(p, size(allReg,[2,3]));
% cluMap = zeros(size(allReg, [2,3]) ); 
% for ii = 1:length(pidxtt)
%     cluMap(pidxtt(ii), pidxfi(ii)) = cluVals(ii); 
% end

% cluMap = load('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject\PPC_TF_mask.mat').cluMap;
% 

% TFconnection = VideoWriter(join(['G:\My Drive\Johnson\MTL_PFC_networkFigs\leadLagAll\' 'TFconnection.avi'],''));
% open(TFconnection); 
% f = figure('visible', false);
% f.Position = [100 100 1000 600];
% for ii = 1:size(allReg,2)
%     disp(ii)
%     for jj = 1:size(allReg,3)
%         cur = reshape(allCors((ii-1)*size(allReg,3)+jj, :), size(allReg, [2,3])); 
%         imagesc(cur')
%         hold on 
%         scatter(ii, jj, 30, 'k', 'filled')
% 
%         yticks([1:2:length(yScale)])
%         yticklabels(round(yScale(1:2:length(yScale))))
%         xticks([1:20:length(xScale)])
%         xticklabels(xScale([1:20:length(xScale)]))
%         set(gca, 'YDir','normal')
%         colorbar
%         caxis([.2, .4])
%         frame = getframe(gcf);
%         writeVideo(TFconnection, frame); 
%     end
% end
% 
% close(TFconnection)

