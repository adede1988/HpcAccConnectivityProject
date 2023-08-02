function [outCluStats] = conditionCluTest(HFB1, pow2, missidx, hitidx, outCluStats, alltim, dstim, chan)
    



  %store the trials X offset X time info for each trial type
        missTemp = zeros(length(missidx), 301, length(dstim)); 
        hitTemp = zeros(length(hitidx), 301, length(dstim)); 
        
        tic
        for offSet = -150:150 %negative means current Channel leads, positive means other channel leads
            
            if offSet<0
                HFB2 = [pow2(abs(offSet)+1:end,:); zeros([abs(offSet), size(pow2,2)] ) ] ;
            elseif offSet>0
                HFB2 = [zeros([abs(offSet), size(pow2,2)] );  pow2(1:end-abs(offSet),:)];
            end
            %sub miss 
         

            missTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-50:x+50, missidx), HFB2(x-50:x+50, missidx)), ...
                                               dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(missidx), length(dstim)] );


            hitTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-50:x+50, hitidx), HFB2(x-50:x+50, hitidx)), ...
                                               dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(hitidx), length(dstim)] );

          

%             for tt = 1:length(missidx)
%                 missTemp(tt, offSet+151, :) = arrayfun(@(x) corr(HFB1(x-50:x+50, missidx(tt)), HFB2(x-50:x+50, missidx(tt))), ...
%                                    dstim+abs(min(alltim))+1 );
%             end
%             %sub hit
%             for tt = 1:length(hitidx)
%                 hitTemp(tt, offSet+151, :) = arrayfun(@(x) corr(HFB1(x-50:x+50, hitidx(tt)), HFB2(x-50:x+50, hitidx(tt))), ...
%                                dstim+abs(min(alltim))+1 );
%             end

    
        end
        toc
        disp('................leadLag calcualted')
        permCount = 1000; 
        obsT = zeros([ size(missTemp, [2,3]),permCount]); 
%         test = reshape(cell2mat(arrayfun(@(x) arrayfun(@(y) ...
%                 corrdiff(mean(hitTemp(:, y, x)), mean(missTemp(:, y, x)), 101, 101), ...
%                 1:301),...
%                 1:size(missTemp,3), 'uniformoutput', false)), size(missTemp, [2,3]));
        nullT = zeros([ size(missTemp, [2,3]),permCount]);
        kVal = min([length(missidx), length(hitidx)]) - 1; 
        for p = 1:permCount
            if mod(p, 100) ==0
                disp(['................permutation: ' num2str(p)])
            end
            % subsample to the lower trial count-1
            h = 1:length(hitidx); 
            m = 1:length(missidx); 
            h = randsample(h, kVal, true); 
            m = randsample(m, kVal, true); 
            obsT(:,:,p) = myArrayT(hitTemp(h, :, :), missTemp(m, :, :),  1);



%             allidx = [h,m+1000];
%             h1 = randsample(1:length(allidx), length(hitidx), false); 
%             m1 = [1:length(allidx)];
%             m1(h1) = []; 
%             %OG miss idx -> new miss or hit idx
%             m2 = allidx(m1); 
%             h2 = allidx(h1); 
% 
%             %misses drawn from hits
%             m1 = m2(m2<1000); 
%             %hits drawn from hits
%             h1 = h2(h2<1000); 
%             %misses drawn from misses
%             m2 = m2(m2>1000) - 1000; 
%             %hits drawn from misses
%             h2 = h2(h2>1000) - 1000; 

            nullT(:,:,p) = myArrayT(hitTemp(h, :, :), missTemp(m, :, :),  2);


%              reshape(cell2mat(arrayfun(@(x) arrayfun(@(y) ...
%                 corrdiff(mean([squeeze(hitTemp(h1, y, x)); squeeze(missTemp(h2, y, x))]),...
%                          mean([squeeze(hitTemp(m1, y, x)); squeeze(missTemp(m2, y, x))]),...
%                          101, 101), 1:301),...
%                 1:size(missTemp,3), 'uniformoutput', false)), size(missTemp, [2,3]));
            

        end
    
        test2 = squeeze(mean(obsT, 3)); 
        
        [h, p, clusterinfo] = cluster_test(test2, nullT); 

        LLtim = -150:150;  
        %get the max LL X time position of positive clusters
        if min([clusterinfo.pos_clusters.p])<.15
            pidx = find([clusterinfo.pos_clusters.p]<.15);
            for cc = 1:length(pidx)

                %stat 1: num points
                outCluStats(chan,1,cc,1) = sum(clusterinfo.pos_clusters(pidx(cc)).inds, 'all'); 
                %stat 2: mean time
                cluTim = dstim(sum(clusterinfo.pos_clusters(pidx(cc)).inds,1)>0);
                outCluStats(chan,1,cc,2) = mean(cluTim); 
                %stat 3: min time
                outCluStats(chan,1,cc,3) = min(cluTim); 
                %stat 4: max time
                outCluStats(chan,1,cc,4) = max(cluTim); 
                %stat 5: mean LL
                cluLL = LLtim(sum(clusterinfo.pos_clusters(pidx(cc)).inds,2)>0);
                outCluStats(chan,1,cc,5) = mean(cluLL); 
                %stat 6: min LL
                outCluStats(chan,1,cc,6) = min(cluLL); 
                %stat 7: max LL
                outCluStats(chan,1,cc,7) = max(cluLL); 
                %stat 8: mean correlation Hit
                outCluStats(chan,1,cc,8) = mean(hitTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds), 'all');
                %stat 9: min correlation Hit
                outCluStats(chan,1,cc,9) = min(hitTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds),[], 'all');
                %stat 10: max correlation Hit
                outCluStats(chan,1,cc,10) = max(hitTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds),[], 'all');
                %stat 11: median correlation Hit
                outCluStats(chan,1,cc,11) = median(hitTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds), 'all');
                %stat 12: mean correlation Miss
                outCluStats(chan,1,cc,12) = mean(missTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds), 'all');
                %stat 13: min correlation Miss
                outCluStats(chan,1,cc,13) = min(missTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds),[], 'all');
                %stat 14: max correlation Miss
                outCluStats(chan,1,cc,14) = max(missTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds),[], 'all');
                %stat 15: median correlation Miss
                outCluStats(chan,1,cc,15) = median(missTemp(:,clusterinfo.pos_clusters(pidx(cc)).inds), 'all');

            end
        end


        %get the max LL X time position of negative clusters
        if min([clusterinfo.neg_clusters.p])<.15
            pidx = find([clusterinfo.neg_clusters.p]<.15);
            for cc = 1:length(pidx)

                %stat 1: num points
                outCluStats(chan,2,cc,1) = sum(clusterinfo.neg_clusters(pidx(cc)).inds, 'all'); 
                %stat 2: mean time
                cluTim = dstim(sum(clusterinfo.neg_clusters(pidx(cc)).inds,1)>0);
                outCluStats(chan,2,cc,2) = mean(cluTim); 
                %stat 3: min time
                outCluStats(chan,2,cc,3) = min(cluTim); 
                %stat 4: max time
                outCluStats(chan,2,cc,4) = max(cluTim); 
                %stat 5: mean LL
                cluLL = LLtim(sum(clusterinfo.neg_clusters(pidx(cc)).inds,2)>0);
                outCluStats(chan,2,cc,5) = mean(cluLL); 
                %stat 6: min LL
                outCluStats(chan,2,cc,6) = min(cluLL); 
                %stat 7: max LL
                outCluStats(chan,2,cc,7) = max(cluLL); 
                %stat 8: mean correlation Hit
                outCluStats(chan,2,cc,8) = mean(hitTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds), 'all');
                %stat 9: min correlation Hit
                outCluStats(chan,2,cc,9) = min(hitTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds),[], 'all');
                %stat 10: max correlation Hit
                outCluStats(chan,2,cc,10) = max(hitTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds),[], 'all');
                %stat 11: median correlation Hit
                outCluStats(chan,2,cc,11) = median(hitTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds), 'all');
                %stat 12: mean correlation Miss
                outCluStats(chan,2,cc,12) = mean(missTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds), 'all');
                %stat 13: min correlation Miss
                outCluStats(chan,2,cc,13) = min(missTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds),[], 'all');
                %stat 14: max correlation Miss
                outCluStats(chan,2,cc,14) = max(missTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds),[], 'all');
                %stat 15: median correlation Miss
                outCluStats(chan,2,cc,15) = median(missTemp(:,clusterinfo.neg_clusters(pidx(cc)).inds), 'all');
                

            end
        end




end













% 
% 
% 
% 
% 
%     chanNum = size(allObs,1);
%     LLtim = -150:150; 
%     %out stats will be chan X chan X pos/neg X clust X stat: 
%     %stat 1: num points
%     %stat 2: mean time
%     %stat 3: min time
%     %stat 4: max time
%     %stat 5: mean LL
%     %stat 6: min LL
%     %stat 7: max LL
% 
%     %Initialize with 30 possible clusters, all valued 999
% 
% 
%     outCluStats = nan(chanNum, chanNum, 2, 30, 7); 
%     for ii = 1:chanNum %check connection in each pair independently
%         for jj = 1:chanNum
%             obs = squeeze(allObs(ii,jj,:,:)); 
%             perms = zeros([size(obs),1000]); 
% 
%             for p = 1:1000
%                 timShuf = randsample(1:size(obs,1), size(obs,1), false); 
%                 llShuf = randsample(1:size(obs,2), size(obs,2), false); 
%                 perms(:,:,p) =   obs(timShuf, llShuf); 
%             end
% 
%             [h, p, clusterinfo] = cluster_test(obs, perms); 
%                     
%             %get the max LL X time position of positive clusters
%             if min([clusterinfo.pos_clusters.p])<.05
%                 pidx = find([clusterinfo.pos_clusters.p]<.05);
%                 for cc = 1:length(pidx)
% 
%                     %stat 1: num points
%                     outCluStats(ii,jj,1,cc,1) = sum(clusterinfo.pos_clusters(pidx(cc)).inds, 'all'); 
%                     %stat 2: mean time
%                     cluTim = tim(sum(clusterinfo.pos_clusters(pidx(cc)).inds,2)>0);
%                     outCluStats(ii,jj,1,cc,2) = mean(cluTim); 
%                     %stat 3: min time
%                     outCluStats(ii,jj,1,cc,3) = min(cluTim); 
%                     %stat 4: max time
%                     outCluStats(ii,jj,1,cc,4) = max(cluTim); 
%                     %stat 5: mean LL
%                     cluLL = LLtim(sum(clusterinfo.pos_clusters(pidx(cc)).inds,1)>0);
%                     outCluStats(ii,jj,1,cc,5) = mean(cluLL); 
%                     %stat 6: min LL
%                     outCluStats(ii,jj,1,cc,6) = min(cluLL); 
%                     %stat 7: max LL
%                     outCluStats(ii,jj,1,cc,7) = max(cluLL); 
%                     
% 
%                 end
%             end
% 
% 
%             %get the max LL X time position of negative clusters
%             if min([clusterinfo.neg_clusters.p])<.05
%                 pidx = find([clusterinfo.neg_clusters.p]<.05);
%                 for cc = 1:length(pidx)
% 
%                     %stat 1: num points
%                     outCluStats(ii,jj,2,cc,1) = sum(clusterinfo.neg_clusters(pidx(cc)).inds, 'all'); 
%                     %stat 2: mean time
%                     cluTim = tim(sum(clusterinfo.neg_clusters(pidx(cc)).inds,2)>0);
%                     outCluStats(ii,jj,2,cc,2) = mean(cluTim); 
%                     %stat 3: min time
%                     outCluStats(ii,jj,2,cc,3) = min(cluTim); 
%                     %stat 4: max time
%                     outCluStats(ii,jj,2,cc,4) = max(cluTim); 
%                     %stat 5: mean LL
%                     cluLL = LLtim(sum(clusterinfo.neg_clusters(pidx(cc)).inds,1)>0);
%                     outCluStats(ii,jj,2,cc,5) = mean(cluLL); 
%                     %stat 6: min LL
%                     outCluStats(ii,jj,2,cc,6) = min(cluLL); 
%                     %stat 7: max LL
%                     outCluStats(ii,jj,2,cc,7) = max(cluLL); 
%                     
% 
%                 end
%             end
% 
% 
%                 
%         end
%     end
% 
%      