function [outCluStats] = conditionCluTest(allObs, tim)
    


    chanNum = size(allObs,1);
    LLtim = -150:150; 
    %out stats will be chan X chan X pos/neg X clust X stat: 
    %stat 1: num points
    %stat 2: mean time
    %stat 3: min time
    %stat 4: max time
    %stat 5: mean LL
    %stat 6: min LL
    %stat 7: max LL

    %Initialize with 30 possible clusters, all valued 999


    outCluStats = nan(chanNum, chanNum, 2, 30, 7); 
    for ii = 1:chanNum %check connection in each pair independently
        for jj = 1:chanNum
            obs = squeeze(allObs(ii,jj,:,:)); 
            perms = zeros([size(obs),1000]); 

            for p = 1:1000
                timShuf = randsample(1:size(obs,1), size(obs,1), false); 
                llShuf = randsample(1:size(obs,2), size(obs,2), false); 
                perms(:,:,p) =   obs(timShuf, llShuf); 
            end

            [h, p, clusterinfo] = cluster_test(obs, perms); 
                    
            %get the max LL X time position of positive clusters
            if min([clusterinfo.pos_clusters.p])<.05
                pidx = find([clusterinfo.pos_clusters.p]<.05);
                for cc = 1:length(pidx)

                    %stat 1: num points
                    outCluStats(ii,jj,1,cc,1) = sum(clusterinfo.pos_clusters(pidx(cc)).inds, 'all'); 
                    %stat 2: mean time
                    cluTim = tim(sum(clusterinfo.pos_clusters(pidx(cc)).inds,2)>0);
                    outCluStats(ii,jj,1,cc,2) = mean(cluTim); 
                    %stat 3: min time
                    outCluStats(ii,jj,1,cc,3) = min(cluTim); 
                    %stat 4: max time
                    outCluStats(ii,jj,1,cc,4) = max(cluTim); 
                    %stat 5: mean LL
                    cluLL = LLtim(sum(clusterinfo.pos_clusters(pidx(cc)).inds,1)>0);
                    outCluStats(ii,jj,1,cc,5) = mean(cluLL); 
                    %stat 6: min LL
                    outCluStats(ii,jj,1,cc,6) = min(cluLL); 
                    %stat 7: max LL
                    outCluStats(ii,jj,1,cc,7) = max(cluLL); 
                    

                end
            end


            %get the max LL X time position of negative clusters
            if min([clusterinfo.neg_clusters.p])<.05
                pidx = find([clusterinfo.neg_clusters.p]<.05);
                for cc = 1:length(pidx)

                    %stat 1: num points
                    outCluStats(ii,jj,2,cc,1) = sum(clusterinfo.neg_clusters(pidx(cc)).inds, 'all'); 
                    %stat 2: mean time
                    cluTim = tim(sum(clusterinfo.neg_clusters(pidx(cc)).inds,2)>0);
                    outCluStats(ii,jj,2,cc,2) = mean(cluTim); 
                    %stat 3: min time
                    outCluStats(ii,jj,2,cc,3) = min(cluTim); 
                    %stat 4: max time
                    outCluStats(ii,jj,2,cc,4) = max(cluTim); 
                    %stat 5: mean LL
                    cluLL = LLtim(sum(clusterinfo.neg_clusters(pidx(cc)).inds,1)>0);
                    outCluStats(ii,jj,2,cc,5) = mean(cluLL); 
                    %stat 6: min LL
                    outCluStats(ii,jj,2,cc,6) = min(cluLL); 
                    %stat 7: max LL
                    outCluStats(ii,jj,2,cc,7) = max(cluLL); 
                    

                end
            end


                
        end
    end

     











end