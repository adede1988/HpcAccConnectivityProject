function [] = makeHFBTimeScatter(dat, perm2, fi, type)

[~, timMax2] = max(abs(mean(perm2.hitVals(:,:,fi),[1,3])));

hold off
chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], dat.hitSub, ...
         num2cell(dat.hitChi), 'UniformOutput',false);
chanIDuni = unique(chanID); 
chanIDuni(perm2.eliminate) = []; 
chanRT = cellfun(@(x) mean(dat.hitLat(ismember(chanID,x))), ...
    chanIDuni); 
                    %  mean(dat.hitRT(ismember(chanID, x))), ...
    
scatter(chanRT, mean(perm2.hitVals(:,timMax2,fi),3))
hold on 
chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], dat.missSub, ...
         num2cell(dat.missChi), 'UniformOutput',false);
chanIDuni = unique(chanID);  
chanRT = cellfun(@(x) mean(dat.missLat(ismember(chanID,x))), ...
    chanIDuni); 
scatter(chanRT, mean(perm2.missVals(:,timMax2,1:23),3))

xlabel('Image time')
ylabel(['HFB locked ' type])

end