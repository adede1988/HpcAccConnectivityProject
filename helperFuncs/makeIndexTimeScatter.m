function [] = makeIndexTimeScatter(dat, perm, perm2, fi)


adjNegs = @(x,y) x - min([min(x), min(y)]);
getScaledIndex = @(x,y) (adjNegs(x,y) - adjNegs(y,x)) ./...
    (adjNegs(x,y) + adjNegs(y,x));

[~, timMax] = max(abs(mean(perm.hitVals(:,:,fi),[1,3])));
[~, timMax2] = max(abs(mean(perm2.hitVals(:,:,fi),[1,3])));

hitIndex = getScaledIndex(mean(perm.hitVals(:,timMax,fi),3), ...
    mean(perm2.hitVals(:,timMax2,fi),3));
missIndex = getScaledIndex(mean(perm.missVals(:,timMax,fi),3), ...
    mean(perm2.missVals(:,timMax2,fi),3));


hold off
chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], dat.hitSub, ...
         num2cell(dat.hitChi), 'UniformOutput',false);
chanIDuni = unique(chanID); 
chanIDuni(perm2.eliminate) = []; 
chanRT = cellfun(@(x) mean(dat.hitLat(ismember(chanID,x))), ...
    chanIDuni); 
                    %  mean(dat.hitRT(ismember(chanID, x))), ...
    
scatter(chanRT, hitIndex)
hold on 
chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], dat.missSub, ...
         num2cell(dat.missChi), 'UniformOutput',false);
chanIDuni = unique(chanID);  
chanRT = cellfun(@(x) mean(dat.missLat(ismember(chanID,x))), ...
    chanIDuni); 
scatter(chanRT, missIndex)
ylim([-1,1])

xlabel('Image time')
ylabel('Image - HFB index')




