function [perm, perm2] = fixMissing(dat, perm, perm2)    
    hitVals = dat.hits; 
    missVals = dat.misses;
    tim = dat.tim; 
  
    hitVals(tim<-450 | tim>3000, :, :) = []; 
    missVals(tim<-450 | tim>3000, :, :) = []; 
    tim(tim<-450 | tim>3000) = []; 
    

  

    chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], dat.hitSub, ...
         num2cell(dat.hitChi), 'UniformOutput',false);

    chanID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], dat.missSub, ...
         num2cell(dat.missChi), 'UniformOutput',false);

    subID = [dat.hitSub; dat.missSub]; 
    sub_uni = unique(subID); 
    chanUni = unique(chanID); 


     hitTrials = hitVals; 
    missTrials = missVals; 
    %chan X time X frequency
    hitVals = zeros([length(chanUni), ...
        size(hitTrials,1), size(hitTrials,3)]); 
    missVals = zeros([length(chanUni), ...
        size(hitTrials,1), size(hitTrials,3)]); 
    sub_out = cell(size(chanUni)); 


    eliminate = []; 

    for chan = 1:length(chanUni)

        hitidx = ismember(chanID, chanUni(chan));
        missidx = ismember(chanID2, chanUni(chan)); 
        
            hitVals(chan, :,:) = mean(hitTrials(:, hitidx,:), 2); 
            missVals(chan, :,:) = mean(missTrials(:,missidx, :), 2); 
            tmp = subID(hitidx); 
            sub_out{chan} = tmp{1}; 
        if sum(hitidx)>0 && sum(missidx)>0
        else
            eliminate = [eliminate, chan]; 
        end
    end

    missVals(eliminate,:,:) = [];
    hitVals(eliminate,:, :) = [];
    perm.missVals = missVals; 
    perm.hitVals = hitVals; 


    %do it for the HFB locked values: 
     hitVals = dat.hits; 
    missVals = dat.misses;
    tim = dat.tim; 
  
%         hitVals(tim<-450 | tim>3000, :, :) = []; 
%         missVals(tim<-450 | tim>3000, :, :) = []; 
%         tim(tim<-450 | tim>3000) = []; 

    w = 20; 

    hitLatIdx = arrayfun(@(x) find(x==tim), dat.hitLat); 
    hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
    hitLatIdx(hitLatIdx<=w) =  w+1; 

    peakHit = arrayfun(@(x) squeeze(hitVals(hitLatIdx(x)-w:...
                                     hitLatIdx(x)+w, x, :)),  [1:length(hitLatIdx)], ...
                                     'uniformoutput', false);
    peakHit = cat(3, peakHit{:});

    missLatIdx = arrayfun(@(x) find(x==tim), dat.missLat); 
    missLatIdx(missLatIdx>length(tim)-w) = length(tim) - w; 
    missLatIdx(missLatIdx<=w) = w+1; 

    peakMiss = arrayfun(@(x) squeeze(missVals(missLatIdx(x)-w:...
                                     missLatIdx(x)+w, x, :)), [1:length(missLatIdx)], ...
                                     'uniformoutput', false);
    peakMiss = cat(3, peakMiss{:});





    %average across trials within channels
    hitTrials = peakHit; 
    missTrials = peakMiss; 
    %chan X time X frequency
    hitVals = zeros([length(chanUni), 41, size(hitTrials,2)]); 
    missVals = zeros([length(chanUni), 41, size(hitTrials,2)]); 
    sub_out = cell(size(chanUni)); 


    eliminate = []; 

    for chan = 1:length(chanUni)

        hitidx = ismember(chanID, chanUni(chan));
        missidx = ismember(chanID2, chanUni(chan)); 
        
            hitVals(chan, :,:) = mean(hitTrials(:, :,hitidx), 3); 
            missVals(chan, :,:) = mean(missTrials(:, :,missidx), 3); 
            tmp = subID(hitidx); 
            sub_out{chan} = tmp{1}; 
        if sum(hitidx)>0 && sum(missidx)>0
        else
            eliminate = [eliminate, chan]; 
        end
    end

    missVals(eliminate, :,:) = []; 
    hitVals(eliminate, :,:) = [];
    perm2.missVals = missVals;
    perm2.hitVals = hitVals; 


end