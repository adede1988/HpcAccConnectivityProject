function [hitVals, missVals] = getSingleTrialHitMiss(HFBdat, stat)


chanID = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.hitSub, ...
         num2cell(HFBdat.hitChi), 'UniformOutput',false);

chanID2 = cellfun(@(x, z) [x,'_', '_', num2str(z)], HFBdat.missSub,...
     num2cell(HFBdat.missChi), 'UniformOutput',false);

subID = [HFBdat.hitSub; HFBdat.missSub]; 
sub_uni = unique(subID); 
chanUni = unique(chanID); 


if stat==1 %image locked

    hitVals = HFBdat.hits'; 
    missVals = HFBdat.misses';
    tim = HFBdat.tim; 
  
    hitVals(:, tim<-450 | tim>3000) = []; 
    missVals(:, tim<-450 | tim>3000) = []; 
    tim(tim<-450 | tim>3000) = []; 


    %average across trials within channels
    hitTrials = hitVals; 
    missTrials = missVals; 
    hitVals = zeros(length(chanUni), length(tim)); 
    missVals = zeros(length(chanUni), length(tim)); 
    sub_out = cell(size(chanUni)); 

    eliminate = zeros(2,length(chanUni)); 

    for chan = 1:length(chanUni)

        hitidx = ismember(chanID, chanUni(chan));
        missidx = ismember(chanID2, chanUni(chan)); 
        
        hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
        missVals(chan, :) = mean(missTrials(missidx, :), 1); 
        tmp = subID(hitidx); 
        sub_out{chan} = tmp{1}; 

        if sum(hitidx)==0 
            eliminate(1,chan) = 1;
        end
        if sum(missidx)==0
            eliminate(2,chan) = 1;
        end
    end
    hitVals(eliminate(1,:)==1,:) = []; 
    missVals(eliminate(2,:)==1,:) = []; 
       

else %HFB locked
    hitVals = HFBdat.hits'; 
    missVals = HFBdat.misses';
    tim = HFBdat.tim; 

    w = 20; 

    hitLatIdx = arrayfun(@(x) find(x==tim), HFBdat.hitLat); 
    hitLatIdx(hitLatIdx>length(tim)-w) = length(tim) - w; 
    hitLatIdx(hitLatIdx<=w) =  w+1; 

    peakHit = arrayfun(@(x) hitVals(x,hitLatIdx(x)-w:...
                                     hitLatIdx(x)+w),...
                                     [1:length(hitLatIdx)], ...
                                     'uniformoutput', false);
    peakHit = cell2mat(peakHit.');

    missLatIdx = arrayfun(@(x) find(x==tim), HFBdat.missLat); 
    missLatIdx(missLatIdx>length(tim)-w) = length(tim) - w; 
    missLatIdx(missLatIdx<=w) = w+1; 

    peakMiss = arrayfun(@(x) missVals(x,missLatIdx(x)-w:...
                                     missLatIdx(x)+w), ...
                                     [1:length(missLatIdx)], ...
                                     'uniformoutput', false);
    peakMiss = cell2mat(peakMiss.');





    %average across trials within channels
    hitTrials = peakHit; 
    missTrials = peakMiss; 
    hitVals = zeros(length(chanUni), 41); 
    missVals = zeros(length(chanUni), 41); 
    sub_out = cell(size(chanUni)); 


    eliminate = zeros(2,length(chanUni)); 

    for chan = 1:length(chanUni)

        hitidx = ismember(chanID, chanUni(chan));
        missidx = ismember(chanID2, chanUni(chan)); 
        
        hitVals(chan, :) = mean(hitTrials(hitidx, :), 1); 
        missVals(chan, :) = mean(missTrials(missidx, :), 1); 
        tmp = subID(hitidx); 
        sub_out{chan} = tmp{1}; 

        if sum(hitidx)==0 
            eliminate(1,chan) = 1;
        end
        if sum(missidx)==0
            eliminate(2,chan) = 1;
        end
    end
    hitVals(eliminate(1,:)==1,:) = []; 
    missVals(eliminate(2,:)==1,:) = []; 



end


end