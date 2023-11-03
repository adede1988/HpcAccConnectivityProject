function LL = getLL(HFB1, pow2, missidx, hitidx, alltim, dstim)

%inputs: 
%HFB1: time X trials (needs to be at 1000hz)
%pow2: time X trials (needs to be at 1000hz)
%missidx: vector of either 1s and 0s coding misses or the indecies
%           themselves
%hitidx: as for missidx
%alltim: the time points at the high sampling rate of the data
%dstim: down sample time, the time points at which measurement is wanted

    %store the trials X offset X time info for each trial type
        missTemp = zeros(length(missidx), 301, length(dstim)); 
        hitTemp = zeros(length(hitidx), 301, length(dstim)); 
        
     
        for offSet = -150:150 %negative means current Channel leads, positive means other channel leads
            
            if offSet<0
                HFB2 = [pow2(abs(offSet)+1:end,:); zeros([abs(offSet), size(pow2,2)] ) ] ;
            elseif offSet>0
                HFB2 = [zeros([abs(offSet), size(pow2,2)] );  pow2(1:end-abs(offSet),:)];
            end
            %sub miss 
         

            missTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-500:x+500, missidx), HFB2(x-500:x+500, missidx)), ...
                                               dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(missidx), length(dstim)] );


            hitTemp(:, offSet+151, :) = reshape(cell2mat( arrayfun(@(x) myArrayCorr(HFB1(x-500:x+500, hitidx), HFB2(x-500:x+500, hitidx)), ...
                                               dstim+abs(min(alltim))+1 , 'uniformoutput', false)), [ length(hitidx), length(dstim)] );

          
 

    
        end


test = squeeze(mean(missTemp));
test2 = squeeze(mean(hitTemp)); 

LL = zeros([2, size(test)]);

%store them as hit = 1; miss = 2
LL(1,:,:) = test2; 
LL(2,:,:) = test; 









end