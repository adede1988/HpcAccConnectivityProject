function [] = HFBleadLagPipeline(chanFiles, idx) 

curSub = load([chanFiles(idx).folder '/' chanFiles(idx).name]).curSub; % go raw if it's not working!

%initialize the output for all channels to all channels 
    %chan X chan X trial X time X offset
    leadLagHit = zeros(length(curSub), length(curSub), size(curSub(1).HFB.subHit,2), size(curSub(1).HFB.subHit,1), 41);
    leadLagMiss = zeros(length(curSub), length(curSub), size(curSub(1).HFB.subMiss,2), size(curSub(1).HFB.subMiss,1), 41);
    leadLagHit_on = zeros(length(curSub), length(curSub), size(curSub(1).HFB.hit_on,2), size(curSub(1).HFB.hit_on,1), 41);
    leadLagMiss_on = zeros(length(curSub), length(curSub), size(curSub(1).HFB.miss_on,2), size(curSub(1).HFB.miss_on,1), 41);
    
tic
    for chan = 1:length(curSub)
        for chan2 = 1:length(curSub)
            for offSet = -20:20 %a time step is 5ms, so -20:20 is -100ms:100ms
                %sub Hit trials 
                HFB1 = curSub(chan).HFB.subHit; 
                HFB2 = curSub(chan2).HFB.subHit; 
                if offSet<0
                    HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
                elseif offSet>0
                    HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
                end
                
                %hard code 30 time step (150ms +- around time point)
                test = reshape(cell2mat(arrayfun(@(x) ... %x is loop on trials
                    arrayfun(@(y) ... %y is loop on time
                        corr(HFB1(y-30:y+30,x), HFB2(y-30:y+30,x)), ...
                    31:size(HFB1,1)-30), ...
                1:size(HFB1,2), 'uniformoutput', false)),[], size(HFB2,2));
                
                leadLagHit(chan, chan2, :, 31:size(test,1)+30, offSet+21) = test'; 
%                 disp([num2str(chan) ' ' num2str(chan2) ' ' num2str(offSet) ' ' num2str(round(toc))])

                 %sub Miss trials 
                HFB1 = curSub(chan).HFB.subMiss; 
                HFB2 = curSub(chan2).HFB.subMiss; 
                if offSet<0
                    HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
                elseif offSet>0
                    HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
                end
                
                %hard code 30 time step (150ms +- around time point)
                test = reshape(cell2mat(arrayfun(@(x) ... %x is loop on trials
                    arrayfun(@(y) ... %y is loop on time
                        corr(HFB1(y-30:y+30,x), HFB2(y-30:y+30,x)), ...
                    31:size(HFB1,1)-30), ...
                1:size(HFB1,2), 'uniformoutput', false)),[], size(HFB2,2));
                
                leadLagMiss(chan, chan2, :, 31:size(test,1)+30, offSet+21) = test'; 

                %hit trials 
                HFB1 = curSub(chan).HFB.hit_on; 
                HFB2 = curSub(chan2).HFB.hit_on; 
                if offSet<0
                    HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
                elseif offSet>0
                    HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
                end
                
                %hard code 30 time step (150ms +- around time point)
                test = reshape(cell2mat(arrayfun(@(x) ... %x is loop on trials
                    arrayfun(@(y) ... %y is loop on time
                        corr(HFB1(y-30:y+30,x), HFB2(y-30:y+30,x)), ...
                    31:size(HFB1,1)-30), ...
                1:size(HFB1,2), 'uniformoutput', false)),[], size(HFB2,2));
                
                leadLagHit_on(chan, chan2, :, 31:size(test,1)+30, offSet+21) = test'; 

                %miss trials 
                HFB1 = curSub(chan).HFB.miss_on; 
                HFB2 = curSub(chan2).HFB.miss_on; 
                if offSet<0
                    HFB2 = [ HFB2(abs(offSet)+1:end, :); zeros([abs(offSet), size(HFB1,2)] )];
                elseif offSet>0
                    HFB2 = [zeros([abs(offSet), size(HFB1,2)] ); HFB2(1:end -abs(offSet), :)];
                end
                
                %hard code 30 time step (150ms +- around time point)
                test = reshape(cell2mat(arrayfun(@(x) ... %x is loop on trials
                    arrayfun(@(y) ... %y is loop on time
                        corr(HFB1(y-30:y+30,x), HFB2(y-30:y+30,x)), ...
                    31:size(HFB1,1)-30), ...
                1:size(HFB1,2), 'uniformoutput', false)),[], size(HFB2,2));
                
                leadLagMiss_on(chan, chan2, :, 31:size(test,1)+30, offSet+21) = test'; 


            end
        end
    end


    curSub(1).leadLagHit = leadLagHit; 
    curSub(1).leadLagMiss = leadLagMiss; 
    curSub(1).leadLagHit_on = leadLagHit_on; 
    curSub(1).leadLagMiss_on = leadLagMiss_on; 

    save([chanFiles(idx).folder '/' chanFiles(idx).name], 'curSub', '-v7.3')
    
end