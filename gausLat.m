function [latency] = gausLat(trials, tim, RT, statChoice)
%statChoice 1 = get the latency
%statChoice 2 = get the peak value



latency = zeros(size(trials,2), 1);

for tt = 1:length(latency)
    if RT(tt)<max(tim)-500
         trial = trials(:,tt); 
                   
        test = gausswin(11);
        test = test ./ sum(test); 
        trial = [zeros(5,1); trial; zeros(5,1)];
        smoothT = conv(trial, test, 'valid'); 
        [maxVal, idx] = max(smoothT(find(tim==0):find(tim>=RT(tt),1))); 
        testTim = tim(find(tim==0):find(tim>=RT(tt),1));
        if statChoice == 1
            latency(tt) =  testTim(idx);
        else
            latency(tt) = maxVal; 
        end
    else
        latency(tt) = -1; 
    end

end

end