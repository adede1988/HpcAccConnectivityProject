function [latency] = gausLat(trials, tim, RT, statChoice, pad)

%trials: timepoints X trials
%tim: timepoints X 1 vector of time stamps
%RT: trials X 1 vector of response time values
%statChoice: allows you to toggle between getting latency v. amplitude
%pad: amount of padding for gaussian window (suggested = 25ms)

%NOTE: this will search for a peak between tim == 0 and tim == RT on each
%trial. If you want a peak between prompt offset and RT, then make a time
%vector with prompt offset at tim==0. If RT is at the same time in each
%trial, then make your RT vector just be that same time repeated for every
%trial. 

%statChoice 1 = get the latency
%statChoice 2 = get the peak value



latency = zeros(size(trials,2), 1);

for tt = 1:length(latency)
    if RT(tt)<max(tim)-500
         trial = trials(:,tt); 
                   
        test = gausswin(pad*2+1);
        test = test ./ sum(test); %normalize to sum to 1 
        padTrial = [zeros(pad,1); trial; zeros(pad,1)];
        smoothT = conv(padTrial, test, 'valid'); 
        [maxVal, idx] = max(smoothT(find(tim==0):find(tim>=RT(tt),1))); 
        testTim = tim(find(tim==0):find(tim>=RT(tt),1));
        if statChoice == 1
            latency(tt) =  testTim(idx);
        else
            idx = idx+find(tim==0)-1; 
            latency(tt) = max(trial(idx:idx)); 
        end
    else
        latency(tt) = -1; 
    end

end

end