function [HFBlockedDat] = getHFBlockedDat(dat, tim, HFBlats, pad)
    %dat: time x trials matrix of input data 
    %tim: time x 1 vector of timestamps that applies to all trials
    %HFBlats: trials x 1 vector of HFB peak latencies
    %pad: scalar how much padding do you want around the HFB latency in
    %terms of time points

    %NOTE: whatever tim vector was used with gausLat to get the HFBlats
    %values should be used here again


    %HFBlockedDat: time X trials matrix of HFB locked data


    HFBidx = arrayfun(@(x) find(tim >= x, 1), HFBlats);
    %hard adjust any HFB latencies that will go outside the range of tim
    HFBidx(HFBidx<pad+1) = pad+1; 
    HFBidx(HFBidx>length(tim)-pad) = length(tim)- pad;
    %get the +- 20 points around the HFB peaks
    HFBlockedDat = arrayfun(@(x, y) dat(x-20:x+20, y), HFBidx', ...
        1:length(HFBidx), 'UniformOutput',false );
    HFBlockedDat = cat(2, HFBlockedDat{:});



end