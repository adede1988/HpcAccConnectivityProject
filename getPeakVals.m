function [peakVals] = getPeakVals(trials, lats, tim)

peakIdx = arrayfun(@(x) find(x==tim), lats , ...
    'UniformOutput', false);
tmp = cell(1,1); 
tmp{1} = -1; 
peakIdx(lats==-1) = tmp; 
peakIdx = cell2mat(peakIdx);

peakVals = zeros(size(peakIdx)); 
for ii = 1:length(peakVals)
    if peakIdx(ii) >-1
        peakVals(ii) = trials(peakIdx(ii), ii); 
    else
        peakVals(ii) = -1; 
    end

end




end