function [offVals] = getOffset(curChan)

offVals = zeros(size(curChan,2), 1); 
LL = -150:150; 
for tt = 1:size(curChan,2)
    cur = curChan(:,tt); 
    cur = cur - min(cur); 
    cur = cur ./ max(cur); 
    offVals(tt) = LL(round(wmean(1:length(LL), cur')));


end





end