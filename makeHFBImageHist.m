function [outDat] = makeHFBImageHist(perm, perm2, fi)


adjNegs = @(x,y) x - min([min(x), min(y)]);
getScaledIndex = @(x,y) (adjNegs(x,y) - adjNegs(y,x)) ./...
    (adjNegs(x,y) + adjNegs(y,x));

[~, timMax] = max(abs(mean(perm.hitVals(:,:,fi),[1,3])));
[~, timMax2] = max(abs(mean(perm2.hitVals(:,:,fi),[1,3])));
x = mean(perm.hitVals(:,timMax,fi),3);
y = mean(perm2.hitVals(:,timMax2,fi),3);
% allMax = max([x;y]); 
% allMin = min([x;y]); 
% 
% 
% degree = 1;
% coefficients = polyfit([allMin, allMax], [allMin, allMax], degree);
% 
% % Generate a linear function using the coefficients
% fitted_y = polyval(coefficients, x);

hold off
hitIndex = getScaledIndex(x, y);
histogram(hitIndex, [-1:.1:1])
hold on 

x = mean(perm.missVals(:,timMax,fi),3);
y = mean(perm2.missVals(:,timMax2,fi),3);
% allMax = max([x;y]); 
% allMin = min([x;y]); 
% 
% 
% degree = 1;
% coefficients = polyfit([allMin, allMax], [allMin, allMax], degree);
% 
% % Generate a linear function using the coefficients
% fitted_y = polyval(coefficients, x);


missIndex = getScaledIndex(x, y);
histogram(missIndex, [-1:.1:1])
xlabel('HFB     Neither    Image')
title('image - HFB index')

outDat = {hitIndex, missIndex}; 




end