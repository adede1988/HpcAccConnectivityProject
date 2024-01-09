function [] = makeHFBImageScatter(perm, perm2, fi, type)

[~, timMax] = max(abs(mean(perm.hitVals(:,:,fi),[1,3])));
[~, timMax2] = max(abs(mean(perm2.hitVals(:,:,fi),[1,3])));

hold off
scatter(mean(perm.hitVals(:,timMax,fi),3), ...
    mean(perm2.hitVals(:,timMax2,fi),3) )
hold on 
scatter(mean(perm.missVals(:,timMax,fi),3), ...
    mean(perm2.missVals(:,timMax2,fi),3) )
xlabel(['image locked ' type])
ylabel(['HFB locked ' type])

maxVal = max([mean(perm.hitVals(:,timMax,fi),3); ...
    mean(perm2.hitVals(:,timMax2,fi),3)]);
minVal = min([mean(perm.hitVals(:,timMax,fi),3); ...
    mean(perm2.hitVals(:,timMax2,fi),3)]);
plot([minVal, maxVal], [minVal, maxVal], ...
    '--', 'color', 'k', 'linewidth', 2)

end