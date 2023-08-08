function [out] = leadLagFigure(inMat, tim, regName, titText)


perms = 1000; 

permVals = zeros([size(inMat), perms]);


for ii = 1:perms
    cur = inMat(:);
    cur = cur(randsample(1:length(cur), length(cur), false)); 
    permvals(:,:, ii) = reshape(cur, size(inMat) ); 

end


 [h, p, clusterinfo] = cluster_test(inMat, permvals); 

p_bin = zeros(size(p)); 
p_bin(p<.05) = 1; %1 = cluster!
out = p_bin; 
phor = diff(p_bin);
pver = diff(p_bin'); 
phor(phor ~= 0) = 1; 
pver(pver ~= 0) = 1; 
phor = [zeros(size(phor,2),1)'; phor]; 
pver = [zeros(size(pver,2),1)'; pver]; 
p_bin = phor + pver'; 
p_bin(p_bin>1) = 1; 

figure
% imagesc(allDat{1}.leadLag.encTim, -150:150,  MTL_all')
imagesc(inMat')
hold on 
for ii = 1:size(p_bin,2)
    cur = p_bin(:,ii); 
    pnt = find(cur==1);
    scatter(pnt, ii, 5, 'red', 'filled', 'square')
end
ylim([0,300])
yticks([0:50:300])
yticklabels([150,100,50,0,50,100,150])
ylabel([ regName ' lags                    ' regName ' leads'])
title(titText)
xlimVals = [0, size(inMat,1)]; 
xlim(xlimVals)
xticks(xlimVals(1):20:xlimVals(2))
xticklabels(tim(xlimVals(1)+1:20:xlimVals(2)))





end