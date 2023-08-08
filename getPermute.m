function [p_bin] = getPermute(inMat)


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
phor = diff(p_bin);
pver = diff(p_bin'); 
phor(phor ~= 0) = 1; 
pver(pver ~= 0) = 1; 
phor = [zeros(size(phor,2),1)'; phor]; 
pver = [zeros(size(pver,2),1)'; pver]; 
p_bin = phor + pver'; 
p_bin(p_bin>1) = 1; 


end