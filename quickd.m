function [d] = quickd(subDat)
T = sum(subDat.retInfo(:,1)==1 | subDat.retInfo(:,1)==2); 
Hr = sum(subDat.retInfo(:,1)==1) / T; 
T = sum(subDat.retInfo(:,1)==3 | subDat.retInfo(:,1)==4); 
F = sum(subDat.retInfo(:,1)==4);
if F == 0
    F = 1; 
end
Fr = F / T; 

d = norminv(Hr) - norminv(Fr); 



end