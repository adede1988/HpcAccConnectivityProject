function [PPCout, ii, latDat] = updatePPCvar(PPCout, PPC, ii, latDat, ...
    inLat1, inLat2, inTF1, inTF2, RT)

n = size(PPC, 3); 
PPCout(ii:ii+n-1, :, :, 1) = permute(PPC, [3,1,2]); 
PPCout(ii:ii+n-1, :, :, 2) = permute(inTF1, [2, 1, 3]);
PPCout(ii:ii+n-1, :, :, 3) = permute(inTF2, [2,1,3]); 

latDat(ii:ii+n-1,1) = inLat1; 
latDat(ii:ii+n-1,2) = inLat2; 
latDat(ii:ii+n-1,3) = RT; 





ii = ii+n; 

end