




enctim = [-1000:3500];
retOtim = [-1000:3000];


    leadLagEncTim = enctim(501:25:end-500);
    leadLagRetTim = retOtim(501:25:end-500); 


    tmp = pairDat.subHitLat2;
    tmp2 = pairDat.LL.subHit_low(tmp>-1,:,:); 
    tmp(tmp==-1) = []; 
    tmp(tmp>=max(pairDat.encTim)-500) = max(pairDat.encTim)-500;
    test = arrayfun(@(x) squeeze(tmp2(x,:, ...
        find(leadLagRetTim>=tmp(x),1)-20 : ...
        find(leadLagRetTim>=tmp(x),1)+20)), ...
        1:length(tmp), 'uniformoutput', false);
    test2 = reshape(cell2mat(test), [61, 41, length(test)]); 

