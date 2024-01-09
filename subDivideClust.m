function clustIDX = subDivideClust(evecCorMat, clustIDX)




[silhouettes, clusti] = getSil(evecCorMat, clustIDX);

%use silhouette values to see if clusters are uniform
%uniform defined as similar silhouette value across the whole cluster
%whether that's good or bad
cIDs = unique(clustIDX); 
silDifs = zeros(length(cIDs),1); 
silMeans = zeros(length(cIDs),1); 
for ci = 1:length(cIDs)
    silMeans(ci) = mean(silhouettes(clusti(:,2)==cIDs(ci)));
    maxSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 90) ; 
    minSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 10) ; 
    silDifs(ci) = maxSil - minSil;
end

for ci = 1:length(cIDs)
    if silDifs(ci)>.3 || silMeans(ci)<.2 %bad uniformity, recluster
        clustCor = evecCorMat(clustIDX == cIDs(ci), ...
                              clustIDX == cIDs(ci)); 
        try
            newVals = DBscanDynamicEpi(clustCor, 3,3,1,0); 
            ogIDX = find(clustIDX == cIDs(ci)); 
            clustIDX(ogIDX) = newVals + 100*ci; 
            clustIDX(clustIDX==(-1+100*ci)) = -1; 
        catch
            clustIDX(clustIDX == cIDs(ci)) = -1; 
        end
    end
end

%fix the values! 
curIDout = 1;
curIDin = clustIDX(1);
for ii = 1:length(clustIDX)
    if clustIDX(ii) ~= -1
        
        if clustIDX(ii) == curIDin; 
            clustIDX(ii) = curIDout; 
        else
            curIDout = curIDout + 1;
            curIDin = clustIDX(ii); 
            clustIDX(ii) = curIDout; 
    
        end
    end
end











end