function [sigFreqs] = GEDbandDiscovery(reg, headFiles, outStatFiles, regions, ...
        phase, outStatFilesPhase, sigFreqs)

%down select files to the current target region and phase of experiment
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, regions{reg}));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, phase));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'stat0'));
HFBperms= outStatFiles(test); 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'stat1'));
Imageperms= outStatFiles(test); 

%down select the phase permutation files
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, regions{reg}));
outStatFilesPhase(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, phase));
outStatFilesPhase(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'stat0'));
HFBpermsPhase= outStatFilesPhase(test); 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'stat1'));
ImagepermsPhase= outStatFilesPhase(test); 

%do the same for the head files, there will only be 1
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, regions{reg}));
headFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, phase));
headFiles(~test) = []; 


%load in the main data
dat = load([headFiles.folder '/' headFiles.name]).statInfo;
perm = load([Imageperms(1).folder '/' Imageperms(1).name]).outDat;

perm2 = load([HFBperms(1).folder '/' HFBperms(1).name]).outDat;

%permutations on TF power took the mean, 
% but I want channel level data back
[perm, perm2] = fixMissing(dat, perm, perm2);


%doing GED-based frequency band discovery Image locked: 
try
[evals,evecs,maps] = deal(zeros(size(perm.hitVals,[1,3])));

for fi = 1:100


hitVals = permute(perm.hitVals, [3,1,2]); 
S = squeeze(hitVals(fi,:,:)) * squeeze(hitVals(fi,:,:))' ./ ...
    prod(size(hitVals,[2,3]));

missVals = permute(perm.missVals, [3,1,2]); 
R = squeeze(missVals(fi,:,:)) * squeeze(missVals(fi,:,:))' ./ ...
    prod(size(missVals,[2,3]));

% regularized R
gamma = .01;
Rr = R*(1-gamma) + eye(size(R,1))*gamma*mean(eig(R));

% global variance normalize (optional; this scales the eigenspectrum)
S = S / (std(S(:))/std(R(:)));

% GED
[W,L] = eig(S,Rr);
[evals(:,fi),sidx] = sort(diag(L),'descend');
W = W(:,sidx);

% store top component map and eigenvector
maps(:,fi) = W(:,1)'*S;
evecs(:,fi) = W(:,1);

end


E = zscore(evecs,[],1)';
evecCorMat = (E*E'/(size(R,1)-1)).^2;

clustIDX = DBscanDynamicEpi(evecCorMat, 3, 3, 1, 0);
changeCheck = true; 
while changeCheck
    newClust = subDivideClust(evecCorMat, clustIDX);
    if sum(newClust - clustIDX) == 0
        changeCheck = false;
    else
        clustIDX = newClust; 
    end
end


[silhouettes, clusti] = getSil(evecCorMat, clustIDX);


cIDs = unique(clustIDX); 
silDifs = zeros(length(cIDs),1); 
silMeans = zeros(length(cIDs),1); 
for ci = 1:length(cIDs)
    silMeans(ci) = mean(silhouettes(clusti(:,2)==cIDs(ci)));
    maxSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 90) ; 
    minSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 10) ; 
    silDifs(ci) = maxSil - minSil;
end

breakVals = find(diff(clustIDX(clustIDX>0))~=0); 
outClusts = clustIDX(clustIDX>0); 
figure('visible', false)
subplot 221
outFrex = dat.frex(clustIDX>0); 
outCorMat = evecCorMat(clustIDX>0, clustIDX>0); 
imagesc(corr(outCorMat))
title('eigenvector correlations')
caxis([.05, .8])
hold on 
for bi = 1:length(breakVals)
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    yline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
yticks(breakVals+.5)
yticklabels(round(outFrex(breakVals)))

subplot(2,2,[3,4])
hold off
[silVals, clusti] = getSil(outCorMat, outClusts);
breakVals2 = [0, breakVals', length(silVals)];

for bi = 1:length(breakVals)

    plot(breakVals2(bi)+1: breakVals2(bi+1),...
        silVals(breakVals2(bi)+1: breakVals2(bi+1)), 'color', 'k')
    hold on 
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
title([regions{reg} ' ' phase ' freq band divisions'])
ylabel('silhouette values')

export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' 'GEDbandDiscovery' '_' regions{reg} '_' phase '.jpg'], '-r300')

outFrexIDX = arrayfun(@(x) find(outFrex(x) == dat.frex), [1:length(outFrex)]);

if strcmp(phase, 'sub')
    sigFreqs(1,outFrexIDX(breakVals),1,1,1) = 1; 
else

    sigFreqs(1,outFrexIDX(breakVals),2,1,1) = 1; 
end
catch
end


%% repeat the whole process for HFB

try

[evals,evecs,maps] = deal(zeros(size(perm2.hitVals,[1,3])));

for fi = 1:100


hitVals = permute(perm2.hitVals, [3,1,2]); 
S = squeeze(hitVals(fi,:,:)) * squeeze(hitVals(fi,:,:))' ./ ...
    prod(size(hitVals,[2,3]));

missVals = permute(perm2.missVals, [3,1,2]); 
R = squeeze(missVals(fi,:,:)) * squeeze(missVals(fi,:,:))' ./ ...
    prod(size(missVals,[2,3]));

% regularized R
gamma = .01;
Rr = R*(1-gamma) + eye(size(R,1))*gamma*mean(eig(R));

% global variance normalize (optional; this scales the eigenspectrum)
S = S / (std(S(:))/std(R(:)));

% GED
[W,L] = eig(S,Rr);
[evals(:,fi),sidx] = sort(diag(L),'descend');
W = W(:,sidx);

% store top component map and eigenvector
maps(:,fi) = W(:,1)'*S;
evecs(:,fi) = W(:,1);

end


E = zscore(evecs,[],1)';
evecCorMat = (E*E'/(size(R,1)-1)).^2;

clustIDX = DBscanDynamicEpi(evecCorMat, 3, 3, 1, 0);
changeCheck = true; 
while changeCheck
    newClust = subDivideClust(evecCorMat, clustIDX);
    if sum(newClust - clustIDX) == 0
        changeCheck = false;
    else
        clustIDX = newClust; 
    end
end


[silhouettes, clusti] = getSil(evecCorMat, clustIDX);


cIDs = unique(clustIDX); 
silDifs = zeros(length(cIDs),1); 
silMeans = zeros(length(cIDs),1); 
for ci = 1:length(cIDs)
    silMeans(ci) = mean(silhouettes(clusti(:,2)==cIDs(ci)));
    maxSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 90) ; 
    minSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 10) ; 
    silDifs(ci) = maxSil - minSil;
end

breakVals = find(diff(clustIDX(clustIDX>0))~=0); 
outClusts = clustIDX(clustIDX>0); 
figure('visible', false)
subplot 221
outFrex = dat.frex(clustIDX>0); 
outCorMat = evecCorMat(clustIDX>0, clustIDX>0); 
imagesc(corr(outCorMat))
title('eigenvector correlations')
caxis([.05, .8])
hold on 
for bi = 1:length(breakVals)
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    yline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
yticks(breakVals+.5)
yticklabels(round(outFrex(breakVals)))

subplot(2,2,[3,4])
hold off
[silVals, clusti] = getSil(outCorMat, outClusts);
breakVals2 = [0, breakVals', length(silVals)];

for bi = 1:length(breakVals)

    plot(breakVals2(bi)+1: breakVals2(bi+1),...
        silVals(breakVals2(bi)+1: breakVals2(bi+1)), 'color', 'k')
    hold on 
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
title([regions{reg} ' ' phase ' freq band divisions'])
ylabel('silhouette values')

export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' 'GEDbandDiscovery' '_HFB_' regions{reg} '_' phase '.jpg'], '-r300')

outFrexIDX = arrayfun(@(x) find(outFrex(x) == dat.frex), [1:length(outFrex)]);

if strcmp(phase, 'sub')
    sigFreqs(1,outFrexIDX(breakVals),1,2,1) = 1; 
else

    sigFreqs(1,outFrexIDX(breakVals),2,2,1) = 1; 
end

catch
end


%% repeat for Image locked phase reset

perm = load([ImagepermsPhase(1).folder '/' ImagepermsPhase(1).name]).outDat;

perm2 = load([HFBpermsPhase(1).folder '/' HFBpermsPhase(1).name]).outDat;


%doing GED-based frequency band discovery Image locked: 
try
[evals,evecs,maps] = deal(zeros(size(perm.hitVals,[1,3])));

for fi = 1:100


hitVals = permute(perm.hitVals, [3,1,2]); 
S = squeeze(hitVals(fi,:,:)) * squeeze(hitVals(fi,:,:))' ./ ...
    prod(size(hitVals,[2,3]));

missVals = permute(perm.missVals, [3,1,2]); 
R = squeeze(missVals(fi,:,:)) * squeeze(missVals(fi,:,:))' ./ ...
    prod(size(missVals,[2,3]));

% regularized R
gamma = .01;
Rr = R*(1-gamma) + eye(size(R,1))*gamma*mean(eig(R));

% global variance normalize (optional; this scales the eigenspectrum)
S = S / (std(S(:))/std(R(:)));

% GED
[W,L] = eig(S,Rr);
[evals(:,fi),sidx] = sort(diag(L),'descend');
W = W(:,sidx);

% store top component map and eigenvector
maps(:,fi) = W(:,1)'*S;
evecs(:,fi) = W(:,1);

end


E = zscore(evecs,[],1)';
evecCorMat = (E*E'/(size(R,1)-1)).^2;

clustIDX = DBscanDynamicEpi(evecCorMat, 3, 3, 1, 0);
changeCheck = true; 
while changeCheck
    newClust = subDivideClust(evecCorMat, clustIDX);
    if sum(newClust - clustIDX) == 0
        changeCheck = false;
    else
        clustIDX = newClust; 
    end
end


[silhouettes, clusti] = getSil(evecCorMat, clustIDX);


cIDs = unique(clustIDX); 
silDifs = zeros(length(cIDs),1); 
silMeans = zeros(length(cIDs),1); 
for ci = 1:length(cIDs)
    silMeans(ci) = mean(silhouettes(clusti(:,2)==cIDs(ci)));
    maxSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 90) ; 
    minSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 10) ; 
    silDifs(ci) = maxSil - minSil;
end

breakVals = find(diff(clustIDX(clustIDX>0))~=0); 
outClusts = clustIDX(clustIDX>0); 
figure('visible', false)
subplot 221
outFrex = dat.frex(clustIDX>0); 
outCorMat = evecCorMat(clustIDX>0, clustIDX>0); 
imagesc(corr(outCorMat))
title('eigenvector correlations')
caxis([.05, .8])
hold on 
for bi = 1:length(breakVals)
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    yline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
yticks(breakVals+.5)
yticklabels(round(outFrex(breakVals)))

subplot(2,2,[3,4])
hold off
[silVals, clusti] = getSil(outCorMat, outClusts);
breakVals2 = [0, breakVals', length(silVals)];

for bi = 1:length(breakVals)

    plot(breakVals2(bi)+1: breakVals2(bi+1),...
        silVals(breakVals2(bi)+1: breakVals2(bi+1)), 'color', 'k')
    hold on 
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
title([regions{reg} ' ' phase ' freq band divisions'])
ylabel('silhouette values')

export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' 'GEDbandDiscovery_phase' '_' regions{reg} '_' phase '.jpg'], '-r300')

outFrexIDX = arrayfun(@(x) find(outFrex(x) == dat.frex), [1:length(outFrex)]);

if strcmp(phase, 'sub')
    sigFreqs(1,outFrexIDX(breakVals),1,1,2) = 1; 
else

    sigFreqs(1,outFrexIDX(breakVals),2,1,2) = 1; 
end
catch
end

%% phase based HFB locked 

try
[evals,evecs,maps] = deal(zeros(size(perm2.hitVals,[1,3])));

for fi = 1:100


hitVals = permute(perm2.hitVals, [3,1,2]); 
S = squeeze(hitVals(fi,:,:)) * squeeze(hitVals(fi,:,:))' ./ ...
    prod(size(hitVals,[2,3]));

missVals = permute(perm2.missVals, [3,1,2]); 
R = squeeze(missVals(fi,:,:)) * squeeze(missVals(fi,:,:))' ./ ...
    prod(size(missVals,[2,3]));

% regularized R
gamma = .01;
Rr = R*(1-gamma) + eye(size(R,1))*gamma*mean(eig(R));

% global variance normalize (optional; this scales the eigenspectrum)
S = S / (std(S(:))/std(R(:)));

% GED
[W,L] = eig(S,Rr);
[evals(:,fi),sidx] = sort(diag(L),'descend');
W = W(:,sidx);

% store top component map and eigenvector
maps(:,fi) = W(:,1)'*S;
evecs(:,fi) = W(:,1);

end


E = zscore(evecs,[],1)';
evecCorMat = (E*E'/(size(R,1)-1)).^2;

clustIDX = DBscanDynamicEpi(evecCorMat, 3, 3, 1, 0);
changeCheck = true; 
while changeCheck
    newClust = subDivideClust(evecCorMat, clustIDX);
    if sum(newClust - clustIDX) == 0
        changeCheck = false;
    else
        clustIDX = newClust; 
    end
end


[silhouettes, clusti] = getSil(evecCorMat, clustIDX);


cIDs = unique(clustIDX); 
silDifs = zeros(length(cIDs),1); 
silMeans = zeros(length(cIDs),1); 
for ci = 1:length(cIDs)
    silMeans(ci) = mean(silhouettes(clusti(:,2)==cIDs(ci)));
    maxSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 90) ; 
    minSil = prctile(silhouettes(clusti(:,2)==cIDs(ci)), 10) ; 
    silDifs(ci) = maxSil - minSil;
end

breakVals = find(diff(clustIDX(clustIDX>0))~=0); 
outClusts = clustIDX(clustIDX>0); 
figure('visible', false)
subplot 221
outFrex = dat.frex(clustIDX>0); 
outCorMat = evecCorMat(clustIDX>0, clustIDX>0); 
imagesc(corr(outCorMat))
title('eigenvector correlations')
caxis([.05, .8])
hold on 
for bi = 1:length(breakVals)
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    yline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
yticks(breakVals+.5)
yticklabels(round(outFrex(breakVals)))

subplot(2,2,[3,4])
hold off
[silVals, clusti] = getSil(outCorMat, outClusts);
breakVals2 = [0, breakVals', length(silVals)];

for bi = 1:length(breakVals)

    plot(breakVals2(bi)+1: breakVals2(bi+1),...
        silVals(breakVals2(bi)+1: breakVals2(bi+1)), 'color', 'k')
    hold on 
    xline(breakVals(bi)+.5, '--', 'linewidth', 2)
    
end
xticks(breakVals+.5)
xticklabels(round(outFrex(breakVals)))
title([regions{reg} ' ' phase ' freq band divisions'])
ylabel('silhouette values')

export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\TF_regional\wavelet'...
    '/' 'GEDbandDiscovery_phase' '_HFB_' regions{reg} '_' phase '.jpg'], '-r300')

outFrexIDX = arrayfun(@(x) find(outFrex(x) == dat.frex), [1:length(outFrex)]);

if strcmp(phase, 'sub')
    sigFreqs(1,outFrexIDX(breakVals),1,2,2) = 1; 
else

    sigFreqs(1,outFrexIDX(breakVals),2,2,2) = 1; 
end
catch 
end





% 
% 
% 
% 
% 
% 
% 
% [~, idxtims] = max(abs(perm.tVals));
% %peak times for each freq
% test = arrayfun(@(x) perm.tVals(idxtims(x),x),1:50); 
% [~, idxfrex] = max(abs(test));
% idxtims = idxtims(idxfrex);
% 
% [~, idxtims2] = max(abs(perm2.tVals));
% %peak times for each freq
% test = arrayfun(@(x) perm2.tVals(idxtims2(x),x),1:50); 
% [~, idxfrex2] = max(abs(test));
% idxtims2 = idxtims2(idxfrex2);
% 
% scatter(perm.hitVals(:,idxtims, idxfrex), ...
%     perm2.hitVals(:, idxtims2, idxfrex2))



end