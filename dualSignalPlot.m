function [] = dualSignalPlot(outStatFiles, regions, outStatFilesPhase,...
    headFiles, reg, phase)

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


[~, idxtims] = max(abs(perm.tVals));
%peak times for each freq
test = arrayfun(@(x) perm.tVals(idxtims(x),x),1:50); 
[~, idxfrex] = max(abs(test));
idxtims = idxtims(idxfrex);

[~, idxtims2] = max(abs(perm2.tVals));
%peak times for each freq
test = arrayfun(@(x) perm2.tVals(idxtims2(x),x),1:50); 
[~, idxfrex2] = max(abs(test));
idxtims2 = idxtims2(idxfrex2);

scatter(perm.hitVals(:,idxtims, idxfrex), ...
    perm2.hitVals(:, idxtims2, idxfrex2))



end