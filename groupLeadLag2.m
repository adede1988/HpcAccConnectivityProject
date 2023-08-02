% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'R:\MSS\Johnson_Lab\dtf8829\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath(genpath([codePre 'mni2atlas']))
addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\fieldtrip-20230118')
ft_defaults
datFolder = [datPre 'CHANDAT']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 




errorChans = cell(length(chanFiles),1); 



clustersHipUp = cell(length(chanFiles),1); 
clustersNotUp = cell(length(chanFiles),1); 
clustersHipDown = cell(length(chanFiles),1); 
clustersNotDown = cell(length(chanFiles),1); 

parfor sub = 1:length(chanFiles)
    sub
    [errorChans{sub}, clustersHipUp{sub}, clustersHipDown{sub}, clustersNotUp{sub}, clustersNotDown{sub}] = getClustScratch(chanFiles, sub); 
end


allUpHip = zeros(177, 301); 
allDownHip = allUpHip; 
allUpNot = allUpHip; 
allDownNot = allUpHip; 
c = 0;
c2 = 0;
c3 = 0; 
c4 = 0; 
for ii = 1:length(clustersHipUp)
    if sum(size(clustersHipUp{ii}))>0

        allUpHip = allUpHip + clustersHipUp{ii};
        c = c+1;
    end

    if sum(size(clustersHipDown{ii}))>0

        allDownHip = allDownHip + clustersHipDown{ii};
        c2 = c2+1;
    end

    if sum(size(clustersNotUp{ii}))>0

        allUpNot = allUpNot + clustersNotUp{ii};
        c3 = c3+1;
    end

    if sum(size(clustersNotDown{ii}))>0

        allDownNot = allDownNot + clustersNotDown{ii};
        c4 = c4+1;
    end

end