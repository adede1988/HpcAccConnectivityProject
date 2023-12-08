%local paths: 


codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
HFBdat = load('R:\MSS\Johnson_Lab\dtf8829\QuestConnect\HFB_KEY_STATS\hip.mat').HFBdat; 
regions = {HFBdat.aggTargs.lab}; 
regions(2) = []; 




%% TF final evaluation 
%it's assumed here that the TF quest pipeline has already been run and that
%the outputs are available 

%TF power files
outStatFiles = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'all.mat'));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'phase'));
outStatFiles(test) = []; 
headFiles = dir([datPre 'TF_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, 'all.mat'));
headFiles(~test) = []; 

%TF phase files
outStatFilesPhase = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'phase'));
outStatFilesPhase(~test) = []; 

%HFB files
outStatFilesHFB = dir([datPre 'HFB_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesHFB.name}, '.mat'));
outStatFilesHFB(~test) = []; 
headFilesHFB = dir([datPre 'HFB_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFilesHFB.name}, '.mat'));
headFilesHFB(~test) = []; 

for reg = 1:length(regions)
    reg
    TFfinalintegrate(reg, headFiles, outStatFiles, regions, ...
        'sub', outStatFilesPhase)
    TFfinalintegrate(reg, headFiles, outStatFiles, regions, ...
        'ret', outStatFilesPhase)

    HFBfinalintegrate(reg, headFilesHFB, outStatFilesHFB, regions, 'sub')
    HFBfinalintegrate(reg, headFilesHFB, outStatFilesHFB, regions, 'ret')


end
