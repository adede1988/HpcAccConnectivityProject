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

outStatFiles = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'all.mat'));
outStatFiles(~test) = []; 
headFiles = dir([datPre 'TF_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, 'all.mat'));
headFiles(~test) = []; 

parfor reg = 2:length(regions)
    reg
    TFfinalintegrate(reg, headFiles, outStatFiles, regions, 'sub')
    TFfinalintegrate(reg, headFiles, outStatFiles, regions, 'ret')


end
