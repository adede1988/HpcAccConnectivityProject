

%github access token: 
%Feb 20: 
%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
%July 3rd: 
% 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


%simulated HPC looping input from shell script

%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%HPC paths: 

codePre = '/projects/p31578/dtf8829/';
datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% initialize 

datFolder = [datPre 'TF_singleTrial']; 
cndFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({cndFiles.name}, 'all.mat'));
cndFiles = cndFiles(test); 

test = readtable([codePre 'HpcAccConnectivityProject/TFstatMaster.csv']);


%% run the pipeline


fileidx = test.filei(start); 
statType = test.stati(start); 
permi = test.permi(start); 


if ~isfile([cndFiles(fileidx).folder '/out/'...
           'stat' num2str(statType) '_' num2str(permi) ...
           '_' cndFiles(fileidx).name])
tic
TFsingleTrialpipeline(cndFiles, fileidx, statType, permi); 
toc

end

