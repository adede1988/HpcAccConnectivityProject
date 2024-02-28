

%github access token: 
%  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
%  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
%  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
% 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
%Feb 20: 
%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



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

datFolder = [datPre 'pairFiles']; 
cndFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({cndFiles.name}, '.mat'));
cndFiles = cndFiles(test); 

test = readtable([codePre 'HpcAccConnectivityProject/ConnectivityStatMaster.csv']);

%% key variables
regions = {'acc', 'dlPFC', 'hip', 'lTemp', 'iTemp', 'mtl', ...
    'pcc', 'pPFC', 'vis'}; 
phase = {'enc' , 'ret'}; 

%% run the pipeline

reg1 = test.reg1(start); 
reg2 = test.reg2(start);
encRet = test.encRet(start); 
statType = test.stati(start); %0 = HFB locked %1 = image locked
permi = test.permi(start); 

idx = cellfun(@(x) length(x)>0, ...
    strfind({cndFiles.name}, [regions{reg1} '_' regions{reg2}]));
cndFiles(~idx) = []; 

tic
connectivitypipeline(cndFiles, reg1, reg2, encRet, statType, permi); 
toc

