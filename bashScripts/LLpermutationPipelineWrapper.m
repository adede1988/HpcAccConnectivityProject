%pipeline wrapper to run permutation cluster significance tests for lead
%lag analysis

%github access token: 
%  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
%  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
% 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

%simulated HPC looping input from shell script
%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% pre = 'R:\MSS\Johnson_Lab\dtf8829\github\';
% preDat = 'R:\MSS\Johnson_Lab\dtf8829\';

%HPC paths: 

pre = '/projects/p31578/dtf8829/';
preDat = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath([pre 'HpcAccConnectivityProject'])
addpath(genpath([pre 'myFrequentUse']))

%% initialize 

sumDat = dir([preDat 'SUMDAT']); 
LLidx = find(cellfun(@(x) length(x)>0, strfind({sumDat.name}, '_LL.mat')));
sumDat = sumDat(LLidx); 

LLpermutationPipeline(sumDat(start))