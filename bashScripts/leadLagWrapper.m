%This is a pipeline wrapper script to analyze individual subject leadlag



%github access token: 
%  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
%  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



%simulated HPC looping input from shell script

%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'R:\MSS\Johnson_Lab\dtf8829\';

%HPC paths: 

codePre = '/projects/p31578/dtf8829/';
datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% initialize 

datFolder = [datPre 'HFBCONDAT']; 
subFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, '.mat'));
subFiles = subFiles(test); 

%run the pipeline! 
tic
leadLagPipeline(subFiles, start)
toc