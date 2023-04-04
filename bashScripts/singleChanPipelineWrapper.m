

%simulated HPC looping input from shell script
%start = 30;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% codePre = 'C:\Users\dtf8829\Documents\GitHub\';
% datPre = 'C:\Users\dtf8829\Documents\QuestConnect\';

%HPC paths: 

codePre = '/projects/p31578/dtf8829/';
datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% initialize 

datFolder = [datPre 'CHANDAT']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 


%% run the pipeline

curChan = chanFiles(start).name; 
subID = split(curChan, '_'); 
subID = subID{2}; 
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, subID)); 
idx = 1:length(chanFiles); 
chanFiles = chanFiles(test);
idx = idx(test); 
curChani = find(idx == start);

disp(['going for ' subID ' ' num2str(curChani)] )

singleChanPipeline(chanFiles, curChani); 