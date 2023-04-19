

%Pipeline wrapper to run permutation cluster significance test across
%condition files

%github access token: 
%  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



%simulated HPC looping input from shell script
%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% pre = 'R:\MSS\Johnson_Lab\dtf8829\';

%HPC paths: 

pre = '/projects/p31578/dtf8829/';

%% set paths

addpath([pre 'github/HpcAccConnectivityProject'])
addpath(genpath([pre 'github/myFrequentUse']))


%% initialize 

datFolder = [pre 'permDat']; 
datFiles = dir(datFolder); 
test = cellfun(@(x) contains(x, '.mat'), {datFiles.name});
datFiles = datFiles(test); 

%% run the pipeline

disp(['going for ' datFiles(start).name] )


conditionClusterTestPipeline(datFiles(start))