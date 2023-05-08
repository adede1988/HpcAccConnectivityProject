

%Pipeline wrapper to run permutation cluster significance test across
%condition files

%github access token: 
%  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
%  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



%simulated HPC looping input from shell script
%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% pre = 'R:\MSS\Johnson_Lab\dtf8829\github\';
% preDat = 'R:\MSS\Johnson_Lab\dtf8829\';

%HPC paths: 

pre = '/projects/p31578/dtf8829/';
preDat = '/projects/p31578/dtf8829/';

%% set paths

addpath([pre 'HpcAccConnectivityProject'])
addpath(genpath([pre 'myFrequentUse']))


%% initialize 

datFolder = [preDat 'permDat']; 
datFiles = dir(datFolder); 
test = cellfun(@(x) contains(x, '.mat'), {datFiles.name});
datFiles = datFiles(test); 

%% run the pipeline

disp(['going for ' datFiles(start).name] )


conditionClusterTestPipeline(datFiles(start))


%% audit pipeline
%NOTE: requires RT error audit in workspace as "errorLog"
%the errorLog is produced by the OG pipelineWrapper script
allErrors = conditionFileAudit(datFiles(1), errorLog); 
ei = length(allErrors)+1; 
for ii = 2:length(datFiles)
    ii

    temp = conditionFileAudit(datFiles(ii), errorLog);
    if isfield(temp, 'subID')
    for jj = 1:length(temp)
        allErrors(ei) = temp(jj); 
        ei = ei+1; 
    end
    end




end





