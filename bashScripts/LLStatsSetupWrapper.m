%% generate leadLag summary files pipeline wrapper


%github access token: 
%  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
%  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
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

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% initialize subject files

datFolder = [datPre 'SUMDAT']; 
subFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, 'leadLag'));
subFiles = subFiles(test); 

%% set the regions we're going for

reg1 = mod(start,10)+1; 
reg2 = ceil(start/10);  

labs = {'ZZZ', 'iTemp', 'acc', 'hip', 'lTemp', 'dlPFC', 'pPFC', 'vis', 'pcc', 'mtl'}; 

disp(['reg1: ' num2str(reg1) ' reg2:' num2str(reg2)])




getSigConnection(subFiles, reg1, reg2, labs, datPre)

