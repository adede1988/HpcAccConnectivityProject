

%github access token: 
%  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
%  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
%  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
% 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
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

% start = 501; 
fileidx = test.filei(start); 
statType = test.stati(start); 
permi = test.permi(start); 


% if ~isfile([cndFiles(fileidx).folder '/out/NEWphase'...
%            'stat' num2str(statType) '_' num2str(permi) ...
%            '_' cndFiles(fileidx).name])



tic
TFphaseTrialpipeline(cndFiles, fileidx, statType, permi); 
toc

% end