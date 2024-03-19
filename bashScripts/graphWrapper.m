

%github access token: 
%  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
%  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
%  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
% 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
%Feb 20: 
%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



%simulated HPC looping input from shell script

%start = 1;
%targLen = 'SHORT'

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
regions = {'acc', 'dlPFC', 'hip', 'lTemp', 'iTemp', 'mtl', ...
    'pcc', 'pPFC', 'vis'}; 
keyRegIdx = [1,2,3,6,8]; 
datFolder = [datPre 'graphAnalysis']; 
connectionFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({connectionFiles.name}, '.mat'));
connectionFiles = connectionFiles(test); 

%make region columns 
splitNames = cellfun(@(x) split(x, '_'), ...
    {connectionFiles.name}, 'uniformOutput', false);
reg1 = cellfun(@(x) x{2}, splitNames, 'UniformOutput',false); 
reg2 = cellfun(@(x) x{3}, splitNames, 'UniformOutput',false); 
lenVal = cellfun(@(x) x{1}, splitNames, 'UniformOutput',false); 
test = cellfun(@(x) ismember( x,regions(keyRegIdx)), reg1); 
test2 =cellfun(@(x) ismember( x,regions(keyRegIdx)), reg2); 
test3 = cellfun(@(x) strcmp(x, targLen), lenVal); 
%cut to the key regions! 
connectionFiles = connectionFiles(test & test2 & test3); 

targLen


%% key variables




phases = {'enc', 'ret'}; 
statName = {'HFB', 'image'};

outFile = load([connectionFiles(start).folder '/' ...
                connectionFiles(start).name]).outFile;
curDat = load([datPre 'Figure3/' outFile.figFile]).outDat; 

phase = outFile.encRet; 
pp = outFile.pp; 
cc = outFile.cc; 
stat = outFile.HFBimage; 


%% run the pipeline
tic
connectionGraph(curDat, regions, phase, stat, cc, pp, datPre);
toc
