

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
tic
singleChanPipeline(chanFiles, curChani, datPre); 
toc

%% put encoding behavioral data onto it
% curSub = 'something'; 
% for ii = 482:length(chanFiles)
%     ii
%     chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
%     if strcmp(curSub, chanDat.subID)
%         chanDat.encInfo = dat.trialinfo; 
%         parsave([chanFiles(ii).folder '/' chanFiles(ii).name], chanDat)
%     else
%         dat = load(fullfile([chanDat.dataDir '/' chanDat.encDatFn])).data; 
%         chanDat.encInfo = dat.trialinfo; 
%         parsave([chanFiles(ii).folder '/' chanFiles(ii).name], chanDat)
%         curSub = chanDat.subID; 
%     end
%    
% end