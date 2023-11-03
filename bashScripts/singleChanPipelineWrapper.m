

%github access token: 
%  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
%  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
% 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
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

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% initialize chanFiles, dif versions for local v. quest running


% chanFiles = load([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat']).chanFiles;
% chanFiles = load([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat']).chanFiles;


% 
datFolder = [datPre 'CHANDAT']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 


% %I need a way of assigning the partner channel 
% partnerChans = struct; 
% pi = 1; 
% tic
% for chan = 1:length(chanFiles)
%    
%     curChan = chanFiles(chan).name; 
%     subID = split(curChan, '_');
% 
%     subID = subID{2}; 
%     
%     subFiles = dir([datPre 'CHANDAT/CHANRAW']);
%     test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
%     subFiles = subFiles(test);
%     
%     curi = find(cellfun(@(x) strcmp(x, chanFiles(chan).name), {subFiles.name}));
%     if curi <length(subFiles)
%     for chan2 = curi+1:length(subFiles)
%         partnerChans(pi).name = chanFiles(chan).name; 
%         partnerChans(pi).folder = [chanFiles(chan).folder '/CHANRAW']; 
%         partnerChans(pi).name2 = subFiles(chan2).name; 
%         partnerChans(pi).folder2 = subFiles(chan2).folder; 
%         pi = pi+1; 
%         if mod(pi,1000)==0
%             disp(pi)
%             toc
%             tic
%         end
%     end
%     end
% 
% end
% 
% save([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat'], 'partnerChans')
% 
% partnerChans = struct; 
% pi = 1; 
% tic
% for chan = 1:length(chanFiles)
%    
%     curChan = chanFiles(chan).name; 
%     subID = split(curChan, '_');
% 
%     subID = subID{2}; 
%     
%     subFiles = dir([datPre 'CHANDAT/CHANRAW']);
%     test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
%     subFiles = subFiles(test);
%     
%     curi = find(cellfun(@(x) strcmp(x, chanFiles(chan).name), {subFiles.name}));
%     if curi <length(subFiles)
%     for chan2 = curi+1:length(subFiles)
%         partnerChans(pi).name = chanFiles(chan).name; 
%         partnerChans(pi).folder = [chanFiles(chan).folder '/CHANRAW']; 
%         partnerChans(pi).name2 = subFiles(chan2).name; 
%         partnerChans(pi).folder2 = subFiles(chan2).folder; 
%         pi = pi+1; 
%         if mod(pi,1000)==0
%             disp(pi)
%             toc
%             tic
%         end
%     end
%     end
% 
% end



%% run the pipeline

curChan = chanFiles(start).name; 
subID = split(curChan, '_'); 
subID = subID{2}; 

subFiles = dir([datPre 'CHANDAT/CHANRAW']);
test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 

subFiles = subFiles(test);



disp(['going for ' subID ' ' chanFiles(start).name] )
tic
singleChanPipeline(chanFiles, start, subFiles, codePre); 
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