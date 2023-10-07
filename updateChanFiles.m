%comparison of finished files vs. remaining files and update


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

%% initialize chanFiles, dif versions for local v. quest running


chanFiles = load([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat']).chanFiles;


finFiles = dir([datPre 'CHANDAT/finished']);
test = cellfun(@(x) length(x)>0, strfind({finFiles.name}, '.mat'));
finFiles = finFiles(test); 
clear test


%for each file in finFiles, find it in chanFiles, remove it, make new
%chanFiles output to be saved for both local and quest running.
%also make sure to eliminate files currently running on quest from planned
%next run
rmidx = zeros(length(finFiles),1); 
for ii = 1:length(finFiles)
    rmidx(ii) = find(cellfun(@(x) strcmp(x, finFiles(ii).name), {chanFiles.name}));
     

end

chanFiles(rmidx) = []; 
save([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat'], 'chanFiles')


chanFiles = load([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat']).chanFiles; 
chanFiles(rmidx) = []; 
save([codePre 'HpcAccConnectivityProject' '/questChanFiles.mat'], 'chanFiles')








