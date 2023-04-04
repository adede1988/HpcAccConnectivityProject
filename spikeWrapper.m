%simulated HPC looping input from shell script
%start = 30;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% codePre = 'C:\Users\dtf8829\Documents\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\DATA\Old Data\iEEG-WM\IR84\raw_data';
saveLoc = 'R:\MSS\Johnson_Lab\smg1656\iEEG-WM\micros\SUMDAT';
%HPC paths: 

codePre = '/projects/p31578/dtf8829/';
datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

% addpath([codePre 'HpcAccConnectivityProject'])
% addpath([codePre 'myFrequentUse'])
% addpath(genpath([codePre 'fieldtrip-20230118']))

%% initialize 


chanFiles = dir(datPre);
test = cellfun(@(x) strcmp(x(1), 'm'), {chanFiles.name} );
chanFiles = chanFiles(test); 


%% run the pipeline


disp(['going for ' subID ' ' num2str(start)] )

singleChanSpikePipeline(chanFiles, start, saveLoc); 