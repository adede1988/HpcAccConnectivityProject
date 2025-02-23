

clear

%github access token: 
%  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
%  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
%July 3rd: 
% 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

%code path
addpath(genpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject'))
addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\myFrequentUse')
addpath(genpath('C:\Users\dtf8829\Documents\MATLAB\fieldtrip-20230118'))

%% initialize the data structures 
prefix = 'R:\';
task = 'MemDev';
datFolder = [prefix 'MSS\Johnson_Lab\DATA\'];
masterSheet = readtable([prefix 'MSS\Johnson_Lab\dtf8829\GitHub\'...
    'HpcAccConnectivityProject\memDevDat2.csv']);
saveFolder = [prefix 'MSS\Johnson_Lab\dtf8829\SUMDAT\'];
chanFolder = [prefix 'MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\CHANRAW'];
allData = getAllDataStruct(datFolder, masterSheet, task);

clear masterSheet task 





%% run pipeline

allErrors = cell(37,1); 


parfor ii = 1:length(allData)

try
allErrors{ii} = readAndSplitPipeline(allData(ii), saveFolder, chanFolder);
catch
allErrors{ii} = "RERUN THIS ONE"; 
end



end



% %% RT error audit 
% 
% errorLog = struct; 
% errorLog(1).subID = 'a'; 
% errorLog(1).numTrials = 0; 
% errorLog(1).trials = [1,2]; 
% 
% for ii = 1:length(allData)
% disp(ii)
% 
% errorLog = auditRTerror(allData(ii), errorLog, ii);
% 
% end






