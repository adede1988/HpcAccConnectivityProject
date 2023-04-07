

clear

%github access token: 
%  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

%code path
addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject')
addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\myFrequentUse')
addpath(genpath('C:\Users\dtf8829\Documents\MATLAB\fieldtrip-20230118'))

%% initialize the data structures 
prefix = 'R:\';
task = 'MemDev';
datFolder = [prefix 'MSS\Johnson_Lab\DATA\'];
masterSheet = readtable([prefix 'MSS\Johnson_Lab\dtf8829\memDevDat.csv']);
saveFolder = [prefix 'MSS\Johnson_Lab\dtf8829\SUMDAT\'];
chanFolder = [prefix 'MSS\Johnson_Lab\dtf8829\CHANDAT\'];
allData = getAllDataStruct(datFolder, masterSheet, task);

clear masterSheet task 




%% bring in the anatomical models
dlPFC = stlread([prefix 'MSS\Johnson_Lab\dtf8829\MNI_based_anat\dlPFC.stl']);
hipp = stlread([prefix 'MSS\Johnson_Lab\dtf8829\MNI_based_anat\bilatHippo.stl']);
phg = stlread([prefix 'MSS\Johnson_Lab\dtf8829\MNI_based_anat\ParahippoCortex.stl']);
acc = stlread([prefix 'MSS\Johnson_Lab\dtf8829\MNI_based_anat\ACC.stl']);

regModels = {dlPFC, hipp, phg, acc}; %following naming convention in function getLabs.m

clear dlPFC hipp phg acc

%% run pipeline

parfor ii = 1:length(allData)
ii

readAndSplitPipeline(allData(ii), prefix, regModels, saveFolder, chanFolder);




end




%% scratch 
%preallocate variables in the struct

%anatomy labels
% allData(1).labels = []; 
% allData(1).labErrors = []; 
% allData(1).elecpos = []; 
% allData(1).roimni = []; 
% allData(1).roiNote = []; 

%irasa output


