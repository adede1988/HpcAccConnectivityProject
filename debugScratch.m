



%tracing through possible missing data

%load all data
load("C:\Users\dtf8829\Documents\QuestConnect\allDat.mat")

%load test file
load('R:\MSS\Johnson_Lab\dtf8829\permDat\hit_on_phg_hip.mat')

input = {allDat.(cndDat.targVar)};
inID = {allDat.subID}; 


bad = cellfun(@(x) isempty(x), input); 

inID(bad) = []; 

outID = {cndDat.metaDat.subID};


%get a subject file 
load('C:\Users\dtf8829\Documents\QuestConnect\SUMDAT\sumDat_LC02.mat')

testConnection = allDat(47).hit_on_phg_hip; 
testConnection = allDat(47).miss_on_phg_hip; 

%data look good at the subject level. They even look good at the allDat
%level. This means something must be wrong in extracting the condition data
%working with code block starting at line 99 of groupLevel.m

find(cellfun(@(x) ~isempty(x),  strfind(varNames, 'hit_on_phg_hip')))

ii = find(cellfun(@(x) ~isempty(x),  strfind(varNames, 'hit_on_phg_hip')))