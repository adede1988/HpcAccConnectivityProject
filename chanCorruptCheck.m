

% checking for corrupt files in CHANDAT
datFolder  = "R:\MSS\Johnson_Lab\dtf8829\CHANDAT";

chanFiles = dir(datFolder); 
chanFiles = chanFiles([chanFiles.isdir]==false); 
chanFiles(1).good = 1; 

parfor ii = 1:length(chanFiles)
    ii
    try 
        temp = load([chanFiles(ii).folder '/' chanFiles(ii).name]);
        chanFiles(ii).good = 1; 
    catch
        chanFiles(ii).good = 0; 
        

    end


end

badFiles = find([chanFiles.good]==0); 
badFiles = 2651;


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



%%
for ii = 1:length(badFiles)

    cur = chanFiles(badFiles(ii)).name;
    subID = split(cur, '_'); 
    ch = split(subID{3}, '.mat'); 
    ch = str2num(ch{1}); 
    subID = subID{2}; 
    subi = find(cellfun(@(x) strcmp(subID, x), {allData.subID})); 

    subDat = allData(subi); 

    %% load in the data

    dat = load([subDat.dataDir '\' subDat.encDatFn]).data;
    dat2 = load([subDat.dataDir '\' subDat.retDatFn]).data;
    
    %make sure there's an elecpos field
    if ~isfield(dat.elec, 'elecpos')
        dat.elec.elecpos = dat.elec.chanpos; 
    end
    
    
    
    % check for corrected_coords.mat file
    %NOTE: THIS MAY NOT BE NEEDED ANYMORE
%     if isfile([subDat.dataDir '\corrected_coords.mat'])
%         dat.elec = load([subDat.dataDir '\corrected_coords.mat']).data_new.elec;
%     end
    
    %get sampling rate
    subDat.fsample = dat.fsample;


    rois = {'dlPFC', 'hip', 'phg', 'acc'}; 

    %% quick stuff grabbing some behavior 
    subDat.use = dat.trialinfo(:,1)==1; %trials where participant paid attention
    subDat.hits = dat.trialinfo(:,2)==1; %trials that resulted in a subsequent hit
    subDat.misses = dat.trialinfo(:,2)==2; %trials that resulted in a subsequent miss
    subDat.retInfo = dat2.trialinfo; %store the trial info for retrieval 


    %% load and assess anatomy

    %get rid of nan electrode locations
    subDat.elecpos = dat.elec.elecpos; 
    subDat.badTrodes = isnan(subDat.elecpos(:,1)); %key for later during further analysis steps
    subDat.elecpos(subDat.badTrodes,:) = [];

    %standardize the labels from the anatomy notes
    [subDat.labels, subDat.labErrors] = getLabs(dat.elec, subDat.elecNotes);

    %convert the labels to ROI specific labels
    roiLab = roiLabel(subDat.labels(:,3));

    %create lookup tables for ROI membership based on labels and MNI space
    %model coregistration 
    [subDat.roimni, subDat.roiNote] = anatomyPlot(regModels, subDat.elecpos, [subDat.labels(:,3)], roiLab, rois, subDat.subID, -1, 0); 

    %% split the data into channel files
    datFolder  = 'R:\MSS\Johnson_Lab\dtf8829\CHANDAT';

    allDatEnc = makeAllDat(dat.trial, dat.time); 
    allDatRetON = makeAllDatRetON(dat2.trial, dat2.time); 
    allDatRetRT = makeAllDatRetRT(dat2.trial, dat2.time, dat2.trialinfo(:,3));
    
    chanDat = subDat; 
    chanDat.enc = squeeze(allDatEnc(ch,:,:)); %data are in trials X time -1000ms:3500ms
    chanDat.retOn = squeeze(allDatRetON(ch,:,:)); %data are in trials X time -1000ms:2000ms
    chanDat.retRT = squeeze(allDatRetRT(ch,:,:)); %data are in trials X time -2000ms:500ms
    chanDat.chi = ch; %note which electrode it is for reference into other structs
    save([datFolder '/' 'chanDat_' chanDat.subID '_' num2str(ch) '.mat'], 'chanDat')

    %% send it to do the singleChanPipeline
    test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, subID)); 
    idx = 1:length(chanFiles); 
    subChans = chanFiles(test);
    idx = idx(test); 
    curChani = find(idx == badFiles(ii));
    singleChanPipeline(subChans, curChani); 

end





