function [subDat] = connectivityPipeline(subDat, prefix, regModels, saveFolder)

%% load in the data

dat = load([subDat.dataDir '\' subDat.encDatFn]).data;

%make sure there's an elecpos field
if ~isfield(dat.elec, 'elecpos')
    dat.elec.elecpos = dat.elec.chanpos; 
end

%get previous work
if isfile([saveFolder 'sumDat_' subDat.subID '.mat'])
    subDat = load([saveFolder 'sumDat_' subDat.subID '.mat']).subDat;
end

% check for corrected_coords.mat file
%NOTE: THIS MAY NOT BE NEEDED ANYMORE
if isfile([subDat.dataDir '\corrected_coords.mat'])
    dat.elec = load([subDat.dataDir '\corrected_coords.mat']).data_new.elec;
end

%% set some general params for frequencies, ROI names
frex = logspace(log10(2),log10(80),100);
numfrex = length(frex); 
stds = linspace(2,5,numfrex);
rois = {'dlPFC', 'hip', 'phg', 'acc'}; 


%% quick stuff grabbing some behavior 
subDat.use = dat.trialinfo(:,1)==1; %trials where participant paid attention
subDat.hits = dat.trialinfo(:,2)==1; %trials that resulted in a subsequent hit
subDat.misses = dat.trialinfo(:,2)==2; %trials that resulted in a subsequent miss


%% load and assess anatomy


if ~isfield(subDat, 'labels')
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
    
    %save out anatomy checking
    save([saveFolder 'sumDat_' subDat.subID '.mat'], 'subDat')
end


%% split the data into channel files

if ~isfield(subDat, 'chanSplit')
    allDat = makeAllDat(dat.trial); 
    for ch = 1:size(allDat,2)


    end

    %record that the chan split has been done
    subDat.chanSplit = 1; 
    save([saveFolder 'sumDat_' subDat.subID '.mat'], 'subDat')
end



%% time frequency 

if ~isfield(subDat, 'tfLow')

    
   
    
    


    subDat.tfLow = getTrialTF(dat.trial, frex, numfrex, stds, dat.fsample);


   
    pow = abs(subDat.tfLow ).^2;
    tic
    powZ = arrayfun(@(x) myZscore(pow(:,:,:,x)), [1:size(pow,4)], 'uniformoutput', false );
    toc
    





end










end








%stuff for getting the IRASA results of the data: 
% 
%    
% 
%     aperiodic = zeros(size(allDat,1), 100); 
%     periodic = zeros(size(allDat,1), 100); 
%     powSpect = zeros(size(allDat,1), 100); 
% %     tic
%     for ch=1:size(allDat,1)
%         [aperiodic(ch,:), periodic(ch,:), powSpect(ch,:)] = IRASA(allDat(ch,:,:)); 
%     end
% %     toc
%     subDat.aperiodic = aperiodic; 
%     subDat.periodic = periodic; 
%     subDat.powSpect = powSpect; 
%     save([saveFolder 'sumDat_' subDat.subID '.mat'], 'subDat')
% end
