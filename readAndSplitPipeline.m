function [] = readAndSplitPipeline(subDat, prefix, regModels, saveFolder, chanFolder)

%get previous work
% if isfile([saveFolder 'sumDat_' subDat.subID '.mat'])
%     subDat = load([saveFolder 'sumDat_' subDat.subID '.mat']).subDat;
% end
% if isfield(subDat, 'chanSplit')
%     return
% end

%% load in the data

dat = load([subDat.dataDir '\' subDat.encDatFn]).data;
dat2 = load([subDat.dataDir '\' subDat.retDatFn]).data;

%make sure there's an elecpos field
if ~isfield(dat.elec, 'elecpos')
    dat.elec.elecpos = dat.elec.chanpos; 
end



% check for corrected_coords.mat file
%NOTE: THIS MAY NOT BE NEEDED ANYMORE
% if isfile([subDat.dataDir '\corrected_coords.mat'])
%     dat.elec = load([subDat.dataDir '\corrected_coords.mat']).data_new.elec;
% end

%get sampling rate
subDat.fsample = dat.fsample;

%% set some general params for frequencies, ROI names
% frex = logspace(log10(2),log10(80),100);
% numfrex = length(frex); 
% stds = linspace(2,5,numfrex);
rois = {'dlPFC', 'hip', 'phg', 'acc'}; 


%% quick stuff grabbing some behavior 
subDat.use = dat.trialinfo(:,1)==1; %trials where participant paid attention
subDat.hits = dat.trialinfo(:,2)==1; %trials that resulted in a subsequent hit
subDat.misses = dat.trialinfo(:,2)==2; %trials that resulted in a subsequent miss
subDat.retInfo = dat2.trialinfo; %store the trial info for retrieval 
subDat.encInfo = dat.trialinfo; 

%check for RT mistake and eliminate trials which don't pass
RT = dat2.trialinfo(:,3); 
trialLengths = cellfun(@(x) size(x,2), dat2.trial);

% trialLengths = zeros(length(dat2.trial),1); 
% for ii = 1:length(trialLengths)
%     trialLengths(ii) = size(dat2.trial{ii},2); 
% end

errorTrials = find(arrayfun(@(x) RT(x)>trialLengths(x)-subDat.fsample, 1:length(RT)));
subDat.retUse = ones(size(RT)); 
subDat.retUse(errorTrials) = 0; 

subDat.retInfo(errorTrials, :) = [];




%% load and assess anatomy


% if ~isfield(subDat, 'labels')
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
%     save([saveFolder 'sumDat_' subDat.subID '.mat'], 'subDat')
% end


%% split the data into channel files

% if ~isfield(subDat, 'chanSplit')
    allDatEnc = makeAllDat(dat.trial, dat.time, dat.fsample); 
    allDatEncRT = makeAllDatEncRT(dat.trial, dat.time, dat.fsample, subDat.encInfo(:,4));

    allDatRetON = makeAllDatRetON(dat2.trial, dat2.time, errorTrials, dat.fsample); 
    allDatRetRT = makeAllDatRetRT(dat2.trial, dat2.time, dat2.trialinfo(:,3), errorTrials, dat.fsample);
    dat.fsample = 1000; 
    for ch = 1:size(allDatEnc,1)
        if ~subDat.badTrodes %only save if it's a good electrode
            chanDat = subDat; 
            chanDat.enc = squeeze(allDatEnc(ch,:,:)); %data are in trials X time -1000ms:3500ms
            chanDat.encRT = squeeze(allDatEncRT(ch,:,:)); %data are in trials X time -2000ms:500ms
            chanDat.retOn = squeeze(allDatRetON(ch,:,:)); %data are in trials X time -1000ms:2000ms
            chanDat.retRT = squeeze(allDatRetRT(ch,:,:)); %data are in trials X time -2000ms:500ms
            chanDat.chi = ch; %note which electrode it is for reference into other structs
            if ch<10
                save([chanFolder '\' 'chanDat_' chanDat.subID '_' '00' num2str(ch) '.mat'], 'chanDat')
            elseif ch<100
                save([chanFolder '\' 'chanDat_' chanDat.subID '_' '0' num2str(ch) '.mat'], 'chanDat')
            else

                save([chanFolder '\' 'chanDat_' chanDat.subID '_' num2str(ch) '.mat'], 'chanDat')
            end

        end

    end

    %record that the chan split has been done
    subDat.chanSplit = 1; 
    save([saveFolder 'sumDat_' subDat.subID '.mat'], 'subDat')
% end





end





% %% time frequency 
% 
% if ~isfield(subDat, 'tfLow')
% 
%     
%    
%     
%     
% 
% 
%     subDat.tfLow = getTrialTF(dat.trial, frex, numfrex, stds, dat.fsample);
% 
% 
%    
%     pow = abs(subDat.tfLow ).^2;
%     tic
%     powZ = arrayfun(@(x) myZscore(pow(:,:,:,x)), [1:size(pow,4)], 'uniformoutput', false );
%     toc
%     
% 
% 
% 
% 
% 
% end
% 
% 
% 
% 







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
