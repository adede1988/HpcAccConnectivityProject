function [errorLog] = auditRTerror(datInfo, errorLog,ii)



dat = load([datInfo.dataDir '\' datInfo.retDatFn]).data;


RT = dat.trialinfo(:,3); 

trialLengths = cellfun(@(x) size(x,2), dat.trial);


errorTrials = find(arrayfun(@(x) RT(x)>trialLengths(x)-1000, 1:length(RT)));


errorLog(ii).subID = datInfo.subID; 
errorLog(ii).numTrials = length(errorTrials);
errorLog(ii).trials = errorTrials; 








end