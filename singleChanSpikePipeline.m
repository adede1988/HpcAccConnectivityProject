function [] = singleChanSpikePipeline(chanFiles, idx, saveFolder)

%% set global params
highpassThresh = 300; 

%% load the data
rawDat = load([chanFiles(idx).folder '/' chanFiles(idx).name]); 

disp(['data loaded: ' num2str(idx)])

%% time frequency decomposition, extract TF summaries for target trial types: 
% subsequent hit / subsequent miss (encoding data)
% hit / miss / CR / FA (retrieval locked to onset data)
% hit / miss / CR / FA (retrieval locked to response data)

if ~isfile([saveFolder '/' 'ChanSummary_' chanFiles(idx).name '.mat'] )
    chanSum = struct; 
    chanSum.fn = chanFiles(idx).name; 
    chanSum.folder = chanFiles(idx).folder; 
    %save your chanSum
else
    %load the previous chanSum
end


if ~isfield(chanSum, 'downsample')
    disp('working on TF')
    

    %DO SOME DOWNSAMPLING 


    save([saveFolder '/' 'ChanSummary_' chanFiles(idx).name '.mat'], 'chanSum'); 

%     disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
else
    disp('downsample already done, exiting')
end


if ~isfield(chanSum, 'spikeWF')
%     disp('working on TF')
    

    %find waveforms


%     save([chanFiles(idx).folder '/' chanFiles(idx).name], 'chanDat'); 

%     disp(['save success: ' chanFiles(idx).folder '/' chanFiles(idx).name])
else
%     disp('downsample already done, exiting')
end





end