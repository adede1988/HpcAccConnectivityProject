

%stitching together at the subject level 


datFolder  = "R:\MSS\Johnson_Lab\dtf8829\CHANDAT";
saveFolder = "C:\Users\dtf8829\Documents\QuestConnect\SUMDAT";

chanFiles = dir(datFolder); 
chanFiles = chanFiles([chanFiles.isdir]==false); 

subIDs = cellfun(@(x) split(x, '_'), {chanFiles.name}, 'uniformoutput', false); 
subIDs = cellfun(@(x) x{2}, subIDs, 'uniformOutput', false);

masterSheet = readtable(['R:\MSS\Johnson_Lab\dtf8829\memDevDat.csv']);
datFolder = ['R:\MSS\Johnson_Lab\DATA\'];
task = 'MemDev';
allDat = getAllDataStruct(datFolder, masterSheet, task);

subIDs_unique = unique({allDat.subID});


parfor sub = 1:length(subIDs_unique)
    try
        if ~isfile(join([saveFolder '/sumDat_' subIDs_unique{sub} '.mat'],''))
    sub
    subFiles = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, subIDs_unique{sub})); 
    subFiles = chanFiles(subFiles); 
    chanNums = cellfun(@(x) split(x, '_'), {subFiles.name}, 'uniformoutput', false); 
    chanNums = cellfun(@(x) split(x{3}, '.mat'), chanNums, 'uniformoutput', false);
    chanNums = cellfun(@(x) str2num(x{1}), chanNums);
    [~, order] = sort(chanNums); 
    

    subDat = load([subFiles(1).folder '/' subFiles(1).name]).chanDat; 
    if isfield(subDat, 'sizeReduce') %check that it's through the pipeline! 
    subDat.chanorder = order; 
    subDat = rmfield(subDat, 'enc');
    subDat = rmfield(subDat, 'retOn'); 
    subDat = rmfield(subDat, 'retRT');
    subDat = rmfield(subDat, 'enctim'); 
    subDat = rmfield(subDat, 'retOtim'); 
    subDat = rmfield(subDat, 'retRtim'); 
    subDat.ISPCout = rmfield(subDat.ISPCout, 'encdi'); 
    subDat.ISPCout = rmfield(subDat.ISPCout, 'ondi'); 
    subDat.ISPCout = rmfield(subDat.ISPCout, 'rtdi'); 

   
        
        %preallocate space for speed
        fields = fieldnames(subDat.TFout);  
        for fi = 1:length(fields)
            curSub = subDat.TFout.(fields{fi});
            curSub = nan([size(curSub), length(subFiles)]);  
            subDat.TFout.(fields{fi}) = curSub; 
        end

        fields = fieldnames(subDat.ISPCout);  
        for fi = 1:length(fields)
            curSub = subDat.ISPCout.(fields{fi});
            curSub = nan([size(curSub), length(subFiles)]);  
            subDat.ISPCout.(fields{fi}) = curSub; 
        end


        %loop on channels to stitch it together! 
        for chan = 1:length(subFiles)
            chanDat = load([subFiles(chan).folder '/' subFiles(chan).name]).chanDat;
            
            if isfield(chanDat, 'sizeReduce') %check that this channel has been done
            chanDat.ISPCout = rmfield(chanDat.ISPCout, 'encdi'); 
            chanDat.ISPCout = rmfield(chanDat.ISPCout, 'ondi'); 
            chanDat.ISPCout = rmfield(chanDat.ISPCout, 'rtdi');  

            %loop over fields in the time frequency output and concatenate
            fields = fieldnames(subDat.TFout);  
            for fi = 1:length(fields)
                subDat.TFout.(fields{fi})(:,:,chan) = chanDat.TFout.(fields{fi});
            end

            %loop over fields in the ISPC output and concatenate
            fields = fieldnames(subDat.ISPCout);
            for fi = 1:length(fields)
                subDat.ISPCout.(fields{fi})(:,:,:,:,chan) = chanDat.ISPCout.(fields{fi}); 
            end

            else
                %need to fill with nans
                disp([subFiles(chan).folder '/' subFiles(chan).name])
            

            end
        end
    end

    

    parsave(join([saveFolder '/sumDat_' subIDs_unique{sub} '.mat'],''), subDat)
        end
    catch
        disp(['subject failure! ' num2str(sub)])
    end

end


