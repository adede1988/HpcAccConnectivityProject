




addpath('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject')

chanFiles = dir("R:\MSS\Johnson_Lab\dtf8829\CHANDAT\CHANRAW");
chanFiles([1:2]) = []; 



outputDat = struct; 
outi = 1; 




for chan = 1:length(chanFiles)

    chan
   if isfield(outputDat, 'subID')
       IDs = unique({outputDat.subID}); 
   else
       IDs = {'nope'}; 
   end

   curOut =   findDataFromChanRaw(chanFiles, chan, IDs); 
    
    if ~isempty(curOut)
%         if length(outputDat)>1
%             IDs = unique(outputDat.subID);
        if ~isfield(outputDat, 'subID') 
            outputDat = curOut; 
            outi = outi+1; 
        else
            outputDat(outi) = curOut; 
            outi = outi+1;
        end
    end





end

save('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject\DatForYes.mat', 'outputDat')

