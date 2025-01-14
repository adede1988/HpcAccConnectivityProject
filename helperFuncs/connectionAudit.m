function [] = connectionAudit(curDat, ...
                regions, phase, stat, datPre, fn)

reg1 = find(cellfun(@(x) strcmp(x, curDat.reg1), regions)); 
reg2 = find(cellfun(@(x) strcmp(x, curDat.reg2), regions));
chanFiles = dir([datPre 'CHANDAT/finished/']);


if isfield(curDat.([phase '_' stat '_clust']), 'pos_clusters')
for cc = 1:length(curDat.([phase '_' stat '_clust']).pos_clusters)
if curDat.([phase '_' stat '_clust']).pos_clusters(cc).p < .05 %sig check
   

   
    for pp = 1:length(curDat.subVals) 
       outFile = struct; 
       outFile.figFile = fn;
       outFile.encRet = phase; 
       outFile.HFBimage = stat; 
       outFile.cc = cc; 
       outFile.pp = pp; 
       subID = curDat.subVals{pp};

       idx = cellfun(@(x) contains(x, subID), ...
        {chanFiles.name}); 
       fileLen = sum(idx); 

        if fileLen < 50

           save([datPre 'graphAnalysis/SHORT_' regions{reg1} '_' regions{reg2} '_'...
               phase '_' stat '_' ...
               num2str(cc) '_' num2str(pp) '.mat'], 'outFile')
        
        elseif fileLen <100 
            
           save([datPre 'graphAnalysis/MED_' regions{reg1} '_' regions{reg2} '_'...
               phase '_' stat '_' ...
               num2str(cc) '_' num2str(pp) '.mat'], 'outFile')

        elseif fileLen <150

           save([datPre 'graphAnalysis/LONG_' regions{reg1} '_' regions{reg2} '_'...
               phase '_' stat '_' ...
               num2str(cc) '_' num2str(pp) '.mat'], 'outFile')

        else

           save([datPre 'graphAnalysis/XL_' regions{reg1} '_' regions{reg2} '_'...
               phase '_' stat '_' ...
               num2str(cc) '_' num2str(pp) '.mat'], 'outFile')
        end
            

    end
end 
end 
end 



end