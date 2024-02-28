%% pairwise connectivity file combine



codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';

%% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
datFolder = [datPre 'CHANDAT/finished']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 







%% loop over channels and get single trial level data for HFB and TF

regions = {'acc', 'dlPFC', 'hip', 'lTemp', 'iTemp', 'mtl', ...
    'pcc', 'pPFC', 'vis'}; 

%store a struct with file path info for all chan pairs 
% available on each connection


parfor ii = 1:length(chanFiles)
 ii

    chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
    lab = chanDat.labels{chanDat.chi,3}; 
    
    T = sum(chanDat.retInfo(:,1)==1 | chanDat.retInfo(:,1)==2); 
    Hr = sum(chanDat.retInfo(:,1)==1) / T; 
    T = sum(chanDat.retInfo(:,1)==3 | chanDat.retInfo(:,1)==4); 
    F = sum(chanDat.retInfo(:,1)==4);
    if F == 0
        acc =  Hr - F/T; 
        F = 1; 
    else
        acc =  Hr - F/T; 
    end
    Fr = F / T; 
    
    d = norminv(Hr) - norminv(Fr); 
    
    if acc>0 && chanDat.age > 13 %memory and age filter
        if sum(cellfun(@(x) strcmp(x, lab), regions)) == 1  &&...
                sum(chanDat.reactiveRes==1)>0 %in a target region and reactive
            
            %get the list of simultaneously recorded channels within
            %subject
            curChan = chanFiles(ii).name; 
            subID = split(curChan, '_'); 
            subID = subID{2}; 
            subFiles = dir([datPre 'CHANDAT/finished']);
            test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
            subFiles = subFiles(test);
            
            reg1 = find(cellfun(@(x) strcmp(x, lab), regions)); 

            %loop over the simultaneous channels to 
            for jj = 1:length(subFiles)
                chanDat2 = load([subFiles(jj).folder '/' subFiles(jj).name]).chanDat; 
                lab = chanDat2.labels{chanDat2.chi,3}; 
                
                T = sum(chanDat2.retInfo(:,1)==1 | chanDat2.retInfo(:,1)==2); 
                Hr = sum(chanDat2.retInfo(:,1)==1) / T; 
                T = sum(chanDat2.retInfo(:,1)==3 | chanDat2.retInfo(:,1)==4); 
                F = sum(chanDat2.retInfo(:,1)==4);
                if F == 0
                    acc =  Hr - F/T; 
                    F = 1; 
                else
                    acc =  Hr - F/T; 
                end
                Fr = F / T; 
                
                d = norminv(Hr) - norminv(Fr); 

                if sum(cellfun(@(x) strcmp(x, lab), regions)) == 1  &&...
                        sum(chanDat2.reactiveRes==1)>0 %in a target region and reactive
                    reg2 = find(cellfun(@(x) strcmp(x, lab), regions));
                    
                    outinfo = struct; 
                    outinfo.setFolder = chanFiles(ii).folder; 
                    outinfo.setName = chanFiles(ii).name; 
                    outinfo.pairFolder = subFiles(jj).folder; 
                    outinfo.pairName = subFiles(jj).name; 

                    parsave(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\pairFiles\'...
                        regions{reg1}, '_' regions{reg2} '_' ...
                        num2str(ii) '_' num2str(jj) '.mat'], ...
                        outinfo)
                  

                end
            



            end






        end

    end


end






