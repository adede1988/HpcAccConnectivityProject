
clear
% gathering all data to do group level summary with ROI and epoch
% specificity

rois = {'dlPFC', 'hip', 'phg', 'acc'}; 
datFolder  = "C:\Users\dtf8829\Documents\QuestConnect\SUMDAT";


subFiles = dir(datFolder); 
subFiles = subFiles([subFiles.isdir]==false); 

%create allDat struct
allDat = struct; 

for sub = 1:length(subFiles)
    sub
    subDat = load([subFiles(sub).folder '/' subFiles(sub).name]).subDat; 

    %get the necessary meta data to track the subject
    allDat(sub).site = subDat.site; 
    allDat(sub).subID = subDat.subID; 
    allDat(sub).age = subDat.age; 
    allDat(sub).sex = subDat.sex{1}; 
    allDat(sub).use = subDat.use; 
    allDat(sub).hits = subDat.hits;
    allDat(sub).misses = subDat.misses; 
    allDat(sub).retInfo = subDat.retInfo; 
    allDat(sub).elecPos = subDat.elecpos; 
    allDat(sub).roimni = subDat.roimni; 
    allDat(sub).roiNote = subDat.roiNote; 
    
    %check that it's through the pipeline
    if isfield(subDat, 'sizeReduce')
       
        %use the note based roi matrix as a default
        if sum(sum(subDat.roiNote))==0
            roiMat = subDat.roimni; 
        else
            roiMat = subDat.roiNote; 
        end
        
        allDat(sub).encepoch = subDat.encepoch; 
        allDat(sub).onepoch = subDat.onepoch; 
        allDat(sub).rtepoch = subDat.rtepoch; 

        %get all the timefrequency stuff
        for rr = 1:size(roiMat,2)
            fields = fieldnames(subDat.TFout); 
            for fi = 1:length(fields)
                cur = subDat.TFout.(fields{fi}); 
                cur = cur(:,:,subDat.chanorder); 
                cur = cur(:,:,roiMat(:,rr)==1); 
                if ~isempty(cur)
                allDat(sub).([fields{fi} '_' rois{rr}]) = mean(cur,3);
                end
            end
        end

        %get all the connectivity stuff
        for rr = 1:size(roiMat,2)
            for rr2 = 1:size(roiMat,2)
                if rr>=rr2
                fields = fieldnames(subDat.ISPCout); 
                for fi = 1:length(fields)
                    cur = subDat.ISPCout.(fields{fi}); 
                    cur = cur(subDat.chanorder,:, :,:,subDat.chanorder); 
                    cur = cur(roiMat(:,rr2)==1,:,:,:,roiMat(:,rr)==1); 
                    if ~isempty(cur)
                    allDat(sub).([fields{fi} '_' rois{rr} '_' rois{rr2}]) = squeeze(mean(cur,[1,5]));
                    end
                end
                end
            end
        end


    end


end




%between area encoding connectivity: 
connectionPlot(allDat, "subHit_acc_hip", "subMiss_acc_hip", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_phg_dlPFC", "subMiss_phg_dlPFC", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_acc_dlPFC", "subMiss_acc_dlPFC", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_phg_hip", "subMiss_phg_hip", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_acc_phg", "subMiss_acc_phg", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_hip_dlPFC", "subMiss_hip_dlPFC", 2, allDat(1).encepoch)

%within area encoding connectivity: 
connectionPlot(allDat, "subHit_acc_acc", "subMiss_acc_acc", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_phg_phg", "subMiss_phg_phg", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_hip_hip", "subMiss_hip_hip", 2, allDat(1).encepoch)

connectionPlot(allDat, "subHit_dlPFC_dlPFC", "subMiss_dlPFC_dlPFC", 2, allDat(1).encepoch)





%between area response locked connectivity: 
connectionPlot(allDat, "hit_rt_acc_hip", "miss_rt_acc_hip", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_phg_dlPFC", "miss_rt_phg_dlPFC", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_acc_dlPFC", "miss_rt_acc_dlPFC", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_phg_hip", "miss_rt_phg_hip", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_acc_phg", "miss_rt_acc_phg", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_hip_dlPFC", "miss_rt_hip_dlPFC", 2, allDat(1).rtepoch)

%within area response locked connectivity: 
connectionPlot(allDat, "hit_rt_acc_acc", "miss_rt_acc_acc", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_phg_phg", "miss_rt_phg_phg", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_hip_hip", "miss_rt_hip_hip", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_dlPFC_dlPFC", "miss_rt_dlPFC_dlPFC", 2, allDat(1).rtepoch)


%between area retrieval image locked connectivity: 
connectionPlot(allDat, "hit_on_acc_hip", "miss_on_acc_hip", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_phg_dlPFC", "miss_on_phg_dlPFC", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_acc_dlPFC", "miss_on_acc_dlPFC", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_phg_hip", "miss_on_phg_hip", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_acc_phg", "miss_on_acc_phg", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_hip_dlPFC", "miss_on_hip_dlPFC", 2, allDat(1).onepoch)

%within area retrieval image locked connectivity: 
connectionPlot(allDat, "hit_on_acc_acc", "miss_on_acc_acc", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_phg_phg", "miss_on_phg_phg", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_hip_hip", "miss_on_hip_hip", 2, allDat(1).onepoch)

connectionPlot(allDat, "hit_on_dlPFC_dlPFC", "miss_on_dlPFC_dlPFC", 2, allDat(1).onepoch)


















connectionPlot(allDat, "subHit_phg_dlPFC", "subMiss_phg_dlPFC", 2, allDat(1).encepoch)

connectionPlot(allDat, "hit_rt_acc_hip", "cr_rt_acc_hip", 2, allDat(1).rtepoch)

connectionPlot(allDat, "hit_rt_phg_dlPFC", "cr_rt_phg_dlPFC", 2, allDat(1).rtepoch)


% 
% 
% figure
% subplot 311
% targField = "subHit_acc_hip"; 
% epoch = allDat(1).encepoch;
% 
% test = find(cellfun(@(x) ~isempty(x), {allDat.(targField)}));
% 
% temp = allDat(test); 
% 
% temp = {temp.(targField)}; 
% out = temp{1}; 
% dimVal = length(size(out)); 
% for ii = 2:length(temp)
%     out = cat(dimVal+1, out, temp{ii} );
% end
% 
% out = mean(out, dimVal+1); 
% 
% imagesc(out(:,:,2)')
% set(gca, 'YDir', 'normal')
% caxis([0,.05])
% colorbar
% xticks([1.5:2:length(epoch)])
% xticklabels(epoch([2:2:end]))
% yticks([20:20:100])
% yticklabels(round(frex([20:20:100])))
% 
% 
% subplot 312
% targField = "subMiss_acc_hip"; 
% 
% 
% test = find(cellfun(@(x) ~isempty(x), {allDat.(targField)}));
% 
% temp = allDat(test); 
% 
% temp = {temp.(targField)}; 
% out = temp{1}; 
% dimVal = length(size(out)); 
% for ii = 2:length(temp)
%     out = cat(dimVal+1, out, temp{ii} );
% end
% 
% out = mean(out, dimVal+1); 
% 
% imagesc(out(:,:,2)')
% set(gca, 'YDir', 'normal')
% caxis([0,.05])
% colorbar
% xticks([1.5:2:length(epoch)])
% xticklabels(epoch([2:2:end]))
% yticks([20:20:100])
% yticklabels(round(frex([20:20:100])))








