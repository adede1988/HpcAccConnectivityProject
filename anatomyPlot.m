function [inArea, noteInArea] = anatomyPlot(regModels, elecpos, labels, roiLab, rois, subID, flipVal, plotIt)






% MNI TEMPLATE BRAIN: 
% https://nist.mni.mcgill.ca/mni-average-brain-305-mri/
%points in a volume: 
% https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume


colors = [[48,155,255];...
          [249,93,167];...
          [186,232,74];...
          [255,253,81]]./255;



inArea = zeros(size(elecpos,1), length(regModels)); 
noteInArea = inArea; 
elecpos(:,2) = elecpos(:,2)*flipVal;

if plotIt==1
    figure
    hold on
end
for ii = 1:length(regModels)
    
    
    TR = regModels{ii}; 
     
    inArea(:,ii) = inpolyhedron(TR.ConnectivityList, TR.Points, elecpos, 'TOL', 5); 
    noteInArea(:,ii) = cellfun(@(x) strcmp(x, rois{ii}), roiLab );
    
    if plotIt == 1
        trimesh(regModels{ii}, 'facecolor', 'none', 'facealpha', .00, 'edgecolor', colors(ii,:), 'linewidth', .01, 'linestyle', ':')
        scatter3(elecpos(inArea(:,ii)==1,1), elecpos(inArea(:,ii)==1,2), elecpos(inArea(:,ii)==1,3), 100, colors(ii,:), 'filled')
        scatter3(elecpos(noteInArea(:,ii)==1,1), elecpos(noteInArea(:,ii)==1,2), elecpos(noteInArea(:,ii)==1,3), 100, 'k', '*')
        
       
        if sum(noteInArea(:,ii)) == 1
        text(elecpos(noteInArea(:,ii)==1,1), elecpos(noteInArea(:,ii)==1,2), elecpos(noteInArea(:,ii)==1,3), {labels{noteInArea(:,ii)==1}})
        elseif sum(noteInArea(:,ii)) > 1
        text(elecpos(noteInArea(:,ii)==1,1), elecpos(noteInArea(:,ii)==1,2), elecpos(noteInArea(:,ii)==1,3), {labels{noteInArea(:,ii)==1}})
        end
        if sum(inArea(:,ii))==1
        text(elecpos(inArea(:,ii)==1,1), elecpos(inArea(:,ii)==1,2), elecpos(inArea(:,ii)==1,3), [labels{inArea(:,ii)==1}])
        elseif sum(inArea(:,ii))>1
        text(elecpos(inArea(:,ii)==1,1), elecpos(inArea(:,ii)==1,2), elecpos(inArea(:,ii)==1,3), {labels{inArea(:,ii)==1}})
        end
    end

end
if plotIt == 1
    view([70,45])
    axis([-70, 70, -50, 50, -50, 80])
    hold on 
    scatter3(elecpos(:,1), elecpos(:,2), elecpos(:,3), 'k')
    % text(subDat.elecpos(:,1), subDat.elecpos(:,2)*-1, subDat.elecpos(:,3), {subDat.labels{:,1}})
    title(subID)
end




end