function [] = makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
    topt, imgopt)

if ~topt
    plotMat(plotMat<0) = 0; 
    plotMat = sqrt(plotMat); 
end
plotMat = plotMat([3,4,1,2,5], [3,4,1,2,5]);
if imgopt
    plotMat = triu(plotMat); 
end
figure('visible', false, 'position', [0,0,600,600])
imagesc(plotMat)
xticks([1:5])
yticks([1:5])
xticklabels(regions(keyRegIdx([3,4,1,2,5])))
yticklabels(regions(keyRegIdx([3,4,1,2,5])))
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
axis square
if topt
    caxis([-4, 4])
else
    caxis([.02, .2])
end
hold on 

pMat = pMat([3,4,1,2,5], [3,4,1,2,5]);
if imgopt
    pMat = tril(ones(size(pMat))) + pMat; 
end
for rowi = 1:5
    for coli = 1:5
        if(pMat(rowi, coli)<.05)
            scatter(coli, rowi, 600, 'red', 'filled')
        end
    end
end

end