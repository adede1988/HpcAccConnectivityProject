function [] = makeConnectionPlot(plotMat, climVals, Xlabs, Ylabs, XtickLocs, YtickLocs)


imagesc(plotMat)
yticks(XtickLocs)
yticklabels(Xlabs)

xticks(YtickLocs)
xticklabels(Ylabs)
caxis(climVals)



end