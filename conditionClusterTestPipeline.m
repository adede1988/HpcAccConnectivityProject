function [] = conditionClusterTestPipeline(fileInfo)

%load in the condition data
cndDat = load([fileInfo.folder '/' fileInfo.name]).cndDat;
frex = logspace(log10(2),log10(80),100);
numTim = size(cndDat.datForPerm,1); 
numFrex = size(cndDat.datForPerm,2); 


%check if it's connectivity data and select variable!
datForPerm = cndDat.datForPerm; 
if length(size(datForPerm))==4
    datForPerm = squeeze(datForPerm(:,:,3,:)); 
end

%the data for permutation testing are already available

    tObs =   reshape(cell2mat(arrayfun(@(x) ...
                    arrayfun(@(y) ...
                        getTstat(datForPerm(x,y,cndDat.condCode==1), datForPerm(x,y,cndDat.condCode==0)), ...
                    1:numFrex),...
                1:numTim, 'uniformoutput', false)), [ numFrex,numTim]);


   perms = zeros([size(tObs),1000]); 

          for p = 1:1000
                curShuf = randsample(cndDat.condCode, length(cndDat.condCode), false); 
                perms(:,:,p) =   reshape(cell2mat(arrayfun(@(x) ...
                            arrayfun(@(y) ...
                                getTstat(datForPerm(x,y,curShuf==1), datForPerm(x,y,curShuf==0)), ...
                            1:numFrex),...
                        1:numTim, 'uniformoutput', false)), [numFrex,numTim]);


          end


    [h, p, clusterinfo] = cluster_test(tObs, perms); 

    cndDat.h = h; 
    cndDat.p = p; 
    cndDat.clusterinfo = clusterinfo; 

    save([fileInfo.folder '/' fileInfo.name], 'cndDat')

%figure out the correct epoch time to use: 
if contains(cndDat.targVar, 'ub')
    epoch = cndDat.encepoch; 
elseif contains(cndDat.targVar, 'rt')
    epoch = cndDat.rtepoch; 
else
    epoch = cndDat.onepoch;
end



    figure('visible', false)
    subplot 231
    plotVals = mean(datForPerm(:,:,cndDat.condCode==1),3)';
    imagesc(plotVals)
    sp1Max = max(max(plotVals));
    sp1Min = min(min(plotVals)); 
    set(gca, 'ydir', 'normal')
    colorbar
    xticks([1.5:2:length(epoch)])
    xticklabels(epoch([2:2:end]))
    yticks([10:10:60])
    yticklabels(round(frex([10:10:60])))
    ylim([1, 60])
    title(cndDat.targVar, 'interpreter', 'none')


    subplot 232
    plotVals2 = mean(cndDat.datForPerm(:,:,cndDat.condCode==0),3)';
    imagesc(plotVals2)
    sp2Max = max(max(plotVals2));
    sp2Min = min(min(plotVals2)); 
    set(gca, 'ydir', 'normal')
    colorbar
    xticks([1.5:2:length(epoch)])
    xticklabels(epoch([2:2:end]))
    yticks([10:10:60])
    yticklabels(round(frex([10:10:60])))
    ylim([1, 60])
    title(cndDat.pairVar, 'interpreter', 'none')

    allMax = max(abs([sp1Max, sp2Max])); 
    allMin = min([sp1Min, sp2Min]); 
    subplot 231
    clim([allMin allMax])
    subplot 232
    clim([allMin allMax])

    subplot 233
    imagesc(plotVals - plotVals2)
    set(gca, 'ydir', 'normal')
    colorbar
    xticks([1.5:2:length(epoch)])
    xticklabels(epoch([2:2:end]))
    yticks([10:10:60])
    yticklabels(round(frex([10:10:60])))
    ylim([1, 60])
    sp3Max = min(min(plotVals-plotVals2)); 
    sp3Min = max(max(plotVals - plotVals2)); 
    allMax = max(abs([sp3Max, sp3Min])); 
    clim([-allMax allMax])
    title('difference')

    subplot 234
    imagesc(tObs)
    set(gca, 'ydir', 'normal')
    colorbar
    xticks([1.5:2:length(epoch)])
    xticklabels(epoch([2:2:end]))
    yticks([10:10:60])
    yticklabels(round(frex([10:10:60])))
    ylim([1, 60])
    sp4Max = min(min(tObs)); 
    sp4Min = max(max(tObs)); 
    allMax = max(abs([sp4Max, sp4Min])); 
    clim([-allMax allMax])
    title('uncorrect t-value')

    subplot 235
    imagesc(h)
    set(gca, 'ydir', 'normal')
    colorbar
    xticks([1.5:2:length(epoch)])
    xticklabels(epoch([2:2:end]))
    yticks([10:10:60])
    yticklabels(round(frex([10:10:60])))
    ylim([1, 60])
    title('cluster significance')

    subplot 236
    imagesc(p)
    set(gca, 'ydir', 'normal')
    colorbar
    xticks([1.5:2:length(epoch)])
    xticklabels(epoch([2:2:end]))
    yticks([10:10:60])
    yticklabels(round(frex([10:10:60])))
    ylim([1, 60])
    title('cluster p-value')
    clim([0 .2])


    set(gcf,'color','w');


    set(gcf, 'Position',  [100, 100, 1500, 800]);

    export_fig(join([fileInfo.folder '/' cndDat.targVar '.jpg'],''), '-r300')
end