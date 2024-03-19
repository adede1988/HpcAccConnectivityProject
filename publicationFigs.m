%% MTL - PFC publication figures 

%description of the general data set 


%% set environment

%figures should not be docked: 
set(0, 'defaultfigurewindowstyle', 'normal')
%local paths: 


codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';
figDat = 'R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\';

% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])

regions = {'acc', 'dlPFC', 'hip', ...
    'lTemp', 'iTemp', 'mtl', 'pcc', 'pPFC', 'vis'}; 




% get all the file paths: 
%it's assumed here that the TF quest pipeline has already been run and that
%the outputs are available 

%TF power files
outStatFiles = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'all.mat'));
outStatFiles(~test) = []; 
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFiles.name}, 'phase'));
outStatFiles(test) = []; 
headFiles = dir([datPre 'TF_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFiles.name}, 'all.mat'));
headFiles(~test) = []; 

%TF phase files
outStatFilesPhase = dir([datPre 'TF_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesPhase.name}, 'NEWphase'));
outStatFilesPhase(~test) = []; 

%HFB files
outStatFilesHFB = dir([datPre 'HFB_singleTrial/out']);
test = cellfun(@(x) length(x)>0, ...
    strfind({outStatFilesHFB.name}, '.mat'));
outStatFilesHFB(~test) = []; 
headFilesHFB = dir([datPre 'HFB_singleTrial']);
test = cellfun(@(x) length(x)>0, ...
    strfind({headFilesHFB.name}, '.mat'));
headFilesHFB(~test) = []; 


%% set graphical parameters: 

hitCol = [88, 113, 206] ./ 256; 
missCol = [240, 172, 99] ./ 256; 
sigCol = [.5,.5,.5]; 
sigAlpha = .3; 
errAlpha = .2; 

regColors = [[204, 153, 204]', ...%ACC: dusky pink
             [145, 162, 80]', ... %dlPFC: light forest green
             [61, 187, 218]', ... %hip light blue
             [145, 162, 80]', ... %iTemp: light forest green
             [145, 162, 80]', ... %lTemp: light forest green
             [214, 77, 97]',...   %parahip: dark pink
             [132, 125, 212]', ...%visual: lilac
             [51, 102, 153]', ... %polarPFC: dark blue
             [184, 78, 83]']' ./ 255;  %PCC: dark red

keyRegIdx = [1,2,3,6,8]; 

b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], ...
    [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)';...
    linspace(255,0,128)']]/255;
b2w2r(129:end, 1) = 1; 
b2w2r(1:128, 3) = 1; 

linWid = 5; 

s = [12, 61, 74]; %sangria
% s = [77, 0, 77]; 
m = [171,189,154]; 
e = [243, 188, 46]; %scarlet
e = [255, 255, 46]; 

s2w2y = [[linspace(s(1),m(1),128)'; linspace(m(1),e(1),128)'], ...
         [linspace(s(2),m(2),128)'; linspace(m(2),e(2),128)'], ...
         [linspace(s(3),m(3),128)'; linspace(m(3),e(3),128)'], ...
         ] / 255;

phaseVals = {'sub', 'ret'}; 

figure
hold on 
for ii = 1:9
    scatter(ii, 2, 100, regColors(ii,:), 'filled')
    text(ii, 2.5, regions{ii})

end
%% Figure 1: 
%the amplitude of the HFB response was compared between subsequent hit and
%miss trials during encoding and between hit and miss trials during
%retrieval. This analysis was carried out using linear mixed effects
%modeling treating channel and subject as random effects and hit/miss as
%the only fixed effect. Analysis was carried out at each time point from
%-450ms before image onset to 2500ms after image onset. Permutation-based
%cluster correction for multiple comparisons was carried out using 5000
%permutations. At encoding, significant differences between subsequent hit 
%subsequent miss trials were revealed in the hippocampus at -450:-375ms
%with hits greater than misses and at 700:775ms again with hits greater
%than misses. A significant difference was also revealed at 

fig1Dat = dir([figDat 'Figure1']); 
fig1Dat(1:2) = []; 

stidx = cellfun(@(x) sum(strfind( x,'singleTrial'))>0, {fig1Dat.name});
singleFig1Dat = fig1Dat(stidx); 
fig1Dat = fig1Dat(~stidx); 
fig1Dat([fig1Dat.isdir]) = []; 

%% basic image locked HFB time courses

%while looping, grab the distribution of effect sizes for later comparison
%to HFB-peak diffs
%effect sizes stored as time X enc/ret X region
%panel 1A
% imgDist = zeros(139, 2, 9); 
parfor ii = 1:length(fig1Dat)
    panelDat = load([fig1Dat(ii).folder '/' fig1Dat(ii).name]).outDat; 

    %quick cleaning for plotting purposes only, all data were used in stats
    panelDat.hitVals(panelDat.eliminate,:) = []; 
    panelDat.missVals(panelDat.eliminate,:) = []; 
%     if strcmp(panelDat.reg, 'acc')
%         panelDat.hitVals(9:10, :) = []; 
%         panelDat.missVals(9:10, :) = []; 
%     end
    panelDat.hitVals(isnan(panelDat.hitVals)) = 0; 
    panelDat.missVals(isnan(panelDat.missVals)) = 0; 

    %get effect sizes
    diffs = panelDat.hitVals - panelDat.missVals; 
    diffs = mean(diffs) ./ std(diffs); 
    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

%     imgDist(:,phasei, regi) = diffs; 
    figure('visible', false, 'position', [0,0,600,600])
    x = panelDat.tim; 
    x(panelDat.p>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
    breakVals = [0, find(diff(x)> 25), length(x)];
    for jj = 1:length(breakVals)-1
        tmpX = x(breakVals(jj)+1:breakVals(jj+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
    
    
    end
    else
        hold on 
    end
    yline(0, 'color', 'k', 'linewidth', linWid)
    xline(0, '--', 'linewidth', linWid, 'color', 'k')
    xline(panelDat.meanHitRT, 'color', hitCol, ...
        'linewidth', linWid, 'linestyle', '--')
    xline(panelDat.meanMissRT, 'color', missCol, ...
        'linewidth', linWid, 'linestyle', '--')
    y = movmean(mean(panelDat.hitVals), 3); 
    x = panelDat.tim; 
    plot(x, y, 'color', hitCol, 'linewidth', 4)
    se = std(panelDat.hitVals) ./ sqrt(size(panelDat.hitVals,1)); 
    y = [y + se, flip(y) - flip(se)]; 
    x = [x, flip(x)]; 
    hold on 
    fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max(y); 
    allMin = min(y); 
    
    y = movmean(mean(panelDat.missVals), 3); 
    x = panelDat.tim; 
    plot(x, y, 'color', missCol, 'linewidth', 4)
    se = std(panelDat.missVals) ./ sqrt(size(panelDat.missVals,1)); 
    y = [y + se, flip(y) - flip(se)]; 
    x = [x, flip(x)]; 
    hold on 
    fill(x, y, missCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max([allMax, max(y)]); 
    allMin = min([allMin, min(y)]); 
    ylim([allMin*1.1, allMax*1.1])
    ylim([-2, 6])
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([-450, 3000])
    export_fig([figDat 'pubFigs/' 'Fig1_' panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')
    
end

%% single trial heatmap figures
% for ii = 1:length(singleFig1Dat)
%hard code for hippocampus at encoding hits only
%panel 1B
    ii = 6;
    panelDat = load([singleFig1Dat(ii).folder '/'...
        singleFig1Dat(ii).name]).outDat; 
    
    figure('visible', true, 'position', [0,0,600,1000])
    [sortedLat, order] = sort(panelDat.hitLat); 
    imagesc(panelDat.tim, [], panelDat.hits(:,order)')
    colormap(s2w2y)
    clim([-10,20])
    colorbar
    hold on 
    xline(0, '--', 'linewidth', linWid)
    scatter(panelDat.hitRT(order), [1:length(panelDat.hitRT)], ...
        20, 'white', 'filled')
%     plot(sortedLat, [1:length(sortedLat)], 'color', regColors(9,:), ...
%         'linestyle', '--', 'linewidth', linWid)
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    yticks([])

    export_fig([figDat 'pubFigs/' 'Fig1_singleTrial' ...
        panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')


% end
%% get the latencies for all regions on a single plot HITS Retrieve
% linew = 5;
retFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])
%extract
for ii = 1:length(singleFig1Dat)
     panelDat = load([singleFig1Dat(ii).folder '/'...
            singleFig1Dat(ii).name]).outDat; 
     
    if strcmp(panelDat.phase, 'sub')

    else
        set(0, 'currentFigure', retFig); 
        if ismember(panelDat.reg, regions(keyRegIdx))
            coli = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
            hold on 
            hitLat = arrayfun(@(x) prctile(panelDat.hitLat ...
                ./ panelDat.hitRT, x), [1:1:100]);
%                 
            plot(hitLat, 1:100, 'color', [regColors(coli,:), .75], ...
                'linewidth', linew)
        
        end
    end

   

end

set(0, 'currentfigure', retFig);
export_fig([figDat 'pubFigs/' 'Fig1_retHitLatencies' ...
    '.jpg'], '-r300')

%% get the latencies for all regions on a single plot HITS Encode
encFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])

%extract
for ii = 1:length(singleFig1Dat)
     panelDat = load([singleFig1Dat(ii).folder '/'...
            singleFig1Dat(ii).name]).outDat; 
     
    if strcmp(panelDat.phase, 'sub')
        set(0, 'currentfigure', encFig);
         %check for target region and plot
        if ismember(panelDat.reg, regions(keyRegIdx))
            coli = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
            hold on 
            hitLat = arrayfun(@(x) prctile(panelDat.hitLat ...
                    ./ panelDat.hitRT, x), [1:1:100]);
%                 
            plot(hitLat, 1:100, 'color', [regColors(coli,:),.75], ...
                'linewidth', linew)
        
        end
    else
       
    end

   

end
set(0, 'currentfigure', encFig);
export_fig([figDat 'pubFigs/' 'Fig1_encHitLatencies' ...
    '.jpg'], '-r300')


%% get the latencies for all regions on a single plot MISSES ENCODE
encFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])

%extract
for ii = 1:length(singleFig1Dat)
     panelDat = load([singleFig1Dat(ii).folder '/'...
            singleFig1Dat(ii).name]).outDat; 
     
    if strcmp(panelDat.phase, 'sub')
        set(0, 'currentfigure', encFig);
         %check for target region and plot
        if ismember(panelDat.reg, regions(keyRegIdx))
            coli = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
            hold on 
            missLat = arrayfun(@(x) prctile(panelDat.missLat ...
                ./ panelDat.missRT,x), [1:1:100]);
%                 , x)
            plot(missLat, 1:100, 'color', [regColors(coli,:), .75], ...
                'linewidth', linew)
        
        end
    else
%         set(0, 'currentFigure', retFig); 
%             if ismember(panelDat.reg, regions(keyRegIdx))
%             coli = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
%             hold on 
%             missLat = arrayfun(@(x) prctile(panelDat.missLat ...
%                 ./ panelDat.missRT, x), [1:1:100]);
%             plot(missLat, 1:100, 'color', [regColors(coli,:), .75], ...
%                 'linewidth', 9)
%         
%         end
    end

   

end
set(0, 'currentfigure', encFig);
export_fig([figDat 'pubFigs/' 'Fig1_encMissLatencies' ...
    '.jpg'], '-r300')
% set(0, 'currentfigure', retFig);
% export_fig([figDat 'pubFigs/' 'Fig1_retMissLatencies' ...
%     '.jpg'], '-r300')

%% get the latencies for all regions on a single plot MISSES RETRIEVE

retFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])
%extract
for ii = 1:length(singleFig1Dat)
     panelDat = load([singleFig1Dat(ii).folder '/'...
            singleFig1Dat(ii).name]).outDat; 
     
    if strcmp(panelDat.phase, 'sub')
  
    else
        set(0, 'currentFigure', retFig); 
            if ismember(panelDat.reg, regions(keyRegIdx))
            coli = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
            hold on 
            missLat = arrayfun(@(x) prctile(panelDat.missLat ...
                ./panelDat.missRT, x), [1:1:100]);
            plot(missLat, 1:100, 'color', [regColors(coli,:), .75], ...
                'linewidth', linew)
        
        end
    end

   

end

set(0, 'currentfigure', retFig);
export_fig([figDat 'pubFigs/' 'Fig1_retMissLatencies' ...
    '.jpg'], '-r300')

%% get the data ready for plotting the mean latency for each region
%need to get all the latency values out
%regions X hit/miss X enc/ret
allLats = cell(9, 2, 2); 
for ii = 1:length(singleFig1Dat)
     panelDat = load([singleFig1Dat(ii).folder '/'...
            singleFig1Dat(ii).name]).outDat; 
    coli = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    if strcmp(panelDat.phase, 'sub')
        allLats{coli,1,1} = panelDat.hitLat ./ panelDat.hitRT; 
        allLats{coli,2,1} = panelDat.missLat ./ panelDat.missRT;; 
    else
        allLats{coli,1,2} = panelDat.hitLat ./ panelDat.hitRT;
        allLats{coli,2,2} = panelDat.missLat ./ panelDat.missRT; 
    end


end


for hm = 1:2
    for er = 1:2
        folderLoc = [singleFig1Dat(1).folder '/Rbox/']; 
        tmp = allLats(keyRegIdx, hm, er); 
        maxLength = max(cellfun(@numel, tmp));
        paddedData = cellfun(@(x) ...
            [x; nan(maxLength - numel(x), 1)], ...
            tmp, 'UniformOutput', false);
        
        % Convert cell array to matrix
        dataMatrix = cat(2, paddedData{:});
        save([folderLoc 'hitMiss_' num2str(hm)...
            '_encRet_' num2str(er) '.mat'], 'dataMatrix')

    end
end


tmp = allLats(keyRegIdx, 1, 1); 
maxLength = max(cellfun(@numel, tmp));
paddedData = cellfun(@(x) ...
    [x; nan(maxLength - numel(x), 1)], ...
    tmp, 'UniformOutput', false);

% Convert cell array to matrix
dataMatrix = cat(2, paddedData{:});
boxplot(dataMatrix, ...
    'Labels', regions(keyRegIdx));
title('Boxplot for Each Cell in the Cell Array');
ylabel('Values');





%% Figure 2: 
%The amplitude difference measured at the time point of the HFB peak.   

fig2Dat = dir([figDat 'Figure2']); 
fig2Dat(1:2) = []; 


%% power spectra at HFB peak

test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TF_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_HFB'));
TF_files_HFB = TF_files(test); 





parfor ii = 1:length(TF_files_HFB)
    try
    panelDat = load([TF_files_HFB(ii).folder '/' ...
        TF_files_HFB(ii).name]).outDat;
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    figure('visible', false, 'position', [0,0,600,600])
    x = 1:length(panelDat.frex); 
%     yline(0, 'color', 'k', 'linewidth', linWid)
    x(panelDat.p_hfb(21,:)>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
    breakVals = [0, find(diff(x)> 25), length(x)];
    for jj = 1:length(breakVals)-1
        tmpX = x(breakVals(jj)+1:breakVals(jj+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
%         yline(0, 'color', 'k', 'linewidth', linWid)
    
    end
    else
        hold on 
    end


%     hitSpect = arrayfun(@(x) find(panelDat.tim>=x, 1), panelDat.hitLat);
%     hitSpect = arrayfun(@(x,y) squeeze(panelDat.hits(x,y,:)), ...
%         hitSpect, [1:length(hitSpect)]', 'uniformoutput', false);
%     hitSpect = cat(2, hitSpect{:});
%     % need to take channel-wise means
%     realID = arrayfun(@(x) [panelDat.hitSub{x} '_'...
%         num2str(panelDat.hitChi(x))], 1:length(panelDat.hitChi), ...
%         'UniformOutput', false);
%     uniIDs = unique(realID); 
%     hitSpect = cellfun(@(x) ...
%         mean(hitSpect(:, ismember( realID,x)), 2), uniIDs, ...
%         'uniformoutput', false);
%     hitSpect = cat(2, hitSpect{:});
    hitSpect = squeeze(panelDat.hits_hfb(:,21,:)); 


    if regi == 9 %outliers in the hippocampus retrieval data, remove for plot
        hitSpect(8:10,:) = []; 
    end
    
%     
% 
%     
%     hitBase = permute(panelDat.hits, [3,1,2]); 
%     hitBase = reshape(hitBase, [100, prod(size(hitBase, [2,3]))]);
% 
%     ybase = median(hitBase, 2); 
% %     plot(ybase, 'color', hitCol, 'linewidth', 4,'linestyle', '--')

    y = mean(hitSpect); 
    plot( y, 'color', hitCol, 'linewidth', 4)
    se = std(hitSpect) ./ sqrt(size(hitSpect,1)); 
    y = [y + se, flip(y) - flip(se)]; 
    x = [[1:100], flip([1:100])]; 
    hold on 
    fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max(y); 
    allMin = min(y);  
    
%     missSpect = arrayfun(@(x) find(panelDat.tim>=x, 1), panelDat.missLat);
%     missSpect = arrayfun(@(x,y) squeeze(panelDat.misses(x,y,:)), ...
%          missSpect, [1:length( missSpect)]', 'uniformoutput', false);
%     missSpect = cat(2,  missSpect{:});
% 
%     realID = arrayfun(@(x) [panelDat.missSub{x} '_'...
%         num2str(panelDat.missChi(x))], 1:length(panelDat.missChi), ...
%         'UniformOutput', false);
%     uniIDs = unique(realID); 
%     missSpect = cellfun(@(x) ...
%         mean(missSpect(:, ismember( realID,x)), 2), uniIDs, ...
%         'uniformoutput', false);
%     missSpect = cat(2, missSpect{:});
    missSpect = squeeze(panelDat.misses_hfb(:,21,:)); 

     if regi == 9 %outliers in the hippocampus retrieval data, remove for plot
        missSpect(8:10,:) = []; 
    end

%     missBase = permute(panelDat.misses, [3,1,2]); 
%     missBase = reshape(missBase, [100, prod(size(missBase, [2,3]))]);
% 
%     ybase = mean(missBase, 2); 
%     plot(ybase, 'color', missCol, 'linewidth', 4,'linestyle', '--')
    y = mean(missSpect);  
    plot( y, 'color', missCol, 'linewidth', 4)
    se = std(missSpect) ./ sqrt(size(missSpect,1));  
    y = [y + se, flip(y) - flip(se)]; 
    x = [[1:100], flip([1:100])]; 
    hold on 
    fill(x, y, missCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max([allMax, max(y)]); 
    allMin = min([allMin, min(y)]); 
    allMin = min([allMin, -2]); 
    ylim([-1, 18])
    xlim([1, 100])
    xticks([1:11:100])
    xticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'Fig2_' panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')
    catch
        ii
    end
end



%% ITPC spectra at HFB peak 

test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TFphase_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_HFB'));
TF_files_HFB = TF_files(test); 





parfor ii = 1:length(TF_files_HFB)
    try
    panelDat = load([TF_files_HFB(ii).folder '/' ...
        TF_files_HFB(ii).name]).outDat;
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    figure('visible', false, 'position', [0,0,600,600])
    x = 1:length(panelDat.frex); 
    yline(0, 'color', 'k', 'linewidth', linWid)
   
    x(panelDat.p_hfb(21,:)>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
    breakVals = [0, find(diff(x)> 25), length(x)];
    for jj = 1:length(breakVals)-1
        tmpX = x(breakVals(jj)+1:breakVals(jj+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
        xline(12, 'color', 'k', 'linewidth',2, 'linestyle', '--')
        xline(33, 'color', 'k', 'linewidth',2, 'linestyle', '--')
    
    end
    else
        hold on 
    end


    hitSpect = squeeze(panelDat.hits_hfb(:,21,:)); 



    


    y = mean(hitSpect); 
    plot( y, 'color', hitCol, 'linewidth', 4)
    se = std(hitSpect) ./ sqrt(size(hitSpect,1)); 
    y = [y + se, flip(y) - flip(se)]; 
    x = [[1:100], flip([1:100])]; 
    hold on 
    fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max(y); 
    allMin = min(y);  
    

    missSpect = squeeze(panelDat.misses_hfb(:,21,:)); 



    y = mean(missSpect);  
    plot( y, 'color', missCol, 'linewidth', 4)
    se = std(missSpect) ./ sqrt(size(missSpect,1));  
    y = [y + se, flip(y) - flip(se)]; 
    x = [[1:100], flip([1:100])]; 
    hold on 
    fill(x, y, missCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max([allMax, max(y)]); 
    allMin = min([allMin, min(y)]); 
    allMin = min([allMin, 0]); 
    ylim([0, 10])
    xlim([1, 100])
    xticks([1:11:100])
    xticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'Fig2_' panelDat.reg '_' ...
        panelDat.phase  '_phase.jpg'], '-r300')

     figure('visible', false, 'position', [0,0,600,600])
     hold off
     polarhistogram(panelDat.hit_angles, 18, ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability');
     hold on 
     polarhistogram(panelDat.miss_angles, 18, ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability');
    set(gcf,'color','w');
        box off;
    ax = gca;
    ax.FontSize = 45; % Set the font size as needed
    ax.FontName = 'Arial'; % Set the font face as needed
    
    % Customize axis labels format (degrees in this case)
    thetaticks([]);
    rticks([]); 
    rlim([0, .2])

    stepSize = 8;
    %3 Hz figure
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
     hold on 
     histogram(panelDat.miss_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
    set(gcf,'color','w');
        box off;
    xticks([-pi:pi/2:pi])
    ax=gca;ax.LineWidth=4;
    yticks([0:1/(stepSize*2):(1/(stepSize))*2])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)
    yline(1/(stepSize*2), 'color', 'k', 'linewidth',2, 'linestyle', '--')

    export_fig([figDat 'pubFigs/' 'Fig2_' panelDat.reg '_' ...
        panelDat.phase  '_phase_circHist.png'], '-transparent', '-r300')

    %6.5 Hz figure
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles_high, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
     hold on 
     histogram(panelDat.miss_angles_high, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
    set(gcf,'color','w');
        box off;
    xticks([-pi:pi/2:pi])
    ax=gca;ax.LineWidth=4;
    yticks([0:1/(stepSize*2):(1/(stepSize))*2])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)
    yline(1/(stepSize*2), 'color', 'k', 'linewidth',2, 'linestyle', '--')

    export_fig([figDat 'pubFigs/' 'Fig2_' panelDat.reg '_' ...
        panelDat.phase  '_phase_circHist_high.png'], '-transparent', '-r300')


    catch
        ii
    end
end

%% power spectra at image lock

test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TF_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_image'));
TF_files_image = TF_files(test); 


%hit, miss, and difference heatmaps 
%the hit and miss heatmaps go into supplemental figure 2
%differences are in main text Figure 3


parfor ii = 1:length(TF_files_image)
   
    panelDat = load([TF_files_image(ii).folder '/' ...
        TF_files_image(ii).name]).outDat;
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    %hit plot
    figure('visible', false, 'position', [0,0,600,600])
    
    imagesc(squeeze(mean(panelDat.hits_image))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    caxis([-2, 2])
    addRedOutline(panelDat.p_image, .05, 'white'); 
    ylim([0,50])
    colorbar
%     colormap(s2w2y)
    xticks([19:20:139])
    xticklabels(panelDat.tim([19:20:139]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'SupFig2_TFHit_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')

    %miss plot
    figure('visible', false, 'position', [0,0,600,600])
    
    imagesc(squeeze(mean(panelDat.misses_image))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    caxis([-2, 2])
    addRedOutline(panelDat.p_image, .05, 'white'); 
    ylim([0,50])
    colorbar
%     colormap(s2w2y)
    xticks([19:20:139])
    xticklabels(panelDat.tim([19:20:139]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'SupFig2_TFMiss_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')


    %dif plot
    figure('visible', false, 'position', [0,0,600,600])
    
    imagesc(squeeze(mean(panelDat.hits_image))' - ...
        squeeze(mean(panelDat.misses_image))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    caxis([-2, 2])
    addRedOutline(panelDat.p_image, .05, 'white'); 
    ylim([0,50])
    colorbar
%     colormap(s2w2y)
    xticks([19:20:139])
    xticklabels(panelDat.tim([19:20:139]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'SupFig2_TFDif_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')

    
    

end


%% ITPC at image lock

test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TFphase_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_image'));
TF_files_image = TF_files(test); 


%hit, miss, and difference heatmaps 
%the hit and miss heatmaps go into supplemental figure 2
%differences are in main text Figure 3


parfor ii = 1:length(TF_files_image)
   
    panelDat = load([TF_files_image(ii).folder '/' ...
        TF_files_image(ii).name]).outDat;
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    %hit plot
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    imagesc(squeeze(mean(panelDat.hits_image))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    caxis([1, 5])
%     if ~ismember(regi, [2,5])
    addRedOutline(panelDat.p_image, .05, 'white'); 
%     end
    ylim([0,50])
    colorbar
%     colormap(s2w2y)
    xticks([19:20:139])
    xticklabels(panelDat.tim([19:20:139]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'SupFig2_TF_phase_Hit_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')

    %miss plot
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    imagesc(squeeze(mean(panelDat.misses_image))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    caxis([1, 5])
%     if ~ismember(regi, [2,5])
    addRedOutline(panelDat.p_image, .05, 'white'); 
%     end
    ylim([0,50])
    colorbar
%     colormap(s2w2y)
    xticks([19:20:139])
    xticklabels(panelDat.tim([19:20:139]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'SupFig2_TF_phase_Miss_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')


    %dif plot
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    imagesc(squeeze(mean(panelDat.hits_image))' - ...
        squeeze(mean(panelDat.misses_image))')
    xline(find(panelDat.tim>=0,1), '--', 'linewidth', linWid, 'color', 'k')
    set(gca, 'ydir', 'normal')
    caxis([-2, 2])
%     if ~ismember(regi, [2,5])
    addRedOutline(panelDat.p_image, .05, 'white'); 
%     end
    ylim([0,50])
    colorbar
%     colormap(s2w2y)
    xticks([19:20:139])
    xticklabels(panelDat.tim([19:20:139]))
    yticks([1:11:100])
    yticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'SupFig2_TF_phase_Dif_' panelDat.reg '_' ...
        panelDat.phase  '.jpg'], '-r300')

    % make plots for the phase distribution! 
      stepSize = 8;
    %3 Hz figure
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
     hold on 
     histogram(panelDat.miss_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
    set(gcf,'color','w');
        box off;
    xticks([-pi:pi/2:pi])
    ax=gca;ax.LineWidth=4;
    yticks([0:1/(stepSize*2):(1/(stepSize))*2])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)
    yline(1/(stepSize*2), 'color', 'k', 'linewidth',2, 'linestyle', '--')

    export_fig([figDat 'pubFigs/' 'Fig2_Image_' panelDat.reg '_' ...
        panelDat.phase  '_phase_circHist.png'], '-transparent', '-r300')

    %6.5 Hz figure
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles_high, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
     hold on 
     histogram(panelDat.miss_angles_high, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4, ...
         'linewidth', 1.2);
    set(gcf,'color','w');
        box off;
    xticks([-pi:pi/2:pi])
    ax=gca;ax.LineWidth=4;
    yticks([0:1/(stepSize*2):(1/(stepSize))*2])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)
    yline(1/(stepSize*2), 'color', 'k', 'linewidth',2, 'linestyle', '--')

    export_fig([figDat 'pubFigs/' 'Fig2_Image_' panelDat.reg '_' ...
        panelDat.phase  '_phase_circHist_high.png'], ...
        '-transparent', '-r300')

end


%% figure for ratio between image and HFB locked ITPC 
test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TFphase_'));
TF_files = fig2Dat(test); 

parfor ii = 1:length(regions)

    curFiles = cellfun(@(x) length(x)>0, ...
        strfind({TF_files.name}, regions{ii})); 
    curFiles = TF_files(curFiles); 
    




end



%% fig 3 
%connectivity 

fig3Dat = dir([figDat 'Figure3']); 
fig3Dat(1:2) = []; 


%% grab all connections: 

%all connections in a reg X reg X time X hit/miss/t/p X enc/ret
allConnections = zeros(9,9,139, 20, 4, 2); 
allConnections2 = zeros(9,9,41,20,4,2); 
for ii = 1:length(fig3Dat)
    curDat = load([fig3Dat(ii).folder '/' fig3Dat(ii).name]).outDat; 
    tim = curDat.enc_image_tim; 
    reg1 = find(cellfun(@(x) strcmp(x, curDat.reg1), regions)); 
    reg2 = find(cellfun(@(x) strcmp(x, curDat.reg2), regions)); 
    %ENCODING: 
    %hit vals
    tmpMat = curDat.enc_image_hitVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 1, 1) = mean(tmpMat,3); 

    tmpMat = curDat.enc_HFB_hitVals; 
    allConnections2(reg1, reg2, :, :, 1, 1) = mean(tmpMat,3); 

    %miss vals
    tmpMat = curDat.enc_image_missVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 2, 1) = mean(tmpMat,3); 

    tmpMat = curDat.enc_HFB_missVals; 
    allConnections2(reg1, reg2, :, :, 2, 1) = mean(tmpMat,3);
    %t vals
    tmpMat = curDat.enc_image_tVal; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 3, 1) = tmpMat; 

    tmpMat = curDat.enc_HFB_tVal; 
    allConnections2(reg1, reg2, :, :, 3, 1) = tmpMat;
    %p vals
    tmpMat = curDat.enc_image_p; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 4, 1) = tmpMat; 

    tmpMat = curDat.enc_HFB_p; 
    allConnections2(reg1, reg2, :, :, 4, 1) = tmpMat;
    %RETRIEVAL: 
    tim = curDat.ret_image_tim; 
    %hit vals
    tmpMat = curDat.ret_image_hitVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 1, 2) = mean(tmpMat,3); 

    tmpMat = curDat.ret_HFB_hitVals; 
    allConnections2(reg1, reg2, :, :, 1, 2) = mean(tmpMat,3);
    %miss vals
    tmpMat = curDat.ret_image_missVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 2, 2) = mean(tmpMat,3); 

    tmpMat = curDat.ret_HFB_missVals; 
    allConnections2(reg1, reg2, :, :, 2, 2) = mean(tmpMat,3); 
    %t vals
    tmpMat = curDat.ret_image_tVal; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 3, 2) = tmpMat; 

    tmpMat = curDat.ret_HFB_tVal; 
    allConnections2(reg1, reg2, :, :, 3, 2) = tmpMat; 
    %p vals
    tmpMat = curDat.ret_image_p; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 4, 2) = tmpMat; 

    tmpMat = curDat.ret_HFB_p; 
    allConnections2(reg1, reg2, :, :, 4, 2) = tmpMat; 

end


%% which frequencies have the most connectivity? and when? 

%connections relative to HFB by frequency: 
pVals = squeeze(allConnections2(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 

figure('visible', false, 'position', [0,0,600,600])
 
yline(0, 'color', 'k', 'linewidth', linWid)
   
encSpect = sum(pVals(:,:,1), 1) ./ 41 ./ 25; 

y = encSpect; 
plot( y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 1) ./ 41 ./ 25;

y =retSpect;  
plot( y, 'color', 'blue', 'linewidth', 4)
ylim([0, .08])
xticks([1:4:20])
xticklabels(round(curDat.frex([1:4:20])))
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
HFBfrexidx =[[2,8]; [3,6]] ; %enc/Ret X first/second ;[[3,8]; [3,8]]
% xline(HFBfrexidx(1,:), 'linestyle', '-.', 'linewidth', ...
%     linWid/2, 'color', 'red')
% xline(HFBfrexidx(2,:), 'linestyle', '--', 'linewidth', ...
%     linWid/2, 'color', 'blue')
export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byFrex_' ...
     '_HFBLocked.jpg'], '-r300')



%connections relative to image by frequency: 
pVals = squeeze(allConnections(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 

figure('visible', false, 'position', [0,0,600,600])
 
yline(0, 'color', 'k', 'linewidth', linWid)
   
encSpect = sum(pVals(:,:,1), 1) ./ 139 ./ 25; 

y = encSpect; 
plot( y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 1) ./ 139 ./ 25;

y =retSpect;  
plot( y, 'color', 'blue', 'linewidth', 4)
ylim([0, .08])
xticks([1:4:20])
xticklabels(round(curDat.frex([1:4:20])))
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
imageFrexidx =[[2, 11]; [3,11]] ; %enc/Ret X first/second [[3, 11]; [3,11]]
% xline(imageFrexidx(1,:), 'linestyle', '-.', 'linewidth', ...
%     linWid/2, 'color', 'red')
% xline(imageFrexidx(2,:), 'linestyle', '--', 'linewidth', ...
%     linWid/2, 'color', 'blue')

export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byFrex_' ...
     '_imageLocked.jpg'], '-r300')

%connections relative to HFB by time: 
pVals = squeeze(allConnections2(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 

figure('visible', false, 'position', [0,0,600,600])
   
encSpect = sum(pVals(:,:,1), 2) ./ 20 ./ 25; 

y = encSpect; 
plot( y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 2) ./ 20 ./ 25;

y =retSpect;  
plot( y, 'color', 'blue', 'linewidth', 4)
ylim([0, .08])
xticks([1:10:41])
xticklabels([-500:250:500])
xline(21, 'linestyle', '--', 'linewidth', linWid, 'color', 'k')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
HFBtimidx =[[8, 25]; [19,23]] ; %enc/Ret X first/second [[8, 25]; [8,25]]
% xline(HFBtimidx(1,:), 'linestyle', '-.', 'linewidth', ...
%     linWid/2, 'color', 'red')
% xline(HFBtimidx(2,:), 'linestyle', '--', 'linewidth', ...
%     linWid/2, 'color', 'blue')

export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byTime' ...
     '_HFBLocked.jpg'], '-r300')

%connections relative to image by time: 
pVals = squeeze(allConnections(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [1,2])); 
x = curDat.enc_image_tim'; 
x(x<-450 | x>3000) = []; 
figure('visible', false, 'position', [0,0,600,600])
 
yline(0, 'color', 'k', 'linewidth', linWid)
   
encSpect = sum(pVals(:,:,1), 2) ./ 20 ./ 25; 

y = encSpect; 
plot(x, y, 'color', 'red', 'linewidth', 4)
hold on   
    
retSpect = sum(pVals(:,:,2), 2) ./ 20 ./ 25;

y =retSpect;  
plot(x, y, 'color', 'blue', 'linewidth', 4)

xline(0, 'linestyle', '--', 'linewidth', linWid, 'color', 'k')
ylim([0, .08])
xlim([-450, 3000])
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
%enc/Ret X 1/2/3 [[0,1050, 2250]; [0,1050,2250]]
imageTimidx =[[-25,1050, 2125]; [50,1050,2525]] ; 
% xline(imageTimidx(1,:), 'linestyle', '-.', ...
%                  'LineWidth', linWid/2, 'color', 'red')
% xline(imageTimidx(2,:), 'linestyle', '--', ...
%                  'LineWidth', linWid/2, 'color', 'blue')

export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byTime' ...
     '_imageLocked.jpg'], '-r300')


%connections relative to HFB by region: 
pVals = squeeze(allConnections2(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [3,4])); 


makeConnectionHeatMap(squeeze(pVals(:,:,1)) ./ (20*41), ones(5), regions,...
                        keyRegIdx, true, false)
caxis([0, .1])
% colorbar
export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byREG_enc' ...
     '_HFBLocked.jpg'], '-r300')

makeConnectionHeatMap(squeeze(pVals(:,:,2)) ./ (20*41), ones(5), regions,...
                        keyRegIdx, true, false)
caxis([0, .1])
% colorbar
export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byREG_ret' ...
     '_HFBLocked.jpg'], '-r300')

%connections relative to image by region: 
pVals = squeeze(allConnections(keyRegIdx, keyRegIdx, :, :, 4, :)); 
pVals = squeeze(sum(pVals<.05, [3,4])); 


makeConnectionHeatMap(squeeze(pVals(:,:,1)) ./ (20*139), ones(5), regions,...
                        keyRegIdx, true, false)
caxis([0, .1])
% colorbar
export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byREG_enc' ...
     '_ImageLocked.jpg'], '-r300')

makeConnectionHeatMap(squeeze(pVals(:,:,2)) ./ (20*139), ones(5), regions,...
                        keyRegIdx, true, false)
caxis([0, .1])
% colorbar
export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byREG_ret' ...
     '_ImageLocked.jpg'], '-r300')

%% schematic plots showing connections 

tim = curDat.enc_image_tim; 
tim(tim<-450 | tim>3000) = []; 
fn = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_image_']; 
fn2 = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_colKey.jpg']; 
fn3 = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_tKey.jpg']; 

%early encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-250 & tim<=500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc1.jpg'], fn2, fn3)

%mid encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=1500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc2.jpg'], fn2, fn3)

%late encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-1500 & tim<=3000, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc3.jpg'], fn2, fn3)

%early retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-250 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret1.jpg'], fn2, fn3)

%mid retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=1500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret2.jpg'], fn2, fn3)

%late retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-1500 & tim<=3000, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret3.jpg'], fn2, fn3)


%% schematics at HFB locked time

tim = [-500:25:500];  
fn = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_HFB_']; 
fn2 = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_colKey.jpg']; 
fn3 = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_tKey.jpg']; 
%early encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, true, [fn 'enc1.jpg'], fn2, fn3)

%mid encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc2.jpg'], fn2, fn3)

%late encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc3.jpg'], fn2, fn3)

%early retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret1.jpg'], fn2, fn3)

%mid retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret2.jpg'], fn2, fn3)

%late retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret3.jpg'], fn2, fn3)

%% same figures but with region colors
tim = [-500:25:500];  
fn = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_HFB_reg_']; 
fn2 = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_colKey.jpg']; 
fn3 = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\pubFigs\'...
        'Fig3_schematic_tKey.jpg']; 
%early encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'enc1.jpg'], fn2, fn3)

%mid encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'enc2.jpg'], fn2, fn3)

%late encoding
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'enc3.jpg'], fn2, fn3)

%early retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-500 & tim<=-200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'ret1.jpg'], fn2, fn3)

%mid retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=-200 & tim<=200, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'ret2.jpg'], fn2, fn3)

%late retrieval
plotMat = squeeze(allConnections2(keyRegIdx, keyRegIdx,...
                       tim>=200 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'reg', regColors, false, [fn 'ret3.jpg'], fn2, fn3)


%% fig 4 
%connectivity ROC analysis

fig3Dat = dir([figDat 'Figure3']); 
fig3Dat(1:2) = []; 




%timing note: early / mid / late 
%HFB: -500:-200 ; -200:200 ; 200:500
%image: -250:500 ; 500:1500 ; 1500:3000

%frequency note: low / mid / high
%ranges: 2-3.2 ; 3.2-5.3 ; 5.3-8.6

%clusters that overlap frequency or time bounds will be split and double
%counted
results = cell(length(fig3Dat),1); 
phases = {'enc', 'ret'}; 
statName = {'HFB', 'image'};

parfor ii = 1:length(fig3Dat)
    
    curDat = load([fig3Dat(ii).folder '/' fig3Dat(ii).name]).outDat; 
    aovDat = table;
    aovDat.subID = repmat("askj", 10000,1); 
    aovDat.chi = zeros(10000,1); 
    aovDat.chi2 = zeros(10000,1); %only applicable to connectivity
    aovDat.HFB_Image = repmat("askj", 10000,1); %HFB / Image timing
    aovDat.statType = repmat("askj", 10000,1); %HFB, TF power, ITPC, connectivity
    aovDat.encRet = repmat("askj", 10000,1); %enc / ret
    aovDat.reg1 = repmat("askj", 10000,1); 
    aovDat.reg2 = repmat("askj", 10000,1); %only applicable for connectivity 
    aovDat.AUC = zeros(10000,1); %area under the curve to differentiate hit/miss
    aovDat.hitVal = zeros(10000,1); %raw stat for hits
    aovDat.missVal = zeros(10000,1); %raw stat for misses
    aovDat.tVal = zeros(10000,1); %t value 
    aovDat.time = zeros(10000,1); % time index
    aovDat.freq = zeros(10000,1); % mean frequency
    ai = 1;
    for pp = 1:length(phases)
        for stat = 1:length(statName)
            tic
            
            [aovDat, ai] = connectionGraph(aovDat, ai, curDat, ...
                            regions, phases{pp}, statName{stat}, datPre);
%             [aovDat, ai] = connectionAUC(aovDat, ai, curDat, ...
%                             regions, phases{pp}, statName{stat}, datPre);

            disp(['file: ' num2str(ii) ' phase: ' ...
                phases{pp} ' stat: ' statName{stat} ' time: ' ...
                num2str(round(toc/60,2)) ])
        end
    end
    
    results{ii} = aovDat; 

end

test = cat(1, results{:});
idx = arrayfun(@(x) strcmp('askj', x), test.subID);
test(idx, :) = []; 
save('R:\MSS\Johnson_Lab\dtf8829\QuestConnect\AUCdat.mat', 'test')

figure
hold off
histogram(test.AUC(arrayfun(@(x) strcmp('lTemp',x), ...
    test.reg1)), [.3:.02:.8], 'normalization', 'probability')
hold on 
histogram(test.AUC, [.3:.02:.8], 'normalization', 'probability')


for pp = 1:2
    for r1 = 1:9
%         for r2 = 1:9
            reg1Screen = arrayfun(@(x) strcmp(regions{r1}, x), test.reg1);
%             reg2Screen = arrayfun(@(x) strcmp(regions{r2}, x), test.reg2); reg2Screen &
            phaseScreen = arrayfun(@(x) strcmp(phases{2}, x), test.encRet);
            cur = test(phaseScreen, :); 
            if size(cur,1) >0
                if length(unique(cur.HFB_Image)) == 2
                    figure
        histogram(cur.AUC(arrayfun(@(x) strcmp('HFB',x), ...
            cur.HFB_Image)), [.3:.02:.8], 'normalization', 'probability')
        hold on 
        histogram(cur.AUC(arrayfun(@(x) strcmp('image',x), ...
            cur.HFB_Image)), [.3:.02:.8], 'normalization', 'probability')
        title([regions{r1} ' ' phases{pp}])
                end
            end
%         end
    end
end
%%






























% %% Make region X region connectivity plots
% %based on the target frequency and time indices found in the line plots
% %above
% tim = curDat.enc_image_tim; 
% tim(tim<-450 | tim>3000 ) = []; 
% for pp = 1:2 %phase index
%     for fi = 1:2 %frequency index
%         %connections locked to image: 
%         tim = curDat.enc_image_tim; 
%         tim(tim<-450 | tim>3000 ) = []; 
%         for ii = 1:3 %three times across image presentation
%             %hits
%             plotMat = allConnections(keyRegIdx, keyRegIdx,...
%                 find(tim==imageTimidx(pp,ii)), ...
%                 imageFrexidx(pp, fi), 1, pp); 
%             pMat = allConnections(keyRegIdx, keyRegIdx,...
%                 find(tim==imageTimidx(pp,ii)), ...
%                 imageFrexidx(pp, fi), 4, pp); 
% 
%             makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
%                 false, true)
%             export_fig([figDat 'pubFigs/' 'Fig3_' ...
%                 'connectionHeatMap_image_' num2str(pp)...
%                 '_frex_' ...
%                 num2str( fi) ...
%                 '_tim_' num2str(ii) ...
%                 '_hits' '.jpg'], '-r300')
% 
%             %misses
%             plotMat = allConnections(keyRegIdx, keyRegIdx,...
%                 find(tim==imageTimidx(pp,ii)), ...
%                 imageFrexidx(pp, fi), 2, pp); 
%             pMat = allConnections(keyRegIdx, keyRegIdx,...
%                 find(tim==imageTimidx(pp,ii)), ...
%                 imageFrexidx(pp, fi), 4, pp); 
% 
%             makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
%                 false, true)
%             export_fig([figDat 'pubFigs/' 'Fig3_' ...
%                 'connectionHeatMap_image_' num2str(pp)...
%                 '_frex_' ...
%                 num2str(fi) ...
%                 '_tim_' num2str(ii) ...
%                 '_misses' '.jpg'], '-r300')
%             %tVals
%             plotMat = allConnections(keyRegIdx, keyRegIdx,...
%                 find(tim==imageTimidx(pp,ii)), ...
%                 imageFrexidx(pp, fi), 3, pp); 
%             pMat = allConnections(keyRegIdx, keyRegIdx,...
%                 find(tim==imageTimidx(pp,ii)), ...
%                 imageFrexidx(pp, fi), 4, pp); 
% 
%             makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
%                 true, true)
%             export_fig([figDat 'pubFigs/' 'Fig3_' ...
%                 'connectionHeatMap_image_' num2str(pp)...
%                 '_frex_' ...
%                 num2str(fi) ...
%                 '_tim_' num2str(ii) ...
%                 '_t' '.jpg'], '-r300')
% 
%         end
%         tim = [-500:25:500]; 
%         for ii = 1:2 %two times across HFB locked time
%              %hits
%             plotMat = allConnections2(keyRegIdx, keyRegIdx,...
%                 HFBtimidx(pp,ii), ...
%                 HFBfrexidx(pp, fi), 1, pp); 
%             pMat = allConnections2(keyRegIdx, keyRegIdx,...
%                 HFBtimidx(pp,ii), ...
%                 HFBfrexidx(pp, fi), 4, pp); 
% 
%             makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
%                 false, false)
%             export_fig([figDat 'pubFigs/' 'Fig3_' ...
%                 'connectionHeatMap_HFB_' num2str(pp)...
%                 '_frex_' ...
%                 num2str(fi) ...
%                 '_tim_' num2str(ii) ...
%                 '_hits' '.jpg'], '-r300')
% 
%             %misses
%             plotMat = allConnections2(keyRegIdx, keyRegIdx,...
%                 HFBtimidx(pp,ii), ...
%                 HFBfrexidx(pp, fi), 2, pp); 
%             pMat = allConnections2(keyRegIdx, keyRegIdx,...
%                 HFBtimidx(pp,ii), ...
%                 HFBfrexidx(pp, fi), 4, pp); 
% 
%             makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
%                 false, false)
%             export_fig([figDat 'pubFigs/' 'Fig3_' ...
%                 'connectionHeatMap_HFB_' num2str(pp)...
%                 '_frex_' ...
%                 num2str( fi) ...
%                 '_tim_' num2str(ii) ...
%                 '_misses' '.jpg'], '-r300')
%             %tVals
%             plotMat = allConnections2(keyRegIdx, keyRegIdx,...
%                 HFBtimidx(pp,ii), ...
%                 HFBfrexidx(pp, fi), 3, pp); 
%             pMat = allConnections2(keyRegIdx, keyRegIdx,...
%                 HFBtimidx(pp,ii), ...
%                 HFBfrexidx(pp, fi), 4, pp); 
% 
%             makeConnectionHeatMap(plotMat, pMat, regions, keyRegIdx,...
%                 true, false)
%             export_fig([figDat 'pubFigs/' 'Fig3_' ...
%                 'connectionHeatMap_HFB_' num2str(pp)...
%                 '_frex_' ...
%                 num2str( fi) ...
%                 '_tim_' num2str(ii) ...
%                 '_t' '.jpg'], '-r300')
% 
% 
% 
% 
%         end
% % %connections relative to HFB by frequency: 
% % pVals = squeeze(allConnections2(keyRegIdx, keyRegIdx, :, :, 4, :)); 
% % pVals = squeeze(sum(pVals<.05, [1,2])); 
% % 
% % figure('visible', false, 'position', [0,0,600,600])
% %  
% % yline(0, 'color', 'k', 'linewidth', linWid)
% %    
% % encSpect = sum(pVals(:,:,1), 1) ./ 41 ./ 25; 
% % 
% % y = encSpect; 
% % plot( y, 'color', 'red', 'linewidth', 4)
% % hold on   
% %     
% % retSpect = sum(pVals(:,:,2), 1) ./ 41 ./ 25;
% % 
% % y =retSpect;  
% % plot( y, 'color', 'blue', 'linewidth', 4)
% % ylim([0, .08])
% % xticks([1:4:20])
% % xticklabels(round(curDat.frex([1:4:20])))
% % set(gcf,'color','w');
% % box off;
% % ax=gca;ax.LineWidth=4;
% % 
% % xline([2,8], 'linestyle', '-.', 'linewidth', linWid/2, 'color', 'red')
% % xline([3,6], 'linestyle', '--', 'linewidth', linWid/2, 'color', 'blue')
% % export_fig([figDat 'pubFigs/' 'Fig3_' 'CountSigConnections_byFrex_' ...
% %      '_HFBLocked.jpg'], '-r300')
% 
% 
%     end
% end
% 
% 
% 
% 
% 
% 
% %%
% 
% 
% % Define the name of the video file
% outputVideoName = ['R:\MSS\Johnson_Lab\dtf8829\publicationFigureData'...
%     '/encoding.avi'];
% 
% % Create a VideoWriter object
% videoObj = VideoWriter(outputVideoName);
% videoObj.FrameRate = 6; % Optional: Specify the frame rate of the video
% 
% % Open the video file
% open(videoObj);
% tim(tim<-450 | tim>3000) = []; 
% fi = 5; 
% for ii = 1:139 % Example loop, adjust according to your needs
%     ii
%     % Generate your figure here
%     figure('visible', false, 'position', [0,0,1200, 400]);
%     subplot 131
%     plotMat = squeeze(allConnections(keyRegIdx,keyRegIdx,ii, fi, 1, 1)); 
%     plotMat(plotMat<0) = 0; 
%     plotMat = sqrt(plotMat); 
%     imagesc(plotMat)
%     xticks([1:9])
%     xticklabels(regions(keyRegIdx))
%     yticks([1:9])
%     yticklabels(regions(keyRegIdx))
%     clim([.05, .2])
%     axis square
%     title(['freq: ' num2str(round(curDat.frex(fi), 1)) ... 
%            ' time: ' num2str(tim(ii))])
%     subplot 132
%     plotMat = squeeze(allConnections(keyRegIdx,keyRegIdx,ii, fi, 2, 1)); 
%     plotMat(plotMat<0) = 0; 
%     plotMat = sqrt(plotMat); 
%     imagesc(plotMat)
%     xticks([1:9])
%     xticklabels(regions(keyRegIdx))
%     yticks([1:9])
%     yticklabels(regions(keyRegIdx))
%     clim([.05, .2])
%     axis square
%     title(['freq: ' num2str(round(curDat.frex(fi), 1)) ... 
%            ' time: ' num2str(tim(ii))])
%     subplot 133
%     tmpMat = squeeze(allConnections(keyRegIdx,keyRegIdx,ii, fi, 3, 1));
%     
%     imagesc(tmpMat)
%     hold on 
%     pMat = squeeze(allConnections(keyRegIdx,keyRegIdx,ii, fi, 4, 1)); 
%     for rri = 1:5
%         for rrj = 1:5
%             if pMat(rri, rrj) < .05
%                 scatter(rri, rrj, 30, 'red', 'filled', 'o')
%             end
%         end
%     end
%     xticks([1:9])
%     xticklabels(regions(keyRegIdx))
%     yticks([1:9])
%     yticklabels(regions(keyRegIdx))
%     clim([-3, 3])
%     axis square
%     title(['freq: ' num2str(round(curDat.frex(fi), 1)) ... 
%            ' time: ' num2str(tim(ii))])
% 
% 
%   
%     
%     % Capture the current figure as a frame
%     frame = getframe(gcf);
%     
%     % Write the frame to the video
%     writeVideo(videoObj, frame);
%     
%     % Close the figure to avoid GUI overload
%     close;
% end
% 
% % Close the video file
% close(videoObj);
% 
% 
% %%
% 
% 
% 
% 
% 
% 
% %hippocampus hard code encoding only 
% dat = load([headFiles(6).folder '/' headFiles(6).name]).statInfo;
% 
% 
% hitSpect = arrayfun(@(x) find(dat.tim>=x, 1), dat.hitLat);
% hitSpect = arrayfun(@(x,y) squeeze(dat.hits(x,y,:)), ...
%     hitSpect, [1:length(hitSpect)]', 'uniformoutput', false);
% hitSpect = cat(2, hitSpect{:});
% 
% missSpect = arrayfun(@(x) find(dat.tim>=x, 1), dat.missLat);
% missSpect = arrayfun(@(x,y) squeeze(dat.misses(x,y,:)), ...
%     missSpect, [1:length(missSpect)]', 'uniformoutput', false);
% missSpect = cat(2, missSpect{:});
% 
% 
% %% basic image locked HFB time courses
% 
% %while looping, grab the distribution of effect sizes for later comparison
% %to HFB-peak diffs
% %effect sizes stored as time X enc/ret X region
% hfbDist = zeros(2, 9); 
% for ii = 1:length(fig2Dat)
%     panelDat = load([fig2Dat(ii).folder '/' fig2Dat(ii).name]).outDat; 
% 
%     %quick cleaning for plotting purposes only, all data were used in stats
%     panelDat.hitVals(panelDat.eliminate,:) = []; 
%     panelDat.missVals(panelDat.eliminate,:) = []; 
%     if strcmp(panelDat.reg, 'acc')
%         panelDat.hitVals(9:10, :) = []; 
%         panelDat.missVals(9:10, :) = []; 
%     end
%     panelDat.hitVals(isnan(panelDat.hitVals)) = 0; 
%     panelDat.missVals(isnan(panelDat.missVals)) = 0; 
% 
%     %get effect sizes
%     diffs = panelDat.hitVals - panelDat.missVals; 
%     diffs = mean(diffs) ./ std(diffs); 
%     regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
%     phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 
% 
%     hfbDist(phasei, regi) = diffs(21); 
% %     figure('visible', false, 'position', [0,0,600,600])
% %     x = panelDat.tim; 
% %     x(panelDat.p>.05) = []; 
% %     if ~isempty(x) %check if we have any sig time points
% %     breakVals = [0, find(diff(x)> 25), length(x)];
% %     for ii = 1:length(breakVals)-1
% %         tmpX = x(breakVals(ii)+1:breakVals(ii+1));
% %         tmpY = ones(length(tmpX),1) * 10000; 
% %         tmpX = [tmpX, flip(tmpX)]; 
% %         tmpY = [tmpY', -tmpY'];
% %         fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
% %         hold on 
% %     
% %     
% %     end
% %     else
% %         hold on 
% %     end
% %     yline(0, 'color', 'k', 'linewidth', linWid)
% %     xline(0, '--', 'linewidth', linWid, 'color', 'k')
% %     xline(panelDat.meanHitRT, 'color', hitCol, ...
% %         'linewidth', linWid, 'linestyle', '--')
% %     xline(panelDat.meanMissRT, 'color', missCol, ...
% %         'linewidth', linWid, 'linestyle', '--')
% %     y = movmean(mean(panelDat.hitVals), 3); 
% %     x = panelDat.tim; 
% %     plot(x, y, 'color', hitCol, 'linewidth', 4)
% %     se = std(panelDat.hitVals) ./ sqrt(size(panelDat.hitVals,1)); 
% %     y = [y + se, flip(y) - flip(se)]; 
% %     x = [x, flip(x)]; 
% %     hold on 
% %     fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
% %     allMax = max(y); 
% %     allMin = min(y); 
% %     
% %     y = movmean(mean(panelDat.missVals), 3); 
% %     x = panelDat.tim; 
% %     plot(x, y, 'color', missCol, 'linewidth', 4)
% %     se = std(panelDat.missVals) ./ sqrt(size(panelDat.missVals,1)); 
% %     y = [y + se, flip(y) - flip(se)]; 
% %     x = [x, flip(x)]; 
% %     hold on 
% %     fill(x, y, missCol, 'facealpha', errAlpha, 'edgealpha', 0)
% %     allMax = max([allMax, max(y)]); 
% %     allMin = min([allMin, min(y)]); 
% %     ylim([allMin*1.1, allMax*1.1])
% %     set(gcf,'color','w');
% %     box off;
% %     ax=gca;ax.LineWidth=4;
% %     xlim([-450, 3000])
%     
% end
% 
% 
% 
% 
















