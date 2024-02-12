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

addpath([codePre 'HpcAccConnectivityProject'])
addpath([codePre 'myFrequentUse'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
HFBdat = load([datPre 'HFB_KEY_STATS\hip.mat']).HFBdat; 
regions = {HFBdat.aggTargs.lab}; 
regions(2) = []; 




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

regColors = [[214, 77, 97]',... %parahip: dark pink
             [145, 162, 80]', ...%dlPFC: light forest green
             [145, 162, 80]', ...%iTemp: light forest green
             [145, 162, 80]', ...%lTemp: light forest green
             [51, 102, 153]', ...%polarPFC: lilac
             [132, 125, 212]', ...%visual: lilac
             [204, 153, 204]', ...%ACC: dark red
             [184, 78, 83]', ...%PCC: dark red
             [61, 187, 218]']' ./ 255; %Hip: dusky blue

keyRegIdx = [1,2,5,7,9]; 

b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], ...
    [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)';...
    linspace(255,0,128)']]/255;
b2w2r(129:end, 1) = 1; 
b2w2r(1:128, 3) = 1; 

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
imgDist = zeros(139, 2, 9); 
for ii = 1:length(fig1Dat)
    panelDat = load([fig1Dat(ii).folder '/' fig1Dat(ii).name]).outDat; 

    %quick cleaning for plotting purposes only, all data were used in stats
    panelDat.hitVals(panelDat.eliminate,:) = []; 
    panelDat.missVals(panelDat.eliminate,:) = []; 
    if strcmp(panelDat.reg, 'acc')
        panelDat.hitVals(9:10, :) = []; 
        panelDat.missVals(9:10, :) = []; 
    end
    panelDat.hitVals(isnan(panelDat.hitVals)) = 0; 
    panelDat.missVals(isnan(panelDat.missVals)) = 0; 

    %get effect sizes
    diffs = panelDat.hitVals - panelDat.missVals; 
    diffs = mean(diffs) ./ std(diffs); 
    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    imgDist(:,phasei, regi) = diffs; 
    figure('visible', false, 'position', [0,0,600,600])
    x = panelDat.tim; 
    x(panelDat.p>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
    breakVals = [0, find(diff(x)> 25), length(x)];
    for ii = 1:length(breakVals)-1
        tmpX = x(breakVals(ii)+1:breakVals(ii+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
    
    
    end
    else
        hold on 
    end
    yline(0, 'color', 'k', 'linewidth', 2)
    xline(0, '--', 'linewidth', 2, 'color', 'k')
    xline(panelDat.meanHitRT, 'color', hitCol, ...
        'linewidth', 2, 'linestyle', '--')
    xline(panelDat.meanMissRT, 'color', missCol, ...
        'linewidth', 2, 'linestyle', '--')
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
    xline(0, '--', 'linewidth', 2)
    scatter(panelDat.hitRT(order), [1:length(panelDat.hitRT)], ...
        20, 'white', 'filled')
%     plot(sortedLat, [1:length(sortedLat)], 'color', regColors(9,:), ...
%         'linestyle', '--', 'linewidth', 2)
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    yticks([])

    export_fig([figDat 'pubFigs/' 'Fig1_singleTrial' ...
        panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')


% end
%% get the latencies for all regions on a single plot HITS Retrieve
linew = 5;
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
%     yline(0, 'color', 'k', 'linewidth', 2)
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
%         yline(0, 'color', 'k', 'linewidth', 2)
    
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
    ylim([allMin*1.1, allMax*1.1])
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
    yline(0, 'color', 'k', 'linewidth', 2)
   
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
    %outliers in the hippocampus retrieval data, remove for plot
     if regi == 9 
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
    allMin = min([allMin, 0]); 
    ylim([allMin*1.1, allMax*1.1])
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
    figure('visible', false, 'position', [0,0,600,600])
    hold off
    histogram(panelDat.hit_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', hitCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4);
     hold on 
     histogram(panelDat.miss_angles, [-pi:pi/stepSize:pi], ...
         'facecolor', missCol, 'edgecolor', 'k',...
         'Normalization', 'probability', 'facealpha', .4);
    set(gcf,'color','w');
        box off;
    xticks([])
    yticks([])
    ylim([0, .32])
    plot(-pi:pi/stepSize:pi, sin(-pi/2:pi/stepSize:3*pi/2) *.09 + .13,...
        'color', 'k', 'linewidth', 3)

    export_fig([figDat 'pubFigs/' 'Fig2_' panelDat.reg '_' ...
        panelDat.phase  '_phase_circHist.png'], '-transparent', '-r300')

    catch
        ii
    end
end




%%






%hippocampus hard code encoding only 
dat = load([headFiles(6).folder '/' headFiles(6).name]).statInfo;


hitSpect = arrayfun(@(x) find(dat.tim>=x, 1), dat.hitLat);
hitSpect = arrayfun(@(x,y) squeeze(dat.hits(x,y,:)), ...
    hitSpect, [1:length(hitSpect)]', 'uniformoutput', false);
hitSpect = cat(2, hitSpect{:});

missSpect = arrayfun(@(x) find(dat.tim>=x, 1), dat.missLat);
missSpect = arrayfun(@(x,y) squeeze(dat.misses(x,y,:)), ...
    missSpect, [1:length(missSpect)]', 'uniformoutput', false);
missSpect = cat(2, missSpect{:});


%% basic image locked HFB time courses

%while looping, grab the distribution of effect sizes for later comparison
%to HFB-peak diffs
%effect sizes stored as time X enc/ret X region
hfbDist = zeros(2, 9); 
for ii = 1:length(fig2Dat)
    panelDat = load([fig2Dat(ii).folder '/' fig2Dat(ii).name]).outDat; 

    %quick cleaning for plotting purposes only, all data were used in stats
    panelDat.hitVals(panelDat.eliminate,:) = []; 
    panelDat.missVals(panelDat.eliminate,:) = []; 
    if strcmp(panelDat.reg, 'acc')
        panelDat.hitVals(9:10, :) = []; 
        panelDat.missVals(9:10, :) = []; 
    end
    panelDat.hitVals(isnan(panelDat.hitVals)) = 0; 
    panelDat.missVals(isnan(panelDat.missVals)) = 0; 

    %get effect sizes
    diffs = panelDat.hitVals - panelDat.missVals; 
    diffs = mean(diffs) ./ std(diffs); 
    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    hfbDist(phasei, regi) = diffs(21); 
%     figure('visible', false, 'position', [0,0,600,600])
%     x = panelDat.tim; 
%     x(panelDat.p>.05) = []; 
%     if ~isempty(x) %check if we have any sig time points
%     breakVals = [0, find(diff(x)> 25), length(x)];
%     for ii = 1:length(breakVals)-1
%         tmpX = x(breakVals(ii)+1:breakVals(ii+1));
%         tmpY = ones(length(tmpX),1) * 10000; 
%         tmpX = [tmpX, flip(tmpX)]; 
%         tmpY = [tmpY', -tmpY'];
%         fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
%         hold on 
%     
%     
%     end
%     else
%         hold on 
%     end
%     yline(0, 'color', 'k', 'linewidth', 2)
%     xline(0, '--', 'linewidth', 2, 'color', 'k')
%     xline(panelDat.meanHitRT, 'color', hitCol, ...
%         'linewidth', 2, 'linestyle', '--')
%     xline(panelDat.meanMissRT, 'color', missCol, ...
%         'linewidth', 2, 'linestyle', '--')
%     y = movmean(mean(panelDat.hitVals), 3); 
%     x = panelDat.tim; 
%     plot(x, y, 'color', hitCol, 'linewidth', 4)
%     se = std(panelDat.hitVals) ./ sqrt(size(panelDat.hitVals,1)); 
%     y = [y + se, flip(y) - flip(se)]; 
%     x = [x, flip(x)]; 
%     hold on 
%     fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
%     allMax = max(y); 
%     allMin = min(y); 
%     
%     y = movmean(mean(panelDat.missVals), 3); 
%     x = panelDat.tim; 
%     plot(x, y, 'color', missCol, 'linewidth', 4)
%     se = std(panelDat.missVals) ./ sqrt(size(panelDat.missVals,1)); 
%     y = [y + se, flip(y) - flip(se)]; 
%     x = [x, flip(x)]; 
%     hold on 
%     fill(x, y, missCol, 'facealpha', errAlpha, 'edgealpha', 0)
%     allMax = max([allMax, max(y)]); 
%     allMin = min([allMin, min(y)]); 
%     ylim([allMin*1.1, allMax*1.1])
%     set(gcf,'color','w');
%     box off;
%     ax=gca;ax.LineWidth=4;
%     xlim([-450, 3000])
%     export_fig([figDat 'pubFigs/' 'Fig1_' panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')
    
end




















