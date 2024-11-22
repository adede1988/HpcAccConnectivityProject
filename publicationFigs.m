%% MTL - PFC publication figures 

%description of the general data set 


%% set environment

%figures should not be docked: 
set(0, 'defaultfigurewindowstyle', 'normal')
%local paths: 


codePre = 'G:\My Drive\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\QuestConnect\';
figDat = 'R:\MSS\Johnson_Lab\dtf8829\publicationFigureData\';

% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath([codePre 'fieldtrip-20230118'])
ft_defaults;
regions = {'acc', 'dlPFC', 'hip', ...
    'lTemp', 'iTemp', 'mtl', 'pcc', 'pPFC', 'vis'}; 

allSig = readtable([codePre 'HpcAccConnectivityProject/demographics.csv']);



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

s = [85,25,86];
m = [186,21,77]; 
e = [249,205,15]; 

purpleYellow = [[linspace(s(1),m(1),128)'; linspace(m(1),e(1),128)'], ...
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

%% methods figure, what is HFB reactive? 
postFigPath = 'G:\My Drive\Johnson\CNS2024/';
test = load([headFilesHFB(12).folder '/' headFilesHFB(12).name]).statInfo;

reactChan = load([datPre 'CHANDAT\finished\chanDat_IR84_034.mat']).chanDat; 
testSub = test.hitChi(cellfun(@(x) strcmp(x, 'IR84'), test.hitSub));
nonReact = load([datPre 'CHANDAT\finished\chanDat_IR84_035.mat']).chanDat;

figure('position', [0,0,600,1200])
allTrials = [reactChan.HFB.subHit, reactChan.HFB.subMiss]; 
imagesc(reactChan.HFB.encMulTim, [], allTrials')
colormap(s2w2y)
clim([0,10])
colorbar
hold on 
xline(0, '--', 'linewidth', linWid, 'color', 'white')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
yticks([])
xlim([-450, 3000])
export_fig([postFigPath 'reactiveMethod_yes_heat' '.jpg'], '-r300')

figure('position', [0,0,600,400])
plot(reactChan.HFB.encMulTim, mean(allTrials,2), 'linewidth', linWid)
xline(0, '--', 'linewidth', linWid)
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
ylim([-5, 5])
yline(1.96, '--', 'linewidth', linWid, 'color', 'red')
yline(-1.96, '--', 'linewidth', linWid, 'color', 'red')
xlim([-450, 3000])
export_fig([postFigPath 'reactiveMethod_yes_mean' '.jpg'], '-r300')

figure('position', [0,0,1000,400])
hold on 
spread = 35;
tim = reactChan.HFB.encMulTim; 
for ii = 9:12
    hitIdx = find(reactChan.hits & reactChan.use); 
    RT = reactChan.encInfo(hitIdx(ii),4);
    plot(tim, allTrials(:,ii)+ii*spread, ...
        'color', purpleYellow((13-9)*60+1,:), 'linewidth', 2)
    scatter(repmat(RT,[50,1]), ii*linspace(spread-1.25,spread+1.25,50), 10, ...
        'red', 'filled')
    trial = allTrials(:,ii);            
    test = gausswin(11);
    test = test ./ sum(test); 
    trial = [zeros(5,1); trial; zeros(5,1)];
    smoothT = conv(trial, test, 'valid'); 
    coli = 14-ii;
    plot(tim, smoothT+ii*spread, ...
        'color', purpleYellow((1-1)*60+1,:), 'linewidth', 2)
    [maxVal, idx] = max(smoothT(find(tim==0):find(tim>=RT,1))); 
    testTim = tim(find(tim==0):find(tim>=RT,1));
    scatter(testTim(idx), smoothT(tim==testTim(idx))+ii*spread,...
            140,'k', 'filled' )
    scatter(testTim(idx), smoothT(tim==testTim(idx))+ii*spread,...
            100,[126, 165, 21]./255, 'filled' )

end
xline(0, '--', 'linewidth', linWid)
set(gcf,'color','w');
box off;
ax=gca; ax.LineWidth=4;
xlim([-450, 3000])
yticks([])
export_fig([figDat 'pubFigs/' 'Fig1_HFBpeakDetection' '.jpg'], '-r300')


figure('position', [0,0,600,1200])
allTrials = [nonReact.HFB.subHit, nonReact.HFB.subMiss]; 
imagesc(nonReact.HFB.encMulTim, [], allTrials')
colormap(s2w2y)
clim([0,10])
colorbar
hold on 
xline(0, '--', 'linewidth', linWid, 'color', 'white')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
yticks([])
xlim([-450, 3000])
export_fig([postFigPath 'reactiveMethod_non_heat' '.jpg'], '-r300')

figure('position', [0,0,600,400])
plot(nonReact.HFB.encMulTim, mean(allTrials,2), 'linewidth', linWid)
xline(0, '--', 'linewidth', linWid)
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
ylim([-5, 5])
yline(1.96, '--', 'linewidth', linWid, 'color', 'red')
yline(-1.96, '--', 'linewidth', linWid, 'color', 'red')
xlim([-450, 3000])
export_fig([postFigPath 'reactiveMethod_non_mean' '.jpg'], '-r300')
    




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


%panel 1A
for ii = 1:length(fig1Dat)
    ii
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
    subIDs = cellfun(@(x) split(x, '_'), panelDat.realID, ...
                         'UniformOutput', false); 
    subIDs = cellfun(@(x) x{1}, subIDs, 'uniformoutput', false); 
    uniIDs = unique(subIDs); 
  
  
    figure('visible', false, 'position', [0,0,600,600])
    x = panelDat.tim; 
    x(panelDat.p>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
    breakVals = [0, find(diff(x)> 25), length(x)];
    
    for jj = 1:length(breakVals)-1
        varName = ['HFB_' panelDat.reg '_' panelDat.phase num2str(jj)];
        allSig.(varName) = nan(length(allSig.subID),1); 
        tmpX = x(breakVals(jj)+1:breakVals(jj+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
        %store raw differences for significant time periods: 
        chanMeans = mean(panelDat.hitVals(:,...
                            breakVals(jj)+1:breakVals(jj+1)) ...
                            - ...
                         panelDat.missVals(:,...
                            breakVals(jj)+1:breakVals(jj+1)), 2);
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
    
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
    export_fig([figDat 'pubFigs/' 'Fig1_' ...
        panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')
    
end






%% single trial heatmap figure and HFB peak
% for ii = 1:length(singleFig1Dat)
%hard code for acc at encoding hits only
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

    for ii = 1:18
    panelDat = load([headFilesHFB(ii).folder '/' ...
        headFilesHFB(ii).name]).statInfo; 

    HFBidx = arrayfun(@(x) find(panelDat.tim >= x, 1), panelDat.hitLat); 
    HFBidx(HFBidx<21) = 21; 
    HFBidx(HFBidx>length(panelDat.tim)-20) = length(panelDat.tim)- 20;
    %get the +- 20 points around the HFB peaks
    HFB_aligned_hits = arrayfun(@(x, y) panelDat.hits(x-20:x+20, y), HFBidx', ...
        1:length(HFBidx), 'UniformOutput',false );
    HFB_aligned_hits = cat(2, HFB_aligned_hits{:});
    %get channel means
    chi = panelDat.hitChi; 
    subIDs = panelDat.hitSub; 
    chanIDs = arrayfun(@(x) [subIDs{x} '_' num2str(chi(x))], 1:length(chi),...
        'uniformoutput', false);
    uniIDs = unique(chanIDs); 
    HFB_aligned_hits = cellfun(@(x)...
        mean(HFB_aligned_hits(:, ismember(chanIDs, x)),2), uniIDs, ...
        'uniformoutput', false);
    HFB_aligned_hits = cat(2, HFB_aligned_hits{:}); 


    HFBidx = arrayfun(@(x) find(panelDat.tim >= x, 1), panelDat.missLat); 
    HFBidx(HFBidx<21) = 21; 
    HFBidx(HFBidx>length(panelDat.tim)-20) = length(panelDat.tim)- 20;
    %get the +- 20 points around the HFB peaks
    HFB_aligned_misses = arrayfun(@(x, y) panelDat.misses(x-20:x+20, y), HFBidx', ...
        1:length(HFBidx), 'UniformOutput',false );
    HFB_aligned_misses = cat(2, HFB_aligned_misses{:});
    %get channel means
    chi = panelDat.missChi; 
    subIDs = panelDat.missSub; 
    chanIDs = arrayfun(@(x) [subIDs{x} '_' num2str(chi(x))], 1:length(chi),...
        'uniformoutput', false);
    uniIDs = unique(chanIDs); 
    HFB_aligned_misses = cellfun(@(x)...
        mean(HFB_aligned_misses(:, ismember(chanIDs, x)),2), uniIDs, ...
        'uniformoutput', false);
    HFB_aligned_misses = cat(2, HFB_aligned_misses{:}); 

    figure('visible', true, 'position', [0,0,600,600])
    x = [-500:25:500]; 
    hold on 
    yline(0, 'color', 'k', 'linewidth', linWid)
    xline(0, '--', 'linewidth', linWid, 'color', 'k')
   
    y = mean(HFB_aligned_hits,2); 
  
    plot(x, y, 'color', hitCol, 'linewidth', 4)
    se = std(HFB_aligned_hits,[],2) ./ sqrt(size(HFB_aligned_hits,2)); 
    y = [y + se; flip(y) - flip(se)]; 
    x = [x, flip(x)]; 
    hold on 
    fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
    allMax = max(y); 
    allMin = min(y); 
    
    y = mean(HFB_aligned_misses,2); 
    x = [-500:25:500]; 
    plot(x, y, 'color', missCol, 'linewidth', 4)
    se = std(HFB_aligned_hits,[],2) ./ sqrt(size(HFB_aligned_hits,2)); 
    y = [y + se; flip(y) - flip(se)]; 
    x = [x, flip(x)]; 
    hold on 
    fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    export_fig([figDat 'pubFigs/' 'Fig1_HFBpeakExample_' ...
        panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')

    end

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


%% basic HFB time courses LOCKED TO HFB PEAK! 


test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'HFB.'));
test2 = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'HFB'));
HFBfiles = fig2Dat(~test & test2); 


%while looping, grab the distribution of effect sizes for later comparison
%to HFB-peak diffs
%effect sizes stored as time X enc/ret X region
%panel 1A
% imgDist = zeros(139, 2, 9); 

parfor ii = 1:length(HFBfiles)
    panelDat = load([HFBfiles(ii).folder '/' HFBfiles(ii).name]).outDat; 

    %quick cleaning for plotting purposes only, all data were used in stats
    panelDat.hitVals(panelDat.eliminate,:) = []; 
    panelDat.missVals(panelDat.eliminate,:) = []; 

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
    ylim([-5, 25])
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
   
    export_fig([figDat 'pubFigs/' 'Fig1_HFBlocked_'...
        panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')
    
end



%% individual peak events


test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TF_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_HFB'));
TF_files_HFB = TF_files(test); 

%grab a .csv of data with the columns: power, frequency, enc/ret, subID,
%chanID
% % % aovDat = table;
% % % aovDat.power = zeros(112200,1); 
% % % aovDat.freq = zeros(112200, 1); 
% % % aovDat.encRet = repmat("askj", 112200,1); 
% % % aovDat.subID = repmat("askj", 112200,1); 
% % % aovDat.reg = repmat("askj", 112200,1); 
% % % ai = 1; 

allSpect = struct; 

parfor ii = 1:length(TF_files_HFB)
    try
    panelDat = load([TF_files_HFB(ii).folder '/' ...
        TF_files_HFB(ii).name]).outDat;
    ii

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals));





  % MEAN PEAK ALIGNED BROADBAND AND HFB BAND DATA

    %HITS
    BroadBandPeaksHITS = zeros([3001, length(panelDat.hitSub)]);
    HFBpeaksHITS = zeros([3001, length(panelDat.hitSub)]); 
    for tt = 1:length(panelDat.hitSub)
       tt
        %load in a single channel of data: 
        if tt==1
            triali = 1; 
            if panelDat.hitChi(tt)<10
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.hitSub{tt} '_00' num2str(panelDat.hitChi(tt)) ...
        '.mat']).chanDat; 
            elseif panelDat.hitChi(tt)<100
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.hitSub{tt} '_0' num2str(panelDat.hitChi(tt)) ...
        '.mat']).chanDat; 
            else
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.hitSub{tt} '_' num2str(panelDat.hitChi(tt)) ...
        '.mat']).chanDat; 
            end
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end

        elseif (panelDat.hitChi(tt-1) ~= panelDat.hitChi(tt)) || ...
                (~strcmp(panelDat.hitSub(tt-1),panelDat.hitSub(tt)))
            triali = 1; 
            if panelDat.hitChi(tt)<10
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.hitSub{tt} '_00' num2str(panelDat.hitChi(tt)) ...
        '.mat']).chanDat; 
            elseif panelDat.hitChi(tt)<100
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.hitSub{tt} '_0' num2str(panelDat.hitChi(tt)) ...
        '.mat']).chanDat; 
            else
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.hitSub{tt} '_' num2str(panelDat.hitChi(tt)) ...
        '.mat']).chanDat; 
            end
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end
        end
        

       
        %HITS
        curTrial = hfbBand(:,hitidx);
        curTrial = curTrial(:, triali); 
        L =length(curTrial); 
        curTrial = mirrorPad(curTrial); 
        phaseTrial = angle(hilbert(curTrial));
        idx = find(tim>=HFBlatHit(triali,1),1); 
        [~, troughOff] = min(phaseTrial(L+idx-10:L+idx+10));
%         troughOff = 10; %% debug check to see how much trough alignment matters
        idx = idx + troughOff - 10; 
        HFBpeaksHITS(:, tt) = curTrial(L+idx-1500: L+idx+1500); 
        curTrial = broadBand(:, hitidx); 
        curTrial = curTrial(:, triali); 
        curTrial = mirrorPad(curTrial); 
        BroadBandPeaksHITS(:, tt) = curTrial(L+idx-1500: L+idx+1500);
%         idx = randsample(HFBlatHit, 1); 
%         idx = find(tim>=idx,1); 
%         BroadBandPeaksHITS(2,:,tt) = curTrial(idx-200:idx+200);
        triali = triali +1; 
    end


    allRatios = []; 

    for nn = 1:1400
        during = mean(abs(hilbert(HFBpeaksHITS(1400:1600,nn))));
        before = mean(abs(hilbert(HFBpeaksHITS(1200:1400,nn))));
        allRatios = [allRatios during/before];

        %plotting individual events: 
        if during/before > 2 && ~exist([figDat 'pubFigs/exampleTrials/' ...
                 'highRatio_'  panelDat.reg '_'...
                 panelDat.phase '_' num2str(nn)  '.jpg'])
            fig = figure('visible', false);
            subplot 211
            plot(BroadBandPeaksHITS(:,nn) - BroadBandPeaksHITS(1501,nn),...
                'color', 'k')
            hold on 
            plot(HFBpeaksHITS(:,nn).*2,...
                'color', [191, 64, 191]./255)
            xlim([0,3001])
            xline(1501)

            title([num2str(round(during/before, 2)) ' during/before power'])

            subplot 212
            plot(BroadBandPeaksHITS(1200:1800,nn) - BroadBandPeaksHITS(1501,nn),...
                'color', 'k')
            hold on 
            plot(HFBpeaksHITS(1200:1800,nn).*2,...
                'color', [191, 64, 191]./255)
            xlim([0,600])
            xline(301)
            saveas(fig, [figDat 'pubFigs/exampleTrials/' ...
                 'highRatio_'  panelDat.reg '_'...
                 panelDat.phase '_' num2str(nn)  '.jpg'])


        end

        if during/before < 1.5 && ~exist([figDat 'pubFigs/exampleTrials/' ...
                 'lowRatio_'  panelDat.reg '_'...
                 panelDat.phase '_' num2str(nn)  '.jpg'])
            fig = figure('visible', false);
            subplot 211
             plot(BroadBandPeaksHITS(:,nn) - BroadBandPeaksHITS(1501,nn),...
                'color', 'k')
             hold on 
            plot(HFBpeaksHITS(:,nn).*2,...
                'color', [191, 64, 191]./255)
            xlim([0,3001])
            xline(1501)

            title([num2str(round(during/before, 2)) ' during/before power'])

            subplot 212
            plot(BroadBandPeaksHITS(1200:1800,nn) - ...
                BroadBandPeaksHITS(1501,nn),'color', 'k')
            hold on 
            plot(HFBpeaksHITS(1200:1800,nn).*2,...
                'color', [191, 64, 191]./255)
            xlim([0,600])
            xline(301)
            saveas(fig, [figDat 'pubFigs/exampleTrials/' ...
                 'lowRatio_'  panelDat.reg '_'...
                 panelDat.phase '_' num2str(nn)  '.jpg'])


        end
        %end plotting individual events
    end





    %MISSES
    BroadBandPeaksMISSES = zeros([401, length(panelDat.missSub)]);
    HFBpeaksMISSES = zeros([401, length(panelDat.missSub)]); 

    for tt = 1:length(panelDat.missSub)
   
        %load in a single channel of data: 
        if tt==1
            triali = 1; 
            if panelDat.missChi(tt)<10
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.missSub{tt} '_00' num2str(panelDat.missChi(tt)) ...
        '.mat']).chanDat; 
            elseif panelDat.missChi(tt)<100
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.missSub{tt} '_0' num2str(panelDat.missChi(tt)) ...
        '.mat']).chanDat; 
            else
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.missSub{tt} '_' num2str(panelDat.missChi(tt)) ...
        '.mat']).chanDat; 
            end
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end
            
        elseif (panelDat.missChi(tt-1) ~= panelDat.missChi(tt)) || ...
                (~strcmp(panelDat.missSub(tt-1),panelDat.missSub(tt)))
            triali = 1; 
            if panelDat.missChi(tt)<10
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.missSub{tt} '_00' num2str(panelDat.missChi(tt)) ...
        '.mat']).chanDat; 
            elseif panelDat.missChi(tt)<100
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.missSub{tt} '_0' num2str(panelDat.missChi(tt)) ...
        '.mat']).chanDat; 
            else
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\' ...
        'chanDat_' panelDat.missSub{tt} '_' num2str(panelDat.missChi(tt)) ...
        '.mat']).chanDat; 
            end
            if strcmp(panelDat.phase, 'sub')
                hfbBand = bandpass(chanDat.enc, [70,150], 1000); 
                broadBand = chanDat.enc;
                tim = chanDat.enctim; 
                hitidx = chanDat.hits & chanDat.use;
                missidx = chanDat.misses & chanDat.use; 
                HFBlatHit = chanDat.HFB_lat.subHit; 
                HFBlatMiss = chanDat.HFB_lat.subMiss; 
            else
                hfbBand = bandpass(chanDat.retOn, [70,150], 1000); 
                broadBand = chanDat.retOn;
                tim = chanDat.retOtim; 
                hitidx = chanDat.retInfo(:,1)==1;
                missidx = chanDat.retInfo(:,1)==2; 
                HFBlatHit = chanDat.HFB_lat.retHit; 
                HFBlatMiss = chanDat.HFB_lat.retMiss; 
            end
        end
        
       
       
        %misses
        curTrial = hfbBand(:,missidx);
        curTrial = curTrial(:, triali); 
        phaseTrial = angle(hilbert(curTrial));
        idx = find(tim>=HFBlatMiss(triali,1),1); 
        [~, troughOff] = min(phaseTrial(idx-10:idx+10));
%         troughOff = 10; %% debug check to see how much trough alignment matters
        idx = idx + troughOff - 10; 
        HFBpeaksMISSES(:, tt) = curTrial(idx-500: idx+500); 
        curTrial = broadBand(:, missidx); 
        curTrial = curTrial(:, triali); 
        BroadBandPeaksMISSES(:, tt) = curTrial(idx-500: idx+500);
        triali = triali +1; 
    end

    
    fig = figure('visible', false, 'position', [0,0,600,600])
    
%     x = 1:401;
    y = mean(BroadBandPeaksHITS,2);
    plot( y, 'color', hitCol, 'linewidth', 4)
    hold on
%     se = std(BroadBandPeaksHITS,[],2) ./ sqrt(size(BroadBandPeaksHITS,2)); 
%     y = [y + se, flip(y) - flip(se)]; 
%     x = [[1:100], flip([1:100])]; 
%     hold on 
%     fill(x, y, hitCol, 'facealpha', errAlpha, 'edgealpha', 0)


    y = mean(BroadBandPeaksMISSES,2);  
    plot( y, 'color', missCol, 'linewidth', 4)
%     se = std(BroadBandPeaksMISSES,[],2) ./ sqrt(size(BroadBandPeaksMISSES,2));  
%     y = [y + se, flip(y) - flip(se)]; 
%     x = [[1:100], flip([1:100])]; 
%     hold on 
%     fill(x, y, missCol, 'facealpha', errAlpha, 'edgealpha', 0)

    xticks([1:50:401])
    xticklabels([-200:50:200])
    xlim([1,401])
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;

    saveas(fig, [figDat 'pubFigs/' 'broadBandHFBpeak_' panelDat.reg '_' panelDat.phase  '.jpg'])

    close(fig)
    catch
        ii
    end
end

%% power spectra at HFB peak

test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TF_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_HFB'));
TF_files_HFB = TF_files(test); 

%grab a .csv of data with the columns: power, frequency, enc/ret, subID,
%chanID
% % % aovDat = table;
% % % aovDat.power = zeros(112200,1); 
% % % aovDat.freq = zeros(112200, 1); 
% % % aovDat.encRet = repmat("askj", 112200,1); 
% % % aovDat.subID = repmat("askj", 112200,1); 
% % % aovDat.reg = repmat("askj", 112200,1); 
% % % ai = 1; 

allSpect = struct; 

for ii = 1:length(TF_files_HFB)
    try
    panelDat = load([TF_files_HFB(ii).folder '/' ...
        TF_files_HFB(ii).name]).outDat;
    ii

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 

    figure('visible', false, 'position', [0,0,600,600])

    imagesc([-500:25:500], [], squeeze(mean(panelDat.hits_hfb,1))')
    yticks([10:10:100])
    yticklabels(round(panelDat.frex([10:10:100])))
    set(gca, 'ydir', 'normal')
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    caxis([0, 10])

%     export_fig([figDat 'pubFigs/' 'SUP2_' panelDat.reg '_HIT_' panelDat.phase  '.jpg'], '-r300')

    figure('visible', false, 'position', [0,0,600,600])

    imagesc([-500:25:500], [], squeeze(mean(panelDat.misses_hfb,1))')
    yticks([10:10:100])
    yticklabels(round(panelDat.frex([10:10:100])))
    set(gca, 'ydir', 'normal')
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    caxis([0, 10])

%     export_fig([figDat 'pubFigs/' 'SUP2_' panelDat.reg '_MISS_' panelDat.phase  '.jpg'], '-r300')


    subIDs = cellfun(@(x) split(x, '_'), panelDat.realID, ...
                         'UniformOutput', false); 
    subIDs = cellfun(@(x) x{1}, subIDs, 'uniformoutput', false); 
    uniIDs = unique(subIDs); 
  


    figure('visible', false, 'position', [0,0,600,600])
    x = 1:length(panelDat.frex); 
%     yline(0, 'color', 'k', 'linewidth', linWid)
    x(panelDat.p_hfb(21,:)>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
    breakVals = [0, find(diff(x)> 25), length(x)];
    for jj = 1:length(breakVals)-1
        varName = ['peakPow_' panelDat.reg '_' panelDat.phase num2str(jj)];
        allSig.(varName) = nan(length(allSig.subID),1); 
        tmpX = x(breakVals(jj)+1:breakVals(jj+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
        %store raw differences for significant time periods: 
        chanMeans = mean(panelDat.hits_hfb(:,21,...
                            breakVals(jj)+1:breakVals(jj+1)) ...
                            - ...
                         panelDat.misses_hfb(:,21,...
                            breakVals(jj)+1:breakVals(jj+1)), 2);
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
    
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
    allMin = min([allMin, -2]); 
    ylim([-1, 18])
    xlim([1, 100])
    xticks([1:11:100])
    xticklabels(round(panelDat.frex([1:11:100])))
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
%     export_fig([figDat 'pubFigs/' 'Fig2_' panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')



    realID = arrayfun(@(x) [panelDat.hitSub{x} '_'...
        num2str(panelDat.hitChi(x))], 1:length(panelDat.hitChi), ...
        'UniformOutput', false);
    uniIDs = unique(realID); 
    comboSpect = (hitSpect + missSpect) /2; 
    
% % %     for subi = 1:size(comboSpect,1)
% % %         aovDat.power(ai:ai+99) = comboSpect(subi,:); 
% % %         aovDat.freq(ai:ai+99) = panelDat.frex; 
% % %         aovDat.encRet(ai:ai+99) = repmat(panelDat.phase, 100,1);
% % %         aovDat.subID(ai:ai+99) = repmat(uniIDs{subi}, 100,1); 
% % %         aovDat.reg(ai:ai+99) = repmat(panelDat.reg, 100,1); 
% % %         ai = ai+100; 
% % %     end

% Get peak frequencies for description of power spectra 
    curSpect = mean([missSpect; hitSpect]);

    [val, idx] = max(curSpect(1:64));
    
    allSpect(ii).reg = panelDat.reg; 
    allSpect(ii).phase = panelDat.phase;
    allSpect(ii).lowPeak = panelDat.frex(idx);
    
    rangeBot = find(diff(movmean(curSpect(idx+4:end), 4))>0, 1) + idx + 2;

    [val, idx] = max(curSpect(rangeBot:77));

    if val>(curSpect(rangeBot)+1) && val>curSpect(77) %edges aren't peaks
        allSpect(ii).highPeak = panelDat.frex(idx+rangeBot-1);
    end

    catch
        ii
    end
end

% % % writetable(aovDat, ...
% % %     ['R:\MSS\Johnson_Lab\dtf8829\GitHub\' ...
% % %     'HpcAccConnectivityProject\HFBpeakPowerSpectra.csv'])


%% ITPC spectra at HFB peak 

test = cellfun(@(x) length(x)>0, ...
    strfind({fig2Dat.name}, 'TFphase_'));
TF_files = fig2Dat(test); 

test = cellfun(@(x) length(x)>0, ...
    strfind({TF_files.name}, '_HFB'));
TF_files_HFB = TF_files(test); 


ITPCindex = struct; 



for ii = 1:length(TF_files_HFB)
    try
    panelDat = load([TF_files_HFB(ii).folder '/' ...
        TF_files_HFB(ii).name]).outDat;
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 
    
    subIDs = cellfun(@(x) split(x, '_'), panelDat.realID, ...
                         'UniformOutput', false); 
    subIDs = cellfun(@(x) x{1}, subIDs, 'uniformoutput', false); 
    uniIDs = unique(subIDs); 
    
    figure('visible', false, 'position', [0,0,600,600])
    x = 1:length(panelDat.frex); 
    yline(0, 'color', 'k', 'linewidth', linWid)
   
    x(panelDat.p_hfb(21,:)>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
        disp([panelDat.reg ' ' panelDat.phase ' ' num2str(panelDat.frex(max(x)))])
    breakVals = [0, find(diff(x)> 1), length(x)];
    for jj = 1:length(breakVals)-1
        tmpX = x(breakVals(jj)+1:breakVals(jj+1));
        tmpY = ones(length(tmpX),1) * 10000; 
        tmpX = [tmpX, flip(tmpX)]; 
        tmpY = [tmpY', -tmpY'];
        fill(tmpX, tmpY, sigCol, 'facealpha', sigAlpha, 'edgealpha', 0)
        hold on 
        xline(12, 'color', 'k', 'linewidth',2, 'linestyle', '--')
        xline(33, 'color', 'k', 'linewidth',2, 'linestyle', '--')
        varName = ['peakITPC_' panelDat.reg '_' panelDat.phase num2str(jj)];
        allSig.(varName) = nan(length(allSig.subID),1);
        %store raw differences for significant time periods: 
        chanMeans = mean(panelDat.hits_hfb(:,21,...
                            breakVals(jj)+1:breakVals(jj+1)) ...
                            - ...
                         panelDat.misses_hfb(:,21,...
                            breakVals(jj)+1:breakVals(jj+1)), 2);
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
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

    ITPCindex(ii).reg = panelDat.reg; 
    ITPCindex(ii).phase = panelDat.phase; 
    [~, fi] = max(abs(mean(hitSpect - missSpect,1)));
    ITPCindex(ii).HFBfreq_hit = panelDat.frex(fi); 
    ITPCindex(ii).HFBfi_hit = fi; 
    ITPCindex(ii).HFBhit = hitSpect(:,fi);

%     [~, fi] = max(abs(mean(missSpect,1)));
    ITPCindex(ii).HFBfreq_miss = panelDat.frex(fi); 
    ITPCindex(ii).HFBfi_miss = fi; 
    ITPCindex(ii).HFBmiss = missSpect(:,fi); 

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


for ii = 1:length(TF_files_image)
   ii
    panelDat = load([TF_files_image(ii).folder '/' ...
        TF_files_image(ii).name]).outDat;
 
    subIDs = cellfun(@(x) split(x, '_'), panelDat.realID, ...
                         'UniformOutput', false); 
    subIDs = cellfun(@(x) x{1}, subIDs, 'uniformoutput', false); 
    uniIDs = unique(subIDs);

    if isfield(panelDat.clusterinfo, 'pos_clusters')
        pos_clusters = panelDat.clusterinfo.pos_clusters; 
        for clu = 1:length(pos_clusters)
            if ~isempty(pos_clusters(clu).p)
            idx = pos_clusters(clu).inds; 
            L = sum(idx, 'all');
            n = size(panelDat.hits_image,1);
            chanMeans = []; 
            for chan = 1:n
                slice = squeeze(panelDat.hits_image(chan,:,:) - ...
                    panelDat.misses_image(chan,:,:));
                chanMeans = [chanMeans, mean(slice(idx))]; 
            end

            varName = ['imagePow_' panelDat.reg '_' ...
                panelDat.phase num2str(clu)];
            allSig.(varName) = nan(length(allSig.subID),1);
       
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
            end
        end


    end

    if isfield(panelDat.clusterinfo, 'neg_clusters')
        neg_clusters = panelDat.clusterinfo.neg_clusters; 
        for clu = 1:length(neg_clusters)
            if ~isempty(neg_clusters(clu).p)
            idx = neg_clusters(clu).inds; 
            L = sum(idx, 'all');
            n = size(panelDat.hits_image,1);

            chanMeans = []; 
            for chan = 1:n
                slice = squeeze(panelDat.hits_image(chan,:,:) - ...
                    panelDat.misses_image(chan,:,:));
                chanMeans = [chanMeans, mean(slice(idx))]; 
            end

            varName = ['imagePow_neg_' panelDat.reg '_' ...
                panelDat.phase num2str(clu)];
            allSig.(varName) = nan(length(allSig.subID),1);
        %store raw differences for significant time periods: 
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
            end
        end


    end



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


for ii = 1:length(TF_files_image)
   
    panelDat = load([TF_files_image(ii).folder '/' ...
        TF_files_image(ii).name]).outDat;
 
    %update allSig
    subIDs = cellfun(@(x) split(x, '_'), panelDat.realID, ...
                         'UniformOutput', false); 
    subIDs = cellfun(@(x) x{1}, subIDs, 'uniformoutput', false); 
    uniIDs = unique(subIDs);

    if isfield(panelDat.clusterinfo, 'pos_clusters')
        pos_clusters = panelDat.clusterinfo.pos_clusters; 
        for clu = 1:length(pos_clusters)
            if ~isempty(pos_clusters(clu).p)
            idx = pos_clusters(clu).inds; 
            L = sum(idx, 'all');
            n = size(panelDat.hits_image,1);

            chanMeans = []; 
            for chan = 1:n
                slice = squeeze(panelDat.hits_image(chan,:,:) - ...
                    panelDat.misses_image(chan,:,:));
                chanMeans = [chanMeans, mean(slice(idx))]; 
            end


           
            varName = ['imageITPC_' panelDat.reg '_' ...
                panelDat.phase num2str(clu)];
            allSig.(varName) = nan(length(allSig.subID),1);
        %store raw differences for significant time periods: 
        
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
            end
        end


    end

    if isfield(panelDat.clusterinfo, 'neg_clusters')
        neg_clusters = panelDat.clusterinfo.neg_clusters; 
        for clu = 1:length(neg_clusters)
            if ~isempty(neg_clusters(clu).p)
            idx = neg_clusters(clu).inds; 
            L = sum(idx, 'all');
            n = size(panelDat.hits_image,1);

            chanMeans = []; 
            for chan = 1:n
                slice = squeeze(panelDat.hits_image(chan,:,:) - ...
                    panelDat.misses_image(chan,:,:));
                chanMeans = [chanMeans, mean(slice(idx))]; 
            end

          
            varName = ['imageITPC_neg_' panelDat.reg '_' ...
                panelDat.phase num2str(clu)];
            allSig.(varName) = nan(length(allSig.subID),1);
        %store raw differences for significant time periods: 
        
        for sub = 1:length(uniIDs)
            %store subject means into allSig for later correlation to
            %memory
            idx = find(strcmp(allSig.reg, panelDat.reg) & ...
                       strcmp(allSig.subID, uniIDs{sub}));
             allSig.(varName)(idx) = ...
                 mean(chanMeans(ismember(subIDs, uniIDs{sub})));

        end
            end
        end


    end

    % end add allSig



    [~, timMax] = max(abs(mean(panelDat.hits_image,[1,3])));
    [~, timMax2] = max(abs(mean(panelDat.misses_image,[1,3])));
    ITPCindex(ii).IMGhit = panelDat.hits_image(:,timMax,ITPCindex(ii).HFBfi_hit);
    ITPCindex(ii).IMGmiss = panelDat.misses_image(:,timMax2,ITPCindex(ii).HFBfi_miss);
    
   
    ITPCindex(ii).IMGtim_hit = panelDat.tim(timMax); 
    ITPCindex(ii).IMGtim_miss = panelDat.tim(timMax2); 

    ITPCindex(ii).idxHIT = getScaledIndex(ITPCindex(ii).HFBhit, ITPCindex(ii).IMGhit);
    ITPCindex(ii).idxMISS = getScaledIndex(ITPCindex(ii).HFBmiss, ITPCindex(ii).IMGmiss);

    subID = cellfun(@(x,y) [x '_' num2str(y)], panelDat.hitSub, num2cell(panelDat.hitChi), 'uniformoutput', false);
    ITPCindex(ii).subID = unique(subID); 


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


%% reformat HFB vs. IMG index into table for .csv export

aovDat = table;
aovDat.index = zeros(1000,1); 
aovDat.freq = zeros(1000, 1); 
aovDat.encRet = repmat("askj", 1000,1); 
aovDat.subID = repmat("askj", 1000,1); 
aovDat.hitMiss = repmat("askj", 1000,1); 
aovDat.reg = repmat("askj", 1000,1); 
ai = 1; 


for ii = 1:length(ITPCindex)
   L = length(ITPCindex(ii).idxHIT); 
   %hits
   aovDat.index(ai:ai+L-1) = ITPCindex(ii).idxHIT; 
   aovDat.freq(ai:ai+L-1) = ones(L,1)*ITPCindex(ii).HFBfreq_hit;
   aovDat.encRet(ai:ai+L-1) = repmat(ITPCindex(ii).phase, L,1);
   aovDat.subID(ai:ai+L-1) = ITPCindex(ii).subID;
   aovDat.hitMiss(ai:ai+L-1) = repmat('hit', L,1);
   aovDat.reg(ai:ai+L-1) = repmat(ITPCindex(ii).reg, L,1);
   ai = ai + L; 
   
   %misses
   aovDat.index(ai:ai+L-1) = ITPCindex(ii).idxMISS; 
   aovDat.freq(ai:ai+L-1) = ones(L,1)*ITPCindex(ii).HFBfreq_hit;
   aovDat.encRet(ai:ai+L-1) = repmat(ITPCindex(ii).phase, L,1);
   aovDat.subID(ai:ai+L-1) = ITPCindex(ii).subID;
   aovDat.hitMiss(ai:ai+L-1) = repmat('miss', L,1);
   aovDat.reg(ai:ai+L-1) = repmat(ITPCindex(ii).reg, L,1);
   ai = ai + L;


end

writetable(aovDat, ...
    ['R:\MSS\Johnson_Lab\dtf8829\GitHub\' ...
    'HpcAccConnectivityProject\HFB_IMG_index.csv'])

%% fig 3 
%connectivity 

fig3Dat = dir([figDat 'Figure3']); 
fig3Dat(1:2) = []; 


%% grab all connections: 
% 
% aovDat = table;
% aovDat.PPC = zeros(10000,1); 
% aovDat.freq = zeros(10000, 1); 
% aovDat.tim = zeros(10000, 1); 
% aovDat.encRet = repmat("askj", 10000,1); 
% aovDat.pairID = repmat("askj", 10000,1); 
% aovDat.hitMiss = repmat("askj", 10000,1); 
% aovDat.reg1 = repmat("askj", 10000,1); 
% aovDat.reg2 = repmat("askj", 10000,1); 
% aovDat.IMG_HFB = repmat("askj", 10000,1); 
% ai = 1; 

%all connections in a reg X reg X time X hit/miss/t/p X enc/ret
allConnections = zeros(9,9,139, 20, 4, 2); 
allConnections2 = zeros(9,9,41,20,4,2); 
for ii = 1:length(fig3Dat)
    curDat = load([fig3Dat(ii).folder '/' fig3Dat(ii).name]).outDat; 
    ii
% % %      %update allSig
% % %     
% % %     subIDs = curDat.subVals; 
% % %     uniIDs = unique(subIDs);
% % %     varCodes = {'enc_image', 'ret_image', 'enc_HFB', 'ret_HFB'};
% % %     direction = {'pos', 'neg'}; 
% % %     for d = 1:2
% % %     for cc = 1:4
% % %     %do encoding image locked
% % %     if isfield(curDat.([varCodes{cc} '_clust']), [direction{d} '_clusters'])
% % %         clust = curDat.([varCodes{cc} '_clust']).([direction{d} '_clusters']); 
% % %         clust(:,[clust.p]>.05) = [];
% % %         for clu = 1:length(clust)
% % %             if ~isempty(clust(clu).p)
% % %             idx = clust(clu).inds; 
% % %             L = sum(idx, 'all');
% % %             n = size(curDat.([varCodes{cc} '_hitVals']),3);
% % % 
% % %             chanMeans = []; 
% % %             chanMeansHit = []; 
% % %             chanMeansMiss = []; 
% % %             for chan = 1:n
% % %                 sliceHit = squeeze(curDat.([varCodes{cc} '_hitVals'])(:,:,chan));
% % %                 sliceMiss= squeeze(curDat.([varCodes{cc} '_missVals'])(:,:,chan));
% % %                 chanMeansHit = [chanMeansHit, mean(sliceHit(idx))];
% % %                 chanMeansMiss = [chanMeansMiss, mean(sliceMiss(idx))]; 
% % %                 slice = sliceHit - sliceMiss; 
% % %                 chanMeans = [chanMeans, mean(slice(idx))]; 
% % %             end
% % % 
% % % 
% % % 
% % % 
% % % % make linked boxplot 
% % %     figure('visible', false, 'position', [0,0,600,400])
% % %     hmSort = [zeros(1,n), ones(1,n)]; 
% % %     b1 = boxchart(zeros(1,n), chanMeansMiss, 'boxfacecolor', missCol);
% % %     hold on 
% % %     b2 = boxchart(ones(1,n), chanMeansHit, 'boxfacecolor', hitCol);
% % %     allChanMeans = [chanMeansMiss, chanMeansHit]; 
% % %     b1.MarkerStyle = 'none'; 
% % %     b2.MarkerStyle = 'none'; 
% % %     b1.LineWidth = 3;  % This changes the outline of the box, but not the whiskers
% % %     b2.LineWidth = 3;  % Same here
% % %     xticks([0,1])
% % %     xticklabels({'Miss', 'Hit'})
% % % %             meanX = round(mean(Xidx(tmp), 'all')); 
% % % %             meanY = round(mean(Yidx(tmp), 'all')); 
% % % %             LLvals = -150:150; 
% % % 
% % %     IQR = prctile(chanMeansMiss, [25,75]);
% % %     Q25 = IQR(1); 
% % %     Q75 = IQR(2); 
% % %     IQR = diff(IQR); 
% % %     tmp = chanMeansMiss(chanMeansMiss>(Q25-1.5*IQR) & chanMeansMiss<(Q75+1.5*IQR));
% % %     % Overlay custom error bars
% % %     % For Miss group
% % %     line([-.15,.15], [max(tmp), max(tmp)],  'LineWidth', 3, 'Color', missCol);
% % %     line([0,0], [Q75, max(tmp)],  'LineWidth', 3, 'Color', missCol);
% % %     
% % %     line([-.15,.15], [min(tmp), min(tmp)],  'LineWidth', 3, 'Color', missCol);
% % %     line([0,0], [Q25, min(tmp)],  'LineWidth', 3, 'Color', missCol);
% % %     
% % %     IQR = prctile(chanMeansHit, [25,75]);
% % %     Q25 = IQR(1); 
% % %     Q75 = IQR(2); 
% % %     IQR = diff(IQR); 
% % %     tmp = chanMeansHit(chanMeansHit>(Q25-1.5*IQR) & chanMeansHit<(Q75+1.5*IQR));
% % %     % Overlay custom error bars
% % %     % For hit group
% % %     line([1-.15,1.15], [max(tmp), max(tmp)],  'LineWidth', 3, 'Color', hitCol);
% % %     line([1,1], [Q75, max(tmp)],  'LineWidth', 3, 'Color', hitCol);
% % %     
% % %     line([1-.15,1.15], [min(tmp), min(tmp)],  'LineWidth', 3, 'Color', hitCol);
% % %     line([1,1], [Q25, min(tmp)],  'LineWidth', 3, 'Color', hitCol);
% % %     PLH = allChanMeans(hmSort==1); 
% % %     PLM = allChanMeans(hmSort==0); 
% % %     randVals = (rand(length(PLH),1)-.5)*.5;
% % %     hold on 
% % %     scatter(randVals, PLM, 10,  'blue')
% % %     scatter(randVals+1, PLH, 10, 'blue')
% % %     
% % %     for ppi = 1:length(PLH)
% % %         plot([0+randVals(ppi),1+randVals(ppi)], [PLM(ppi),PLH(ppi)], 'color', [.2,.1,.5,.2], 'linewidth', .5 )
% % %         
% % % 
% % %     end
% % % 
% % %     %make subject means
% % %     subHits = zeros(length(uniIDs),1); 
% % %     subMisses = zeros(length(uniIDs),1); 
% % %     uniqueSubs = uniIDs; 
% % %     regSubs = subIDs; 
% % %     for sub = 1:length(subHits)
% % %         subidx = cellfun(@(x) strcmp(x, uniqueSubs{sub}), regSubs); 
% % %         subHits(sub) = mean(PLH(subidx));
% % %         subMisses(sub) = mean(PLM(subidx));
% % %         plot([0,1], [subMisses(sub), subHits(sub)], 'color', 'k')
% % % 
% % %     end
% % %     scatter(ones(length(uniqueSubs),1), subHits, 35,  'red', 'filled')
% % %     scatter(zeros(length(uniqueSubs),1), subMisses, 35,  'red', 'filled')
% % %     xlim([-.5, 1.5])
% % % %     ylim([-.05, min([.3; max([PLH, PLM])+.05])])
% % % 
% % % 
% % % 
% % % set(gcf,'color','w');
% % % box off;
% % % ax=gca;ax.LineWidth=4;
% % % export_fig([figDat 'pubFigs/' 'ppc_' direction{d} '_' varCodes{cc} '_' ...
% % %                 curDat.reg1 '_' curDat.reg2 '_' num2str(clu) '.jpg'], '-r300')
% % % 
% % % 
% % % % end linked boxplot
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % %          
% % %             varName = ['ppc_' direction{d} '_' varCodes{cc} '_' ...
% % %                 curDat.reg1 '_' curDat.reg2 '_' num2str(clu)];
% % %             altVar = ['ppc_' direction{d} '_' varCodes{cc} '_' ...
% % %                 curDat.reg2 '_' curDat.reg1 '_' num2str(clu)];
% % %             allSig.(varName) = nan(length(allSig.subID),1);
% % %         %store raw differences for significant time periods: 
% % %         
% % %         for sub = 1:length(uniIDs)
% % %             %store subject means into allSig for later correlation to
% % %             %memory
% % %             if cc<3 %image connections are symmetric but HFB are not! 
% % %                 idx = find((strcmp(allSig.reg, curDat.reg1) & ...
% % %                            strcmp(allSig.subID, uniIDs{sub})) |...
% % %                            (strcmp(allSig.reg, curDat.reg2) & ...
% % %                            strcmp(allSig.subID, uniIDs{sub})));
% % %             else
% % %                 
% % %                 idx = find((strcmp(allSig.reg, curDat.reg1) & ...
% % %                            strcmp(allSig.subID, uniIDs{sub})));
% % %             end
% % %              allSig.(varName)(idx) = ...
% % %                  mean(chanMeans(ismember(subIDs, uniIDs{sub})));
% % % 
% % %         end
% % %             end
% % %         end
% % % 
% % % 
% % %     end
% % % 
% % %     end
% % %     end
% % %     %end update allSig

    tim = curDat.enc_image_tim; 
    reg1 = find(cellfun(@(x) strcmp(x, curDat.reg1), regions)); 
    reg2 = find(cellfun(@(x) strcmp(x, curDat.reg2), regions)); 
    %ENCODING: 
    %hit vals
    tmpMat = curDat.enc_image_hitVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 1, 1) = mean(tmpMat,3); 

    %hit image figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_image_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    pMat = curDat.enc_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    figure('visible', false, 'position', [0,0,600,400])
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_image_hit_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    tmpMat = curDat.enc_HFB_hitVals; 
    allConnections2(reg1, reg2, :, :, 1, 1) = mean(tmpMat,3); 

    %hit HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_HFB_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    pMat = curDat.enc_HFB_p;
    ss = 1; 
    figure('visible', false, 'position', [0,0,600,400])
    makeConnectivityHeatMap2(inMat, pMat, ...
                            ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_HFB_hit_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    %miss vals
    tmpMat = curDat.enc_image_missVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 2, 1) = mean(tmpMat,3); 


    %miss image figure
     %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_image_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.enc_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_image_miss_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    tmpMat = curDat.enc_HFB_missVals; 
    allConnections2(reg1, reg2, :, :, 2, 1) = mean(tmpMat,3);

    %miss HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.enc_HFB_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000);  
    pMat = curDat.enc_HFB_p;
    ss = 1; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_HFB_miss_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    %t vals
    tmpMat = curDat.enc_image_tVal; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 3, 1) = tmpMat; 

    %image t values
    %ss: 1 = HFB, 0 = image
    inMat = curDat.enc_image_tVal;  
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.enc_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_image_tVal_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')



    tmpMat = curDat.enc_HFB_tVal; 
    allConnections2(reg1, reg2, :, :, 3, 1) = tmpMat;

    %HFB t values 
     %ss: 1 = HFB, 0 = image
    inMat = curDat.enc_HFB_tVal;  
    pMat = curDat.enc_HFB_p; 
    ss = 1; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_HFB_tVal_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %p vals
    tmpMat = curDat.enc_image_p; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 4, 1) = tmpMat; 

    tmpMat = curDat.enc_HFB_p; 
    allConnections2(reg1, reg2, :, :, 4, 1) = tmpMat;


    %update the aovDat for export IMAGE data
%     pmask = squeeze(allConnections(reg1, reg2, :,:,4,1));
%     if sum(pmask<.05, 'all') > 0 %only collecting significant connections
%         %hits
%         hitMat = curDat.enc_image_hitVals; 
%         hitMat(tim<-450 | tim>3000, :, :) = []; 
%         hitMat = reshape(hitMat, [prod(size(hitMat,[1,2])), size(hitMat,3)]);
%         %misses
%         missMat = curDat.enc_image_missVals; 
%         missMat(tim<-450 | tim>3000, :, :) = []; 
%         missMat = reshape(missMat, [prod(size(missMat,[1,2])), size(missMat,3)]);
%         %get the clusters of significance
%         pmask = pmask<.05; 
%         [L, num] = bwlabel(pmask, 8);
%         L = L(:);
%         freqMat = repmat(curDat.frex, size(pmask,1),1);
%         freqMat = freqMat(:); 
%         timMat = repmat(tim(tim>=-450 & tim<=3000)', 1,size(pmask,2));
%         timMat = timMat(:); 
%         for mi = 1:num %loop over the clusters of significance
%             N = size(hitMat, 2); 
%            
%             aovDat.PPC(ai:ai+N-1) =  mean(hitMat(L==mi, :),1); 
%             aovDat.freq(ai:ai+N-1) = mean(freqMat(L==mi)); 
%             aovDat.tim(ai:ai+N-1) = mean(timMat(L==mi)); 
%             aovDat.encRet(ai:ai+N-1) = repmat("enc", N,1); 
%             aovDat.pairID(ai:ai+N-1) = curDat.pairIDs; 
%             aovDat.hitMiss(ai:ai+N-1) = repmat("hit", N,1); 
%             aovDat.reg1(ai:ai+N-1) = repmat(curDat.reg1, N,1); 
%             aovDat.reg2(ai:ai+N-1) = repmat(curDat.reg2, N,1); 
%             aovDat.IMG_HFB(ai:ai+N-1) = repmat("IMG", N,1); 
% 
%         end
% 
%     end

    %RETRIEVAL: 
    tim = curDat.ret_image_tim; 
    %hit vals
    tmpMat = curDat.ret_image_hitVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 1, 2) = mean(tmpMat,3); 

    %hit image figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_image_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    pMat = curDat.ret_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    % % % figure('visible', false, 'position', [0,0,600,400])
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_ret_image_hit_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')


    tmpMat = curDat.ret_HFB_hitVals; 
    allConnections2(reg1, reg2, :, :, 1, 2) = mean(tmpMat,3);

    %hit HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_HFB_hitVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    pMat = curDat.ret_HFB_p;
    ss = 1; 
    % % % figure('visible', false, 'position', [0,0,600,400])
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_ret_HFB_hit_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %miss vals
    tmpMat = curDat.ret_image_missVals; 
    tmpMat(tim<-450 | tim>3000, :, :) = []; 
    allConnections(reg1, reg2, :, :, 2, 2) = mean(tmpMat,3); 


    %miss image figure
     %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_image_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.ret_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_ret_image_miss_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    tmpMat = curDat.ret_HFB_missVals; 
    allConnections2(reg1, reg2, :, :, 2, 2) = mean(tmpMat,3); 

    %miss HFB figure
    %ss: 1 = HFB, 0 = image
    inMat = squeeze(mean(curDat.ret_HFB_missVals,3)); 
    inMat(inMat<0) = 0; 
    inMat = sqrt(inMat); 
    pltim = tim(tim>=-450 & tim<=3000);  
    pMat = curDat.ret_HFB_p;
    ss = 1; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_ret_HFB_miss_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %t vals
    tmpMat = curDat.ret_image_tVal; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 3, 2) = tmpMat; 

     %image t values
    %ss: 1 = HFB, 0 = image
    inMat = curDat.ret_image_tVal;  
    pltim = tim(tim>=-450 & tim<=3000); 
    inMat(tim<-450 | tim>3000, :) = []; 
    pMat = curDat.ret_image_p; 
    pMat(tim<-450 | tim>3000, :) = []; 
    ss = 0; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_ret_image_tVal_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    tmpMat = curDat.ret_HFB_tVal; 
    allConnections2(reg1, reg2, :, :, 3, 2) = tmpMat; 

     %HFB t values 
     %ss: 1 = HFB, 0 = image
    inMat = curDat.ret_HFB_tVal;  
    pMat = curDat.ret_HFB_p; 
    ss = 1; 
    % % % figure('visible', false, 'position', [0,0,600,400]);
    % % % makeConnectivityHeatMap2(inMat, pMat, ...
    % % %                         ss, curDat.frex, pltim)
    % % % export_fig([figDat 'pubFigs/' 'sup3_' 'connectionHeatmap_ret_HFB_tVal_'...
    % % %     regions{reg1} '_' regions{reg2} '.jpg'], '-r300')

    %p vals
    tmpMat = curDat.ret_image_p; 
    tmpMat(tim<-450 | tim>3000, :) = []; 
    allConnections(reg1, reg2, :, :, 4, 2) = tmpMat; 

    tmpMat = curDat.ret_HFB_p; 
    allConnections2(reg1, reg2, :, :, 4, 2) = tmpMat; 

end

% writetable(allSig, ...
%     [codePre 'HpcAccConnectivityProject\allSig.csv'])

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
                       tim>=500 & tim<=1500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc2.jpg'], fn2, fn3)

%late encoding
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=1500 & tim<=2500, :, :, 1));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'enc3.jpg'], fn2, fn3)

%early retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=-250 & tim<=500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret1.jpg'], fn2, fn3)

%mid retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=500 & tim<=1500, :, :, 2));
plotConnectionSchematic(plotMat, curDat.frex, "tVal", ...
    regions, keyRegIdx, 'freq', regColors, false, [fn 'ret2.jpg'], fn2, fn3)

%late retrieval
plotMat = squeeze(allConnections(keyRegIdx, keyRegIdx,...
                       tim>=1500 & tim<=2500, :, :, 2));
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



%% Fig 4
%graph based connectivity metrics getting a .csv for R analysis

graphDat = dir([datPre 'graphAnalysis/out']);
L = length(graphDat); 
graphDat(1:2) = [];
ii = 1;
cur = load([graphDat(ii).folder '/' graphDat(ii).name]).outDat;
graphDat(ii).reg1 = cur.reg1; 
graphDat(ii).reg2 = cur.reg2; 
graphDat(ii).timeSet = cur.HFB_Image;
graphDat(ii).time = cur.time; 
graphDat(ii).freq = cur.freq; 
graphDat(ii).encRet = cur.encRet; 
graphDat(ii).hitBC = cur.hitBC(cur.chi); 
graphDat(ii).missBC = cur.missBC(cur.chi); 
graphDat(ii).hitST = cur.hitST(cur.chi); 
graphDat(ii).missST = cur.missST(cur.chi); 
for ii = 1:L
    tic
    cur = load([graphDat(ii).folder '/' graphDat(ii).name]).outDat;

    cur.hitSparce = cur.hitMat; 
    cur.missSparce = cur.missMat; 
    cur.hitSparce(cur.hitMat<.1) = 0; 
    cur.missSparce(cur.missMat<.1) = 0; 
    cur.hitSparce(cur.hitSparce>0) = 1; 
    cur.missSparce(cur.missSparce>0) = 1; 

    hitPos = cur.hitMat; 
    hitPos(hitPos<0) = 0; 
    missPos = cur.missMat; 
    missPos(missPos<0) = 0; 

    
    graphDat(ii).reg1 = cur.reg1; 
    graphDat(ii).reg2 = cur.reg2; 
    graphDat(ii).subID = cur.subID; 
    graphDat(ii).chi = cur.chi; 
    graphDat(ii).chi2 = cur.chi2; 
    graphDat(ii).timeSet = cur.HFB_Image;
    graphDat(ii).time = cur.time; 
    graphDat(ii).freq = cur.freq; 
    graphDat(ii).encRet = cur.encRet; 
    graphDat(ii).hitBC = cur.hitBC(cur.chi); 
    graphDat(ii).missBC = cur.missBC(cur.chi); 
    graphDat(ii).hitST = cur.hitST(cur.chi); 
    graphDat(ii).missST = cur.missST(cur.chi); 
    graphDat(ii).hitSaturation = sum(cur.hitMat(cur.chi,:)>.1) / ...
                                    size(cur.hitMat,1); 
    graphDat(ii).missSaturation = sum(cur.missMat(cur.chi,:)>.1) / ...
                                    size(cur.hitMat,1);
    graphDat(ii).hitEff = efficiency_bin(cur.hitSparce); 
    graphDat(ii).missEff = efficiency_bin(cur.missSparce); 
    graphDat(ii).hitChar = charpath(distance_wei(1./hitPos)); 
    graphDat(ii).missChar = charpath(distance_wei(1./missPos)); 

    if strcmp(cur.HFB_Image, 'image') %also put the flip in there for image connections
        graphDat(ii+L).reg1 = cur.reg2; 
        graphDat(ii+L).reg2 = cur.reg1; 
        graphDat(ii+L).subID = cur.subID; 
        graphDat(ii+L).chi = cur.chi; 
        graphDat(ii+L).chi2 = cur.chi2; 
        graphDat(ii+L).timeSet = cur.HFB_Image;
        graphDat(ii+L).time = cur.time; 
        graphDat(ii+L).freq = cur.freq; 
        graphDat(ii+L).encRet = cur.encRet; 
        graphDat(ii+L).hitBC = cur.hitBC(cur.chi); 
        graphDat(ii+L).missBC = cur.missBC(cur.chi); 
        graphDat(ii+L).hitST = cur.hitST(cur.chi); 
        graphDat(ii+L).missST = cur.missST(cur.chi); 
        graphDat(ii+L).hitSaturation = sum(cur.hitMat(cur.chi,:)>.1) / ...
                                        size(cur.hitMat,1); 
        graphDat(ii+L).missSaturation = sum(cur.missMat(cur.chi,:)>.1) / ...
                                        size(cur.hitMat,1);
        graphDat(ii+L).hitEff = efficiency_bin(cur.hitSparce); 
        graphDat(ii+L).missEff = efficiency_bin(cur.missSparce); 
        graphDat(ii+L).hitChar = charpath(distance_wei(1./hitPos)); 
        graphDat(ii+L).missChar = charpath(distance_wei(1./missPos)); 


    end

toc
end

test = cellfun(@(x) isempty(x), {graphDat.reg1});
graphDat(test) = []; 

t = struct2table(graphDat);

% Save the table as a CSV file
writetable(t, [figDat 'graphDat.csv']);




%% Fig 4, make an example figure showing how graph measures are calcualted

pPFC_enc_HFB_files = {'pPFC_acc_enc_HFB_1_1.mat',...
                      'pPFC_acc_enc_HFB_1_2.mat',...
                      'pPFC_acc_enc_HFB_1_3.mat',...
                      'pPFC_acc_enc_HFB_1_4.mat',...
                      'pPFC_acc_enc_HFB_1_5.mat',...
                      'pPFC_acc_enc_HFB_1_6.mat',...
                      'pPFC_acc_enc_HFB_1_7.mat',...
                      'pPFC_acc_enc_HFB_1_8.mat',...
                      'pPFC_acc_enc_HFB_1_9.mat',...
                      'pPFC_acc_enc_HFB_1_10.mat',...
                      'pPFC_acc_enc_HFB_1_11.mat'};

for ii = 1:length(pPFC_enc_HFB_files)

exampDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\graphAnalysis\out\'...
                 pPFC_enc_HFB_files{ii}]).outDat;
hipModel =  stlread('R:\MSS\Johnson_Lab\dtf8829\GitHub\HpcAccConnectivityProject\bilatHippo.stl');
%get the 3d locations of the electrodes: 
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\chanDat_'...
                exampDat.subID '_001.mat']).chanDat;
exampDat.chanpos = chanDat.chanpos; 
set(0, 'defaultfigurewindowstyle', 'normal')
missName = ['Fig4_miss_' num2str(ii) '_pPFC_acc_' exampDat.subID];
hitName = ['Fig4_hit_' num2str(ii) '_pPFC_acc_' exampDat.subID]; 

ftpath = [codePre 'fieldtrip-20230118']; 
yeo7 = ft_read_atlas(...
    fullfile(ftpath, ...
    'template/atlas/yeo/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'));
hold on 
test = yeo7.tissue>0; 
indices = zeros(1003785, 3); 
vals = zeros(1003785, 1); 
indi = 1; 
for ii = 1:256
    for jj = 1:256
        for kk = 1:256
            if test(ii,jj,kk)
                indices(indi, :) = [ii*-1+127,...
                                    kk*-1+145,...
                                    jj*-1+147]; 
                vals(indi) = yeo7.tissue(ii, jj, kk); 
                indi = indi + 1; 
            end
        end
    end
end
figure('visible', false, 'position', [50,50,1000,1000])
test = randsample(1:indi, 10000, false); 
hold off
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')
chi = exampDat.chi; 
chi2 = exampDat.chi2; 
scatter3(exampDat.chanpos(chi,1), ...
         -exampDat.chanpos(chi,2), ...
         exampDat.chanpos(chi,3),'red', 'filled')
scatter3(exampDat.chanpos(chi2,1), ...
         -exampDat.chanpos(chi2,2), ...
         exampDat.chanpos(chi2,3),'green', 'filled')
plot3([exampDat.chanpos(chi,1), exampDat.chanpos(chi2,1)], ...
      [-exampDat.chanpos(chi,2), -exampDat.chanpos(chi2,2)], ...
      [exampDat.chanpos(chi,3), exampDat.chanpos(chi2,3)], ...
      'color', 'red', 'linewidth', 3)






view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 

hitBin = exampDat.hitMat>.10; 
missBin = exampDat.missMat>.10; 

for ii = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if hitBin(ii,jj)
            plot3([exampDat.chanpos(ii,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(ii,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(ii,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end

export_fig([figDat 'pubFigs/' hitName '.jpg'], '-r300')


figure('visible', false, 'position', [50,50,1000,1000])
hold off
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')
chi = exampDat.chi; 
chi2 = exampDat.chi2; 
scatter3(exampDat.chanpos(chi,1), ...
         -exampDat.chanpos(chi,2), ...
         exampDat.chanpos(chi,3),'red', 'filled')
scatter3(exampDat.chanpos(chi2,1), ...
         -exampDat.chanpos(chi2,2), ...
         exampDat.chanpos(chi2,3),'green', 'filled')
plot3([exampDat.chanpos(chi,1), exampDat.chanpos(chi2,1)], ...
      [-exampDat.chanpos(chi,2), -exampDat.chanpos(chi2,2)], ...
      [exampDat.chanpos(chi,3), exampDat.chanpos(chi2,3)], ...
      'color', 'red', 'linewidth', 3)

for ii = 1:size(missBin,1)
    for jj = 1:size(missBin,1)
        if missBin(ii,jj)
            plot3([exampDat.chanpos(ii,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(ii,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(ii,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])

set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 

export_fig([figDat 'pubFigs/' missName '.jpg'], '-r300')

end

%% EXAMPLE 2! HIP PHG

hip_enc_HFB_files = {'hip_acc_ret_HFB_2_1.mat',...
                      'hip_acc_ret_HFB_2_2.mat',...
                      'hip_acc_ret_HFB_2_3.mat',...
                      'hip_acc_ret_HFB_2_5.mat',...
                      'hip_acc_ret_HFB_2_7.mat',...
                      'hip_acc_ret_HFB_2_9.mat',...
                      'hip_acc_ret_HFB_2_11.mat',...
                      'hip_acc_ret_HFB_2_13.mat',...
                      'hip_acc_ret_HFB_2_15.mat',...
                      'hip_acc_ret_HFB_2_17.mat',...
                      'hip_acc_ret_HFB_2_19.mat',...
                      'hip_acc_ret_HFB_2_21.mat',...
                      'hip_acc_ret_HFB_2_23.mat',...
                      'hip_acc_ret_HFB_2_25.mat',...
                      'hip_acc_ret_HFB_2_26.mat',...
                      'hip_acc_ret_HFB_2_27.mat',... %using an image example! 
                      'hip_acc_ret_HFB_3_1.mat'};

for ii = 1:length(hip_enc_HFB_files)
% ii = 12 is used in figure
exampDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\graphAnalysis\out\'...
                 hip_enc_HFB_files{ii}]).outDat;

missName = ['Fig4_miss_' num2str(ii) '_hip_acc_ret_HFB_' exampDat.subID];
hitName = ['Fig4_hit_' num2str(ii) '_hip_acc_ret_HFB_' exampDat.subID]; 

hipModel =  stlread('G:\My Drive\GitHub\HpcAccConnectivityProject\bilatHippo.stl');
%get the 3d locations of the electrodes: 
chanDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished\chanDat_'...
                exampDat.subID '_001.mat']).chanDat;
exampDat.chanpos = chanDat.chanpos; 
set(0, 'defaultfigurewindowstyle', 'normal')


ftpath = [codePre 'fieldtrip-20230118']; 
yeo7 = ft_read_atlas(...
    fullfile(ftpath, ...
    'template/atlas/yeo/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'));
hold on 
test = yeo7.tissue>0; 
indices = zeros(1003785, 3); 
vals = zeros(1003785, 1); 
indi = 1; 
for nn = 1:256
    for jj = 1:256
        for kk = 1:256
            if test(nn,jj,kk)
                indices(indi, :) = [nn*-1+127,...
                                    kk*-1+145,...
                                    jj*-1+147]; 
                vals(indi) = yeo7.tissue(nn, jj, kk); 
                indi = indi + 1; 
            end
        end
    end
end
figure('visible', false, 'position', [50,50,1000,1000])
test = randsample(1:indi, 10000, false); 
hold off
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')
chi = exampDat.chi; 
chi2 = exampDat.chi2; 
scatter3(exampDat.chanpos(chi,1), ...
         -exampDat.chanpos(chi,2), ...
         exampDat.chanpos(chi,3),'red', 'filled')
scatter3(exampDat.chanpos(chi2,1), ...
         -exampDat.chanpos(chi2,2), ...
         exampDat.chanpos(chi2,3),'green', 'filled')
plot3([exampDat.chanpos(chi,1), exampDat.chanpos(chi2,1)], ...
      [-exampDat.chanpos(chi,2), -exampDat.chanpos(chi2,2)], ...
      [exampDat.chanpos(chi,3), exampDat.chanpos(chi2,3)], ...
      'color', 'red', 'linewidth', 3)






view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 

hitBin = exampDat.hitMat>.10; 
missBin = exampDat.missMat>.10; 

for nn = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if hitBin(nn,jj)
            plot3([exampDat.chanpos(nn,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(nn,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(nn,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end

export_fig([figDat 'pubFigs/' hitName '.jpg'], '-r300')


figure('visible', false, 'position', [50,50,1000,1000])
hold off
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(exampDat.chanpos(:,1), ...
         -exampDat.chanpos(:,2), ...
         exampDat.chanpos(:,3),'blue', 'filled')
chi = exampDat.chi; 
chi2 = exampDat.chi2; 
scatter3(exampDat.chanpos(chi,1), ...
         -exampDat.chanpos(chi,2), ...
         exampDat.chanpos(chi,3),'red', 'filled')
scatter3(exampDat.chanpos(chi2,1), ...
         -exampDat.chanpos(chi2,2), ...
         exampDat.chanpos(chi2,3),'green', 'filled')
plot3([exampDat.chanpos(chi,1), exampDat.chanpos(chi2,1)], ...
      [-exampDat.chanpos(chi,2), -exampDat.chanpos(chi2,2)], ...
      [exampDat.chanpos(chi,3), exampDat.chanpos(chi2,3)], ...
      'color', 'red', 'linewidth', 3)

for nn = 1:size(missBin,1)
    for jj = 1:size(missBin,1)
        if missBin(nn,jj)
            plot3([exampDat.chanpos(nn,1), exampDat.chanpos(jj,1)], ...
                  [-exampDat.chanpos(nn,2), -exampDat.chanpos(jj,2)], ...
                  [exampDat.chanpos(nn,3), exampDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])

set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 

export_fig([figDat 'pubFigs/' missName '.jpg'], '-r300')

end

%% example hip - PHG connections summary 


ii = 12;
exampDat = load(['R:\MSS\Johnson_Lab\dtf8829\QuestConnect\graphAnalysis\out\'...
                 hip_enc_HFB_files{ii}]).outDat;
hitVal = triu(exampDat.hitMat,1); 
missVal = triu(exampDat.missMat,1);
hitVal = hitVal(:); 
missVal = missVal(:); 
hitVal(hitVal==0) = []; 
missVal(missVal==0) = []; 
hitValSIG = hitVal; 
missValSIG = missVal; 
% figure
% histogram(hitVal, [min([hitVal; missVal]):.01:1], 'facecolor', hitCol)
% hold on 
% histogram(missVal, [min([hitVal; missVal]):.01:1], 'facecolor', missCol)
% hitVal(hitVal<0) = 0;  
% missVal(missVal<0) = 0; 

figure
histogram((hitVal - missVal))  %./ (hitVal + missVal))
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
xline(0, 'linewidth', 3, 'linestyle', '--')


%% EXAMPLE 2 IMAGE LOCKED
subID = 'SLCH010'; 
hipModel =  stlread('G:\My Drive\GitHub\HpcAccConnectivityProject\bilatHippo.stl');
ftpath = [codePre 'fieldtrip-20230118']; 
yeo7 = ft_read_atlas(...
    fullfile(ftpath, ...
    'template/atlas/yeo/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'));

test = yeo7.tissue>0; 
indices = zeros(1003785, 3); 
vals = zeros(1003785, 1); 
indi = 1; 
for ii = 1:256
    for jj = 1:256
        for kk = 1:256
            if test(ii,jj,kk)
                indices(indi, :) = [ii*-1+127,...
                                    kk*-1+145,...
                                    jj*-1+147]; 
                vals(indi) = yeo7.tissue(ii, jj, kk); 
                indi = indi + 1; 
            end
        end
    end
end

test = randsample(1:indi, 10000, false); 




subChans = dir('R:\MSS\Johnson_Lab\dtf8829\QuestConnect\CHANDAT\finished'); 
subi = cellfun(@(x) ~contains(x, subID), {subChans.name});
subChans(subi,:) = []; 

imHit = zeros(length(subChans));
imMiss = imHit; 

parfor ii = 1:length(subChans)
    ii
    chanDat = load([subChans(ii).folder '/' subChans(ii).name]).chanDat;
    imHit(ii,:) = mean(chanDat.ISPC.hit_on(:, chanDat.HFB.onMulTim>=0 &...
                           chanDat.HFB.onMulTim<=1500, 5:11, 2), [2,3]);
    imMiss(ii,:) = mean(chanDat.ISPC.miss_on(:, chanDat.HFB.onMulTim>=0 &...
                           chanDat.HFB.onMulTim<=1500, 5:11, 2), [2,3]);


end
chanDat = load([subChans(1).folder '/' subChans(1).name]).chanDat;
imHit = imHit - eye(size(imHit)); 
imMiss = imMiss - eye(size(imMiss)); 

figure('position', [50,50,1000,1000])
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(chanDat.chanpos(:,1), ...
         -chanDat.chanpos(:,2), ...
         chanDat.chanpos(:,3),'blue', 'filled')


view([250,-50, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 


hitBin = imHit>.1; 
missBin = imMiss>.1; 

for ii = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if hitBin(ii,jj)
            plot3([chanDat.chanpos(ii,1), chanDat.chanpos(jj,1)], ...
                  [-chanDat.chanpos(ii,2), -chanDat.chanpos(jj,2)], ...
                  [chanDat.chanpos(ii,3), chanDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
% export_fig([figDat 'pubFigs/' 'Fig4_' 'ExampleConnectionsHIT_image.jpg'], '-r300')


figure('position', [50,50,1000,1000])
scatter3(indices(test,1), indices(test,2), indices(test,3), 1, 'k', 'MarkerEdgeAlpha', .3)
hold on 
trimesh(hipModel, 'facecolor', 'none', 'facealpha', .00, ...
    'edgecolor', regColors(3,:), 'linewidth', .01, 'linestyle', ':')
scatter3(chanDat.chanpos(:,1), ...
         -chanDat.chanpos(:,2), ...
         chanDat.chanpos(:,3),'blue', 'filled')
% chi = chanDat.chi; 
% chi2 = chanDat.chi2; 
% scatter3(chanDat.chanpos(chi,1), ...
%          -chanDat.chanpos(chi,2), ...
%          chanDat.chanpos(chi,3),'red', 'filled')
% scatter3(chanDat.chanpos(chi2,1), ...
%          -chanDat.chanpos(chi2,2), ...
%          chanDat.chanpos(chi2,3),'green', 'filled')
% plot3([chanDat.chanpos(chi,1), chanDat.chanpos(chi2,1)], ...
%       [-chanDat.chanpos(chi,2), -chanDat.chanpos(chi2,2)], ...
%       [chanDat.chanpos(chi,3), chanDat.chanpos(chi2,3)], ...
%       'color', 'red', 'linewidth', 3)

view([250,-50, 140])
axis([-100, 100, -90, 90, -80, 80])
set(gcf,'color','w');
box off;
set(gca, 'color', 'none');
axis off; 


hitBin = imHit>.1; 
missBin = imMiss>0.1; 

for ii = 1:size(hitBin,1)
    for jj = 1:size(hitBin,1)
        if missBin(ii,jj)
            plot3([chanDat.chanpos(ii,1), chanDat.chanpos(jj,1)], ...
                  [-chanDat.chanpos(ii,2), -chanDat.chanpos(jj,2)], ...
                  [chanDat.chanpos(ii,3), chanDat.chanpos(jj,3)], ...
                  'color', 'k')
        end

    end
end
view([250,-130, 140])
axis([-100, 100, -90, 90, -80, 80])
% export_fig([figDat 'pubFigs/' 'Fig4_' 'ExampleConnectionsMISS_image.jpg'], '-r300')


%% example image locked full graph connection 

hitVal = triu(imHit,1); 
missVal = triu(imMiss,1);
hitVal = hitVal(:); 
missVal = missVal(:); 
hitVal(hitVal==0) = []; 
missVal(missVal==0) = []; 
% hitVal(hitVal<0) = 0; 
% missVal(missVal<0) = 0; 


missValz = zscore(missVal); 
hitValz = zscore(hitVal);
missValSIGz = zscore(missValSIG); 
hitValSIGz = zscore(hitValSIG); 


figure('position', [50,50,500,500])
scatter(missVal, hitVal, 'filled', 'markeredgealpha', .1, ...
    'markerfacealpha', .1,'markeredgecolor', '#66023C',...
    'markerfacecolor', '#66023C')
hold on 
plot([-.1, .8], [-.1, .8], 'color', 'k', 'linewidth', 2)
ylim([-.1, .8])
xlim([-.1, .8])
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
export_fig([figDat 'pubFigs/' 'Fig4_' 'hitMiss_overallScatter.jpg'], '-r300')
corr(hitVal, missVal)

figure('position', [50,50,500,500])
scatter(missValSIG, hitValSIG, 'filled', 'markeredgealpha', .1, ...
    'markerfacealpha', .1,'markeredgecolor', '#66023C',...
    'markerfacecolor', '#66023C')
hold on 
plot([-.1, .8], [-.1, .8], 'color', 'k', 'linewidth', 2)
ylim([-.1, .8])
xlim([-.1, .8])
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
export_fig([figDat 'pubFigs/' 'Fig4_' 'hitMiss_SIGScatter.jpg'], '-r300')
corr(hitValSIG, missValSIG)
% 
% figure
% scatter(missValz, missValSIGz, 'markeredgealpha', .3)
% hold on 
% plot([-.1, .8], [-.1, .8])
% figure
% scatter(hitValz, hitValSIGz, 'markeredgealpha', .1)
% hold on 
% plot([-.1, .8], [-.1, .8])

% figure
% histogram(hitValSIG - hitVal, [-1.6:.025:1.7], 'facecolor', hitCol, ...
%     'normalization', 'probability')
% hold on 
% histogram(missValSIG - missVal, [-1.6:.025:1.7], 'facecolor', missCol, ...
%     'normalization', 'probability')

figure
histogram(hitVal, [-.1:.01:.2], 'facecolor', hitCol)
hold on 
histogram(missVal, [-.1:.01:.2], 'facecolor', missCol)
xline(0, 'linewidth', 3, 'linestyle', '--')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;

export_fig([figDat 'pubFigs/' 'Fig4_' 'longTimeHist.jpg'], '-r300')

figure
histogram(hitVal, [-.1:.02:.7], 'facecolor', hitCol)
hold on 
histogram(missVal, [-.1:.02:.7], 'facecolor', missCol)
xline(0, 'linewidth', 3, 'linestyle', '--')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
xlim([.025, .7])
ylim([0,10])
export_fig([figDat 'pubFigs/' 'Fig4_' 'longTimeHist2.jpg'], '-r300')

figure
histogram(hitVal, [-.1:.01:.7], 'facecolor', hitCol)
hold on 
histogram(missVal, [-.1:.01:.7], 'facecolor', missCol)
xline(0, 'linewidth', 3, 'linestyle', '--')
set(gcf,'color','w');
box off;
ax=gca;ax.LineWidth=4;
xlim([-.025, .025])
export_fig([figDat 'pubFigs/' 'Fig4_' 'longTimeHist3.jpg'], '-r300')


hitVal = hitVal - min([hitVal; missVal]); 
missVal = missVal - min([missVal; hitVal]);

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
















