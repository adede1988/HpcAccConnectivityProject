%% MTL - PFC publication figures 

%description of the general data set 


%% set environment

%figures should not be docked: 
set(0, 'defaultfigurewindowstyle', 'normal')
%local paths: 

%edit here for local file paths: 
codePre = 'G:\My Drive\GitHub\';
datPre = 'R:\MSS\Johnson_Lab\dtf8829\Pubdat\';
figDat = datPre;
figSavePath = 'R:\MSS\Johnson_Lab\dtf8829\Pubdat\FiguresOut\';

% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
% addpath([codePre 'myFrequentUse'])
% addpath([codePre 'subNetworkDynamics'])
addpath([codePre 'myFrequentUse/export_fig_repo'])
% addpath([codePre 'fieldtrip-20230118'])
% ft_defaults;
regions = {'acc', 'dlPFC', 'hip', 'mtl',  'pPFC'}; 




%% set graphical parameters: 

hitCol = [88, 113, 206] ./ 256; 
missCol = [240, 172, 99] ./ 256; 
sigCol = [.5,.5,.5]; 
sigAlpha = .3; 
errAlpha = .2; 

regColors = [[204, 153, 204]', ...%ACC: dusky pink
             [145, 162, 80]', ... %dlPFC: light forest green
             [61, 187, 218]', ... %hip light blue
             [214, 77, 97]',...   %parahip: dark pink
             [51, 102, 153]']' ./ 255; %polarPFC: dark blue
            

keyRegIdx = [1,2,3,4,5]; 

b2w2r = [[linspace(0,255,128)'; linspace(255,0,128)'], ...
    [linspace(0,255,128)'; linspace(255,0,128)'], [linspace(0,255,128)';...
    linspace(255,0,128)']]/255;
b2w2r(129:end, 1) = 1; 
b2w2r(1:128, 3) = 1; 

linWid = 5; 

%color anchors for custom color map
s = [12, 61, 74]; 
m = [171,189,154]; 
e = [255, 255, 46]; 

s2w2y = [[linspace(s(1),m(1),128)'; linspace(m(1),e(1),128)'], ...
         [linspace(s(2),m(2),128)'; linspace(m(2),e(2),128)'], ...
         [linspace(s(3),m(3),128)'; linspace(m(3),e(3),128)'], ...
         ] / 255;

%color anchors for custom color map
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
for ii = 1:5
    scatter(ii, 2, 100, regColors(ii,:), 'filled')
    text(ii, 2.5, regions{ii})

end




%% Supplemental Figure 1A: image locked HFB time courses



fig1Dat = dir([datPre 'HFB_image']); 
fig1Dat(1:2) = []; 

for ii = 1:length(fig1Dat)
    ii
    panelDat = load([fig1Dat(ii).folder '/' fig1Dat(ii).name]).outDat; 

    %remove rare instances where participant had no miss trials
    panelDat.hitVals(panelDat.eliminate,:) = []; 
    panelDat.missVals(panelDat.eliminate,:) = []; 
    
    %nan values at first and last time point, zero out
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
    export_fig([figSavePath 'supFig1a_' ...
        panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')
    
end















