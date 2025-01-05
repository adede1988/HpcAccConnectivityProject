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



%% Supplemental Figure 3 ABC power spectra at HFB peak


TF_files_HFB = dir([datPre 'TF_HFB/']); 
TF_files_HFB(1:2) = []; 



for ii = 1:length(TF_files_HFB)
   
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

    export_fig([figSavePath 'supFig3AB_' panelDat.reg '_HIT_' panelDat.phase  '.jpg'], '-r300')

    figure('visible', false, 'position', [0,0,600,600])

    imagesc([-500:25:500], [], squeeze(mean(panelDat.misses_hfb,1))')
    yticks([10:10:100])
    yticklabels(round(panelDat.frex([10:10:100])))
    set(gca, 'ydir', 'normal')
    set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    caxis([0, 10])

    export_fig([figSavePath 'supFig3AB_' panelDat.reg '_MISS_' panelDat.phase  '.jpg'], '-r300')


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
    export_fig([figSavePath 'supFig3C_' panelDat.reg '_' panelDat.phase  '.jpg'], '-r300')



% Get peak frequencies for description of power spectra 
    curSpect = mean([missSpect; hitSpect]);
    
    [val, idx] = max(curSpect(1:64));
    
    disp(['summary for ' panelDat.reg ' ' panelDat.phase])
    disp(['low peak: ' num2str(round(panelDat.frex(idx),2))])

    rangeBot = find(diff(movmean(curSpect(idx+4:end), 4))>0, 1) + idx + 2;

    [val, idx] = max(curSpect(rangeBot:77));

    if val>(curSpect(rangeBot)+1) && val>curSpect(77) %edges aren't peaks
        disp(['high peak: ' num2str(round(panelDat.frex(idx+rangeBot-1),2))])
    end

  
end



















