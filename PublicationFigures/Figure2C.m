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


 
%% Fig 2C: latencies for all regions and all trials 
%HITS Retrieve
linew = 5;

singleTrialDat = dir([datPre 'HFB_singleTrial']); 
singleTrialDat(1:2) = []; 
retFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])
%extract
for ii = 1:length(singleTrialDat)
     panelDat = load([singleTrialDat(ii).folder '/'...
            singleTrialDat(ii).name]).outDat; 
     
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
export_fig([figSavePath 'Fig2C_retHitLatencies' ...
    '.jpg'], '-r300')

% HITS Encode
encFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])

%extract
for ii = 1:length(singleTrialDat)
     panelDat = load([singleTrialDat(ii).folder '/'...
            singleTrialDat(ii).name]).outDat; 
     
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
export_fig([figSavePath 'Fig2C_encHitLatencies' ...
    '.jpg'], '-r300')


% MISSES ENCODE
encFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])

%extract
for ii = 1:length(singleTrialDat)
     panelDat = load([singleTrialDat(ii).folder '/'...
            singleTrialDat(ii).name]).outDat; 
     
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

    end

   

end
set(0, 'currentfigure', encFig);
export_fig([figSavePath 'Fig2C__encMissLatencies' ...
    '.jpg'], '-r300')


% MISSES RETRIEVE

retFig = figure('visible', true, 'position', [0,0,600,600]);
set(gca, 'ydir', 'reverse')
 set(gcf,'color','w');
    box off;
    ax=gca;ax.LineWidth=4;
    xlim([0,1])
%extract
for ii = 1:length(singleTrialDat)
     panelDat = load([singleTrialDat(ii).folder '/'...
            singleTrialDat(ii).name]).outDat; 
     
    if strcmp(panelDat.phase, 'ret')
  
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
export_fig([figSavePath 'Fig2C_retMissLatencies' ...
    '.jpg'], '-r300')















