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


%% Figure 3A and Supplemental Figure 4A ITPC spectra at HFB peak 

TF_files_HFB = dir([datPre 'TFphase_HFB/']); 
TF_files_HFB(1:2) = []; 


for ii = 1:length(TF_files_HFB)

    panelDat = load([TF_files_HFB(ii).folder '/' ...
        TF_files_HFB(ii).name]).outDat;
 

    regi = find(cellfun(@(x) strcmp(x, panelDat.reg), regions)); 
    phasei = find(cellfun(@(x) strcmp(x, panelDat.phase), phaseVals)); 
    

    
    figure('visible', false, 'position', [0,0,600,600])
    x = 1:length(panelDat.frex); 
    yline(0, 'color', 'k', 'linewidth', linWid)
   
    x(panelDat.p_hfb(21,:)>.05) = []; 
    if ~isempty(x) %check if we have any sig time points
        disp(['Max significantly different frequency for ' ...
            panelDat.reg ' ' panelDat.phase ' is: ' ...
            num2str(panelDat.frex(max(x)))])
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
    export_fig([figSavePath 'Fig3A_' panelDat.reg '_' ...
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

    export_fig([figSavePath 'supFig4A_' panelDat.reg '_' ...
        panelDat.phase  '_phase_Hist.png'], '-transparent', '-r300')

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

    export_fig([figSavePath 'supFig4A_' panelDat.reg '_' ...
        panelDat.phase  '_phase_Hist_high.png'], '-transparent', '-r300')


end
















